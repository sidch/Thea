//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Princeton University
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#include "JointBoost.hpp"
#include "../Math.hpp"
#include "../Serializable.hpp"
#include "../UnorderedSet.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>

// #define JOINT_BOOST_TEST_ALL_THRESHOLDS

namespace Thea {
namespace Algorithms {

std::string
JointBoost::SharedStump::toString() const
{
  std::ostringstream oss;
  oss << "Stump[f: " << f << ", n: " << n << ']';
  return oss.str();
}

JointBoost::Options::Options()
: min_boosting_rounds(-1),
  max_boosting_rounds(-1),
  min_fractional_error_reduction(-1),
  feature_sampling_fraction(-1),
  max_thresholds_fraction(-1),
  force_exhaustive(false),
  force_greedy(false),
  verbose(true)
{
}

JointBoost::JointBoost(long num_classes_, long num_features_, Options const & options_)
: num_classes(num_classes_),
  num_features(num_features_),
  options(options_)
{
  alwaysAssertM(num_classes_ >= 2, "JointBoost: At least two classes required");
  alwaysAssertM(num_features_ >= 1, "JointBoost: At least one feature required");
}

JointBoost::JointBoost(std::string const & path)
: num_classes(0), num_features(0)
{
  if (!load(path))
    throw Error("JointBoost: Could not load classifier from '" + path + '\'');
}

JointBoost::~JointBoost()
{
  clear();
}

void
JointBoost::clear()
{
  weights.resize(0, 0);
  stumps.clear();
}

std::string
JointBoost::getClassName(long i) const
{
  if (class_names.empty())
  {
    std::ostringstream rout;
    rout << i;
    return rout.str();
  }
  else
    return class_names[(array_size_t)i];
}

long
JointBoost::train(TrainingData const & training_data_, TrainingData const * validation_data_)
{
  alwaysAssertM(training_data_.numExamples() > 0, "JointBoost: No training examples provided");
  alwaysAssertM(num_features == training_data_.numFeatures(), "JointBoost: Unexpected number of features in training data");

  // Create a temporary reference to the training data so we don't have to pass it to every function
  training_data = &training_data_;

  // Set values of parameters
  long min_rounds = options.min_boosting_rounds;
  long max_rounds = options.max_boosting_rounds;

  if (min_rounds <= 0)
    min_rounds = 1;

  if (max_rounds < 0)
    max_rounds = num_classes;

  if (max_rounds < min_rounds)
    max_rounds = min_rounds;

  double min_fractional_error_reduction = options.min_fractional_error_reduction;
  if (min_fractional_error_reduction < 0)
    min_fractional_error_reduction = 0.001;

  feature_sampling_fraction = options.feature_sampling_fraction;
  if (feature_sampling_fraction <= 0)
    feature_sampling_fraction = 1.0 / max_rounds;

  max_thresholds_fraction = options.max_thresholds_fraction;
  if (max_thresholds_fraction <= 0)
    max_thresholds_fraction = 1.0;  // test all features as thresholds by default

  // Get class names from training data, if available
  training_data->getClassNames(class_names);
  alwaysAssertM(class_names.empty() || (long)class_names.size() == num_classes,
                "JointBoost: Incorrect number of class names specified");

  // Get classes for training data
  TheaArray<long> classes;
  training_data->getClasses(classes);

  // Initialize the weights
  TheaArray<double> training_weights;
  training_data->getWeights(training_weights);
  weights.resize(num_classes, (long)classes.size());
  if (!training_weights.empty())
  {
    for (array_size_t i = 0; i < training_weights.size(); ++i)
      for (long c = 0; c < num_classes; ++c)
        weights(c, (long)i) = training_weights[i];
  }
  else
    weights.fill(1);

  if (options.verbose)
  {
    THEA_CONSOLE << "JointBoost: Starting training (" << num_classes << " classes, " << num_features << " features, "
                 << training_data_.numExamples() << " examples)";
    THEA_CONSOLE << "JointBoost:     -- min_rounds = " << min_rounds;
    THEA_CONSOLE << "JointBoost:     -- max_rounds = " << max_rounds;
    THEA_CONSOLE << "JointBoost:     -- min_fractional_error_reduction = " << min_fractional_error_reduction;
    THEA_CONSOLE << "JointBoost:     -- feature_sampling_fraction = " << feature_sampling_fraction;
    THEA_CONSOLE << "JointBoost:     -- max_thresholds_fraction = " << max_thresholds_fraction;
    THEA_CONSOLE << "JointBoost:     -- force_exhaustive = " << options.force_exhaustive;
    THEA_CONSOLE << "JointBoost:     -- force_greedy = " << options.force_greedy;
    THEA_CONSOLE << "JointBoost:     -- verbose = " << options.verbose;
  }

  // Do several rounds of boosting, adding a new stump (weak learner) in each round
  stumps.clear();
  double error = -1;
  double validation_error = -1;
  for (long round = 0; round < max_rounds; ++round)
  {
    if (options.verbose) THEA_CONSOLE << "JointBoost: Training round " << round;

    // Create the stump for this round
    SharedStump::Ptr stump(new SharedStump);

    // Compute the k values for the stump
    stump->k.resize((array_size_t)num_classes, 0);
    for (long c = 0; c < num_classes; ++c)
    {
      double k_numer = 0, k_denom = 0;
      for (array_size_t i = 0; i < classes.size(); ++i)
      {
        double w = weights(c, (long)i);
        int z = (classes[i] == c ? +1 : -1);

        k_numer += (w * z);
        k_denom += w;
      }

      stump->k[(array_size_t)c] = (k_denom != 0 ? k_numer / k_denom : 0);
    }

    // Optimize this stump
    double round_err = optimizeStump(*stump, classes);
    double round_validation_err = -1;

    if (options.verbose)
      THEA_CONSOLE << "JointBoost:     Error in round " << round << " = " << round_err << " with " << stump->toString();

    if (round_err <= 0)
      break;

    // Check if we should stop at this point
    if (round >= min_rounds)
    {
      // Check if error reduction is below threshold
      if (error >= 0)
      {
        double err_reduction = error - round_err;
        if (err_reduction < min_fractional_error_reduction * error)
        {
          if (options.verbose)
            THEA_CONSOLE << "JointBoost: Stopping training, fractional error reduction less than threshold";

          break;
        }
      }

      // Check if error has increased on validation set
      if (validation_data_)
      {
        round_validation_err = computeValidationError(*validation_data_, stump);

        if (validation_error >= 0 && round_validation_err > validation_error)
        {
          if (options.verbose)
            THEA_CONSOLE << "JointBoost: Stopping training, validation error increased from " << validation_error << " to "
                         << round_validation_err;

          break;
        }
      }
    }

    // We have new error values
    error = round_err;
    validation_error = round_validation_err;

    // Add the stump to the strong classifier
    stumps.push_back(stump);

    // Update weights
    TheaArray<double> stump_features;
    training_data->getFeature(stump->f, stump_features);

    for (long c = 0; c < num_classes; ++c)
      for (array_size_t i = 0; i < classes.size(); ++i)
      {
        // THEA_CONSOLE << "      weights[c = " << c << ", i = " << i << "] = " << weights(c, (long)i);

        int z = (classes[i] == c ? +1 : -1);
        double h = (*stump)(stump_features[i], c);

        // THEA_CONSOLE << "      h[i = " << i << ", c = " << c << "] = " << h << ", z = " << z;

        weights(c, (long)i) *= std::exp(-z * h);

        // THEA_CONSOLE << "      weights[c = " << c << ", i = " << i << "] = " << weights(c, (long)i);
      }
  }

  training_data = NULL;

  if (options.verbose)
  {
    if (validation_data_)
    {
      if (validation_error < 0)
        validation_error = computeValidationError(*validation_data_);

      THEA_CONSOLE << "JointBoost: Completed training, added " << stumps.size() << " stump(s) with final error " << error
                   << ", validation error " << validation_error;
    }
    else
      THEA_CONSOLE << "JointBoost: Completed training, added " << stumps.size() << " stump(s) with final error " << error;
  }

  return (long)stumps.size();
}

double
JointBoost::optimizeStump(SharedStump & stump, TheaArray<long> const & stump_classes)
{
  if (num_features <= 0)
  {
    THEA_WARNING << "JointBoost:     Can't optimize stump with " << num_features << " features";
    return -1;
  }

  // Select a random fraction of the features
  TheaArray<int32> candidate_features;
  if (feature_sampling_fraction >= 1)
  {
    for (long f = 0; f < num_features; ++f)
      candidate_features.push_back(f);
  }
  else
  {
    int32 num_candidate_features = std::max((int32)Math::round(feature_sampling_fraction * num_features), (int32)1);
    candidate_features.resize((array_size_t)num_candidate_features);
    Random::common().integers(0, (int32)num_features - 1, num_candidate_features, &candidate_features[0]);
  }

  if (options.verbose && (long)candidate_features.size() < num_features)
    THEA_CONSOLE << "JointBoost:     Optimizing over " << candidate_features.size() << " randomly selected feature(s)";

  SharedStump test = stump;
  double min_err = -1;
  for (array_size_t f = 0; f < candidate_features.size(); ++f)
  {
    test.f = candidate_features[f];

    // Get feature values
    TheaArray<double> stump_features;
    training_data->getFeature(test.f, stump_features);

    // Optimize the stump via an exhaustive or a greedy search over subsets of classes
    double err = 0;
    if ((num_classes <= 5 && !options.force_greedy) || options.force_exhaustive)
      err = optimizeStumpExhaustive(test, stump_features, stump_classes);
    else
      err = optimizeStumpGreedy(test, stump_features, stump_classes);

    // THEA_CONSOLE << "JointBoost:     Error when splitting on feature " << test.f << " = " << err;

    // Did we improve the classification?
    if (err >= 0 && (min_err < 0 || err < min_err))
    {
      min_err = err;
      stump = test;
    }
  }

  return min_err;
}

double
JointBoost::optimizeStumpExhaustive(SharedStump & stump, TheaArray<double> const & stump_features,
                                    TheaArray<long> const & stump_classes)
{
  // Loop over all possible subsets of classes
  SharedStump test = stump;
  unsigned long num_subsets = (1L << num_classes);  // assume no more classes than an unsigned long can hold
  double min_err = -1;
  double cum_num_thresholds = 0;
  for (unsigned long subset = 1; subset + 1 < num_subsets; ++subset)
  {
    test.n = SharingSet((SharingSet::size_type)num_classes, subset);

    // THEA_CONSOLE << "JointBoost:         Testing stump " << test.toString();

    long num_thresholds = 0;
    double err = fitStump(test, stump_features, stump_classes, &num_thresholds);
    if (err >= 0 && (min_err < 0 || err < min_err))
    {
      stump = test;
      min_err = err;
    }

    cum_num_thresholds += num_thresholds;
  }

  if (options.verbose)
  {
    THEA_CONSOLE << "JointBoost:     Generated " << cum_num_thresholds / (num_subsets - 2)
                 << " candidate feature threshold(s) on average";
  }

  return min_err;
}

double
JointBoost::optimizeStumpGreedy(SharedStump & stump, TheaArray<double> const & stump_features,
                                TheaArray<long> const & stump_classes)
{
  TheaArray<SharedStump> candidate_stumps;
  TheaArray<double> candidate_errors;

  SharedStump test = stump;
  SharingSet current_n((SharingSet::size_type)num_classes, 0);  // initially empty
  double cum_num_thresholds = 0;
  long num_fitted_stumps = 0;
  for (long c = 0; c < num_classes - 1; ++c)
  {
    // Add the single new class to the current set that jointly gives the minimum error
    SharedStump best_stump;
    double min_err = -1;
    bool min_found = false;
    for (long c = 0; c < num_classes - 1; ++c)
    {
      if (current_n[(SharingSet::size_type)c])
        continue;

      test.n = current_n;
      test.n[(SharingSet::size_type)c] = true;

      // THEA_CONSOLE << "JointBoost:         Testing stump " << test.toString();

      long num_thresholds = 0;
      double err = fitStump(test, stump_features, stump_classes, &num_thresholds);
      if (err >= 0 && (min_err < 0 || err < min_err))
      {
        best_stump = test;
        min_err = err;
        min_found = true;
      }

      cum_num_thresholds += num_thresholds;
      num_fitted_stumps++;
    }

    if (!min_found)  // this shouldn't normally happen, it means every attempt to fit gave an error
      break;

    candidate_stumps.push_back(best_stump);
    candidate_errors.push_back(min_err);

    // THEA_CONSOLE << "JointBoost:         Greedy search candidate " << best_stump.toString() << " with error " << min_err;

    current_n = best_stump.n;
  }

  double min_err = -1;
  for (array_size_t i = 0; i < candidate_stumps.size(); ++i)
  {
    if (min_err < 0 || candidate_errors[i] < min_err)
    {
      stump = candidate_stumps[i];
      min_err = candidate_errors[i];
    }
  }

  if (options.verbose)
  {
    THEA_CONSOLE << "JointBoost:     Generated " << cum_num_thresholds / num_fitted_stumps
                 << " candidate feature threshold(s) on average";
  }

  return min_err;
}

namespace JointBoostInternal {

struct IndexedComparator
{
  IndexedComparator(double const * values_) : values(values_) {}
  bool operator()(array_size_t i, array_size_t j) const { return values[i] < values[j]; }

  double const * values;
};

void
sortIndexed(TheaArray<double> const & values, TheaArray<array_size_t> & sorted_indices)
{
  sorted_indices.resize(values.size());
  for (array_size_t i = 0; i < sorted_indices.size(); ++i)
    sorted_indices[i] = i;

  std::sort(sorted_indices.begin(), sorted_indices.end(), IndexedComparator(&values[0]));
}

void
getClassificationAccuracy(double split_value, SharingSet const & pos_classes, TheaArray<double> const & features,
                          TheaArray<long> const & classes, TheaArray<array_size_t> const & sorted_indices,
                          TheaArray<int> & accuracy)
{
  accuracy.resize(features.size());
  for (array_size_t i = 0; i < features.size(); ++i)
  {
    array_size_t index = sorted_indices[i];

    long c = classes[index];
    double f = features[index];
    bool is_positive = pos_classes[(SharingSet::size_type)c];
    bool correct = ((is_positive && f > split_value) || (!is_positive && f <= split_value));
    accuracy[index] = (correct ? +1 : -1);
  }
}

double
splitQuality(Matrix<double> const & weights, TheaArray<long> const & classes, TheaArray<int> const & accuracy)
{
  // We'll define the split quality as
  //
  //   sum { weights[classes[i], i] * accuracy[i] }
  //
  // where accuracy[i] is +1 if example i has a positive class and its feature above the split index, or has a negative class
  // and is below the split index, and -1 otherwise.
  //
  // As mentioned in getCandidateThresholds(), the absolute value of the split quality is relevant, not its sign.

  double quality = 0;
  for (array_size_t i = 0; i < accuracy.size(); ++i)
    quality += (accuracy[i] * weights(classes[i], (long)i));

  return quality;
}

void
getCandidateThresholds(SharingSet const & pos_classes, Matrix<double> const & weights, TheaArray<double> const & features,
                       TheaArray<long> const & classes, long max_thresholds, TheaArray<double> & thresholds)
{
  // A good threshold separates the positive classes from the negative classes. In other words, we want either many positive
  // examples > theta and many negative examples <= theta, or vice versa. Note that the positive examples need *not* be mostly
  // on the "positive side of" (greater than) theta -- we just require theta to separate the positive and negative examples
  // as well as possible. (Hence, the absolute value of splitQuality() is relevant, not its sign.)

  // Sort the examples by feature value
  TheaArray<array_size_t> sorted_indices;
  sortIndexed(features, sorted_indices);

  // Get the minimum and maximum feature values
  double lo = features[0], hi = features[0];
  for (array_size_t i = 1; i < features.size(); ++i)
  {
    if (features[i] < lo) lo = features[i];
    if (features[i] > hi) hi = features[i];
  }

  // Choose a tolerance in comparing features
  double tolerance = 1.0e-10 * (hi - lo);

  // Generate samples and adjust their position to be better cuts
  typedef TheaUnorderedSet<array_size_t> IndexSet;
  IndexSet threshold_set;
  for (long t = 0; t < max_thresholds; ++t)
  {
    // Start with a random seed
    array_size_t index = (array_size_t)Random::common().integer(0, (int32)sorted_indices.size() - 1);
    array_size_t mapped_index = sorted_indices[index];
    double threshold = features[mapped_index];

    // Measure the classification accuracy (+/-1) per example, and the overall quality of the split
    TheaArray<int> accuracy;
    getClassificationAccuracy(threshold, pos_classes, features, classes, sorted_indices, accuracy);
    double quality = splitQuality(weights, classes, accuracy);

    // Iteratively try to improve it by moving one step to the right or left
    for (array_size_t round = 0; round < sorted_indices.size(); ++round)
    {
      // Get the quality change in each direction
      double offset_qualities[2] = { 0, 0 };
      bool equal_quality[2] = { false, false };
      for (int offset = -1; offset <= 1; offset += 2)
      {
        if ((offset < 0 && index == 0) || (offset > 0 && index == sorted_indices.size() - 1))
          continue;

        array_size_t offset_index = (array_size_t)(index + offset);
        array_size_t mapped_offset_index = sorted_indices[offset_index];

        // THEA_CONSOLE << "index = " << index << ", mapped_offset_index = " << mapped_offset_index;

        double offset_quality = quality;

        if (Math::fuzzyEq(features[mapped_offset_index], features[mapped_index], tolerance))
        {
          equal_quality[(offset + 1) >> 1] = true;
        }
        else
        {
          if (offset < 0)  // move left
          {
            // The current example's classification accuracy is inverted. The offset example remains unchanged.
            long c = classes[mapped_index];
            offset_quality -= 2 * (accuracy[mapped_index] * weights(c, (long)mapped_index));
          }
          else // move right
          {
            // The offset example's classification accuracy is inverted. The current example remains unchanged.
            long c = classes[mapped_offset_index];
            offset_quality -= 2 * (accuracy[mapped_offset_index] * weights(c, (long)mapped_offset_index));
          }
        }

        offset_qualities[(offset + 1) >> 1] = offset_quality;
      }

      // Which is the better direction to move?
      int best_offset = (std::fabs(offset_qualities[0]) > std::fabs(offset_qualities[1]) ? 0 : 1);
      if (std::fabs(offset_qualities[best_offset]) > std::fabs(quality)
       || (equal_quality[best_offset] && Random::common().uniform01() < 0.5))  // if equal, move with 50% probability
      {
        // If move left, current example's accuracy is inverted
        if (best_offset == 0)
          accuracy[mapped_index] = -accuracy[mapped_index];

        // Move to index + offset
        index = (array_size_t)(index + 2 * best_offset - 1);
        mapped_index = sorted_indices[index];
        quality = offset_qualities[best_offset];

        // If move right, offset example's accuracy is inverted
        if (best_offset == 1)
          accuracy[mapped_index] = -accuracy[mapped_index];
      }
      else if (!equal_quality[best_offset])
        break;
    }

    threshold_set.insert(mapped_index);
  }

  thresholds.clear();
  for (IndexSet::const_iterator ti = threshold_set.begin(); ti != threshold_set.end(); ++ti)
    thresholds.push_back(features[*ti]);
}

} // namespace JointBoostInternal

double
JointBoost::fitStump(SharedStump & stump, TheaArray<double> const & stump_features, TheaArray<long> const & stump_classes,
                     long * num_generated_thresholds)
{
  using namespace JointBoostInternal;

  // Find good places to cut the features (candidate values of the threshold theta)
  TheaArray<double> thresholds;

#ifdef JOINT_BOOST_TEST_ALL_THRESHOLDS
  thresholds = stump_features;
#else
  long max_thresholds = (long)std::ceil(max_thresholds_fraction * stump_features.size());
  if (max_thresholds < (long)stump_features.size())
    getCandidateThresholds(stump.n, weights, stump_features, stump_classes, max_thresholds, thresholds);
  else
    thresholds = stump_features;
#endif

  if (num_generated_thresholds)
    *num_generated_thresholds = (long)thresholds.size();

  // THEA_CONSOLE << "JointBoost: Generated " << thresholds.size() << " candidate threshold(s) for " << stump.toString();
  // for (array_size_t t = 0; t < thresholds.size(); ++t)
  //   THEA_CONSOLE << "JointBoost:     -- threshold[" << t << "] = " << thresholds[t];

  // Precalculate the error term depending on stump.k -- it does not depend on threshold selection
  double err_k = 0;
  for (long c = 0; c < num_classes; ++c)
  {
    if (!stump.n[(SharingSet::size_type)c])
    {
      for (array_size_t i = 0; i < stump_classes.size(); ++i)
      {
        double w = weights(c, (long)i);
        int z = (stump_classes[i] == c ? +1 : -1);

        err_k += w * Math::square(z - stump.k[(array_size_t)c]);  // k is precalculated outside this function
      }
    }
  }

  // THEA_CONSOLE << "err_k = " << err_k;

  // Fit the stump parameters to the observed data
  long best_threshold = -1;
  double min_err = -1;
  for (array_size_t t = 0; t < thresholds.size(); ++t)
  {
    double theta = thresholds[t];
    double a_numer = 0, a_denom = 0, b_numer = 0, b_denom = 0;
    for (long c = 0; c < num_classes; ++c)
    {
      if (stump.n[(SharingSet::size_type)c])
      {
        for (array_size_t i = 0; i < stump_classes.size(); ++i)
        {
          double w = weights(c, (long)i);
          int z = (stump_classes[i] == c ? +1 : -1);

          if (stump_features[i] > theta)
          {
            a_numer += w * z;
            a_denom += w;
          }
          else
          {
            b_numer += w * z;
            b_denom += w;
          }
        }
      }
    }

    double a = (a_denom != 0 ? a_numer / a_denom : 0);
    double b = (b_denom != 0 ? b_numer / b_denom : 0);

    double err = (1 - Math::square(a)) * a_denom + (1 - Math::square(b)) * b_denom + err_k;

    if (best_threshold < 0 || err < min_err)
    {
      best_threshold = (long)t;
      min_err = err;

      stump.a = a;
      stump.b = b;
      stump.theta = theta;
    }
  }

  return min_err;
}

long
JointBoost::predict(double const * features, double * class_probabilities) const
{
  long best_class = -1;
  double best_H = -1;
  double sum_probs = 0;

  for (long c = 0; c < num_classes; ++c)
  {
    double H = 0;
    for (array_size_t m = 0; m < stumps.size(); ++m)
    {
      SharedStump const & stump = *stumps[m];
      double h = stump(features[stump.f], c);
      H += h;
    }

    if (class_probabilities)
    {
      double p = std::exp(H);  // softmax transform, see TextonBoost paper and Friedman/Hastie/Tibshirani
      class_probabilities[c] = p;
      sum_probs += p;
    }

    if (best_class < 0 || H > best_H)
    {
      best_class = c;
      best_H = H;
    }
  }

  // Normalize probabilities
  if (class_probabilities && sum_probs > 0)
  {
    for (long c = 0; c < num_classes; ++c)
      class_probabilities[c] /= sum_probs;
  }

  return best_class;
}

double
JointBoost::computeValidationError(TrainingData const & validation_data_, SharedStump::Ptr new_stump)
{
  alwaysAssertM(validation_data_.numFeatures() == num_features, "JointBoost: Validation set has different number of features");

  if (new_stump)
    stumps.push_back(new_stump);

  Matrix<double, MatrixLayout::ROW_MAJOR> validation_features(validation_data_.numExamples(), num_features);
  TheaArray<double> feat;
  for (long i = 0; i < num_features; ++i)
  {
    validation_data_.getFeature(i, feat);
    validation_features.setColumn(i, &feat[0]);
  }

  TheaArray<long> validation_classes;
  TheaArray<double> validation_weights;

  validation_data_.getClasses(validation_classes);
  validation_data_.getWeights(validation_weights);

  alwaysAssertM((long)validation_classes.size() == validation_features.numRows(),
                "JointBoost: Ground truth classes for validation data do not match number of examples");

  double err = 0;
  for (long i = 0; i < validation_features.numRows(); ++i)
  {
    long predicted = predict(&validation_features.getMutable(i, 0));
    if (predicted != validation_classes[i])
      err += (validation_weights.empty() ? 1.0 : validation_weights[i]);
  }

  if (new_stump)
    stumps.pop_back();

  return err;
}

bool
JointBoost::Options::load(std::string const & path)
{
  try
  {
    // FIXME: This is currently needed because TextInputStream crashes if the file cannot be read
    // UPDATE: Was this fixed by importing the G3D I/O classes into Thea?
    {
      std::ifstream in(path.c_str());
      if (!in)
      {
        THEA_ERROR << "JointBoost: Could not load options from input file '" << path << '\'';
        return false;
      }
    }

    TextInputStream in(path, Serializable::configReadSettings());
    return deserialize(in);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "JointBoost: Could not load options from input file '%s'", path.c_str())
}

bool
JointBoost::Options::save(std::string const & path) const
{
  try
  {
    // FIXME: This is currently needed because TextOutputStream may crash (?) if the file cannot be written (going by the
    // corresponding crash for TextInputStream)
    // UPDATE: Was this fixed by importing the G3D I/O classes into Thea?
    {
      std::ofstream out(path.c_str());
      if (!out)
      {
        THEA_ERROR << "JointBoost: Could not save options to output file '" << path << '\'';
        return false;
      }
    }

    TextOutputStream out(path, Serializable::configWriteSettings());
    if (!serialize(out))
      return false;

    out.commit();
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "JointBoost: Could not save options to output file '%s'", path.c_str())

  return true;
}

bool
JointBoost::load(std::string const & path)
{
  std::ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "JointBoost: Could not load classifier from file '" << path << '\'';
    return false;
  }

  if (!options.deserialize(in))
    return false;

  if (!deserialize(in))
    return false;

  return true;
}

bool
JointBoost::save(std::string const & path) const
{
  std::ofstream out(path.c_str());
  if (!out)
  {
    THEA_ERROR << "JointBoost: Could not save classifier to file '" << path << '\'';
    return false;
  }

  if (!options.serialize(out))
    return false;

  out << std::endl;  // break up the output a bit

  if (!serialize(out))
    return false;

  return true;
}

bool
JointBoost::Options::deserialize(std::istream & in)
{
  in >> min_boosting_rounds
     >> max_boosting_rounds
     >> min_fractional_error_reduction
     >> feature_sampling_fraction
     >> max_thresholds_fraction
     >> force_exhaustive
     >> force_greedy
     >> verbose;

  if (!in)
  {
    THEA_ERROR << "JointBoost: Error reading options";
    return false;
  }
  else
    return true;
}

bool
JointBoost::Options::deserialize(TextInputStream & input)
{
  *this = defaults();

  while (input.hasMore())
  {
    std::string field = input.readSymbol();
    input.readSymbol("=");

    if (field == "min_boosting_rounds")
      min_boosting_rounds = (long)input.readNumber();
    else if (field == "max_boosting_rounds")
      max_boosting_rounds = (long)input.readNumber();
    else if (field == "min_fractional_error_reduction")
      min_fractional_error_reduction = input.readNumber();
    else if (field == "feature_sampling_fraction")
      feature_sampling_fraction = input.readNumber();
    else if (field == "max_thresholds_fraction")
      max_thresholds_fraction = input.readNumber();
    else if (field == "force_exhaustive")
      force_exhaustive = input.readBoolean();
    else if (field == "force_greedy")
      force_greedy = input.readBoolean();
    else if (field == "verbose")
      verbose = input.readBoolean();
  }

  return true;
}

bool
JointBoost::Options::serialize(std::ostream & out) const
{
  out << min_boosting_rounds << '\n'
      << max_boosting_rounds << '\n'
      << min_fractional_error_reduction << '\n'
      << feature_sampling_fraction << '\n'
      << max_thresholds_fraction << '\n'
      << force_exhaustive << '\n'
      << force_greedy << '\n'
      << verbose << std::endl;

  if (!out)
  {
    THEA_ERROR << "JointBoost: Error writing options";
    return false;
  }
  else
    return true;
}

bool
JointBoost::Options::serialize(TextOutputStream & output) const
{
  output.printf("min_boosting_rounds = %ld\n", min_boosting_rounds);
  output.printf("max_boosting_rounds = %ld\n", max_boosting_rounds);
  output.printf("min_fractional_error_reduction = %lg\n", min_fractional_error_reduction);
  output.printf("feature_sampling_fraction = %lg\n", feature_sampling_fraction);
  output.printf("max_thresholds_fraction = %lg\n", max_thresholds_fraction);
  output.printf("force_exhaustive = %s\n", (force_exhaustive ? "true" : "false"));
  output.printf("force_greedy = %s\n", (force_greedy ? "true" : "false"));
  output.printf("verbose = %s\n", (verbose ? "true" : "false"));

  return true;
}

bool
JointBoost::SharedStump::deserialize(std::istream & in)
{
  if (!(in >> f >> n >> a >> b >> theta))
    return false;

  long num_k = 0;
  in >> num_k;
  if (!in || num_k < 0)
    return false;

  k.resize((std::size_t)num_k);
  for (std::size_t i = 0; i < k.size(); ++i)
    in >> k[i];

  return !in.fail();
}

bool
JointBoost::SharedStump::serialize(std::ostream & out) const
{
  out << f << '\n'
      << n << '\n'
      << a << '\n'
      << b << '\n'
      << theta << '\n';

  out << k.size();
  for (std::size_t i = 0; i < k.size(); ++i)
    out << ' ' << k[i];

  out << std::endl;

  return !out.fail();
}

bool
JointBoost::deserialize(std::istream & in)
{
  if (!(in >> num_classes >> num_features))
  {
    THEA_ERROR << "JointBoost: Could not read class and feature counts";
    return false;
  }

  if (num_classes < 2)
  {
    THEA_ERROR << "JointBoost: Classifier has invalid number of classes";
    return false;
  }

  if (num_features < 1)
  {
    THEA_ERROR << "JointBoost: Classifier has invalid number of features";
    return false;
  }

  class_names.resize((array_size_t)num_classes);
  for (array_size_t i = 0; i < class_names.size(); )
  {
    if (!std::getline(in, class_names[i]))
    {
      THEA_ERROR << "JointBoost: Could not read name of class " << i;
      return false;
    }

    class_names[i] = trimWhitespace(class_names[i]);

    if (class_names[i].empty())
      continue;
    else
      i++;
  }

  long num_stumps = 0;
  if (!(in >> num_stumps))
  {
    THEA_ERROR << "JointBoost: Could not read number of stumps";
    return false;
  }

  if (num_stumps < 0)
  {
    THEA_ERROR << "JointBoost: Invalid number of stumps (" << num_stumps << ')';
    return false;
  }

  stumps.resize((std::size_t)num_stumps);
  for (std::size_t i = 0; i < stumps.size(); ++i)
  {
    stumps[i] = SharedStump::Ptr(new SharedStump);
    if (!stumps[i]->deserialize(in))
    {
      THEA_ERROR << "JointBoost: Could not read stump " << i;
      return false;
    }
  }

  return true;
}

bool
JointBoost::serialize(std::ostream & out) const
{
  out << num_classes << '\n'
      << num_features << '\n' << std::endl;

  if (class_names.empty())
  {
    for (int i = 0; i < num_classes; ++i)
      out << i << std::endl;
  }
  else
  {
    for (array_size_t i = 0; i < class_names.size(); ++i)
      out << class_names[i] << std::endl;
  }

  out << std::endl;

  out << stumps.size() << std::endl;
  for (std::size_t i = 0; i < stumps.size(); ++i)
    if (!stumps[i]->serialize(out))
      return false;

  return !out.fail();
}

void
JointBoost::dumpToConsole() const
{
  THEA_CONSOLE << "JointBoost: Classifier has " << stumps.size() << " stump(s)";

  for (array_size_t m = 0; m < stumps.size(); ++m)
  {
    SharedStump const & stump = *stumps[m];

    THEA_CONSOLE << "JointBoost:     -- Stump " << m << ':';
    THEA_CONSOLE << "JointBoost:         -- f      =  " << stump.f;
    THEA_CONSOLE << "JointBoost:         -- n      =  " << stump.n;
    THEA_CONSOLE << "JointBoost:         -- a      =  " << stump.a;
    THEA_CONSOLE << "JointBoost:         -- b      =  " << stump.b;
    THEA_CONSOLE << "JointBoost:         -- theta  =  " << stump.theta;

    std::ostringstream oss; oss << "[ ";
    for (array_size_t i = 0; i < stump.k.size(); ++i)
    {
      if (i > 0) oss << ", ";
      oss << stump.k[i];
    }
    oss << " ]";

    THEA_CONSOLE << "JointBoost:         -- k      =  " << oss.str();
  }
}

} // namespace Algorithms
} // namespace Thea
