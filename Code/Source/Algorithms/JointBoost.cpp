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
#include "../UnorderedSet.hpp"
#include <algorithm>
#include <sstream>

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
  force_greedy(false)
{
}

JointBoost::JointBoost(long num_classes_, long num_features_, Options const & options_)
: num_classes(num_classes_),
  num_features(num_features_),
  options(options_)
{
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

void
JointBoost::train(TrainingData const & training_data_)
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

  // Get classes for training data
  TheaArray<long> stump_classes;
  training_data->getClasses(stump_classes);

  // Initialize the weights
  weights.resize(num_classes, (long)stump_classes.size());
  weights.fill(1);

  THEA_CONSOLE << "JointBoost: Starting training";
  THEA_CONSOLE << "JointBoost:     -- min_rounds = " << min_rounds;
  THEA_CONSOLE << "JointBoost:     -- max_rounds = " << max_rounds;
  THEA_CONSOLE << "JointBoost:     -- min_fractional_error_reduction = " << min_fractional_error_reduction;
  THEA_CONSOLE << "JointBoost:     -- feature_sampling_fraction = " << feature_sampling_fraction;
  THEA_CONSOLE << "JointBoost:     -- force_exhaustive = " << options.force_exhaustive;
  THEA_CONSOLE << "JointBoost:     -- force_greedy = " << options.force_greedy;

  // Do several rounds of boosting, adding a new stump (weak learner) in each round
  stumps.clear();
  double error = -1;
  for (long round = 0; round < max_rounds; ++round)
  {
    THEA_CONSOLE << "JointBoost: Training round " << round;

    // Create the stump for this round
    SharedStump::Ptr stump(new SharedStump);

    // Compute the k values for the stump
    stump->k.resize((array_size_t)num_classes, 0);
    for (long c = 0; c < num_classes; ++c)
    {
      double k_numer = 0, k_denom = 0;
      for (array_size_t i = 0; i < stump_classes.size(); ++i)
      {
        double w = weights((long)i, stump_classes[i]);
        int z = (stump_classes[i] == c ? +1 : -1);

        k_numer += w * z;
        k_denom += w;
      }

      stump->k[(array_size_t)c] = (k_denom != 0 ? k_numer / k_denom : 0);
    }

    // Optimize this stump
    double round_err = optimizeStump(*stump, stump_classes);

    THEA_CONSOLE << "JointBoost:     Error in round " << round << " = " << round_err << " with " << stump->toString();

    if (round_err >= 0)
    {
      // Check if we should stop at this point
      if (round >= min_rounds && error >= 0)
      {
        double err_reduction = error - round_err;
        if (err_reduction < min_fractional_error_reduction * error)
          break;
      }

      // We have a new error value
      error = round_err;

      // Add the stump to the strong classifier
      stumps.push_back(stump);

      // Update weights
      TheaArray<double> stump_features;
      training_data->getFeature(stump->f, stump_features);

      for (long c = 0; c < num_classes; ++c)
        for (array_size_t i = 0; i < stump_classes.size(); ++i)
        {
          // THEA_CONSOLE << "      weights[c = " << c << ", i = " << i << "] = " << weights(c, (long)i);

          int z = (stump_classes[i] == c ? +1 : -1);
          double h = (*stump)(stump_features[i], c);

          // THEA_CONSOLE << "      h[i = " << i << ", c = " << c << "] = " << h << ", z = " << z;

          weights(c, (long)i) *= std::exp(-z * h);
        }
    }
  }

  training_data = NULL;

  THEA_CONSOLE << "JointBoost: Completed training, added " << stumps.size() << " stump(s) with final error " << error;
}

double
JointBoost::optimizeStump(SharedStump & stump, TheaArray<long> const & stump_classes)
{
  SharedStump test = stump;
  double min_err = -1;
  for (long i = 0; i < num_features; ++i)
  {
    if (feature_sampling_fraction < 1)
      if (Math::rand01() > feature_sampling_fraction)
        continue;

    test.f = i;

    // Get feature values
    TheaArray<double> stump_features;
    training_data->getFeature(test.f, stump_features);

    // Optimize the stump via an exhaustive or a greedy search over subsets of classes
    double err = 0;
    if ((num_classes <= 5 && !options.force_greedy) || options.force_exhaustive)
      err = optimizeStumpExhaustive(test, stump_features, stump_classes);
    else
      err = optimizeStumpGreedy(test, stump_features, stump_classes);

    THEA_CONSOLE << "JointBoost:     Error when splitting on feature " << i << " = " << err;

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
  unsigned long num_subsets = (1L << num_classes);
  double min_err = -1;
  for (unsigned long subset = 1; subset + 1 < num_subsets; ++subset)
  {
    test.n = SharingSet(num_classes, subset);
    double err = fitStump(test, stump_features, stump_classes);
    if (err >= 0 && (min_err < 0 || err < min_err))
    {
      min_err = err;
      stump = test;
    }
  }

  return min_err;
}

double
JointBoost::optimizeStumpGreedy(SharedStump & stump, TheaArray<double> const & stump_features,
                                TheaArray<long> const & stump_classes)
{
  // TODO
  return -1;
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
  // and is below the split index , and -1 otherwise

  double quality = 0;
  for (array_size_t i = 0; i < accuracy.size(); ++i)
    quality += (accuracy[i] * weights(classes[i], (long)i));

  return quality;
}

void
getCandidateThresholds(SharingSet const & pos_classes, Matrix<double> const & weights, TheaArray<double> const & features,
                       TheaArray<long> const & classes, long max_thresholds, TheaArray<double> & thresholds)
{
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
    array_size_t index = (array_size_t)Math::randIntegerInRange(0, (long)sorted_indices.size() - 1);
    array_size_t mapped_index = sorted_indices[index];
    double feature = features[mapped_index];

    // Measure the classification accuracy (+/-1) per example, and the overall quality of the split
    TheaArray<int> accuracy;
    getClassificationAccuracy(feature, pos_classes, features, classes, sorted_indices, accuracy);
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
       || (equal_quality[best_offset] && Math::rand01() < 0.5))  // if equal, move with 50% probability
      {
        // Move to index + offset
        index = (array_size_t)(index + 2 * best_offset - 1);
        mapped_index = sorted_indices[index];
        quality = offset_qualities[best_offset];
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
JointBoost::fitStump(SharedStump & stump, TheaArray<double> const & stump_features, TheaArray<long> const & stump_classes)
{
  using namespace JointBoostInternal;

  // Find good places to cut the features (candidate values of the threshold theta)
  TheaArray<double> thresholds;
  long max_thresholds = (long)std::ceil(max_thresholds_fraction * stump_features.size());
  getCandidateThresholds(stump.n, weights, stump_features, stump_classes, max_thresholds, thresholds);

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
        double w = weights(c, stump_classes[i]);
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
          double w = weights((long)i, stump_classes[i]);
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
