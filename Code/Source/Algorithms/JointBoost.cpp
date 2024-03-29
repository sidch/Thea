//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2012
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

namespace JointBoostInternal {

void
write(SharingSet const & s, std::ostream & out)
{
  for (auto b : s)
    out << (b ? 1 : 0);
}

bool
read(SharingSet & s, std::istream & in)
{
  for (size_t i = 0; i < s.size(); ++i)
  {
    auto c = in.get();
    if (c == in.widen('0'))
      s[i] = false;
    else if (c == in.widen('1'))
      s[i] = true;
    else
    {
      THEA_ERROR << "JointBoost: Could not read bitfield with " << s.size() << " expected elements";
      return false;
    }
  }

  return true;
}

} // namespace JointBoostInternal

std::string
JointBoost::SharedStump::toString() const
{
  std::ostringstream oss;
  oss << "Stump[f: " << f << ", n: ";
  JointBoostInternal::write(n, oss);
  oss << ']';

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

JointBoost::JointBoost(intx num_classes_, intx num_features_, Options const & options_)
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
JointBoost::getClassName(intx i) const
{
  if (class_names.empty())
  {
    std::ostringstream rout;
    rout << i;
    return rout.str();
  }
  else
    return class_names[(size_t)i];
}

intx
JointBoost::train(TrainingData const & training_data_, TrainingData const * validation_data_)
{
  alwaysAssertM(training_data_.numExamples() > 0, "JointBoost: No training examples provided");
  alwaysAssertM(num_features == training_data_.numFeatures(), "JointBoost: Unexpected number of features in training data");

  // Create a temporary reference to the training data so we don't have to pass it to every function
  training_data = &training_data_;

  // Set values of parameters
  intx min_rounds = options.min_boosting_rounds;
  intx max_rounds = options.max_boosting_rounds;

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
  alwaysAssertM(class_names.empty() || (intx)class_names.size() == num_classes,
                "JointBoost: Incorrect number of class names specified");

  // Get classes for training data
  Array<intx> classes;
  training_data->getClasses(classes);

  // Initialize the weights
  Array<double> training_weights;
  training_data->getWeights(training_weights);
  weights.resize(num_classes, (intx)classes.size());
  if (!training_weights.empty())
  {
    for (size_t i = 0; i < training_weights.size(); ++i)
      for (intx c = 0; c < num_classes; ++c)
        weights(c, (intx)i) = training_weights[i];
  }
  else
    weights.fill(1);

  if (options.verbose)
  {
    THEA_CONSOLE << "JointBoost: Starting training (" << num_classes << " classes, " << num_features << " features, "
                 << training_data->numExamples() << " examples)";
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
  for (intx round = 0; round < max_rounds; ++round)
  {
    if (options.verbose) THEA_CONSOLE << "JointBoost: Training round " << round;

    // Create the stump for this round
    auto stump = std::make_shared<SharedStump>(num_classes);

    // Compute the k values for the stump
    for (intx c = 0; c < num_classes; ++c)
    {
      double k_numer = 0, k_denom = 0;
      for (size_t i = 0; i < classes.size(); ++i)
      {
        double w = weights(c, (intx)i);
        int z = (classes[i] == c ? +1 : -1);

        k_numer += (w * z);
        k_denom += w;
      }

      stump->k[(size_t)c] = (k_denom != 0 ? k_numer / k_denom : 0);
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
        round_validation_err = computeValidationError(validation_data_, stump);

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
    Array<double> stump_features;
    training_data->getFeature(stump->f, stump_features);

    for (intx c = 0; c < num_classes; ++c)
      for (size_t i = 0; i < classes.size(); ++i)
      {
        // THEA_CONSOLE << "      weights[c = " << c << ", i = " << i << "] = " << weights(c, (intx)i);

        int z = (classes[i] == c ? +1 : -1);
        double h = (*stump)(stump_features[i], c);

        // THEA_CONSOLE << "      h[i = " << i << ", c = " << c << "] = " << h << ", z = " << z;

        weights(c, (intx)i) *= std::exp(-z * h);

        // THEA_CONSOLE << "      weights[c = " << c << ", i = " << i << "] = " << weights(c, (intx)i);
      }
  }

  training_data = nullptr;

  if (options.verbose)
  {
    if (validation_data_)
    {
      if (validation_error < 0)
        validation_error = computeValidationError(validation_data_);

      THEA_CONSOLE << "JointBoost: Completed training, added " << stumps.size() << " stump(s) with final error " << error
                   << ", validation error " << validation_error;
    }
    else
      THEA_CONSOLE << "JointBoost: Completed training, added " << stumps.size() << " stump(s) with final error " << error;
  }

  return (intx)stumps.size();
}

double
JointBoost::optimizeStump(SharedStump & stump, Array<intx> const & stump_classes)
{
  if (num_features <= 0)
  {
    THEA_WARNING << "JointBoost:     Can't optimize stump with " << num_features << " features";
    return -1;
  }

  // Select a random fraction of the features
  Array<int32> candidate_features;
  if (feature_sampling_fraction >= 1)
  {
    for (intx f = 0; f < num_features; ++f)
      candidate_features.push_back(f);
  }
  else
  {
    int32 num_candidate_features = std::max((int32)std::round(feature_sampling_fraction * num_features), (int32)1);
    candidate_features.resize((size_t)num_candidate_features);
    Random::common().integers(0, (int32)num_features - 1, num_candidate_features, &candidate_features[0]);
  }

  if (options.verbose && (intx)candidate_features.size() < num_features)
    THEA_CONSOLE << "JointBoost:     Optimizing over " << candidate_features.size() << " randomly selected feature(s)";

  SharedStump test = stump;
  double min_err = -1;
  for (size_t f = 0; f < candidate_features.size(); ++f)
  {
    test.f = candidate_features[f];

    // Get feature values
    Array<double> stump_features;
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
JointBoost::optimizeStumpExhaustive(SharedStump & stump, Array<double> const & stump_features,
                                    Array<intx> const & stump_classes)
{
  // Loop over all possible subsets of classes
  alwaysAssertM((size_t)num_classes < sizeof(uintx), "JointBoost: Class count exceeds bitcount of largest unsigned integer");
  uintx num_subsets = ((uintx)1 << num_classes);
  SharedStump test = stump;
  double min_err = -1;
  double cum_num_thresholds = 0;
  for (uintx subset = 1; subset + 1 < num_subsets; ++subset)
  {
    for (intx b = 0; b < num_classes; ++b)
      test.n[(size_t)b] = (bool)((subset >> b) & 1);

    // THEA_CONSOLE << "JointBoost:         Testing stump " << test.toString();

    intx num_thresholds = 0;
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
JointBoost::optimizeStumpGreedy(SharedStump & stump, Array<double> const & stump_features, Array<intx> const & stump_classes)
{
  Array<SharedStump> candidate_stumps;
  Array<double> candidate_errors;

  SharedStump test = stump;
  SharingSet current_n((size_t)num_classes, false);  // initially empty
  double cum_num_thresholds = 0;
  intx num_fitted_stumps = 0;
  for (intx c = 0; c < num_classes - 1; ++c)
  {
    // Add the single new class to the current set that jointly gives the minimum error
    SharedStump best_stump;
    double min_err = -1;
    bool min_found = false;
    for (intx c = 0; c < num_classes - 1; ++c)
    {
      if (current_n[(size_t)c])
        continue;

      test.n = current_n;
      test.n[(size_t)c] = true;

      // THEA_CONSOLE << "JointBoost:         Testing stump " << test.toString();

      intx num_thresholds = 0;
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
  for (size_t i = 0; i < candidate_stumps.size(); ++i)
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
  bool operator()(size_t i, size_t j) const { return values[i] < values[j]; }

  double const * values;
};

void
sortIndexed(Array<double> const & values, Array<size_t> & sorted_indices)
{
  sorted_indices.resize(values.size());
  for (size_t i = 0; i < sorted_indices.size(); ++i)
    sorted_indices[i] = i;

  std::sort(sorted_indices.begin(), sorted_indices.end(), IndexedComparator(&values[0]));
}

void
getClassificationAccuracy(double split_value, SharingSet const & pos_classes, Array<double> const & features,
                          Array<intx> const & classes, Array<size_t> const & sorted_indices, Array<int> & accuracy)
{
  accuracy.resize(features.size());
  for (size_t i = 0; i < features.size(); ++i)
  {
    size_t index = sorted_indices[i];

    intx c = classes[index];
    double f = features[index];
    bool is_positive = pos_classes[(size_t)c];
    bool correct = ((is_positive && f > split_value) || (!is_positive && f <= split_value));
    accuracy[index] = (correct ? +1 : -1);
  }
}

double
splitQuality(MatrixX<double> const & weights, Array<intx> const & classes, Array<int> const & accuracy)
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
  for (size_t i = 0; i < accuracy.size(); ++i)
    quality += (accuracy[i] * weights(classes[i], (intx)i));

  return quality;
}

void
getCandidateThresholds(SharingSet const & pos_classes, MatrixX<double> const & weights, Array<double> const & features,
                       Array<intx> const & classes, intx max_thresholds, Array<double> & thresholds)
{
  // A good threshold separates the positive classes from the negative classes. In other words, we want either many positive
  // examples > theta and many negative examples <= theta, or vice versa. Note that the positive examples need *not* be mostly
  // on the "positive side of" (greater than) theta -- we just require theta to separate the positive and negative examples
  // as well as possible. (Hence, the absolute value of splitQuality() is relevant, not its sign.)

  // Sort the examples by feature value
  Array<size_t> sorted_indices;
  sortIndexed(features, sorted_indices);

  // Get the minimum and maximum feature values
  double lo = features[0], hi = features[0];
  for (size_t i = 1; i < features.size(); ++i)
  {
    if (features[i] < lo) lo = features[i];
    if (features[i] > hi) hi = features[i];
  }

  // Choose a tolerance in comparing features
  double tolerance = 1.0e-10 * (hi - lo);

  // Generate samples and adjust their position to be better cuts
  typedef UnorderedSet<size_t> IndexSet;
  IndexSet threshold_set;
  for (intx t = 0; t < max_thresholds; ++t)
  {
    // Start with a random seed
    size_t index = (size_t)Random::common().integer(0, (int32)sorted_indices.size() - 1);
    size_t mapped_index = sorted_indices[index];
    double threshold = features[mapped_index];

    // Measure the classification accuracy (+/-1) per example, and the overall quality of the split
    Array<int> accuracy;
    getClassificationAccuracy(threshold, pos_classes, features, classes, sorted_indices, accuracy);
    double quality = splitQuality(weights, classes, accuracy);

    // Iteratively try to improve it by moving one step to the right or left
    for (size_t round = 0; round < sorted_indices.size(); ++round)
    {
      // Get the quality change in each direction
      double offset_qualities[2] = { 0, 0 };
      bool equal_quality[2] = { false, false };
      for (int offset = -1; offset <= 1; offset += 2)
      {
        if ((offset < 0 && index == 0) || (offset > 0 && index == sorted_indices.size() - 1))
          continue;

        size_t offset_index = (size_t)(index + offset);
        size_t mapped_offset_index = sorted_indices[offset_index];

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
            intx c = classes[mapped_index];
            offset_quality -= 2 * (accuracy[mapped_index] * weights(c, (intx)mapped_index));
          }
          else // move right
          {
            // The offset example's classification accuracy is inverted. The current example remains unchanged.
            intx c = classes[mapped_offset_index];
            offset_quality -= 2 * (accuracy[mapped_offset_index] * weights(c, (intx)mapped_offset_index));
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
        index = (size_t)(index + 2 * best_offset - 1);
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
JointBoost::fitStump(SharedStump & stump, Array<double> const & stump_features, Array<intx> const & stump_classes,
                     intx * num_generated_thresholds)
{
  using namespace JointBoostInternal;

  // Find good places to cut the features (candidate values of the threshold theta)
  Array<double> thresholds;

#ifdef JOINT_BOOST_TEST_ALL_THRESHOLDS
  thresholds = stump_features;
#else
  intx max_thresholds = (intx)std::ceil(max_thresholds_fraction * stump_features.size());
  if (max_thresholds < (intx)stump_features.size())
    getCandidateThresholds(stump.n, weights, stump_features, stump_classes, max_thresholds, thresholds);
  else
    thresholds = stump_features;
#endif

  if (num_generated_thresholds)
    *num_generated_thresholds = (intx)thresholds.size();

  // THEA_CONSOLE << "JointBoost: Generated " << thresholds.size() << " candidate threshold(s) for " << stump.toString();
  // for (size_t t = 0; t < thresholds.size(); ++t)
  //   THEA_CONSOLE << "JointBoost:     -- threshold[" << t << "] = " << thresholds[t];

  // Precalculate the error term depending on stump.k -- it does not depend on threshold selection
  double err_k = 0;
  for (intx c = 0; c < num_classes; ++c)
  {
    if (!stump.n[(size_t)c])
    {
      for (size_t i = 0; i < stump_classes.size(); ++i)
      {
        double w = weights(c, (intx)i);
        int z = (stump_classes[i] == c ? +1 : -1);

        err_k += w * Math::square(z - stump.k[(size_t)c]);  // k is precalculated outside this function
      }
    }
  }

  // THEA_CONSOLE << "err_k = " << err_k;

  // Fit the stump parameters to the observed data
  intx best_threshold = -1;
  double min_err = -1;
  for (size_t t = 0; t < thresholds.size(); ++t)
  {
    double theta = thresholds[t];
    double a_numer = 0, a_denom = 0, b_numer = 0, b_denom = 0;
    for (intx c = 0; c < num_classes; ++c)
    {
      if (stump.n[(size_t)c])
      {
        for (size_t i = 0; i < stump_classes.size(); ++i)
        {
          double w = weights(c, (intx)i);
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
      best_threshold = (intx)t;
      min_err = err;

      stump.a = a;
      stump.b = b;
      stump.theta = theta;
    }
  }

  return min_err;
}

intx
JointBoost::predict(double const * features, double * class_probabilities) const
{
  intx best_class = -1;
  double best_H = -1;
  double sum_probs = 0;

  for (intx c = 0; c < num_classes; ++c)
  {
    double H = 0;
    for (size_t m = 0; m < stumps.size(); ++m)
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
    for (intx c = 0; c < num_classes; ++c)
      class_probabilities[c] /= sum_probs;
  }

  return best_class;
}

double
JointBoost::computeValidationError(TrainingData const * validation_data_, SharedStump::Ptr new_stump)
{
  alwaysAssertM(validation_data_, "JointBoost: No validation set provided");
  alwaysAssertM(validation_data_->numFeatures() == num_features, "JointBoost: Validation set has different number of features");

  if (new_stump)
    stumps.push_back(new_stump);

  MatrixX<double, MatrixLayout::ROW_MAJOR> validation_features(validation_data_->numExamples(), num_features);
  Array<double> feat;
  for (intx i = 0; i < num_features; ++i)
  {
    validation_data_->getFeature(i, feat);
    validation_features.col(i) = VectorXd::Map(&feat[0], (intx)feat.size());
  }

  Array<intx> validation_classes;
  Array<double> validation_weights;

  validation_data_->getClasses(validation_classes);
  validation_data_->getWeights(validation_weights);

  alwaysAssertM((intx)validation_classes.size() == validation_features.rows(),
                "JointBoost: Ground truth classes for validation data do not match number of examples");

  double err = 0;
  for (intx i = 0; i < validation_features.rows(); ++i)
  {
    intx predicted = predict(validation_features.row(i).data());
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
    return read(in);
  }
  THEA_CATCH(return false;, ERROR, "JointBoost: Could not load options from input file '%s'", path.c_str())
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
    if (!write(out))
      return false;

    out.commit();
  }
  THEA_CATCH(return false;, ERROR, "JointBoost: Could not save options to output file '%s'", path.c_str())

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

  if (!options.read(in))
    return false;

  if (!read(in))
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

  if (!options.write(out))
    return false;

  out << std::endl;  // break up the output a bit

  if (!write(out))
    return false;

  return true;
}

bool
JointBoost::Options::read(std::istream & in)
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
JointBoost::Options::read(TextInputStream & input)
{
  *this = defaults();

  while (input.hasMore())
  {
    std::string field = input.readSymbol();
    input.readSymbol("=");

    if (field == "min_boosting_rounds")
      min_boosting_rounds = (intx)input.readNumber();
    else if (field == "max_boosting_rounds")
      max_boosting_rounds = (intx)input.readNumber();
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
JointBoost::Options::write(std::ostream & out) const
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
JointBoost::Options::write(TextOutputStream & output) const
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
JointBoost::SharedStump::read(std::istream & in)
{
  if (!(in >> f)) { return false; }
  if (!JointBoostInternal::read(n, in)) { return false; }
  if (!(in >> a >> b >> theta)) { return false; }
  for (auto & kv : k) { if (!(in >> kv)) return false; }

  return true;
}

bool
JointBoost::SharedStump::write(std::ostream & out) const
{
  out << f << '\n';
  JointBoostInternal::write(n, out); out << '\n';
  out << a << '\n'
      << b << '\n'
      << theta << '\n'
      << stringJoin(k.begin(), k.end(), ' ') << std::endl;

  return !out.fail();
}

bool
JointBoost::read(std::istream & in)
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

  class_names.resize((size_t)num_classes);
  for (size_t i = 0; i < class_names.size(); )
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

  intx num_stumps = 0;
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

  stumps.resize((size_t)num_stumps);
  for (size_t i = 0; i < stumps.size(); ++i)
  {
    stumps[i] = std::make_shared<SharedStump>(num_classes);
    if (!stumps[i]->read(in))
    {
      THEA_ERROR << "JointBoost: Could not read stump " << i;
      return false;
    }
  }

  return true;
}

bool
JointBoost::write(std::ostream & out) const
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
    for (size_t i = 0; i < class_names.size(); ++i)
      out << class_names[i] << std::endl;
  }

  out << std::endl;

  out << stumps.size() << std::endl;
  for (size_t i = 0; i < stumps.size(); ++i)
    if (!stumps[i]->write(out))
      return false;

  return !out.fail();
}

void
JointBoost::dumpToConsole() const
{
  THEA_CONSOLE << "JointBoost: Classifier has " << stumps.size() << " stump(s)";

  for (size_t m = 0; m < stumps.size(); ++m)
  {
    SharedStump const & stump = *stumps[m];

    THEA_CONSOLE << "JointBoost:     -- Stump " << m << ':';
    THEA_CONSOLE << "JointBoost:         -- f      =  " << stump.f;

    std::ostringstream oss; JointBoostInternal::write(stump.n, oss);
    THEA_CONSOLE << "JointBoost:         -- n      =  " << oss.str();

    THEA_CONSOLE << "JointBoost:         -- a      =  " << stump.a;
    THEA_CONSOLE << "JointBoost:         -- b      =  " << stump.b;
    THEA_CONSOLE << "JointBoost:         -- theta  =  " << stump.theta;
    THEA_CONSOLE << "JointBoost:         -- k      =  [" << stringJoin(stump.k.begin(), stump.k.end(), ',') << ']';
  }
}

} // namespace Algorithms
} // namespace Thea
