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
#include "LogisticRegression.hpp"
#include "../Math.hpp"
#include <sstream>

namespace Thea {
namespace Algorithms {

std::string
JointBoost::SharedStump::toString() const
{
  std::ostringstream oss;
  oss << "Stump[f: " << feature_index << ", n: " << n << ']';
  return oss.str();
}

JointBoost::Options::Options()
: feature_sampling_fraction(-1)
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
JointBoost::train(TrainingData const & training_data_, long min_rounds, long max_rounds, double min_error_reduction)
{
  alwaysAssertM(training_data_.numExamples() > 0, "JointBoost: No training examples provided");
  alwaysAssertM(num_features == training_data_.numFeatures(), "JointBoost: Unexpected number of features in training data");
  alwaysAssertM(min_rounds > 0, "JointBoost: Training requires at least one round of boosting");

  // Create a temporary reference to the training data so we don't have to pass it to every function
  training_data = &training_data_;

  if (max_rounds < min_rounds)
    max_rounds = min_rounds;

  if (min_error_reduction < 0)
  {
    // TODO: Choose default value
  }

  double sampling_fraction = options.feature_sampling_fraction;
  if (sampling_fraction <= 0) sampling_fraction = 1.0 / max_rounds;

  // Get classes for training data
  TheaArray<long> stump_classes;
  training_data->getClasses(stump_classes);

  // Initialize the weights
  weights.resize(num_classes, (long)stump_classes.size());
  weights.fill(1);

  // Do several rounds of boosting, adding a new stump (weak learner) in each round
  stumps.clear();
  double error = std::numeric_limits<double>::max();
  for (long round = 0; round < max_rounds; ++round)
  {
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
    double round_err = optimizeStump(*stump, stump_classes, sampling_fraction);

    if (round_err >= 0)
    {
      // Check if we should stop at this point
      if (round >= min_rounds)
      {
        double err_reduction = error - round_err;
        if (err_reduction < min_error_reduction)
          break;
      }

      // Add the stump to the strong classifier
      stumps.push_back(stump);

      // Update weights
      TheaArray<double> stump_features;
      training_data->getSingleFeature(stump->feature_index, stump_features);

      for (long c = 0; c < num_classes; ++c)
        for (array_size_t i = 0; i < stump_classes.size(); ++i)
        {
          int z = (stump_classes[i] == c ? +1 : -1);
          double h = (*stump)(stump_features[i], c);

          weights(c, (long)i) *= std::exp(-z * h);
        }
    }
  }

  training_data = NULL;
}

double
JointBoost::optimizeStump(SharedStump & stump, TheaArray<long> const & stump_classes, double feature_sampling_fraction)
{
  SharedStump test = stump;
  double min_err = -1;
  for (long i = 0; i < num_features; ++i)
  {
    if (options.feature_sampling_fraction < 1)
      if (Math::rand01() > options.feature_sampling_fraction)
        continue;

    test.feature_index = i;

    // Get feature values
    TheaArray<double> stump_features;
    training_data->getSingleFeature(stump.feature_index, stump_features);

    // Optimize the stump via an exhaustive or a greedy search over subsets of classes
    double err = 0;
    if ((num_classes <= 5 && !options.force_greedy) || options.force_exhaustive)
      err = optimizeStumpExhaustive(test, stump_features, stump_classes);
    else
      err = optimizeStumpGreedy(test, stump_features, stump_classes);

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
  for (unsigned long subset = 0; subset < num_subsets; ++subset)
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

double
JointBoost::fitStump(SharedStump & stump, TheaArray<double> const & stump_features, TheaArray<long> const & stump_classes)
{
  // Find the best place to cut the features (value of the threshold theta)
  LogisticRegression logreg(1);
  for (array_size_t i = 0; i < stump_features.size(); ++i)
  {
    long class_index = stump_classes[i];
    if (stump.n[(SharingSet::size_type)class_index])
      logreg.addObservation(&stump_features[i], 1);
    else
      logreg.addObservation(&stump_features[i], 0);
  }

  if (!logreg.solve())
  {
    THEA_WARNING << "JointBoost: Could not find splitting threshold for " << stump.toString() << " via logistic regression";
    return -1;
  }

  // Select the pivot of the best-fit logistic curve as the threshold value
  double theta = logreg.getSolution()[0];

  // Fit the stump parameters to the observed data
  double a_numer = 0, a_denom = 0, b_numer = 0, b_denom = 0;
  double err = 0;
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
    else
    {
      for (array_size_t i = 0; i < stump_classes.size(); ++i)
      {
        double w = weights((long)i, stump_classes[i]);
        int z = (stump_classes[i] == c ? +1 : -1);

        err += w * Math::square(z - stump.k[(array_size_t)c]);  // k is precalculated outside this function
      }
    }
  }

  stump.a = (a_denom != 0 ? a_numer / a_denom : 0);
  stump.b = (b_denom != 0 ? b_numer / b_denom : 0);

  err += (1 - Math::square(stump.a)) * a_denom + (1 - Math::square(stump.b)) * b_denom;

  return err;
}

} // namespace Algorithms
} // namespace Thea
