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

#ifndef __Thea_Algorithms_JointBoost_hpp__
#define __Thea_Algorithms_JointBoost_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Matrix.hpp"
#include <boost/dynamic_bitset.hpp>

namespace Thea {
namespace Algorithms {

namespace JointBoostInternal {

typedef boost::dynamic_bitset<> SharingSet;

} // namespace JointBoostInternal

/**
 * A boosting classifier using shared features. Based on
 *
 * A. Torralba, K. Murphy and W. Freeman, "Sharing features: efficient boosting procedures for multiclass object detection",
 * Proc. CVPR, 2004.
 *
 * To use the class, implement an appropriate subclass of TrainingData, call train(), and then call classify().
 */
class THEA_API JointBoost
{
  public:
    THEA_DEF_POINTER_TYPES(JointBoost, shared_ptr, weak_ptr)

    /** Interface for accessing training data. */
    class TrainingData
    {
      public:
        THEA_DEF_POINTER_TYPES(TrainingData, shared_ptr, weak_ptr)

        /** Destructor. */
        virtual ~TrainingData() {}

        /** Get the number of training examples. */
        virtual long numExamples() const = 0;

        /** Get the number of features per example. */
        virtual long numFeatures() const = 0;

        /**
         * Get the values of a particular feature for all training examples. \a feature_index must be in the range
         * 0... numFeatures() - 1.
         */
        virtual void getFeature(long feature_index, TheaArray<double> & values) const = 0;

        /** Get the class of each training example. */
        virtual void getClasses(TheaArray<long> & classes) const = 0;

    }; // class TrainingData

    /**
     * %Options for the classifier. In most cases, passing a negative value for a normally non-negative parameter sets the
     * default value for that parameter.
     */
    struct Options
    {
      long min_boosting_rounds;  /**< Minimum number of boosting rounds that must be performed even if the error reduction
                                      between successive rounds is below the threshold min_fractional_error_reduction. */
      long max_boosting_rounds;  ///< Maximum number of boosting rounds. This also limits the maximum number of stumps added.
      double min_fractional_error_reduction;  /**< Minimum error reduction required for boosting rounds to continue beyond
                                                   min_boosting_rounds. */
      double feature_sampling_fraction;  ///< Fraction of features sampled in a round.
      double max_thresholds_fraction;  /**< Set the maximum number of candidate thresholds generated, as a fraction of the
                                            number of features. */
      bool force_exhaustive;  ///< Force exhaustive O(2^C) optimization over all possible subsets of classes.
      bool force_greedy;  ///< Force greedy O(C^2) optimization over subsets of classes.

      /** Constructor. Sets default options. */
      Options();

      /**
       * Set the minimum number of boosting rounds that must be performed even if the error reduction between successive rounds
       * is below the threshold min_fractional_error_reduction.
       */
      Options & setMinBoostingRounds(long rounds) { min_boosting_rounds = rounds; return *this; }

      /** Maximum number of boosting rounds. This also limits the maximum number of stumps added. */
      Options & setMaxBoostingRounds(long rounds) { max_boosting_rounds = rounds; return *this; }

      /** Minimum error reduction required for boosting rounds to continue beyond min_boosting_rounds. */
      Options & setMinFractionalErrorReduction(double frac) { min_fractional_error_reduction = frac; return *this; }

      /** Set the fraction of features sampled in a round. */
      Options & setFeatureSamplingFraction(double frac) { feature_sampling_fraction = frac; return *this; }

      /** Set the maximum number of candidate thresholds generated, as a fraction of the number of features. */
      Options & setMaxThresholdsFraction(double frac) { max_thresholds_fraction = frac; return *this; }

      /** Set if exhaustive O(2^C) optimization over all possible subsets of classes will be forced or not. */
      Options & setForceExhaustive(bool value) { force_exhaustive = value; return *this; }

      /** Set if greedy O(C^2) optimization over subsets of classes will be forced or not. */
      Options & setForceGreedy(bool value) { force_greedy = value; return *this; }

      /** Get the set of default options. */
      static Options const & defaults() { static Options const def; return def; }

    }; // class Options;

    /**
     * Constructor.
     *
     * @param num_classes_ Number of classes for classification. The classes are numbered 0 to \a num_classes - 1.
     * @param num_features_ Number of features per object.
     * @param options_ Addition options controlling the behaviour of the classifier.
     */
    JointBoost(long num_classes_, long num_features_, Options const & options_ = Options::defaults());

    /** Destructor. */
    ~JointBoost();

    /** Reset the classifier to the initial state. */
    void clear();

    /** Get the number of classes into which objects may fall. The classes are numbered 0 to numClasses() - 1. */
    long numClasses() const { return num_classes; }

    /** Get the number of features for an object. */
    long numFeatures() const { return num_features; }

    /**
     * Train the strong classifier by several rounds of boosting.
     *
     * @param training_data_ Data used for training the classifier.
     * @param min_rounds The minimum number of boosting rounds that must be performed even if the error reduction between
     *   successive rounds is below the threshold \a min_error_reduction.
     * @param max_rounds The maximum number of boosting rounds. This also limits the maximum number of stumps added.
     * @param min_error_reduction The minimum (subtractive) reduction in the classification error needed to continue iteration
     *   beyond \a min_rounds rounds. A negative argument chooses a default value.
     */
    void train(TrainingData const & training_data_);

    /**
     * Predict the most likely class for an object with a given set of features, optionally also returning the probability of
     * belonging to each class.
     *
     * @param features The features of the object. Must contain numFeatures() values.
     * @param class_probabilities If non-null, used to return the probability of belonging to each class. Must be preallocated
     *   to numClasses() elements.
     *
     * @return The index of the most likely class of the object.
     */
    long predict(double const * features, double * class_probabilities = NULL) const;

    /** Print debugging information about this classifier to the console. */
    void dumpToConsole() const;

  private:
    /** A subset of classes encoded as a set of bits, one per class in global set. */
    typedef JointBoostInternal::SharingSet SharingSet;

    /** Decision stump (weak learner). */
    struct SharedStump
    {
      THEA_DEF_POINTER_TYPES(SharedStump, shared_ptr, weak_ptr)

      long f;               ///< The index of the feature for this stump.
      SharingSet n;         ///< Indices of positive classes for this stump.
      double a;             ///< Regression weight a.
      double b;             ///< Regression weight b.
      double theta;         ///< Cut threshold.
      TheaArray<double> k;  ///< Constants for classes not in the sharing set.

      /** Evaluate the stump for a given feature and class. */
      double operator()(double feature_value, long class_index) const
      {
        if (n[(SharingSet::size_type)class_index])
          return feature_value > theta ? a : b;
        else
          return k[(array_size_t)class_index];
      }

      /** Get a printable description of the stump. */
      std::string toString() const;

    }; // struct SharedStump

    /** Optimize a stump, for the current set of weights, by searching over features and subsets of classes. */
    double optimizeStump(SharedStump & stump, TheaArray<long> const & stump_classes);

    /** Optimize a stump via exhaustive O(2^C) search over all subsets of classes. */
    double optimizeStumpExhaustive(SharedStump & stump, TheaArray<double> const & stump_features,
                                   TheaArray<long> const & stump_classes);

    /** Optimize a stump via greedy O(C^2) search over subsets of classes. */
    double optimizeStumpGreedy(SharedStump & stump, TheaArray<double> const & stump_features,
                               TheaArray<long> const & stump_classes);

    /** Fit stump parameters a, b and theta, for a particular subset of classes. */
    double fitStump(SharedStump & stump, TheaArray<double> const & stump_features, TheaArray<long> const & stump_classes);

    long num_classes;   ///< Number of object classes.
    long num_features;  ///< Number of features per object.
    Options options;    ///< Additional options.

    TrainingData const * training_data;  ///< Cached handle to training data, valid only during training.
    double feature_sampling_fraction;  ///< Cached fraction of features test per round, valid only during training.
    double max_thresholds_fraction;  /**< Cached number of candidate thresholds, expressed as a fraction of the number of
                                          features, valid only during training. */

    Matrix<double> weights;  ///< Weights indexed by (class index, object_index).
    TheaArray<SharedStump::Ptr> stumps;  ///< Current set of selected decision stumps.

}; // class JointBoost

} // namespace Algorithms
} // namespace Thea

#endif
