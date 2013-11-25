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
#include "../IOStream.hpp"
#include "../Matrix.hpp"
#include <boost/dynamic_bitset.hpp>
#include <iostream>

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
 * To use the class, implement an appropriate subclass of TrainingData, call train(), and then call predict().
 */
class THEA_API JointBoost
{
  public:
    THEA_DEF_POINTER_TYPES(JointBoost, shared_ptr, weak_ptr)

    /** Interface for accessing training data. */
    class THEA_API TrainingData
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

        /** Get the weight of each training example. If an empty array is returned, all weights are set to 1. */
        virtual void getWeights(TheaArray<double> & weights) const { weights.clear(); }

        /**
         * Get the name of each class. If it returns an empty vector, the numeric index of the class will be used as the label
         * name.
         */
        virtual void getClassNames(TheaArray<std::string> & names) const { names.clear(); };

    }; // class TrainingData

    /**
     * %Options for the classifier. In most cases, passing a negative value for a normally non-negative parameter auto-selects a
     * suitable value for that parameter.
     */
    class THEA_API Options
    {
      public:
        /** Constructor. Sets default options. */
        Options();

        /**
         * Set the minimum number of boosting rounds that must be performed even if the error reduction between successive
         * rounds is below the threshold min_fractional_error_reduction (default -1).
         */
        Options & setMinBoostingRounds(long rounds) { min_boosting_rounds = rounds; return *this; }

        /** Maximum number of boosting rounds (default -1). This also limits the maximum number of stumps added. */
        Options & setMaxBoostingRounds(long rounds) { max_boosting_rounds = rounds; return *this; }

        /** Minimum error reduction required for boosting rounds to continue beyond min_boosting_rounds (default -1). */
        Options & setMinFractionalErrorReduction(double frac) { min_fractional_error_reduction = frac; return *this; }

        /** Set the fraction of features sampled in a round (default -1). */
        Options & setFeatureSamplingFraction(double frac) { feature_sampling_fraction = frac; return *this; }

        /** Set the maximum number of candidate thresholds generated, as a fraction of the number of features (default -1). */
        Options & setMaxThresholdsFraction(double frac) { max_thresholds_fraction = frac; return *this; }

        /**
         * Set if exhaustive O(2^C) optimization over all possible subsets of classes will be forced or not (default false).
         * Currently, this requires that there be no more classes than the number of bits in an unsigned long, minus 1. If you
         * have a large number of classes, use the greedy optimization instead. The latter is turned on by default when there
         * are too many classes, or if you call setForceGreedy(true).
         */
        Options & setForceExhaustive(bool value) { force_exhaustive = value; return *this; }

        /** Set if greedy O(C^2) optimization over subsets of classes will be forced or not (default false). */
        Options & setForceGreedy(bool value) { force_greedy = value; return *this; }

        /** Set whether progress information will be printed to the console or not (default true). */
        Options & setVerbose(bool value) { verbose = value; return *this; }

        /** Load classifier options from a disk file. */
        bool load(std::string const & path);

        /** Save classifier options to a disk file. */
        bool save(std::string const & path) const;

        /** Load classifier options from a standard input stream. */
        bool deserialize(std::istream & in);

        /** Load classifier options from a text input stream. */
        bool deserialize(TextInputStream & in);

        /** Save classifier options to a standard output stream. */
        bool serialize(std::ostream & out) const;

        /** Save classifier options to a text output stream. */
        bool serialize(TextOutputStream & out) const;

        /** Get the set of default options. */
        static Options const & defaults() { static Options const def; return def; }

      private:
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
        bool verbose;  ///< Print progress information to the console.

        friend class JointBoost;

    }; // class Options

    /**
     * Constructor.
     *
     * @param num_classes_ Number of classes for classification. The classes are numbered 0 to \a num_classes - 1.
     * @param num_features_ Number of features per object.
     * @param options_ Additional options controlling the behaviour of the classifier.
     */
    JointBoost(long num_classes_, long num_features_, Options const & options_ = Options::defaults());

    /** Construct a JointBoost classifier by loading it from a file. */
    JointBoost(std::string const & path);

    /** Destructor. */
    ~JointBoost();

    /** Reset the classifier to the initial state. */
    void clear();

    /** Get the number of classes into which objects may fall. The classes are numbered 0 to numClasses() - 1. */
    long numClasses() const { return num_classes; }

    /** Get the number of features for an object. */
    long numFeatures() const { return num_features; }

    /** Get the name of a class. */
    std::string getClassName(long i) const;

    /** Get the current options for the classifier. */
    Options const & getOptions() const { return options; }

    /**
     * Train the strong classifier by several rounds of boosting.
     *
     * @param training_data_ Data used for training the classifier.
     * @param validation_data_ [Optional] Data used for cross-validation. Training is stopped if the error on the validation set
     *   starts increasing.
     *
     * @return The number of stumps added during the training process.
     */
    long train(TrainingData const & training_data_, TrainingData const * validation_data_ = NULL);

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

    /** Load the classifier from a disk file. */
    bool load(std::string const & path);

    /** Save the trained classifier to disk. */
    bool save(std::string const & path) const;

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

      /** Load the stump from an input stream. */
      bool deserialize(std::istream & in);

      /** Save the stump to an output stream. */
      bool serialize(std::ostream & out) const;

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
    double fitStump(SharedStump & stump, TheaArray<double> const & stump_features, TheaArray<long> const & stump_classes,
                    long * num_generated_thresholds = NULL);

    /**
     * Compute the prediction error (number of misclassified examples) on a validation set.
     *
     * @param validation_data_ Validation data set.
     * @param new_stump If non-null, added to the current set of stumps before measuring error (and removed afterwards).
     */
    double computeValidationError(TrainingData const & validation_data_, SharedStump::Ptr new_stump = SharedStump::Ptr());

    /** Load the classifier from an input stream. */
    bool deserialize(std::istream & in);

    /** Save the trained classifier to an output stream. */
    bool serialize(std::ostream & out) const;

    long num_classes;                    ///< Number of object classes.
    long num_features;                   ///< Number of features per object.
    TheaArray<std::string> class_names;  ///< Name of each class, if specified in training data.
    Options options;                     ///< Additional options.

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
