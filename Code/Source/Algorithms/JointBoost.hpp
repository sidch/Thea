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

#ifndef __Thea_Algorithms_JointBoost_hpp__
#define __Thea_Algorithms_JointBoost_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Iostream.hpp"
#include "../MatVec.hpp"
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
    THEA_DECL_SMART_POINTERS(JointBoost)

    /** Interface for accessing training data. */
    class THEA_API TrainingData
    {
      public:
        THEA_DECL_SMART_POINTERS(TrainingData)

        /** Destructor. */
        virtual ~TrainingData() {}

        /** Get the number of training examples. */
        virtual intx numExamples() const = 0;

        /** Get the number of features per example. */
        virtual intx numFeatures() const = 0;

        /**
         * Get the values of a particular feature for all training examples. \a feature_index must be in the range
         * 0... numFeatures() - 1.
         */
        virtual void getFeature(intx feature_index, Array<double> & values) const = 0;

        /** Get the class of each training example. */
        virtual void getClasses(Array<intx> & classes) const = 0;

        /** Get the weight of each training example. If an empty array is returned, all weights are set to 1. */
        virtual void getWeights(Array<double> & weights) const { weights.clear(); }

        /**
         * Get the name of each class. If it returns an empty vector, the numeric index of the class will be used as the label
         * name.
         */
        virtual void getClassNames(Array<std::string> & names) const { names.clear(); };

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
        Options & setMinBoostingRounds(intx rounds) { min_boosting_rounds = rounds; return *this; }

        /** Maximum number of boosting rounds (default -1). This also limits the maximum number of stumps added. */
        Options & setMaxBoostingRounds(intx rounds) { max_boosting_rounds = rounds; return *this; }

        /** Minimum error reduction required for boosting rounds to continue beyond min_boosting_rounds (default -1). */
        Options & setMinFractionalErrorReduction(double frac) { min_fractional_error_reduction = frac; return *this; }

        /** Set the fraction of features sampled in a round (default -1). */
        Options & setFeatureSamplingFraction(double frac) { feature_sampling_fraction = frac; return *this; }

        /** Set the maximum number of candidate thresholds generated, as a fraction of the number of features (default -1). */
        Options & setMaxThresholdsFraction(double frac) { max_thresholds_fraction = frac; return *this; }

        /**
         * Set if exhaustive O(2^C) optimization over all possible subsets of classes will be forced or not (default false).
         * Currently, this requires that there be no more classes than the number of bits in an uintx, minus 1. If you have a
         * large number of classes, use the greedy optimization instead. The latter is turned on by default when there are too
         * many classes, or if you call setForceGreedy(true).
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
        bool read(std::istream & in);

        /** Load classifier options from a text input stream. */
        bool read(TextInputStream & in);

        /** Save classifier options to a standard output stream. */
        bool write(std::ostream & out) const;

        /** Save classifier options to a text output stream. */
        bool write(TextOutputStream & out) const;

        /** Get the set of default options. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        intx min_boosting_rounds;  /**< Minimum number of boosting rounds that must be performed even if the error reduction
                                        between successive rounds is below the threshold min_fractional_error_reduction. */
        intx max_boosting_rounds;  ///< Maximum number of boosting rounds. This also limits the maximum number of stumps added.
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
    JointBoost(intx num_classes_, intx num_features_, Options const & options_ = Options::defaults());

    /** Construct a JointBoost classifier by loading it from a file. */
    JointBoost(std::string const & path);

    /** Destructor. */
    ~JointBoost();

    /** Reset the classifier to the initial state. */
    void clear();

    /** Get the number of classes into which objects may fall. The classes are numbered 0 to numClasses() - 1. */
    intx numClasses() const { return num_classes; }

    /** Get the number of features for an object. */
    intx numFeatures() const { return num_features; }

    /** Get the name of a class. */
    std::string getClassName(intx i) const;

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
    intx train(TrainingData const & training_data_, TrainingData const * validation_data_ = nullptr);

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
    intx predict(double const * features, double * class_probabilities = nullptr) const;

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
      THEA_DECL_SMART_POINTERS(SharedStump)

      intx f;               ///< The index of the feature for this stump.
      SharingSet n;         ///< Indices of positive classes for this stump.
      double a;             ///< Regression weight a.
      double b;             ///< Regression weight b.
      double theta;         ///< Cut threshold.
      Array<double> k;      ///< Constants for classes not in the sharing set.

      /** Evaluate the stump for a given feature and class. */
      double operator()(double feature_value, intx class_index) const
      {
        if (n[(SharingSet::size_type)class_index])
          return feature_value > theta ? a : b;
        else
          return k[(size_t)class_index];
      }

      /** Get a printable description of the stump. */
      std::string toString() const;

      /** Load the stump from an input stream. */
      bool read(std::istream & in);

      /** Save the stump to an output stream. */
      bool write(std::ostream & out) const;

    }; // struct SharedStump

    /** Optimize a stump, for the current set of weights, by searching over features and subsets of classes. */
    double optimizeStump(SharedStump & stump, Array<intx> const & stump_classes);

    /** Optimize a stump via exhaustive O(2^C) search over all subsets of classes. */
    double optimizeStumpExhaustive(SharedStump & stump, Array<double> const & stump_features,
                                   Array<intx> const & stump_classes);

    /** Optimize a stump via greedy O(C^2) search over subsets of classes. */
    double optimizeStumpGreedy(SharedStump & stump, Array<double> const & stump_features, Array<intx> const & stump_classes);

    /** Fit stump parameters a, b and theta, for a particular subset of classes. */
    double fitStump(SharedStump & stump, Array<double> const & stump_features, Array<intx> const & stump_classes,
                    intx * num_generated_thresholds = nullptr);

    /**
     * Compute the prediction error (number of misclassified examples) on a validation set.
     *
     * @param validation_data_ Validation data set.
     * @param new_stump If non-null, added to the current set of stumps before measuring error (and removed afterwards).
     */
    double computeValidationError(TrainingData const * validation_data_, SharedStump::Ptr new_stump = SharedStump::Ptr());

    /** Load the classifier from an input stream. */
    bool read(std::istream & in);

    /** Save the trained classifier to an output stream. */
    bool write(std::ostream & out) const;

    intx num_classes;                ///< Number of object classes.
    intx num_features;               ///< Number of features per object.
    Array<std::string> class_names;  ///< Name of each class, if specified in training data.
    Options options;                 ///< Additional options.

    TrainingData const * training_data;  ///< Cached handle to training data, valid only during training.
    double feature_sampling_fraction;  ///< Cached fraction of features test per round, valid only during training.
    double max_thresholds_fraction;  /**< Cached number of candidate thresholds, expressed as a fraction of the number of
                                          features, valid only during training. */

    MatrixX<double> weights;  ///< Weights indexed by (class index, object_index).
    Array<SharedStump::Ptr> stumps;  ///< Current set of selected decision stumps.

}; // class JointBoost

} // namespace Algorithms
} // namespace Thea

#endif
