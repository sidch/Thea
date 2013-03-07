//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_HoughForest_hpp__
#define __Thea_Algorithms_HoughForest_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include <iostream>

namespace Thea {
namespace Algorithms {

namespace HoughForestInternal {

// Forward declaration of Hough tree class
class HoughTree;

} // namespace HoughForestInternal

/**
 * An implementation of class-specific Hough forests. Based on
 *
 * J. Gall and V. Lempitsky, "Class-Speciﬁc Hough Forests for Object Detection", Pro. CVPR, 2009.
 *
 * To use the class, implement an appropriate subclass of TrainingData, call train(), and then call vote().
 */
class THEA_API HoughForest
{
  public:
    THEA_DEF_POINTER_TYPES(HoughForest, shared_ptr, weak_ptr)

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

        /** Get the number of parameters (dimensions) of the Hough space for a particular class. */
        virtual long numVoteParameters(long class_index) const = 0;

        /**
         * Get the values of a particular feature for a subset of training examples. \a feature_index must be in the range
         * 0... numFeatures() - 1.
         *
         * @param feature_index The index of the requested feature.
         * @param num_selected_examples The size of the selected subset.
         * @param selected_examples The indices of the selected subset.
         * @param values Used to return the feature values (assumed to be pre-allocated to \a num_selected_examples elements).
         */
        virtual void getFeatures(long feature_index, long num_selected_examples, long const * selected_examples,
                                 double * values) const = 0;

        /**
         * Get the classes of a subset of training examples.
         *
         * @param num_selected_examples The size of the selected subset.
         * @param selected_examples The indices of the selected subset.
         * @param classes Used to return the classes (assumed to be pre-allocated to \a num_selected_examples elements).
         */
        virtual void getClasses(long num_selected_examples, long const * selected_examples, long * classes) const = 0;

        /**
         * Get the parameters of a Hough vote by a particular example for its parent object.
         *
         * @param example_index The index of the voting example.
         * @param params Used to return the parameters of the Hough vote, assumed to be preallocated to the appropriate number
         *   of dimensions (see numVoteParameters()).
         *
         * @return True if a valid vote was found, else false.
         */
        virtual bool getSelfVote(long example_index, double * params) const = 0;

    }; // class TrainingData

    /**
     * %Options for the forest. In most cases, passing a negative value for a normally non-negative parameter sets the default
     * value for that parameter.
     */
    class Options
    {
      public:
        /** Constructor. Sets default options. */
        Options();

        /** Set the maximum depth of the tree (default -1). */
        Options & setMaxDepth(long value) { max_depth = value; return *this; }

        /**
         * Set the maximum number of elements in a leaf node of the tree, unless the maximum depth has been reached (default
         * -1).
         */
        Options & setMaxLeafElements(long value) { max_leaf_elements = value; return *this; }

        /** Set the maximum number of features to consider for splitting per iteration (default -1). */
        Options & setMaxCandidateFeatures(long value) { max_candidate_features = value; return *this; }

        /** Set the maximum number of randomly selected thresholds to consider for splitting along a feature (default -1). */
        Options & setMaxCandidateThresholds(long value) { max_candidate_thresholds = value; return *this; }

        /** Set the minimum class uncertainty required to split a node by class uncertainty (default -1). */
        Options & setMinClassUncertainty(double value) { min_class_uncertainty = value; return *this; }

        /** Set whether progress information will be printed to the console or not (default true). */
        Options & setVerbose(bool value) { verbose = value; return *this; }

        /** Load forest options from a disk file. */
        bool load(std::string const & path);

        /** Save forest options to a disk file. */
        bool save(std::string const & path) const;

        /** Get the set of default options. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        /** Load forest options from an input stream. */
        bool load(std::istream & in);

        /** Save forest options to an output stream. */
        bool save(std::ostream & out) const;

        long max_depth;                   ///< Maximum depth of tree.
        long max_leaf_elements;           ///< Maximum number of elements in leaf node, unless the maximum depth is reached.
        long max_candidate_features;      ///< Maximum number of features to consider for splitting per iteration.
        long max_candidate_thresholds;    ///< Maximum number of randomly selected thresholds for splitting along a feature.
        double min_class_uncertainty;     ///< Minimum class uncertainty required to split a node by class uncertainty.
        bool verbose;                     ///< Print progress information to the console.

        friend class HoughForest;
        friend class HoughForestInternal::HoughTree;

    }; // class Options

    /** Interface for a callback that is called for each Hough vote. */
    class VoteCallback
    {
      public:
        /** Destructor. */
        virtual ~VoteCallback() {}

        /**
         * Function called for each Hough vote. Specialize this to suit your needs.
         *
         * @param target_class The class for which the Hough vote is cast (NOT the predicted class of the voting element).
         * @param num_vote_params Number of parameters (dimensions) of Hough voting space for the target class. Matches the
         *   value used to construct the Hough forest, provided again here for convenience.
         * @param vote_params The array of Hough vote parameters, with \a num_vote_params values.
         * @param weight The weight of the vote.
         */
        virtual bool operator()(long target_class, long num_vote_params, double const * vote_params, double weight) = 0;
    };

    /**
     * Constructor.
     *
     * @param num_classes_ Number of classes for classification. The classes are numbered 0 to \a num_classes - 1.
     * @param num_features_ Number of features per object.
     * @param num_vote_params_ Number of parameters (dimensions) of Hough space per class.
     * @param options_ Addition options controlling the behaviour of the forest.
     */
    HoughForest(long num_classes_, long num_features_, long const * num_vote_params_,
                Options const & options_ = Options::defaults());

    /** Destructor. */
    ~HoughForest();

    /** Reset the forest to the initial state. */
    void clear();

    /** Get the number of classes into which objects may fall. The classes are numbered 0 to numClasses() - 1. */
    long numClasses() const { return num_classes; }

    /** Get the number of features for an object. */
    long numFeatures() const { return num_features; }

    /** Get the number of parameters (dimensions) of the Hough voting space for a given class. */
    long numVoteParameters(long class_index) const { return num_vote_params[(array_size_t)class_index]; }

    /**
     * Train the Hough forest.
     *
     * @param num_trees Number of trees in the forest.
     * @param training_data_ Data used for training the forest.
     */
    void train(long num_trees, TrainingData const & training_data_);

    /**
     * Sample the Hough votes for an object with a given set of features.
     *
     * @param features The features of the object. Must contain numFeatures() values.
     * @param max_votes Maximum number of votes to sample. The votes are drawn evenly from the trees in the forest.
     * @param callback Called once for every vote.
     *
     * @return The index of the most likely class of the object.
     */
    void vote(double const * features, long max_votes, VoteCallback & callback) const;

    /** Load the forest from a disk file. */
    bool load(std::string const & path);

    /** Save the trained forest to disk. */
    bool save(std::string const & path) const;

    /** Print debugging information about this forest to the console. */
    void dumpToConsole() const;

  private:
    friend class HoughForestInternal::HoughTree;

    typedef HoughForestInternal::HoughTree Tree;  ///< Hough tree class.
    typedef shared_ptr<Tree> TreePtr;  ///< Shared pointer to a Hough tree.

    /** Automatically choose suitable values for unspecified options. */
    void autoSelectUnspecifiedOptions(Options & options_, TrainingData const * training_data_) const;

    /** Load the forest from an input stream. */
    bool load(std::istream & in);

    /** Save the trained forest to an output stream. */
    bool save(std::ostream & out) const;

    long num_classes;                 ///< Number of object classes.
    long num_features;                ///< Number of features per object.
    TheaArray<long> num_vote_params;  ///< Number of Hough parameters per class.
    Options options;                  ///< Additional options.

    TheaArray<TreePtr> trees;  ///< Set of Hough trees in the forest.

}; // class HoughForest

} // namespace Algorithms
} // namespace Thea

#endif
