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

#include "HoughForest.hpp"
#include "../Math.hpp"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <stack>

namespace Thea {
namespace Algorithms {

namespace HoughForestInternal {

static long const BACKGROUND_CLASS = 0;

struct Gaussian
{
  double k;
  TheaArray<Real> mean;
  Matrix<double> invcov;

}; // struct Gaussian

class Gaussian1D
{
  private:
    double mean;
    double variance;
    double k;  // normalizing constant

  public:
    Gaussian1D() : mean(0), variance(0), k(0) {}
    Gaussian1D(double mean_, double variance_) : mean(mean_), variance(variance_), k(std::sqrt(2 * Math::pi() * variance)) {}

    double prob(double x)
    {
      if (variance > 0)
        return k * (std::exp(-Math::square(x - mean) / variance));
      else
        return 1.0e-20;
    }

    double probFast(double x)
    {
      if (variance > 0)
        return k * (Math::fastMinusExp((double)(Math::square(x - mean) / variance)));
      else
        return 1.0e-20;
    }

}; // class Gaussian1D

// A node of a Hough tree.
class HoughNode
{
  public:
    HoughNode(long depth_ = 0) : depth(depth_), split_feature(-1), split_value(0), left(NULL), right(NULL) {}

    long depth;
    long split_feature;
    double split_value;
    HoughNode * left;
    HoughNode * right;
    TheaArray<long> elems;
    Gaussian1D feature_distribution;

    ~HoughNode()
    {
      clear();
    }

    void clear()
    {
      delete left;
      delete right;
      left = right = NULL;
    }

    bool isLeaf() const
    {
      return !left && !right;
    }

    void serializeSubtree(BinaryOutputStream & out) const
    {
      out.writeInt64(depth);
      out.writeInt64(split_feature);
      out.writeFloat64(split_value);
    }

    void deserializeSubtree(BinaryInputStream & in)
    {
      clear();

      depth = (long)in.readInt64();
      split_feature = (long)in.readInt64();
      split_value = in.readFloat64();

      if (in.readInt8())
      {
        left = new HoughNode(depth + 1);
        left->deserializeSubtree(in);
      }

      if (in.readInt8())
      {
        right = new HoughNode(depth + 1);
        right->deserializeSubtree(in);
      }
    }

}; // class HoughNode

// A single tree in a Hough forest.
class HoughTree
{
  private:
    enum MeasureMode
    {
      CLASS_UNCERTAINTY,
      VOTE_UNCERTAINTY
    };

  public:
    typedef HoughNode Node;
    typedef HoughForest::TrainingData TrainingData;
    typedef HoughForest::Options Options;
    typedef HoughForest::VoteCallback VoteCallback;

    // Constructor.
    HoughTree(HoughForest * parent_, long num_classes_, long num_features_, TheaArray<long> num_vote_params_,
              Options const & options_)
    : parent(parent_),
      num_classes(num_classes_),
      num_features(num_features_),
      num_vote_params(num_vote_params_),
      max_vote_params(*std::max_element(num_vote_params_.begin(), num_vote_params_.end())),
      options(options_),
      root(NULL)
    {
      alwaysAssertM(parent_, "HoughTree: Can't create tree without parent forest");
    }

    // Destructor.
    ~HoughTree()
    {
      clear();
    }

    // Reset to an empty tree.
    void clear()
    {
      delete root;
      root = NULL;
    }

    // Train the tree.
    void train(TrainingData const & training_data)
    {
      clear();

      if (training_data.numExamples() <= 0)
      {
        THEA_CONSOLE << "HoughForest:  - No training examples provided, tree is empty";
        return;
      }

      root = new Node(0);
      root->elems.resize((array_size_t)training_data.numExamples());
      for (long i = 0; i < training_data.numExamples(); ++i)
        root->elems[(array_size_t)i] = i;

      long num_nodes = 1;
      long max_depth = 0;

      // Depth-first expansion of the tree
      std::stack<Node *> nodes;
      nodes.push(root);

      while (!nodes.empty())
      {
        Node * node = nodes.top();
        nodes.pop();

        // Split the node
        if (node->depth < options.max_depth
         && (long)node->elems.size() > options.max_leaf_elements)
        {
          long split_feature;
          double split_value;
          MeasureMode measure_mode;
          if (optimizeTest(node, training_data, split_feature, split_value, measure_mode))
          {
            // Set up the node to split the given feature at the given value
            node->split_feature = split_feature;
            node->split_value = split_value;

            if (options.verbose >= 2)
              THEA_CONSOLE << "HoughForest:    - Node test splits feature " << split_feature << " at value " << split_value
                           << " to reduce " << (measure_mode == CLASS_UNCERTAINTY ? "class" : "vote") << " uncertainty";

            // Create left (< feature_value) and right (>= feature_value) children
            node->left   =  new Node(node->depth + 1);
            node->right  =  new Node(node->depth + 1);

            // Update stats
            num_nodes += 2;
            if (node->left->depth > max_depth)
              max_depth = node->left->depth;

            // Divide the elements into left and right subtrees, depending on their feature values
            TheaArray<double> features;
            getNodeFeatures(node, split_feature, training_data, features);

            TheaArray<double> left_features, right_features;
            for (array_size_t i = 0; i < node->elems.size(); ++i)
            {
              if (features[i] < split_value)
              {
                node->left->elems.push_back(node->elems[i]);
                left_features.push_back(features[i]);
              }
              else
              {
                node->right->elems.push_back(node->elems[i]);
                right_features.push_back(features[i]);
              }
            }

            if (options.verbose >= 2)
              THEA_CONSOLE << "HoughForest:    - Left node has " << node->left->elems.size() << " elements, right node has "
                           << node->right->elems.size() << " elements";

            // Estimation distribution of this feature within each child
            estimateFeatureDistribution(left_features,   node->left->feature_distribution);
            estimateFeatureDistribution(right_features,  node->right->feature_distribution);

            // Free up memory at this node
            node->elems.clear();

            // Push children onto stack
            nodes.push(node->left);
            nodes.push(node->right);
          }
        }
      }

      THEA_CONSOLE << "HoughForest:  - Trained tree with " << num_nodes << " node(s) and " << max_depth + 1 << " level(s) from "
                   << training_data.numExamples() << " example(s)";
    }

    // Cast a single vote for a query point with given features. */
    bool singleVoteSelf(double const * features, VoteCallback & callback) const
    {
      Node * curr = root;
      while (curr)
      {
        if (curr->isLeaf())
        {
          if (curr->elems.empty())
            return false;

          if (options.verbose >= 2)
            THEA_CONSOLE << "HoughForest: Reached leaf at depth " << curr->depth << ", casting vote";

          long index = Math::randIntegerInRange(0, (long)curr->elems.size() - 1);
          parent->singleSelfVoteByLookup(index, 1.0, callback);
          break;
        }
        else
        {
          // Make a probabilistic choice for which child to step into, for smooth vote distributions
          double feat = features[curr->split_feature];
          bool done = false;
          if (options.probabilistic_sampling)
          {
            double p_left   =  curr->left->feature_distribution.probFast(feat);
            double p_right  =  curr->right->feature_distribution.probFast(feat);

            double p_sum = p_left + p_right;
            if (p_sum > 0)
            {
              double coin_toss = Math::rand01();
              if (coin_toss < p_left / p_sum)
              {
                if (options.verbose >= 2)
                  THEA_CONSOLE << "HoughForest: Traversing left on feature " << curr->split_feature;

                curr = curr->left;
              }
              else
              {
                if (options.verbose >= 2)
                  THEA_CONSOLE << "HoughForest: Traversing right on feature " << curr->split_feature;

                curr = curr->right;
              }

              done = true;
            }
          }

          if (!done)  // make the traditional hard choice
          {
            if (feat < curr->split_value)
            {
              if (options.verbose >= 2)
                THEA_CONSOLE << "HoughForest: Traversing left on feature " << curr->split_feature;

              curr = curr->left;
            }
            else
            {
              if (options.verbose >= 2)
                THEA_CONSOLE << "HoughForest: Traversing right on feature " << curr->split_feature;

              curr = curr->right;
            }
          }
        }
      }

      return true;
    }

    void serializeNodes(BinaryOutputStream & out) const
    {
      out.setEndian(Endianness::LITTLE);
      if (root)
      {
        out.writeInt8(1);
        root->serializeSubtree(out);
      }
      else
        out.writeInt8(0);
    }

    void deserializeNodes(BinaryInputStream & in)
    {
      clear();

      in.setEndian(Endianness::LITTLE);
      if (in.readInt8())
      {
        root = new Node(0);
        root->deserializeSubtree(in);
      }
    }

  private:
    // Find a suitable splitting decision for a node, based on minimizing classification/regression uncertainty.
    bool optimizeTest(Node const * node, TrainingData const & training_data, long & split_feature, double & split_value,
                      MeasureMode & measure_mode) const
    {
      if (node->elems.size() < 4)  // doesn't make sense splitting smaller sets, and the quadrant index math will fail anyway
        return false;

      measure_mode = (std::rand() % 2 == 0 ? CLASS_UNCERTAINTY : VOTE_UNCERTAINTY);
      if (measure_mode == CLASS_UNCERTAINTY)
      {
        double uncertainty = measureUncertainty(node->elems, training_data, measure_mode);
        if (uncertainty < options.min_class_uncertainty)
          measure_mode = VOTE_UNCERTAINTY;
      }

      TheaArray<double> features, reordered_features;
      double min_uncertainty = -1;
      split_feature = -1;
      split_value = 0;

      long max_feat_iters    =  std::min(options.max_candidate_features,   num_features);
      long max_thresh_iters  =  std::min(options.max_candidate_thresholds, (long)node->elems.size());

      for (long feat_iter = 0; feat_iter < max_feat_iters; ++feat_iter)
      {
        // Generate a random feature
        long test_feature = (max_feat_iters < num_features ? std::rand() % num_features : feat_iter);

        // Find a suitable split threshold
        getNodeFeatures(node, test_feature, training_data, features);

        for (long thresh_iter = 0; thresh_iter < max_thresh_iters; ++thresh_iter)
        {
          double test_threshold = 0;
          if (max_thresh_iters < (long)features.size())
          {
            // Generate a splitting value in the middle half (second and third quadrants) in the sorted order
            array_size_t index = features.size() / 4 + (std::rand() % (features.size() / 2));

            if (options.verbose >= 3)
              THEA_CONSOLE << "HoughForest:      - Testing split index " << index << " for feature " << test_feature;

            // Operate on a scratch array so the original features array retains its ordering and can be passed to
            // measureUncertaintyAfterSplit()
            if (thresh_iter == 0)
              reordered_features = features;

            std::nth_element(reordered_features.begin(), reordered_features.begin() + index, reordered_features.end());
            test_threshold = reordered_features[index];
          }
          else
          {
            // We're considering all possible features as thresholds so no need to sort
            test_threshold = features[(array_size_t)thresh_iter];
          }

          // Measure the uncertainty after the split
          bool valid_split;
          double uncertainty = measureUncertaintyAfterSplit(node->elems, features, test_threshold, training_data, measure_mode,
                                                            valid_split);
          if (valid_split && (split_feature < 0 || uncertainty < min_uncertainty))
          {
            split_feature = test_feature;
            split_value = test_threshold;
            min_uncertainty = uncertainty;

            if (options.verbose >= 3)
              THEA_CONSOLE << "HoughForest:      - Updated split feature " << split_feature << ", split value " << split_value
                           << ", uncertainty " << min_uncertainty;
          }
        }
      }

      return (split_feature >= 0);
    }

    // Get the value of a single feature for each element in a list.
    void getNodeFeatures(Node const * node, long feature_index, TrainingData const & training_data,
                         TheaArray<double> & features) const
    {
      features.resize(node->elems.size());
      training_data.getFeatures(feature_index, (long)node->elems.size(), &node->elems[0], &features[0]);
    }

    // Measure the classification/regression uncertainty of a collection training examples.
    double measureUncertainty(TheaArray<long> const & elems, TrainingData const & training_data, MeasureMode measure_mode) const
    {
      if (measure_mode == CLASS_UNCERTAINTY)
      {
        // The uncertainty is defined to be the Shannon entropy of the labeling:
        //   Entropy = -\sum_{c \in Classes} p(c) log p(c)

        TheaArray<long> classes(elems.size());
        training_data.getClasses((long)elems.size(), &elems[0], &classes[0]);

        TheaArray<long> class_freq((array_size_t)num_classes, 0);
        for (array_size_t i = 0; i < classes.size(); ++i)
          class_freq[(array_size_t)classes[i]]++;

        double entropy = 0;
        for (array_size_t i = 0; i < class_freq.size(); ++i)
          if (class_freq[i] > 0)
          {
            double prob = class_freq[i] / (double)elems.size();
            entropy -= prob * std::log(prob);
          }

        return entropy;
      }
      else
      {
        // The uncertainty is defined to be the average variance in Hough votes for the same object. Context votes (elements of
        // one object voting for parameters of another object) are ignored.

        // First get the class of each training example
        TheaArray<long> classes(elems.size());
        training_data.getClasses((long)elems.size(), &elems[0], &classes[0]);

        // Now measure the vote variance per class
        Matrix<double> sum_votes(num_classes, max_vote_params, 0.0);       // to measure "square of mean"
        TheaArray<double> sum_vote_sqlen((array_size_t)num_classes, 0.0);  // to measure "mean of squares"
        TheaArray<long> class_freq((array_size_t)num_classes, 0);          // only count elements that have valid Hough votes
                                                                           // for their own classes
        TheaArray<double> vote((array_size_t)max_vote_params);
        for (array_size_t i = 0; i < elems.size(); ++i)
        {
          if (classes[i] == BACKGROUND_CLASS)  // ignore background class, assumed to have index 0
            continue;

          training_data.getSelfVote(elems[i], &vote[0]);

          // Add up:
          //   - the vote vector
          //   - the square(d length) of the vote
          for (long j = 0; j < num_vote_params[(array_size_t)classes[i]]; ++j)
          {
            double v = vote[(array_size_t)j];
            sum_votes(classes[i], j) += v;
            sum_vote_sqlen[(array_size_t)classes[i]] += (v * v);
          }

          class_freq[(array_size_t)classes[i]]++;
        }

        // (Possible) FIXME: The following sum assumes variances for different classes (with possibly different numbers of Hough
        // parameters) are directly comparable.

        double tot_var = 0;
        for (long i = 0; i < num_classes; ++i)
        {
          if (i == BACKGROUND_CLASS)
            continue;

          long num_class_members = class_freq[(array_size_t)i];
          if (num_class_members > 0)
          {
            double sqlen_of_mean = 0;
            for (long j = 0; j < num_vote_params[(array_size_t)i]; ++j)
            {
              double v = sum_votes(i, j);
              sqlen_of_mean += (v * v);
            }
            sqlen_of_mean /= (num_class_members * num_class_members);

            // Apply: variance = mean of squares - square of mean
            tot_var += (sum_vote_sqlen[i] / num_class_members - sqlen_of_mean);
          }
        }

        return tot_var / (num_classes - 1);  // ignore background class
      }
    }

    // Measure the classification/regression uncertainty after splitting elements along a feature.
    double measureUncertaintyAfterSplit(TheaArray<long> const & elems, TheaArray<double> const & features, double split_value,
                                        TrainingData const & training_data, MeasureMode measure_mode, bool & valid_split) const
    {
      static Real const MAX_ASYMMETRY = 7;

      valid_split = true;

      if (elems.empty())
      {
        valid_split = false;
        return 0;
      }

      TheaArray<long> left_elems, right_elems;
      for (array_size_t i = 0; i < elems.size(); ++i)
      {
        if (features[i] < split_value)
          left_elems.push_back(elems[i]);
        else
          right_elems.push_back(elems[i]);
      }

      // If the split is too asymmetrical, it is not valid
      if (left_elems.size()  > MAX_ASYMMETRY * right_elems.size()
       || right_elems.size() > MAX_ASYMMETRY * left_elems.size())
      {
        valid_split = false;
        return 0;
      }

      if (options.verbose >= 3)
        THEA_CONSOLE << "HoughForest:        - Split uncertainty test: left node has " << left_elems.size()
                     << " elements, right node has " << right_elems.size() << " elements";

      return measureUncertainty(left_elems,  training_data, measure_mode)
           + measureUncertainty(right_elems, training_data, measure_mode);
    }

    // Estimate the distribution of a set of 1D feature values
    void estimateFeatureDistribution(TheaArray<double> const & features, Gaussian1D & distrib)
    {
      alwaysAssertM(!features.empty(), "HoughForest: Can't model distribution of empty set");

      double sum = 0, sum_squares = 0;
      for (array_size_t i = 0; i < features.size(); ++i)
      {
        sum += features[i];
        sum_squares += (features[i] * features[i]);
      }

      double mean = sum / features.size();
      double mean_squares = sum_squares / features.size();
      double variance = mean_squares - mean * mean;

      distrib = Gaussian1D(mean, variance);
    }

    void estimateClassGaussians(TheaArray<long> const & elems, TrainingData const & training_data,
                                TheaArray<Gaussian> & gaussians)
    {
      // NOTE: Incomplete

      gaussians.resize((array_size_t)num_classes);

      TheaArray<long> classes(elems.size());
      training_data.getClasses((long)elems.size(), &elems[0], &classes[0]);

      Matrix<double, MatrixLayout::ROW_MAJOR> all_features(num_features, (long)elems.size());
      for (long i = 0; i < num_features; ++i)
        training_data.getFeatures(i, (long)elems.size(), &elems[0], &all_features(i, 0));

      for (long i = 0; i < num_classes; ++i)
        estimateClassGaussian(i);
    }

    void estimateClassGaussian(long c)
    {
      // NOTE: Incomplete
    }

    HoughForest * parent;
    long num_classes;
    long num_features;
    TheaArray<long> num_vote_params;
    long max_vote_params;
    Options options;
    Node * root;

}; // class HoughTree

double
uncertaintyForDominantFraction(double dominant_fraction, long num_classes)
{
  double rem_fraction = (1.0 - dominant_fraction) / num_classes;

  return -(dominant_fraction * std::log(dominant_fraction)
         + (num_classes - 1) * (rem_fraction * std::log(rem_fraction)));
}

} // namespace HoughForestInternal

HoughForest::Options::Options()
: max_depth(-1),
  max_leaf_elements(-1),
  max_candidate_features(-1),
  max_candidate_thresholds(-1),
  min_class_uncertainty(-1),
  max_dominant_fraction(-1),
  probabilistic_sampling(true),
  verbose(1)
{
}

HoughForest::Options &
HoughForest::Options::setMinClassUncertainty(double value)
{
  min_class_uncertainty = value;
  max_dominant_fraction = -1;
  return *this;
}

HoughForest::Options &
HoughForest::Options::setMaxDominantFraction(double value)
{
  min_class_uncertainty = -1;
  max_dominant_fraction = value;
  return *this;
}

bool
HoughForest::Options::load(std::string const & path)
{
  try
  {
    TextInputStream in(path, configReadSettings());
    deserialize(in);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not load options from input file '%s'", path.c_str())

  return true;
}

bool
HoughForest::Options::save(std::string const & path) const
{
  try
  {
    TextOutputStream out(path, configWriteSettings());
    serialize(out);
    out.commit();
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not save options to output file '%s'", path.c_str())

  return true;
}

void
HoughForest::Options::deserialize(BinaryInputStream & input, Codec const & codec)
{
  max_depth                 =  (long)input.readInt64();
  max_leaf_elements         =  (long)input.readInt64();
  max_candidate_features    =  (long)input.readInt64();
  max_candidate_thresholds  =  (long)input.readInt64();
  min_class_uncertainty     =  (double)input.readFloat64();
  max_dominant_fraction     =  (double)input.readFloat64();
  probabilistic_sampling    =  (input.readInt8() != 0);
  verbose                   =  (int)input.readInt32();
}

void
HoughForest::Options::deserialize(TextInputStream & input, Codec const & codec)
{
  *this = defaults();

  while (input.hasMore())
  {
    std::string field = input.readSymbol();
    input.readSymbol("=");

    if (field == "max_depth")
      max_depth = (long)input.readNumber();
    else if (field == "max_leaf_elements")
      max_leaf_elements = (long)input.readNumber();
    else if (field == "max_candidate_features")
      max_candidate_features = (long)input.readNumber();
    else if (field == "max_candidate_thresholds")
      max_candidate_thresholds = (long)input.readNumber();
    else if (field == "min_class_uncertainty")
      min_class_uncertainty = input.readNumber();
    else if (field == "max_dominant_fraction")
      max_dominant_fraction = input.readNumber();
    else if (field == "probabilistic_sampling")
      probabilistic_sampling = input.readBoolean();
    else if (field == "verbose")
      verbose = (int)input.readNumber();
  }
}

void
HoughForest::Options::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  output.writeInt64(max_depth);
  output.writeInt64(max_leaf_elements);
  output.writeInt64(max_candidate_features);
  output.writeInt64(max_candidate_thresholds);
  output.writeFloat64(min_class_uncertainty);
  output.writeFloat64(max_dominant_fraction);
  output.writeInt8(probabilistic_sampling ? 1 : 0);
  output.writeInt32(verbose);
}

void
HoughForest::Options::serialize(TextOutputStream & output, Codec const & codec) const
{
  output.printf("max_depth = %ld\n", max_depth);
  output.printf("max_leaf_elements = %ld\n", max_leaf_elements);
  output.printf("max_candidate_features = %ld\n", max_candidate_features);
  output.printf("max_candidate_thresholds = %ld\n", max_candidate_thresholds);
  output.printf("min_class_uncertainty = %lf\n", min_class_uncertainty);
  output.printf("max_dominant_fraction = %lf\n", max_dominant_fraction);
  output.printf("probabilistic_sampling = %s\n", (probabilistic_sampling ? "true" : "false"));
  output.printf("verbose = %d\n", verbose);
}

HoughForest::HoughForest(long num_classes_, long num_features_, long const * num_vote_params_, Options const & options_)
: num_classes(num_classes_),
  num_features(num_features_),
  num_vote_params(num_vote_params_, num_vote_params_ + num_classes_),
  max_vote_params(*std::max_element(num_vote_params_, num_vote_params_ + num_classes_)),
  options(options_)
{
  alwaysAssertM(num_classes_ >= 2, "HoughForest: Can't create Hough forest for less than 2 classes");
  alwaysAssertM(num_features_ >= 1, "HoughForest: Can't create Hough forest with no features");

  for (long i = 0; i < num_classes_; ++i)
  {
    alwaysAssertM(num_vote_params_[i] >= 1,
                  format("HoughForest: Can't create Hough forest with no Hough vote parameters for class %ld", i));
  }
}

HoughForest::~HoughForest()
{
  clear();
}

void
HoughForest::clear()
{
  trees.clear();

  all_classes.clear();
  all_features.resize(0, 0);
  all_self_votes.resize(0, 0);
}

void
HoughForest::autoSelectUnspecifiedOptions(Options & opts_, TrainingData const * training_data_) const
{
  if (opts_.max_leaf_elements < 0)
    opts_.max_leaf_elements = 10;

  if (opts_.max_depth < 0)
  {
    if (training_data_)
      opts_.max_depth = Math::binaryTreeDepth(training_data_->numExamples(), opts_.max_leaf_elements);
    else
      opts_.max_depth = 10;
  }

  if (opts_.max_candidate_features < 0)
    opts_.max_candidate_features = (long)std::ceil(num_features / 10.0);

  if (opts_.max_candidate_thresholds < 0)
  {
    if (training_data_)
      opts_.max_candidate_thresholds = std::max(5L, (long)std::ceil(training_data_->numExamples() / 1000.0));
    else
      opts_.max_candidate_thresholds = 5;
  }

  if (opts_.min_class_uncertainty < 0)  // we don't need to check max_dominant_fraction, it's not directly used anywhere
  {
    // Split to improve classification only if no class accounts for more than 95% of the examples
    if (opts_.max_dominant_fraction < 0) opts_.max_dominant_fraction = 0.95;
    opts_.min_class_uncertainty = HoughForestInternal::uncertaintyForDominantFraction(opts_.max_dominant_fraction, num_classes);
  }
}

void
HoughForest::train(long num_trees, TrainingData const & training_data_)
{
  alwaysAssertM(num_trees >= 0, "HoughForest: Forest cannot have a negative number of trees");
  alwaysAssertM(training_data_.numFeatures() == numFeatures(), "HoughForest: Training data has different number of features");
  alwaysAssertM(training_data_.numClasses() == numClasses(), "HoughForest: Training data has different number of classes");

  for (long i = 0; i < num_classes; ++i)
  {
    alwaysAssertM(training_data_.numVoteParameters(i) == numVoteParameters(i),
                  format("HoughForest: Training data has different number of vote parameters for class %ld", i));
  }

  THEA_CONSOLE << "HoughForest: Training forest with " << num_trees << " tree(s)";

  G3D::Stopwatch overall_timer;
  overall_timer.tick();

  // Auto-select appropriate values for unspecified options
  Options full_opts = options;
  autoSelectUnspecifiedOptions(full_opts, &training_data_);

  clear();
  trees.resize((array_size_t)num_trees);

  G3D::Stopwatch timer;
  for (array_size_t i = 0; i < trees.size(); ++i)
  {
    TreePtr tree(new Tree(this, num_classes, num_features, num_vote_params, full_opts));

    timer.tick();
      tree->train(training_data_);
    timer.tock();

    trees[i] = tree;

    if (options.verbose)
      THEA_CONSOLE << "HoughForest:  - Trained tree " << i << " in " << timer.elapsedTime() << 's';
  }

  timer.tick();
    cacheTrainingData(training_data_);
  timer.tock();

  if (options.verbose)
    THEA_CONSOLE << "HoughForest:  - Cached training data in " << timer.elapsedTime() << 's';

  overall_timer.tock();
  THEA_CONSOLE << "HoughForest: Training completed in " << overall_timer.elapsedTime() << 's';
}

long
HoughForest::voteSelf(double const * features, long num_votes, VoteCallback & callback) const
{
  if (trees.empty())
    return 0;

  long votes_cast = 0;
  for (long i = 0; i < num_votes; ++i)
  {
    long tree_index = Math::randIntegerInRange(0, (long)trees.size() - 1);
    if (trees[(array_size_t)tree_index]->singleVoteSelf(features, callback))
      votes_cast++;
  }

  return votes_cast;
}

void
HoughForest::singleSelfVoteByLookup(long index, double weight, VoteCallback & callback) const
{
  long c = all_classes[(array_size_t)index];
  long nv = num_vote_params[(array_size_t)c];
  callback(c, nv, all_self_votes.data() + index * all_self_votes.numColumns(), weight);
}

void
HoughForest::cacheTrainingData(TrainingData const & training_data)
{
  long num_examples = training_data.numExamples();

  // Cache classes
  {
    all_classes.resize((array_size_t)num_examples);
    training_data.getClasses(&all_classes[0]);
  }

  // Cache features, one example per column
  {
    all_features.resize(num_features, num_examples);
    TheaArray<double> features(num_examples);
    for (long i = 0; i < num_features; ++i)
    {
      training_data.getFeatures(i, &features[0]);
      all_features.setRow(i, &features[0]);
    }
  }

  // Cache self-votes, one example per row
  {
    all_self_votes.resize(num_examples, max_vote_params);
    for (long i = 0; i < num_examples; ++i)
    {
      double * row_start = all_self_votes.data() + i * all_self_votes.numColumns();
      training_data.getSelfVote(i, row_start);
    }
  }
}

bool
HoughForest::load(std::string const & path)
{
  try
  {
    BinaryInputStream in(path, Endianness::LITTLE);
    deserialize(in);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not load from input file '%s'", path.c_str())

  return true;
}

bool
HoughForest::save(std::string const & path) const
{
  try
  {
    BinaryOutputStream out(path, Endianness::LITTLE);
    if (!out.ok())
    {
      THEA_ERROR << "HoughForest: Could not save to output file '" << path << '\'';
    }

    serialize(out);
    out.commit();
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not save to output file '%s'", path.c_str())

  return true;
}

void
HoughForest::deserialize(BinaryInputStream & input, Codec const & codec)
{
  options.deserialize(input);

  input.setEndian(Endianness::LITTLE);

  // Write nodes in
}

void
HoughForest::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  options.serialize(output);

  output.setEndian(Endianness::LITTLE);
}

void
HoughForest::dumpToConsole() const
{
  THEA_CONSOLE << "HoughForest: Forest has " << trees.size() << " tree(s)";
}

} // namespace Algorithms
} // namespace Thea
