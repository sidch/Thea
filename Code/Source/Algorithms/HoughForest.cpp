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

// A node of a Hough tree.
struct HoughNode
{
  HoughNode(long depth_ = 0) : depth(depth_), split_feature(-1), split_value(0), left(NULL), right(NULL) {}

  long depth;
  long split_feature;
  double split_value;
  HoughNode * left;
  HoughNode * right;
  TheaArray<long> elems;

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

}; // struct HoughNode

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

            if (options.verbose)
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

            for (array_size_t i = 0; i < node->elems.size(); ++i)
            {
              if (features[i] < split_feature)
                node->left->elems.push_back(node->elems[i]);
              else
                node->right->elems.push_back(node->elems[i]);
            }

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

    // Cast votes for a query point with given features. Returns the number of votes cast. */
    long voteSelf(double const * features, long approx_max_votes, VoteCallback & callback) const
    {
      Node * curr = root;
      while (curr)
      {
        if (curr->isLeaf())
        {
          if (curr->elems.empty())
            return 0;

          double accept_prob = approx_max_votes / (double)curr->elems.size();
          long num_votes = 0;
          for (array_size_t i = 0; i < curr->elems.size(); ++i)
            if (approx_max_votes < 0 || Math::rand01() < accept_prob)
            {
              parent->castSelfVoteByLookup(curr->elems[i], 1.0, callback);
              num_votes++;
            }

          return num_votes;
        }
        else
        {
          if (features[curr->split_feature] < curr->split_value)
            curr = curr->left;
          else
            curr = curr->right;
        }
      }

      return 0;
    }

  private:
    // Find a suitable splitting decision for a node, based on minimizing classification/regression uncertainty.
    bool optimizeTest(Node const * node, TrainingData const & training_data, long & split_feature, double & split_value,
                      MeasureMode & measure_mode) const
    {
      if (node->elems.empty())
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
          double uncertainty = measureUncertaintyAfterSplit(node->elems, features, test_threshold, training_data, measure_mode);
          if (split_feature < 0 || uncertainty < min_uncertainty)
          {
            split_feature = test_feature;
            split_value = test_threshold;
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
    double measureUncertaintyAfterSplit(TheaArray<long> const & elems, TheaArray<double> features, double split_value,
                                        TrainingData const & training_data, MeasureMode measure_mode) const
    {
      TheaArray<long> left_elems, right_elems;
      for (array_size_t i = 0; i < elems.size(); ++i)
      {
        if (features[i] < split_value)
          left_elems.push_back(elems[i]);
        else
          right_elems.push_back(elems[i]);
      }

      return measureUncertainty(left_elems,  training_data, measure_mode)
           + measureUncertainty(right_elems, training_data, measure_mode);
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
  verbose(false)
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
  std::ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "HoughForest: Could not load options from file '" << path << '\'';
    return false;
  }

  return load(in);
}

bool
HoughForest::Options::save(std::string const & path) const
{
  std::ofstream out(path.c_str());
  if (!out)
  {
    THEA_ERROR << "HoughForest: Could not save options to file '" << path << '\'';
    return false;
  }

  return save(out);
}

bool
HoughForest::Options::load(std::istream & in)
{
  in >> verbose;

  if (!in)
  {
    THEA_ERROR << "HoughForest: Error loading options";
    return false;
  }
  else
    return true;
}

bool
HoughForest::Options::save(std::ostream & out) const
{
  out << verbose << std::endl;

  if (!out)
  {
    THEA_ERROR << "HoughForest: Error saving options";
    return false;
  }
  else
    return true;
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
      THEA_CONSOLE << "HoughForest:  - Trained tree " << i << " in " << timer.elapsedTime() << "s";
  }

  overall_timer.tock();
  THEA_CONSOLE << "HoughForest: Training completed in " << overall_timer.elapsedTime() << "s";
}

void
HoughForest::voteSelf(double const * features, long approx_max_votes_per_tree, VoteCallback & callback) const
{
  for (array_size_t i = 0; i < trees.size(); ++i)
    trees[i]->voteSelf(features, approx_max_votes_per_tree, callback);
}

void
HoughForest::castSelfVoteByLookup(long index, double weight, VoteCallback & callback) const
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

  // Cache features
  {
    all_features.resize(num_features, num_examples);
    TheaArray<double> features(num_examples);
    for (long i = 0; i < num_features; ++i)
    {
      training_data.getFeatures(i, &features[0]);
      all_features.setRow(i, &features[0]);
    }
  }

  // Cache self-votes
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
  std::ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "HoughForest: Could not load from file '" << path << '\'';
    return false;
  }

  if (!options.load(in))
    return false;

  if (!load(in))
    return false;

  return true;
}

bool
HoughForest::save(std::string const & path) const
{
  std::ofstream out(path.c_str());
  if (!out)
  {
    THEA_ERROR << "HoughForest: Could not save to file '" << path << '\'';
    return false;
  }

  if (!options.save(out))
    return false;

  out << std::endl;  // break up the output a bit

  if (!save(out))
    return false;

  return true;
}

bool
HoughForest::load(std::istream & in)
{
  return false;
}

bool
HoughForest::save(std::ostream & out) const
{
  return false;
}

void
HoughForest::dumpToConsole() const
{
  THEA_CONSOLE << "HoughForest: Forest has " << trees.size() << " tree(s)";
}

} // namespace Algorithms
} // namespace Thea
