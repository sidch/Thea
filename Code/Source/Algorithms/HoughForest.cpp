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
// First version: 2013
//
//============================================================================

#include "HoughForest.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "../Stopwatch.hpp"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <stack>

#define THEA_HOUGH_SYMMETRIC_VARIANCE

namespace Thea {
namespace Algorithms {

namespace HoughForestInternal {

static intx const BACKGROUND_CLASS = 0;

// A 1-dimensional Gaussian distribution.
class Gaussian1D
{
  private:
    double mean;
    double variance;
    double k;  // normalizing constant

    static double kFromVariance(double var) { return 1.0 / std::sqrt(2 * Math::pi() * var); }

  public:
    Gaussian1D() : mean(0), variance(0), k(0) {}
    Gaussian1D(double mean_, double variance_) : mean(mean_), variance(variance_), k(kFromVariance(variance_)) {}

    double prob(double x) const
    {
      if (variance > 0)
        return k * (std::exp(-Math::square(x - mean) / (2 * variance)));
      else
        return 1.0e-20;
    }

    double probFast(double x) const
    {
      if (variance > 0)
        return k * (Math::fastMinusExp((double)(Math::square(x - mean) / (2 * variance))));
      else
        return 1.0e-20;
    }

    double getMean() const { return mean; }
    double getVariance() const { return variance; }

    void setMean(double value) { mean = value; }
    void setVariance(double value)
    {
      variance = value;
      k = kFromVariance(value);
    }

    void read(BinaryInputStream & in)
    {
      mean = in.readFloat64();
      variance = in.readFloat64();
      k = kFromVariance(variance);
    }

    void write(BinaryOutputStream & out) const
    {
      out.writeFloat64(mean);
      out.writeFloat64(variance);
    }

}; // class Gaussian1D

// An n-dimensional Gaussian distribution.
class Gaussian
{
  private:
    VectorX<double> mean;
    MatrixX<double> inv_cov;

    void resize(intx ndims)
    {
      mean.resize(ndims);
      inv_cov.resize(ndims, ndims);
    }

  public:
    Gaussian() {}
    Gaussian(intx ndims) { resize(ndims); }

    double unnormalizedProb(double const * x) const
    {
      double xCx = 0;
      for (intx i = 0; i < mean.size(); ++i)
      {
        double prod = 0;
        for (intx j = 0; j < mean.size(); ++j)
          prod += (inv_cov(i, j) * (mean[j] - x[j]));

        xCx += (mean[i] - x[i]) * prod;
      }

      return std::exp(-0.5 * xCx);
    }

    double unnormalizedProbFast(double const * x) const
    {
      double xCx = 0;
      for (intx i = 0; i < mean.size(); ++i)
      {
        double prod = 0;
        for (intx j = 0; j < mean.size(); ++j)
          prod += (inv_cov(i, j) * (mean[j] - x[j]));

        xCx += (mean[i] - x[i]) * prod;
      }

      return Math::fastMinusExp(0.5 * xCx);
    }

    intx dims() const { return (intx)mean.size(); }

    VectorX<double> const & getMean() const { return mean; }
    MatrixX<double> const & getInvCov() const { return inv_cov; }

    void setMean(VectorX<double> const & value) { mean = value; }
    void setCov(MatrixX<double> const & cov)
    {
      inv_cov = cov.inverse();  // FIXME: Add invertibility check via methods of Eigen::FullPivLU
    }

    void read(BinaryInputStream & in)
    {
      intx ndims = (intx)in.readInt64();
      resize(ndims);

      for (intx i = 0; i < mean.size(); ++i)
        mean[i] = in.readFloat64();

      for (intx i = 0; i < inv_cov.rows(); ++i)
        for (intx j = 0; j < inv_cov.cols(); ++j)
          inv_cov(i, j) = in.readFloat64();
    }

    void write(BinaryOutputStream & out) const
    {
      out.writeInt64(dims());

      for (intx i = 0; i < mean.size(); ++i)
        out.writeFloat64(mean[i]);

      for (intx i = 0; i < inv_cov.rows(); ++i)
        for (intx j = 0; j < inv_cov.cols(); ++j)
          out.writeFloat64(inv_cov(i, j));
    }

}; // struct Gaussian

// A node of a Hough tree.
class HoughNode
{
  public:
    intx depth;
    intx split_feature;
    double split_value;
    Array<intx> elems;  // Only in leaf nodes. Sorted by class for efficient sampling from a particular class.
    Array<intx> class_cum_freq;  // i'th entry is number of elements with class ID <= i
    Array<Gaussian1D> class_feature_distrib;
    HoughNode * left;
    HoughNode * right;

    HoughNode(intx depth_ = 0)
    : depth(depth_), split_feature(-1), split_value(0), left(nullptr), right(nullptr) {}

    ~HoughNode()
    {
      clear();
    }

    void clear()
    {
      delete left;
      delete right;
      left = right = nullptr;
    }

    bool isLeaf() const
    {
      return !left && !right;
    }

    intx numElements() const { return class_cum_freq.empty() ? 0 : class_cum_freq.back(); }

    intx getClassFrequency(intx class_id) const
    {
      return class_id == 0 ? class_cum_freq[0]
                           : class_cum_freq[(size_t)class_id] - class_cum_freq[(size_t)class_id - 1];
    }

    void readSubtree(BinaryInputStream & in)
    {
      clear();

      depth = (intx)in.readInt64();
      split_feature = (intx)in.readInt64();
      split_value = in.readFloat64();

      int64 num_elems = in.readInt64();
      elems.resize((size_t)num_elems);
      for (size_t i = 0; i < elems.size(); ++i)
        elems[i] = (intx)in.readInt64();

      int64 num_classes = in.readInt64();

      class_cum_freq.resize((size_t)num_classes);
      for (size_t i = 0; i < class_cum_freq.size(); ++i)
        class_cum_freq[i] = (intx)in.readInt64();

      class_feature_distrib.resize((size_t)num_classes);
      for (size_t i = 0; i < class_feature_distrib.size(); ++i)
        class_feature_distrib[i].read(in);

      if (in.readInt8())
      {
        left = new HoughNode(depth + 1);
        left->readSubtree(in);
      }

      if (in.readInt8())
      {
        right = new HoughNode(depth + 1);
        right->readSubtree(in);
      }
    }

    void writeSubtree(BinaryOutputStream & out) const
    {
      out.writeInt64(depth);
      out.writeInt64(split_feature);
      out.writeFloat64(split_value);

      out.writeInt64((int64)elems.size());
      for (size_t i = 0; i < elems.size(); ++i)
        out.writeInt64(elems[i]);

      alwaysAssertM(class_cum_freq.size() == class_feature_distrib.size(),
                    "HoughForest: class_cum_freq.size() != class_feature_distrib.size()");

      out.writeInt64((int64)class_cum_freq.size());

      for (size_t i = 0; i < class_cum_freq.size(); ++i)
        out.writeInt64(class_cum_freq[i]);

      for (size_t i = 0; i < class_feature_distrib.size(); ++i)
        class_feature_distrib[i].write(out);

      if (left)
      {
        out.writeInt8(1);
        left->writeSubtree(out);
      }
      else
        out.writeInt8(0);

      if (right)
      {
        out.writeInt8(1);
        right->writeSubtree(out);
      }
      else
        out.writeInt8(0);
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
    HoughTree(HoughForest * parent_, intx num_classes_, intx num_features_, Array<intx> num_vote_params_,
              Options const & options_)
    : parent(parent_),
      num_classes(num_classes_),
      num_features(num_features_),
      num_vote_params(num_vote_params_),
      max_vote_params(*std::max_element(num_vote_params_.begin(), num_vote_params_.end())),
      options(options_),
      root(nullptr)
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
      root = nullptr;
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

      // Set up root node
      root = new Node(0);
      root->elems.resize((size_t)training_data.numExamples());
      for (intx i = 0; i < training_data.numExamples(); ++i)
        root->elems[(size_t)i] = i;

      measureClassCumulativeFrequencies(root->elems, training_data, root->class_cum_freq);
      root->class_feature_distrib.resize(root->class_cum_freq.size());

      intx num_nodes = 1;
      intx max_depth = 0;

      // Depth-first expansion of the tree
      std::stack<Node *> nodes;
      nodes.push(root);

      while (!nodes.empty())
      {
        Node * node = nodes.top();
        nodes.pop();

        // Split the node
        if (node->depth < options.max_depth
         && (intx)node->elems.size() > options.max_leaf_elements)
        {
          intx split_feature;
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
            Array<double> features;
            getNodeFeatures(node, split_feature, training_data, features);

            Array<double> left_features, right_features;
            for (size_t i = 0; i < node->elems.size(); ++i)
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

            // Measure frequency of each class
            measureClassCumulativeFrequencies(node->left->elems,  training_data, node->left->class_cum_freq);
            measureClassCumulativeFrequencies(node->right->elems, training_data, node->right->class_cum_freq);

            // Estimate distribution of this feature for each class within each child
            estimateClassFeatureDistributions(left_features, node->left->elems, training_data,
                                              node->left->class_feature_distrib);
            estimateClassFeatureDistributions(right_features, node->right->elems, training_data,
                                              node->right->class_feature_distrib);

            if (options.verbose >= 2)
              THEA_CONSOLE << "HoughForest:    - Left node has " << node->left->numElements() << " elements, right node has "
                           << node->right->numElements() << " elements";

#ifdef THEA_HOUGH_SYMMETRIC_VARIANCE
            // Override the data variances to have equal fuzziness on the left and right sides. This is to compensate for
            // situations when all the feature values on one or both sides are identical (e.g. all zero), so the variance is
            // zero.
            for (size_t i = 0; i < node->left->class_feature_distrib.size(); ++i)
            {
              double sep = node->left->class_feature_distrib[i].getMean() - node->right->class_feature_distrib[i].getMean();
              double var = Math::square(sep / 4);
              node->left->class_feature_distrib[i].setVariance(var);
              node->right->class_feature_distrib[i].setVariance(var);
            }
#endif

            // Free up memory at this node
            node->elems.clear();

            // Push children onto stack
            nodes.push(node->left);
            nodes.push(node->right);
          }
        }

        // Leaf nodes contain lists of elements which must be sorted by class for efficient sampling from a particular class
        if (node->isLeaf())
          sortByClass(node->elems, training_data);
      }

      THEA_CONSOLE << "HoughForest:  - Trained tree with " << num_nodes << " node(s) and " << max_depth + 1 << " level(s) from "
                   << training_data.numExamples() << " example(s)";
    }

    // Cast a single vote for a query point with given features. */
    bool singleVoteSelf(intx query_class, double const * features, VoteCallback & callback) const
    {
      Node * curr = root;
      while (curr)
      {
        if (curr->isLeaf())
        {
          // Elements are sorted by class, so we can access the range of elements for the class directly
          intx start_index = (query_class == 0 ? 0 : curr->class_cum_freq[(size_t)query_class - 1]);
          intx end_index = curr->class_cum_freq[(size_t)query_class] - 1;
          if (end_index < start_index)
            return false;

          intx index = Random::common().integer(start_index, end_index);
          double weight = curr->getClassFrequency(query_class) / (double)curr->numElements();

          if (options.verbose >= 2)
          {
            intx actual_index = curr->elems[(size_t)index];
            THEA_CONSOLE << "HoughForest: Reached leaf with " << end_index - start_index + 1 << " element(s) of class "
                         << query_class << " at depth " << curr->depth << ", casting vote for element " << actual_index
                         << " with weight " << weight;

            std::ostringstream oss;
            oss << "HoughForest:  - Element features = [";
            for (intx i = 0; i < num_features; ++i)
            {
              if (i > 0) oss << ", ";
              oss << parent->all_features(actual_index, i);
            }
            oss << ']';
            THEA_CONSOLE << oss.str();
          }

          parent->singleSelfVoteByLookup(curr->elems[(size_t)index], weight, callback);
          break;
        }
        else
        {
          double feat = features[curr->split_feature];
          bool done = false;

          if (options.verbose >= 3)
            THEA_CONSOLE << "HoughForest: Testing feature " << curr->split_feature << ", value " << feat
                         << " against split value " << curr->split_value << " at depth " << curr->depth;

          //===================================================================================================================
          // If either the left or the right subtree has zero samples from the query class, we can just visit the other subtree
          //===================================================================================================================
          intx left_freq   =  curr->left->getClassFrequency(query_class);
          intx right_freq  =  curr->right->getClassFrequency(query_class);
          if (left_freq > 0)
          {
            if (right_freq <= 0)
            {
              curr = curr->left;
              done = true;
            }
          }
          else if (right_freq > 0)
          {
            curr = curr->right;
            done = true;
          }
          else  // this subtree has no samples from the query class
            return false;

          //===================================================================================================================
          // Make a probabilistic choice for which child to step into, for smooth vote distributions. Similar to Gaussian
          // kd-trees [Adams et al. 2009] but the left and right probabilities are evaluated a little differently.
          //===================================================================================================================
          if (!done && options.probabilistic_sampling)
          {
            double p_left   =  curr->left ->class_feature_distrib[(size_t)query_class].probFast(feat);
            double p_right  =  curr->right->class_feature_distrib[(size_t)query_class].probFast(feat);

            p_left   *=  left_freq;
            p_right  *=  right_freq;

            double p_sum = p_left + p_right;
            if (p_sum > 0)
            {
              double coin_toss = Random::common().uniform01();
              if (coin_toss * p_sum < p_left)
              {
                if (options.verbose >= 3)
                  THEA_CONSOLE << "HoughForest: Traversing left with probabilistic sampling on feature " << curr->split_feature;

                curr = curr->left;
              }
              else
              {
                if (options.verbose >= 3)
                  THEA_CONSOLE << "HoughForest: Traversing right with probabilistic sampling on feature "
                               << curr->split_feature;

                curr = curr->right;
              }

              done = true;
            }
          }

          //===================================================================================================================
          // Make the traditional hard choice
          //===================================================================================================================
          if (!done)
          {
            if (feat < curr->split_value)
            {
              if (options.verbose >= 3)
                THEA_CONSOLE << "HoughForest: Traversing left with hard decision on feature " << curr->split_feature;

              curr = curr->left;
            }
            else
            {
              if (options.verbose >= 3)
                THEA_CONSOLE << "HoughForest: Traversing right with hard decision on feature " << curr->split_feature;

              curr = curr->right;
            }
          }
        }
      }

      return true;
    }

    void readNodes(BinaryInputStream & in)
    {
      clear();

      if (in.readInt8())
      {
        root = new Node(0);
        root->readSubtree(in);
      }
    }

    void writeNodes(BinaryOutputStream & out) const
    {
      if (root)
      {
        out.writeInt8(1);
        root->writeSubtree(out);
      }
      else
        out.writeInt8(0);
    }

    void setVerbose(int level)
    {
      options.setVerbose(level);
    }

  private:
    // Find a suitable splitting decision for a node, based on minimizing classification/regression uncertainty.
    bool optimizeTest(Node const * node, TrainingData const & training_data, intx & split_feature, double & split_value,
                      MeasureMode & measure_mode) const
    {
      if (node->elems.size() < 4)  // doesn't make sense splitting smaller sets, and the quadrant index math will fail anyway
        return false;

      if (num_classes <= 1)
        measure_mode = VOTE_UNCERTAINTY;
      else
      {
        measure_mode = (Random::common().coinToss() ? CLASS_UNCERTAINTY : VOTE_UNCERTAINTY);
        if (measure_mode == CLASS_UNCERTAINTY)
        {
          double uncertainty = measureUncertainty(node->elems, training_data, measure_mode);
          if (uncertainty < options.min_class_uncertainty)
            measure_mode = VOTE_UNCERTAINTY;
        }
      }

      Array<double> features, reordered_features;
      double min_uncertainty = -1;
      split_feature = -1;
      split_value = 0;

      intx max_feat_iters    =  std::min((options.num_feature_expansions + 1) * options.max_candidate_features, num_features);
      intx max_thresh_iters  =  std::min(options.max_candidate_thresholds, (intx)node->elems.size());

      for (intx feat_iter = 0; feat_iter < max_feat_iters; ++feat_iter)
      {
        // Generate a random feature
        // TODO: Rewrite this using randIntegersInRange to avoid repeating features
        intx test_feature = (max_feat_iters < num_features ? Random::common().integer() % num_features : feat_iter);

        // Find a suitable split threshold
        getNodeFeatures(node, test_feature, training_data, features);

        for (intx thresh_iter = 0; thresh_iter < max_thresh_iters; ++thresh_iter)
        {
          double test_threshold = 0;
          if (max_thresh_iters < (intx)features.size())
          {
            // Generate a splitting value in the middle half (second and third quadrants) in the sorted order
            size_t index = features.size() / 4 + (Random::common().integer() % (features.size() / 2));

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
            test_threshold = features[(size_t)thresh_iter];
          }

          // Measure the uncertainty after the split
          bool valid_split;
          double uncertainty = measureUncertaintyAfterSplit(node->elems, features, test_threshold, training_data, measure_mode,
                                                            valid_split);
          if (valid_split && (split_feature < 0 || uncertainty < min_uncertainty))
          {
            split_feature = test_feature;
            split_value = refineSplitValue(test_threshold, features);
            min_uncertainty = uncertainty;

            if (options.verbose >= 3)
              THEA_CONSOLE << "HoughForest:      - Updated split feature " << split_feature << ", split value " << split_value
                           << ", uncertainty " << min_uncertainty;
          }
        }

        // If we've reached the end of one expansion and we've found a valid splitting feature, we can stop
        if (split_feature >= 0 && (feat_iter + 1) % options.max_candidate_features == 0)
          break;
      }

      return (split_feature >= 0);
    }

    // Get the value of a single feature for each element in a list.
    void getNodeFeatures(Node const * node, intx feature_index, TrainingData const & training_data,
                         Array<double> & features) const
    {
      features.resize(node->elems.size());
      training_data.getFeatures(feature_index, (intx)node->elems.size(), &node->elems[0], &features[0]);
    }

    // Measure the frequency of occurrence of each class in a set of elements.
    void measureClassFrequencies(Array<intx> const & elems, TrainingData const & training_data,
                                 Array<intx> & class_freq) const
    {
      class_freq.resize((size_t)num_classes);
      std::fill(class_freq.begin(), class_freq.end(), 0);

      Array<intx> classes(elems.size());
      training_data.getClasses((intx)elems.size(), &elems[0], &classes[0]);

      for (size_t i = 0; i < classes.size(); ++i)
        class_freq[(size_t)classes[i]]++;
    }

    // Measure the cumulative frequency of occurrence of classes, in order of ID, in a set of elements.
    void measureClassCumulativeFrequencies(Array<intx> const & elems, TrainingData const & training_data,
                                           Array<intx> & class_cum_freq) const
    {
      measureClassFrequencies(elems, training_data, class_cum_freq);

      for (size_t i = 1; i < class_cum_freq.size(); ++i)
        class_cum_freq[i] = class_cum_freq[i - 1] + class_cum_freq[i];
    }

    // Sort a set of elements by their class ID's.
    void sortByClass(Array<intx> & elems, TrainingData const & training_data) const
    {
      Array<intx> classes(elems.size());
      training_data.getClasses((intx)elems.size(), &elems[0], &classes[0]);

      // Quadratic-time sort but we expect this to be called only by leaf nodes with few elements
      for (size_t i = 0; i < elems.size(); ++i)
        for (size_t j = i + 1; j < elems.size(); ++j)
          if (classes[i] > classes[j])
          {
            std::swap(elems[i], elems[j]);
            std::swap(classes[i], classes[j]);
          }
    }

    // Measure the classification/regression uncertainty of a collection training examples.
    double measureUncertainty(Array<intx> const & elems, TrainingData const & training_data, MeasureMode measure_mode) const
    {
      if (measure_mode == CLASS_UNCERTAINTY)
      {
        // The uncertainty is defined to be the Shannon entropy of the labeling:
        //   Entropy = -\sum_{c \in Classes} p(c) log p(c)

        Array<intx> class_freq;
        measureClassFrequencies(elems, training_data, class_freq);

        double entropy = 0;
        for (size_t i = 0; i < class_freq.size(); ++i)
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
        Array<intx> classes(elems.size());
        training_data.getClasses((intx)elems.size(), &elems[0], &classes[0]);

        // Now measure the vote variance per class
        MatrixX<double> sum_votes(num_classes, max_vote_params); sum_votes.setZero();  // to measure "square of mean"
        Array<double> sum_vote_sqlen((size_t)num_classes, 0.0);                        // to measure "mean of squares"
        Array<intx> class_freq((size_t)num_classes, 0);                                // only count elements that have valid
                                                                                       // Hough votes for their own classes
        Array<double> vote((size_t)max_vote_params);
        for (size_t i = 0; i < elems.size(); ++i)
        {
          if (classes[i] == BACKGROUND_CLASS)  // ignore background class, assumed to have index 0
            continue;

          training_data.getSelfVote(elems[i], &vote[0]);

          // Add up:
          //   - the vote vector
          //   - the square(d length) of the vote
          for (intx j = 0; j < num_vote_params[(size_t)classes[i]]; ++j)
          {
            double v = vote[(size_t)j];
            sum_votes(classes[i], j) += v;
            sum_vote_sqlen[(size_t)classes[i]] += (v * v);
          }

          class_freq[(size_t)classes[i]]++;
        }

        // (Possible) FIXME: The following sum assumes variances for different classes (with possibly different numbers of Hough
        // parameters) are directly comparable.

        double tot_var = 0;
        for (intx i = 0; i < num_classes; ++i)
        {
          if (i == BACKGROUND_CLASS)
            continue;

          intx num_class_members = class_freq[(size_t)i];
          if (num_class_members > 0)
          {
            double sqlen_of_mean = 0;
            for (intx j = 0; j < num_vote_params[(size_t)i]; ++j)
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
    double measureUncertaintyAfterSplit(Array<intx> const & elems, Array<double> const & features, double split_value,
                                        TrainingData const & training_data, MeasureMode measure_mode, bool & valid_split) const
    {
      static Real const MAX_ASYMMETRY = 7;

      valid_split = true;

      if (elems.empty())
      {
        valid_split = false;
        return 0;
      }

      Array<intx> left_elems, right_elems;
      for (size_t i = 0; i < elems.size(); ++i)
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

    // Adjust the split value to bisect the closest pair from the left and right sides.
    double refineSplitValue(double split_value, Array<double> const & features) const
    {
      double left_max = 0, right_min = 0;
      bool found_left = false, found_right = false;
      for (size_t i = 0; i < features.size(); ++i)
      {
        double feat = features[i];
        if (feat < split_value)
        {
          if (!found_left)
          {
            left_max = feat;
            found_left = true;
          }
          else if (feat > left_max)
            left_max = feat;
        }
        else
        {
          if (!found_right)
          {
            right_min = feat;
            found_right = true;
          }
          else if (feat < right_min)
            right_min = feat;
        }
      }

      if (found_left && found_right)
        return 0.5 * (left_max + right_min);
      else
        return split_value;
    }

    // Estimate the distribution of a set of 1D feature values
    void estimateClassFeatureDistributions(Array<double> const & features, Array<intx> const & elems,
                                           TrainingData const & training_data, Array<Gaussian1D> & class_feature_distrib)
    {
      alwaysAssertM(!elems.empty(), "HoughForest: Can't model distribution of empty set");

      Array<intx> classes(elems.size());
      training_data.getClasses((intx)elems.size(), &elems[0], &classes[0]);

      Array<double> sum((size_t)num_classes, 0.0);
      Array<double> sum_squares((size_t)num_classes, 0.0);
      Array<intx>   class_freq((size_t)num_classes, 0);

      for (size_t i = 0; i < elems.size(); ++i)
      {
        size_t c = (size_t)classes[i];
        sum        [c] += features[i];
        sum_squares[c] += (features[i] * features[i]);
        class_freq [c]++;
      }

      class_feature_distrib.resize(class_freq.size());
      for (size_t i = 0; i < class_freq.size(); ++i)
      {
        intx n = class_freq[i];
        if (n > 0)
        {
          double mean = sum[i] / n;
          double mean_squares = sum_squares[i] / n;
          double variance = mean_squares - mean * mean;
          class_feature_distrib[i] = Gaussian1D(mean, variance);
        }
        else
          class_feature_distrib[i] = Gaussian1D(0, 0);
      }
    }

    void estimateClassGaussians(Array<intx> const & elems, TrainingData const & training_data, Array<Gaussian> & gaussians)
    {
      // FIXME: Incomplete

      gaussians.resize((size_t)num_classes);

      Array<intx> classes(elems.size());
      training_data.getClasses((intx)elems.size(), &elems[0], &classes[0]);

      MatrixX<double, MatrixLayout::ROW_MAJOR> all_features(num_features, (intx)elems.size());
      for (intx i = 0; i < num_features; ++i)
        training_data.getFeatures(i, (intx)elems.size(), &elems[0], &all_features(i, 0));

      for (intx i = 0; i < num_classes; ++i)
        estimateClassGaussian(i);
    }

    void estimateClassGaussian(intx c)
    {
      // FIXME: Incomplete
    }

    HoughForest * parent;
    intx num_classes;
    intx num_features;
    Array<intx> num_vote_params;
    intx max_vote_params;
    Options options;
    Node * root;

}; // class HoughTree

double
uncertaintyForDominantFraction(double dominant_fraction, intx num_classes)
{
  double rem_fraction = (1.0 - dominant_fraction) / num_classes;
  if (rem_fraction < 1.0e-10)
    return 0;

  return -(dominant_fraction * std::log(dominant_fraction)
         + (num_classes - 1) * (rem_fraction * std::log(rem_fraction)));
}

} // namespace HoughForestInternal

HoughForest::Options::Options()
: max_depth(-1),
  max_leaf_elements(-1),
  max_candidate_features(-1),
  num_feature_expansions(-1),
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
    // FIXME: This is currently needed because TextInputStream crashes if the file cannot be read
    {
      std::ifstream in(path.c_str());
      if (!in)
      {
        THEA_ERROR << "HoughForest: Could not load options from input file '" << path << '\'';
        return false;
      }
    }

    TextInputStream in(path, configReadSettings());
    read(in);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not load options from input file '%s'", path.c_str())

  return true;
}

bool
HoughForest::Options::save(std::string const & path) const
{
  try
  {
    // FIXME: This is currently needed because TextOutputStream may crash (?) if the file cannot be written (going by the
    // corresponding crash for TextInputStream)
    {
      std::ofstream out(path.c_str());
      if (!out)
      {
        THEA_ERROR << "HoughForest: Could not save options to output file '" << path << '\'';
        return false;
      }
    }

    TextOutputStream out(path, configWriteSettings());
    write(out);
    out.commit();
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not save options to output file '%s'", path.c_str())

  return true;
}

void
HoughForest::Options::read(BinaryInputStream & input, Codec const & codec, bool read_block_header)
{
  (void)read_block_header;  // ignored

  { BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);

    max_depth                 =  (intx)input.readInt64();
    max_leaf_elements         =  (intx)input.readInt64();
    max_candidate_features    =  (intx)input.readInt64();
    num_feature_expansions    =  (intx)input.readInt64();
    max_candidate_thresholds  =  (intx)input.readInt64();
    min_class_uncertainty     =  (double)input.readFloat64();
    max_dominant_fraction     =  (double)input.readFloat64();
    probabilistic_sampling    =  (input.readInt8() != 0);
    verbose                   =  (int)input.readInt32();
  }
}

void
HoughForest::Options::read(TextInputStream & input, Codec const & codec)
{
  *this = defaults();

  while (input.hasMore())
  {
    std::string field = input.readSymbol();
    input.readSymbol("=");

    if (field == "max_depth")
      max_depth = (intx)input.readNumber();
    else if (field == "max_leaf_elements")
      max_leaf_elements = (intx)input.readNumber();
    else if (field == "max_candidate_features")
      max_candidate_features = (intx)input.readNumber();
    else if (field == "num_feature_expansions")
      num_feature_expansions = (intx)input.readNumber();
    else if (field == "max_candidate_thresholds")
      max_candidate_thresholds = (intx)input.readNumber();
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
HoughForest::Options::write(BinaryOutputStream & output, Codec const & codec, bool write_block_header) const
{
  (void)write_block_header;  // ignored

  { BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);

    output.writeInt64(max_depth);
    output.writeInt64(max_leaf_elements);
    output.writeInt64(max_candidate_features);
    output.writeInt64(num_feature_expansions);
    output.writeInt64(max_candidate_thresholds);
    output.writeFloat64(min_class_uncertainty);
    output.writeFloat64(max_dominant_fraction);
    output.writeInt8(probabilistic_sampling ? 1 : 0);
    output.writeInt32(verbose);
  }
}

void
HoughForest::Options::write(TextOutputStream & output, Codec const & codec) const
{
  output.printf("max_depth = %ld\n", max_depth);
  output.printf("max_leaf_elements = %ld\n", max_leaf_elements);
  output.printf("max_candidate_features = %ld\n", max_candidate_features);
  output.printf("num_feature_expansions = %ld\n", num_feature_expansions);
  output.printf("max_candidate_thresholds = %ld\n", max_candidate_thresholds);
  output.printf("min_class_uncertainty = %lg\n", min_class_uncertainty);
  output.printf("max_dominant_fraction = %lg\n", max_dominant_fraction);
  output.printf("probabilistic_sampling = %s\n", (probabilistic_sampling ? "true" : "false"));
  output.printf("verbose = %d\n", verbose);
}

HoughForest::HoughForest(intx num_classes_, intx num_features_, intx const * num_vote_params_, Options const & options_)
: num_classes(num_classes_),
  num_features(num_features_),
  num_vote_params(num_vote_params_, num_vote_params_ + num_classes_),
  max_vote_params(*std::max_element(num_vote_params_, num_vote_params_ + num_classes_)),
  options(options_)
{
  alwaysAssertM(num_classes_ >= 2, "HoughForest: Can't create Hough forest for less than 2 classes");
  alwaysAssertM(num_features_ >= 1, "HoughForest: Can't create Hough forest with no features");

  for (intx i = 0; i < num_classes_; ++i)
  {
    alwaysAssertM(num_vote_params_[i] >= 1,
                  format("HoughForest: Can't create Hough forest with no Hough vote parameters for class %ld", i));
  }
}

HoughForest::HoughForest(std::string const & path)
: num_classes(0),
  num_features(0),
  max_vote_params(0)
{
  if (!load(path))
    throw Error("HoughForest: Could not load Hough forest from '" + path + '\'');
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
HoughForest::autoSelectUnspecifiedOptions(Options & opts_, TrainingData const & training_data) const
{
  if (opts_.max_leaf_elements < 0)
    opts_.max_leaf_elements = 10;

  if (opts_.max_depth < 0)
    opts_.max_depth = 3 * Math::binaryTreeDepth(training_data.numExamples(), opts_.max_leaf_elements);

  if (opts_.max_candidate_features < 0)
  {
    // opts_.max_candidate_features = (intx)std::ceil(num_features / 10.0);
    // opts_.max_candidate_features = num_features;
    opts_.max_candidate_features = (intx)std::ceil(std::sqrt((double)num_features));
  }

  if (opts_.num_feature_expansions < 0)
  {
    opts_.num_feature_expansions = 2;
  }

  if (opts_.max_candidate_thresholds < 0)
    opts_.max_candidate_thresholds = std::max(5L, (intx)std::ceil(training_data.numExamples() / 1000.0));

  if (opts_.min_class_uncertainty < 0)  // we don't need to check max_dominant_fraction, it's not directly used anywhere
  {
    // Split to improve classification only if no class accounts for more than 95% of the examples
    if (opts_.max_dominant_fraction < 0) opts_.max_dominant_fraction = 0.95;
    opts_.min_class_uncertainty = HoughForestInternal::uncertaintyForDominantFraction(opts_.max_dominant_fraction, num_classes);
  }
}

void
HoughForest::train(intx num_trees, TrainingData const & training_data)
{
  alwaysAssertM(num_trees >= 0, "HoughForest: Forest cannot have a negative number of trees");
  alwaysAssertM(training_data.numFeatures() == numFeatures(), "HoughForest: Training data has different number of features");
  alwaysAssertM(training_data.numClasses() == numClasses(), "HoughForest: Training data has different number of classes");

  for (intx i = 0; i < num_classes; ++i)
  {
    alwaysAssertM(training_data.numVoteParameters(i) == numVoteParameters(i),
                  format("HoughForest: Training data has different number of vote parameters for class %ld", i));
  }

  THEA_CONSOLE << "HoughForest: Training forest with " << num_trees << " tree(s)";

  Stopwatch overall_timer;
  overall_timer.tick();

  // Auto-select appropriate values for unspecified options
  Options full_opts = options;
  autoSelectUnspecifiedOptions(full_opts, training_data);

  if (full_opts.verbose)
  {
    TextOutputStream opts_out;
    full_opts.write(opts_out);
    THEA_CONSOLE << "HoughForest: Options:\n" << opts_out.commitToString();
  }

  clear();
  trees.resize((size_t)num_trees);

  Stopwatch timer;
  for (size_t i = 0; i < trees.size(); ++i)
  {
    TreePtr tree(new Tree(this, num_classes, num_features, num_vote_params, full_opts));

    timer.tick();
      tree->train(training_data);
    timer.tock();

    trees[i] = tree;

    if (options.verbose)
      THEA_CONSOLE << "HoughForest:  - Trained tree " << i << " in " << timer.elapsedTime() << 's';
  }

  timer.tick();
    cacheTrainingData(training_data);
  timer.tock();

  if (options.verbose)
    THEA_CONSOLE << "HoughForest:  - Cached training data in " << timer.elapsedTime() << 's';

  overall_timer.tock();
  THEA_CONSOLE << "HoughForest: Training completed in " << overall_timer.elapsedTime() << 's';
}

intx
HoughForest::voteSelf(intx query_class, double const * features, intx num_votes, VoteCallback & callback) const
{
  if (trees.empty())
    return 0;

  intx votes_cast = 0;
  for (intx i = 0; i < num_votes; ++i)
  {
    intx tree_index = Random::common().integer(0, (intx)trees.size() - 1);
    if (trees[(size_t)tree_index]->singleVoteSelf(query_class, features, callback))
      votes_cast++;
  }

  return votes_cast;
}

void
HoughForest::singleSelfVoteByLookup(intx index, double weight, VoteCallback & callback) const
{
  intx c = all_classes[(size_t)index];
  intx nv = num_vote_params[(size_t)c];
  callback(Vote(c,
                nv,
                all_self_votes.data() + index * all_self_votes.cols(),
                weight,
                index,
                num_features,
                all_features.data() + index * num_features));
}

void
HoughForest::cacheTrainingData(TrainingData const & training_data)
{
  intx num_examples = training_data.numExamples();

  // Cache classes
  {
    all_classes.resize((size_t)num_examples);
    training_data.getClasses(&all_classes[0]);
  }

  // Cache features, one example per row
  {
    all_features.resize(num_examples, num_features);
    VectorXd features(num_examples);
    for (intx i = 0; i < num_features; ++i)
    {
      training_data.getFeatures(i, &features[0]);  // gets feature i for all training examples
      all_features.col(i) = features;
    }
  }

  // Cache self-votes, one example per row
  {
    all_self_votes.resize(num_examples, max_vote_params);
    for (intx i = 0; i < num_examples; ++i)
    {
      double * row_start = all_self_votes.data() + i * all_self_votes.cols();
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
    read(in);
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

    write(out);
    out.commit();
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "HoughForest: Could not save to output file '%s'", path.c_str())

  return true;
}

void
HoughForest::read(BinaryInputStream & input, Codec const & codec, bool read_block_header)
{
  using namespace HoughForestInternal;

  clear();

  if (read_block_header)
    input.skip(Codec::BlockHeader::SERIALIZED_LENGTH);  // header is not needed

  options.read(input);

  { BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);

    num_classes = (intx)input.readInt64();
    num_features = (intx)input.readInt64();

    num_vote_params.resize((size_t)num_classes);
    intx max_vote_params = 0;
    for (size_t i = 0; i < num_vote_params.size(); ++i)
    {
      num_vote_params[i] = (intx)input.readInt64();
      if (i == 0 || num_vote_params[i] > max_vote_params)
        max_vote_params = num_vote_params[i];
    }

    intx num_trees = (intx)input.readInt64();
    trees.resize(num_trees);
    for (size_t i = 0; i < trees.size(); ++i)
    {
      trees[i] = TreePtr(new HoughTree(this, num_classes, num_features, num_vote_params, options));
      trees[i]->readNodes(input);
    }

    int64 num_examples = input.readInt64();
    all_classes.resize((size_t)num_examples);
    for (size_t i = 0; i < all_classes.size(); ++i)
      all_classes[i] = (intx)input.readInt64();

    intx nrows = (intx)input.readInt64();
    intx ncols = (intx)input.readInt64();
    all_features.resize(nrows, ncols);
    for (intx i = 0; i < nrows; ++i)
      for (intx j = 0; j < ncols; ++j)
        all_features(i, j) = input.readFloat64();

    nrows = (intx)input.readInt64();
    ncols = (intx)input.readInt64();
    all_self_votes.resize(nrows, ncols);
    for (intx i = 0; i < nrows; ++i)
      for (intx j = 0; j < ncols; ++j)
        all_self_votes(i, j) = input.readFloat64();
  }
}

void
HoughForest::write(BinaryOutputStream & output, Codec const & codec, bool write_block_header) const
{
  using namespace HoughForestInternal;

  Codec::BlockHeader header("HoughF");
  if (write_block_header)
    header.markAndSkip(output);

  options.write(output);

  { BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);

    output.writeInt64(num_classes);
    output.writeInt64(num_features);
    for (size_t i = 0; i < num_vote_params.size(); ++i)
      output.writeInt64(num_vote_params[i]);

    output.writeInt64((int64)trees.size());
    for (size_t i = 0; i < trees.size(); ++i)
      trees[i]->writeNodes(output);

    output.writeInt64((int64)all_classes.size());
    for (size_t i = 0; i < all_classes.size(); ++i)
      output.writeInt64(all_classes[i]);

    output.writeInt64(all_features.rows());
    output.writeInt64(all_features.cols());
    for (intx i = 0; i < all_features.rows(); ++i)
      for (intx j = 0; j < all_features.cols(); ++j)
        output.writeFloat64(all_features(i, j));

    output.writeInt64(all_self_votes.rows());
    output.writeInt64(all_self_votes.cols());
    for (intx i = 0; i < all_self_votes.rows(); ++i)
      for (intx j = 0; j < all_self_votes.cols(); ++j)
        output.writeFloat64(all_self_votes(i, j));
  }

  if (write_block_header)
    header.calcAndWrite(output);
}

void
HoughForest::dumpToConsole() const
{
  THEA_CONSOLE << "HoughForest: Forest has " << trees.size() << " tree(s)";
}

void
HoughForest::setVerbose(int level)
{
  options.setVerbose(level);

  for (size_t i = 0; i < trees.size(); ++i)
    trees[i]->setVerbose(level);
}

} // namespace Algorithms
} // namespace Thea
