#include "../Common.hpp"
#include "../Algorithms/JointBoost.hpp"
#include "../StringAlg.hpp"
#include "../UnorderedMap.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define JB_FAST
// #define JB_VERBOSE

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testJointBoostFile(string const & path);

int
main(int argc, char * argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " <training-file>" << endl;
    return 0;
  }

  try
  {
    if (!testJointBoostFile(argv[1]))
      return -1;
  }
  catch(...)
  {
    cerr << "An error occurred" << endl;
    return -1;
  }

  return 0;
}

class ExampleSet : public JointBoost::TrainingData
{
  public:
    typedef std::shared_ptr<ExampleSet> Ptr;

    ExampleSet(intx nfeatures) : features((intx)0, nfeatures) {}

    template <typename DVecT>
    void addExample(DVecT const & example_features, intx example_class)
    {
      addExample(&example_features[0], example_class);
    }

    void addExample(double const * example_features, intx example_class)
    {
      features.conservativeResize(features.rows() + 1, Eigen::NoChange);
      features.row(features.rows() - 1) = RowVectorXd::Map(example_features, features.cols());
      classes.push_back(example_class);
    }

    intx numExamples() const { return features.rows(); }
    intx numFeatures() const { return features.cols(); }

    void getFeature(intx feature_index, Array<double> & values) const
    {
      values.resize((size_t)features.rows());
      VectorXd::Map(&values[0], features.rows()) = features.col(feature_index);
    }

    void getClasses(Array<intx> & classes_) const
    {
      classes_ = classes;
    }

    void getExampleFeatures(intx example, double * f) const
    {
      RowVectorXd::Map(f, features.cols()) = features.row(example);
    }

    intx getExampleClass(intx example) const
    {
      return classes[(size_t)example];
    }

  private:
    MatrixX<double, MatrixLayout::ROW_MAJOR> features;
    vector<intx> classes;
};

template <typename ArrayT>
string
arrayToString(ArrayT const & arr)
{
  ostringstream oss; oss << '[';
  for (intx i = 0; i < (intx)arr.size(); ++i)
  {
    if (i > 0) oss << ", ";
    oss << arr[i];
  }
  oss << ']';

  return oss.str();
}

bool
test(JointBoost const & jb, ExampleSet const & test_set, bool get_class_probs = false)
{
  vector<double> class_probabilities((size_t)jb.numClasses());
  vector<double> f((size_t)jb.numFeatures());
  intx num_correct = 0;
  for (intx i = 0; i < test_set.numExamples(); ++i)
  {
    test_set.getExampleFeatures(i, &f[0]);
    intx c = test_set.getExampleClass(i);

    intx pc = jb.predict(&f[0], get_class_probs ? &class_probabilities[0] : nullptr);
    if (pc == c)
      num_correct++;

#ifdef JB_VERBOSE
    cout << "Features[" << i << "] = " << arrayToString(f) << ", actual class = " << c << ", predicted class = " << pc << endl;
    if (get_class_probs)
    {
      for (size_t c = 0; c < class_probabilities.size(); ++c)
        cout << "    Probability of class " << c << " = " << class_probabilities[c] << endl;

      cout << endl;
    }
#endif
  }

#ifndef JB_VERBOSE
  cout << endl;
#endif

  cout << "Accuracy on " << test_set.numExamples() << " test cases = " << 100.0 * num_correct / (float)test_set.numExamples()
       << '%' << endl << endl;

  return true;
}

bool
testJointBoostFile(string const & path)
{
  ifstream in(path.c_str());
  if (!in)
  {
    cerr << "Couldn't open file " << path << endl;
    return false;
  }

  cout << "======================================" << endl;
  cout << "Reading features from file" << endl;
  cout << "======================================" << endl << endl;

  ExampleSet::Ptr all_training;
  ExampleSet::Ptr training_subset;
  ExampleSet::Ptr holdout_subset;

  typedef UnorderedMap<string, intx> LabelIndexMap;
  LabelIndexMap labels;

  vector<double> features;
  double feature;
  string label;

  string line;
  while (getline(in, line))
  {
    line = trimWhitespace(line);
    if (line.empty())
      continue;

    vector<string> fields;
    stringSplit(line, ",\t", fields, /* skip_empty_fields = */ false);

    if (fields.size() < 2)
    {
      cerr << "Data has too few features per example" << endl;
      return false;
    }

    intx nfeat = (intx)fields.size() - 1;
    if (!all_training)
    {
      cout << "Data has " << nfeat << " features per example" << endl;
      all_training = ExampleSet::Ptr(new ExampleSet(nfeat));
    }
    else
    {
      if (nfeat != all_training->numFeatures())
      {
        cout << "Inconsistent number of features for example " << all_training->numExamples() << endl;
        return false;
      }
    }

    features.clear();
    for (intx i = 0; i < nfeat; ++i)
    {
      istringstream iss(fields[(size_t)i]);
      iss >> feature;
      features.push_back(feature);
    }

    label = fields.back();

    LabelIndexMap::const_iterator existing_label = labels.find(label);
    intx index;
    if (existing_label == labels.end())
    {
      index = (intx)labels.size();
      labels[label] = index;
      cout << "Added class with label '" << label << "' and index " << index << endl;
    }
    else
      index = existing_label->second;

    all_training->addExample(features, index);

    if (rand() % 2)  // holdout half the input set for testing
    {
      if (!training_subset) training_subset = ExampleSet::Ptr(new ExampleSet(nfeat));
      training_subset->addExample(features, index);
    }
    else
    {
      if (!holdout_subset) holdout_subset = ExampleSet::Ptr(new ExampleSet(nfeat));
      holdout_subset->addExample(features, index);
    }
  }

  if (!all_training)
  {
    cout << "Could not read any lines from file" << endl;
    return false;
  }

  intx num_classes = (intx)labels.size();

  cout << "Read " << all_training->numExamples() << " examples from " << num_classes << " classes from file" << endl;

  JointBoost::Options opts;

#ifdef JB_FAST
  opts.setMinBoostingRounds(num_classes)
      .setMaxBoostingRounds(4 * num_classes)
      .setMinFractionalErrorReduction(-1)
      .setFeatureSamplingFraction(3.0 / all_training->numFeatures())
      .setMaxThresholdsFraction(0.25)
      .setForceGreedy(true)
      .setVerbose(false);

  // Options for bupa
  // opts.setMinBoostingRounds(10 * num_classes)
  //     .setMaxBoostingRounds(40 * num_classes)
  //     .setMinFractionalErrorReduction(0.00001)
  //     .setFeatureSamplingFraction(1)
  //     .setMaxThresholdsFraction(0.25)
  //     .setForceGreedy(true)
  //     .setVerbose(false);

  // Options for pendigits
  // opts.setMinBoostingRounds(min(num_classes, 3L))
  //     .setMaxBoostingRounds(4 * num_classes)
  //     .setMinFractionalErrorReduction(0.0000001)
  //     .setFeatureSamplingFraction(0.25)
  //     .setMaxThresholdsFraction(0.001)
  //     .setForceGreedy(true)
  //     .setVerbose(true);
#else
  opts.setMinBoostingRounds(10 * num_classes)
      .setMaxBoostingRounds(40 * num_classes)
      .setMinFractionalErrorReduction(0.00001)
      .setFeatureSamplingFraction(1)
      .setMaxThresholdsFraction(1)
      .setForceExhaustive(true)
      .setVerbose(true);
#endif

  JointBoost jb(num_classes, all_training->numFeatures(), opts);

  // Self-testing
  {
    cout << endl;
    cout << "======================================" << endl;
    cout << "Self-testing" << endl;
    cout << "======================================" << endl << endl;

    jb.train(*all_training);
    jb.dumpToConsole();

    if (!test(jb, *all_training, true))
      return false;
  }

  jb.clear();

  // Holdout-testing
  if (training_subset && holdout_subset)
  {
    cout << endl;
    cout << "======================================" << endl;
    cout << "Holdout-testing" << endl;
    cout << "======================================" << endl << endl;

    jb.train(*training_subset);
    jb.dumpToConsole();

    if (!test(jb, *holdout_subset, true))
      return false;
  }
  else
  {
    cerr << "Not enough examples to do holdout testing" << endl;
    return false;
  }

  return true;
}
