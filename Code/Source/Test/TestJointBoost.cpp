#include "../Algorithms/JointBoost.hpp"
#include "../Array.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../Matrix.hpp"
#include "../UnorderedMap.hpp"
#include "../VectorN.hpp"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testJointBoost();
bool testJointBoostFile(string const & path);

int
main(int argc, char * argv[])
{
  try
  {
    if (!testJointBoost()) return -1;

    THEA_CONSOLE << ""; THEA_CONSOLE << "";
    if (argc > 1)
    {
      if (!testJointBoostFile(argv[1])) return -1;
    }
    else
      THEA_CONSOLE << "Not testing JointBoost with data from file";
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}

class TrainingData : public JointBoost::TrainingData
{
  public:
    THEA_DEF_POINTER_TYPES(TrainingData, shared_ptr, weak_ptr)

    TrainingData(long nfeatures) : features(0, nfeatures) {}

    template <typename DVecT>
    void addExample(DVecT const & example_features, long example_class)
    {
      addExample(&example_features[0], example_class);
    }

    void addExample(double const * example_features, long example_class)
    {
      features.appendRow();
      features.setRow(features.numRows() - 1, example_features);
      classes.push_back(example_class);
    }

    long numExamples() const { return features.numRows(); }
    long numFeatures() const { return features.numColumns(); }

    void getFeature(long feature_index, TheaArray<double> & values) const
    {
      values.resize((array_size_t)numExamples());
      features.getColumn(feature_index, &values[0]);
    }

    void getClasses(TheaArray<long> & classes_) const
    {
      classes_ = classes;
    }

    void getExampleFeatures(long example, double * f) const
    {
      features.getRow(example, f);
    }

    long getExampleClass(long example) const
    {
      return classes[(array_size_t)example];
    }

  private:
    Matrix<double, MatrixLayout::ROW_MAJOR> features;
    TheaArray<long> classes;
};

template <typename ArrayT>
string
arrayToString(ArrayT const & arr)
{
  ostringstream oss; oss << '[';
  for (long i = 0; i < (long)arr.size(); ++i)
  {
    if (i > 0) oss << ", ";
    oss << arr[i];
  }
  oss << ']';

  return oss.str();
}

bool
selfTest(JointBoost const & jb, TrainingData const & td, bool get_class_probs = false)
{
  THEA_CONSOLE << "";
  THEA_CONSOLE << "======================================";
  THEA_CONSOLE << "Self-testing";
  THEA_CONSOLE << "======================================";

  TheaArray<double> class_probabilities((array_size_t)jb.numClasses());
  TheaArray<double> f((array_size_t)jb.numFeatures());
  long num_correct = 0;
  for (long i = 0; i < td.numExamples(); ++i)
  {
    td.getExampleFeatures(i, &f[0]);
    long c = td.getExampleClass(i);

    long pc = jb.predict(&f[0], get_class_probs ? &class_probabilities[0] : NULL);
    if (pc == c)
      num_correct++;

    THEA_CONSOLE << "Features[" << i << "] = " << arrayToString(f) << ", actual class = " << c << ", predicted class = " << pc;
    if (get_class_probs)
    {
      for (array_size_t c = 0; c < class_probabilities.size(); ++c)
        THEA_CONSOLE << "    Probability of class " << c << " = " << class_probabilities[c];

      THEA_CONSOLE << "";
    }
  }

  THEA_CONSOLE << "Accuracy = " << 100.0 * num_correct / (float)td.numExamples() << '%';

  return true;
}

bool
testJointBoost()
{
  THEA_CONSOLE << "======================================";
  THEA_CONSOLE << "Training";
  THEA_CONSOLE << "======================================";

  static int const NFEATURES = 2;
  static long const PTS_PER_RANGE = 10;

  typedef VectorN<NFEATURES, double> Vec;
  typedef AxisAlignedBoxN<NFEATURES, double> AABB;

  AABB ranges[] = { AABB(Vec(0, 0), Vec(1, 1)),
                    AABB(Vec(1, 1), Vec(4, 2)),
                    AABB(Vec(2, 0), Vec(3, 1))
                  };

  TrainingData td(NFEATURES);

  long num_classes = (long)(sizeof(ranges) / sizeof(ranges[0]));

  for (long r = 0; r < num_classes; ++r)
  {
    Vec ext = ranges[r].getExtent();
    for (long i = 0; i < PTS_PER_RANGE; ++i)
    {
      Vec p = ranges[r].getLow() + Vec(Math::rand01() * ext.x(), Math::rand01() * ext.y());
      td.addExample(p, r);
    }
  }

  JointBoost::Options opts;
  opts.setMinBoostingRounds(1)
      .setMaxBoostingRounds(3 * num_classes)
      .setMinFractionalErrorReduction(-1)
      .setFeatureSamplingFraction(1)
      .setForceExhaustive(true);

  JointBoost jb(num_classes, NFEATURES, opts);
  jb.train(td);
  jb.dumpToConsole();

  return selfTest(jb, td, true);
}

bool
testJointBoostFile(string const & path)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_WARNING << "Couldn't open file " << path;
    return false;
  }

  THEA_CONSOLE << "======================================";
  THEA_CONSOLE << "Training from file";
  THEA_CONSOLE << "======================================";

  TrainingData::Ptr td;

  typedef TheaUnorderedMap<string, long> LabelIndexMap;
  LabelIndexMap labels;

  TheaArray<double> features;
  double feature;
  string label;

  string line;
  while (getline(in, line))
  {
    boost::trim(line);
    if (line.empty())
      continue;

    vector<string> fields;
    boost::split(fields, line, boost::is_any_of(",\t"), boost::token_compress_on);

    if (fields.size() < 2)
    {
      THEA_ERROR << "Data has too few features per example";
      return false;
    }

    long nfeat = (long)fields.size() - 1;
    if (!td)
    {
      THEA_CONSOLE << "Data has " << nfeat << " features per example";
      td = TrainingData::Ptr(new TrainingData(nfeat));
    }
    else
    {
      if (nfeat != td->numFeatures())
      {
        THEA_ERROR << "Inconsistent number of features for example " << td->numExamples();
        return false;
      }
    }

    features.clear();
    for (long i = 0; i < nfeat; ++i)
    {
      istringstream iss(fields[(array_size_t)i]);
      iss >> feature;
      features.push_back(feature);
    }

    label = fields.back();

    LabelIndexMap::const_iterator existing_label = labels.find(label);
    long index;
    if (existing_label == labels.end())
    {
      index = (long)labels.size();
      labels[label] = index;
      THEA_CONSOLE << "Added class with label '" << label << "' and index " << index;
    }
    else
      index = existing_label->second;

    td->addExample(features, index);
  }

  if (!td)
  {
    THEA_ERROR << "Could not read any lines from file";
    return false;
  }

  long num_classes = (long)labels.size();

  THEA_CONSOLE << "Read " << td->numExamples() << " examples from " << num_classes << " classes from file";

  JointBoost::Options opts;
  opts.setMinBoostingRounds(num_classes)
      .setMaxBoostingRounds(4 * num_classes)
      .setMinFractionalErrorReduction(-1)
      .setFeatureSamplingFraction(3.0 / td->numFeatures())
      .setForceGreedy(true);

  JointBoost jb(num_classes, td->numFeatures(), opts);
  jb.train(*td);
  jb.dumpToConsole();

  return selfTest(jb, *td, true);
}
