#include "../Algorithms/JointBoost.hpp"
#include "../Array.hpp"
#include "../Matrix.hpp"
#include "../VectorN.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testJointBoost();

int
main(int argc, char * argv[])
{
  try
  {
    if (!testJointBoost()) return -1;
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}

static int const NCLASSES = 2;
static int const NFEATURES = 3;
typedef VectorN<NFEATURES, double> DVec;

class TrainingData : public JointBoost::TrainingData
{
  public:
    TrainingData() : features(0, NFEATURES) {}

    void addExample(double const * example_features, long example_class)
    {
      features.appendRow();
      features.setRow(features.numRows() - 1, example_features);
      classes.push_back(example_class);
    }

    void addExample(DVec const & example_features, long example_class)
    {
      features.appendRow();
      features.setRow(features.numRows() - 1, &example_features[0]);
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

    DVec getExampleFeatures(long example)
    {
      DVec v; features.getRow(example, &v[0]);
      return v;
    }

    long getExampleClass(long example)
    {
      return classes[(array_size_t)example];
    }

  private:
    Matrix<double, MatrixLayout::ROW_MAJOR> features;
    TheaArray<long> classes;
};

bool
testJointBoost()
{
  THEA_CONSOLE << "======================================";
  THEA_CONSOLE << "Training";
  THEA_CONSOLE << "======================================";

  TrainingData td;
  td.addExample(DVec(0, 0, 0) * Math::randInRange(0.99, 1.01), 0);
  td.addExample(DVec(0, 0, 1) * Math::randInRange(0.99, 1.01), 1);
  td.addExample(DVec(0, 1, 0) * Math::randInRange(0.99, 1.01), 0);
  td.addExample(DVec(0, 1, 1) * Math::randInRange(0.99, 1.01), 1);
  td.addExample(DVec(1, 0, 0) * Math::randInRange(0.99, 1.01), 0);
  td.addExample(DVec(1, 0, 1) * Math::randInRange(0.99, 1.01), 1);
  td.addExample(DVec(1, 1, 0) * Math::randInRange(0.99, 1.01), 0);
  td.addExample(DVec(1, 1, 1) * Math::randInRange(0.99, 1.01), 1);

  JointBoost::Options opts;
  opts.setMinBoostingRounds(1)
      .setMaxBoostingRounds(3)
      .setMinFractionalErrorReduction(-1)
      .setFeatureSamplingFraction(1)
      .setForceExhaustive(true);

  JointBoost jb(NCLASSES, NFEATURES, opts);
  jb.train(td);

  jb.dumpToConsole();

  THEA_CONSOLE << "";
  THEA_CONSOLE << "======================================";
  THEA_CONSOLE << "Testing";
  THEA_CONSOLE << "======================================";

  TheaArray<double> class_probabilities((array_size_t)NCLASSES);
  for (long i = 0; i < td.numExamples(); ++i)
  {
    DVec f = td.getExampleFeatures(i);
    long c = td.getExampleClass(i);

    long pc = jb.predict(&f[0], &class_probabilities[0]);

    THEA_CONSOLE << "Features = " << f.toString() << ", actual class = " << c << ", predicted class = " << pc;
    for (array_size_t c = 0; c < class_probabilities.size(); ++c)
      THEA_CONSOLE << "    Probability of class " << c << " = " << class_probabilities[c];

    THEA_CONSOLE << "";
  }

  return true;
}
