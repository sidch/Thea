#include "../Algorithms/BagOfWords.hpp"
#include "../MatrixWrapper.hpp"
#include "../MatVec.hpp"
#include "../Random.hpp"
#include "../EnumClass.hpp"

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testBagOfWords();

int
main(int argc, char * argv[])
{
  if (!testBagOfWords()) return -1;

  return 0;
}

bool
testBagOfWords()
{
  static intx const NUM_POINT_FEATURES = 100;
  static intx const NUM_WORDS = 5;
  static intx const NUM_TRAINING_POINTS = 100;
  static intx const NUM_TEST_POINTS = 5000;

  Random rng(0x01234567);

  MatrixXd training_points(NUM_TRAINING_POINTS, NUM_POINT_FEATURES);
  for (intx i = 0; i < NUM_TRAINING_POINTS; ++i)
    for (intx j = 0; j < NUM_POINT_FEATURES; ++j)
      training_points(i, j) = rng.uniform01();

  BagOfWords bow;
  bow.train(NUM_WORDS, MatrixWrapper<MatrixXd>(&training_points));

  MatrixXd test_points(NUM_TEST_POINTS, NUM_POINT_FEATURES);
  for (intx i = 0; i < NUM_TEST_POINTS; ++i)
    for (intx j = 0; j < NUM_POINT_FEATURES; ++j)
      test_points(i, j) = rng.uniform01();

  double histogram[NUM_WORDS];
  bow.computeWordFrequencies(MatrixWrapper<MatrixXd>(&test_points), &histogram[0]);

  THEA_CONSOLE << "Word frequencies = [";

  for (intx i = 0; i < NUM_WORDS; ++i)
    THEA_CONSOLE << "      " << histogram[i];

  THEA_CONSOLE << "  ]";

  return true;
}
