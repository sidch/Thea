#include "../Algorithms/BagOfWords.hpp"
#include "../Matrix.hpp"
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
  static long const NUM_POINT_FEATURES = 100;
  static long const NUM_WORDS = 5;
  static long const NUM_TRAINING_POINTS = 100;
  static long const NUM_TEST_POINTS = 5000;

  Random rng(0x01234567);

  Matrix<double> training_points(NUM_TRAINING_POINTS, NUM_POINT_FEATURES);
  for (long i = 0; i < NUM_TRAINING_POINTS; ++i)
    for (long j = 0; j < NUM_POINT_FEATURES; ++j)
      training_points(i, j) = rng.uniform01();

  BagOfWords bow;
  bow.train(NUM_WORDS, training_points);

  Matrix<double> test_points(NUM_TEST_POINTS, NUM_POINT_FEATURES);
  for (long i = 0; i < NUM_TEST_POINTS; ++i)
    for (long j = 0; j < NUM_POINT_FEATURES; ++j)
      test_points(i, j) = rng.uniform01();

  double histogram[NUM_WORDS];
  bow.computeWordFrequencies(test_points, &histogram[0]);

  THEA_CONSOLE << "Word frequencies = [";

  for (long i = 0; i < NUM_WORDS; ++i)
    THEA_CONSOLE << "      " << histogram[i];

  THEA_CONSOLE << "  ]";

  return true;
}
