#include "../Common.hpp"
#include "../Algorithms/PyramidMatch.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

int
main(int argc, char * argv[])
{
  try
  {
    PyramidMatch::test();
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  // Hooray, all tests passed
  cout << "PyramidMatch: Test completed" << endl;
  return 0;
}
