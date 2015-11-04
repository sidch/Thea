#include "../Algorithms/LinearScalarFunction.hpp"
#include "../Algorithms/NumericalOptimizer.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

class MyOpt : public NumericalOptimizer
{
  public:
    MyOpt() : NamedObject("MyOpt") {}

    bool minimize(ScalarFunction const & objective, double const * hint = NULL, Options const & options = Options())
    {
      return true;
    }
};

int
main(int argc, char * argv[])
{
  return 0;
}
