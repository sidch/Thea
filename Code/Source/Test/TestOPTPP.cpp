#include "../Algorithms/ScalarFunction.hpp"
#include "../Algorithms/NumericalOptimizer.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

class MyOpt : public NumericalOptimizer, public virtual NamedObject
{
  public:
    MyOpt() : NamedObject("MyOpt") {}

    int8 minimize(ScalarFunction const & objective, float64 const * hint = NULL, AbstractOptions const * options = NULL)
    {
      return true;
    }
};

int
main(int argc, char * argv[])
{
  return 0;
}
