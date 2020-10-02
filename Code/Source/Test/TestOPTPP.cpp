#include "../Algorithms/IScalarFunction.hpp"
#include "../Algorithms/INumericalOptimizer.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

class MyOpt : public INumericalOptimizer, public virtual NamedObject
{
  public:
    MyOpt() : NamedObject("MyOpt") {}

    int8 THEA_ICALL
    minimize(IScalarFunction const & objective, float64 const * hint = nullptr, IOptions const * options = nullptr)
    {
      return true;
    }
};

int
main(int argc, char * argv[])
{
  return 0;
}
