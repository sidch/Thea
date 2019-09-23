#include "../Algorithms/IteratorModifiers.hpp"
#include <vector>

using namespace std;
using namespace Thea;
using namespace Algorithms;

int
main(int argc, char * argv[])
{
  int a = 1, b = 2, c = 3;
  vector<int> ivec; ivec.push_back(a);
  vector<int *> ipvec; ipvec.push_back(&b);
  vector<int const *> ipcvec; ipcvec.push_back(&c);

  auto i1 = makeRefIterator(ivec.begin());
  auto i2 = makeRefIterator(ipvec.begin());
  auto i3 = makeRefIterator(ipcvec.begin());

  auto i4 = makePtrIterator(ivec.begin());
  auto i5 = makePtrIterator(ipvec.begin());
  auto i6 = makePtrIterator(ipcvec.begin());

  auto i7 = makeRefIterator(&ivec[0]);
  auto i8 = makeRefIterator(&ipvec[0]);
  auto i9 = makeRefIterator((int const * const *)&ipcvec[0]);

  auto i10 = makePtrIterator(&ivec[0]);
  auto i11 = makePtrIterator(&ipvec[0]);
  auto i12 = makePtrIterator((int const * const *)&ipcvec[0]);

  // All should print "1 2 3"
  THEA_CONSOLE << *i1 << " " << *i2 << " " << *i3;
  THEA_CONSOLE << **i4 << " " << **i5 << " " << **i6;
  THEA_CONSOLE << *i7 << " " << *i8 << " " << *i9;
  THEA_CONSOLE << **i10 << " " << **i11 << " " << **i12;

  return 0;
}
