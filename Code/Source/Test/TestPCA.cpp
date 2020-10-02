#include "../Algorithms/PcaN.hpp"
#include "../Algorithms/SparsePcaN.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../Common.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testPCA();

int
main(int argc, char * argv[])
{
  if (!testPCA()) return -1;

  return 0;
}

bool
testPCA()
{
  static intx const DATA_SIZE = 10;
  Vector3 data[DATA_SIZE];
  cout << "Data vectors:" << endl;
  for (intx i = 0; i < DATA_SIZE; ++i)
  {
    data[i] = Vector3(rand() / (Real)RAND_MAX,
                      rand() / (Real)RAND_MAX,
                      rand() / (Real)RAND_MAX);

    cout << "  " << data[i] << endl;
  }

  //==========================================================================================================================
  // PCA
  //==========================================================================================================================

  Vector3 axes[3], centroid;
  Real vars[3];
  PcaN<Vector3, 3>::compute(&data[0], &data[DATA_SIZE], vars, axes, &centroid);

  cout << "\nPrincipal axes = " << endl;
  for (intx i = 0; i < 3; ++i)
    cout << "  " << axes[i] << " (variance = " << vars[i] << ')' << endl;
  cout << "Centroid = " << centroid << endl;

  for (intx i = 0; i < 3; ++i)
  {
    intx j = (i + 1) % 3;
    cout << i << " x " << j << " = " << toString(axes[i].cross(axes[j])) << endl;
  }

  //==========================================================================================================================
  // Sparse PCA
  //==========================================================================================================================

  SparsePcaN<Vector3, 3>::compute(&data[0], &data[DATA_SIZE], vars, axes, &centroid);

  cout << "\nSparse principal axes = " << endl;
  for (intx i = 0; i < 3; ++i)
    cout << "  " << axes[i] << " (variance = " << vars[i] << ')' << endl;
  cout << "Centroid = " << centroid << endl;

  for (intx i = 0; i < 3; ++i)
  {
    intx j = (i + 1) % 3;
    cout << i << " x " << j << " = " << toString(axes[i].cross(axes[j])) << endl;
  }

  return true;
}
