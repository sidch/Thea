#include "../Algorithms/PCA_N.hpp"
#include "../Algorithms/SparsePCA_N.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../Common.hpp"
#include "../Math.hpp"
#include "../Vector3.hpp"
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
  static long const DATA_SIZE = 10;
  Vector3 data[DATA_SIZE];
  cout << "Data vectors:" << endl;
  for (long i = 0; i < DATA_SIZE; ++i)
  {
    data[i] = Vector3(rand() / (Real)RAND_MAX,
                      rand() / (Real)RAND_MAX,
                      rand() / (Real)RAND_MAX);

    cout << "  " << data[i].toString() << endl;
  }

  //==========================================================================================================================
  // PCA
  //==========================================================================================================================

  Vector3 axes[3], centroid;
  Real vars[3];
  PCA_N<Vector3, 3>::compute(&data[0], &data[DATA_SIZE], vars, axes, &centroid);

  cout << "\nPrincipal axes = " << endl;
  for (long i = 0; i < 3; ++i)
    cout << "  " << axes[i].toString() << " (variance = " << vars[i] << ')' << endl;
  cout << "Centroid = " << centroid.toString() << endl;

  for (long i = 0; i < 3; ++i)
  {
    long j = (i + 1) % 3;
    cout << i << " x " << j << " = " << axes[i].cross(axes[j]).toString() << endl;
  }

  //==========================================================================================================================
  // Sparse PCA
  //==========================================================================================================================

  SparsePCA_N<Vector3, 3>::compute(&data[0], &data[DATA_SIZE], vars, axes, &centroid);

  cout << "\nSparse principal axes = " << endl;
  for (long i = 0; i < 3; ++i)
    cout << "  " << axes[i].toString() << " (variance = " << vars[i] << ')' << endl;
  cout << "Centroid = " << centroid.toString() << endl;

  for (long i = 0; i < 3; ++i)
  {
    long j = (i + 1) % 3;
    cout << i << " x " << j << " = " << axes[i].cross(axes[j]).toString() << endl;
  }

  return true;
}
