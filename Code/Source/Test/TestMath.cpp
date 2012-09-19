#include "../Algorithms/SVD.hpp"
#include "../Algorithms/LinearLeastSquares.hpp"
#include "../Algorithms/LogisticRegression.hpp"
#include "../Array.hpp"
#include "../Matrix.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testSVD();
bool testLogisticRegression();

int
main(int argc, char * argv[])
{
  try
  {
    if (!testSVD()) return -1;
    if (!testLogisticRegression()) return -1;
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}

void
printMatrix(Matrix<float> const & m)
{
  for (long r = 0; r < m.numRows(); ++r)
  {
    cout << " [";

    for (long c = 0; c < m.numColumns(); ++c)
      cout << '\t' << m(r, c);

    cout << "\t]" << endl;
  }
}

bool
testSVD()
{
  //====================================================================
  // SVD
  //====================================================================

  Matrix<float> a(5, 4, 0.0f), u, v;
  TheaArray<float> d;

  a(0, 0) = 1;
  a(1, 3) = 4;
  a(2, 2) = 3;
  a(4, 0) = 2;

  // Make sure things work with more columns than rows
  a = a.transpose();

  if (!SVD::compute(a, u, d, v))
  {
    THEA_ERROR << "Could not compute SVD";
    return false;
  }

  cout.precision(3);

  cout << "U: " << u.numRows() << " x " << u.numColumns() << endl;
  cout << "D: " << d.size() << " elements" << endl;
  cout << "V: " << v.numRows() << " x " << v.numColumns() << endl;

  // Print A
  cout << "\nA =" << endl;
  printMatrix(a);

  // Print U
  cout << "\nU =" << endl;
  printMatrix(u);

  // Print D
  cout << "\nD = [\t";
  for (array_size_t i = 0; i < d.size(); ++i)
  {
    if (i > 0) cout << '\t';
    cout << d[i];
  }
  cout << " \t]" << endl;

  // Print V
  cout << "\nV =" << endl;
  printMatrix(v);

  // Print U * D * V
  cout << "\nU * D * V =" << endl;
  printMatrix(u * Matrix<float>::fromDiagonal(d) * v);

  // Compute the pseudo-inverse of A via SVD
  Matrix<float> svd_inv;
  SVD::pseudoInverse(a, svd_inv);
  cout << "\nSVD pseudo-inverse =" << endl;
  printMatrix(svd_inv);

  // Compute the pseudo-inverse of A directly
  // cout << "\nDirect pseudo-inverse =" << endl;
  // printMatrix((a.transpose() * a).inverse() * a.transpose());

  //====================================================================
  // Linear least squares
  //====================================================================

  // First, unconstrained (expected result [-3.5, 1.4])
  LinearLeastSquares lsq(2);
  double coeffs[2] = { -1, 0 };
  coeffs[1] = 1; lsq.addObjective(coeffs,  6);
  coeffs[1] = 2; lsq.addObjective(coeffs,  5);
  coeffs[1] = 3; lsq.addObjective(coeffs,  7);
  coeffs[1] = 4; lsq.addObjective(coeffs, 10);

  cout << endl;
  if (lsq.solve())
    cout << "Unconstrained linear least squares solution: [" << lsq.getSolution()[0] << ", " << lsq.getSolution()[1] << "]"
         << endl;
  else
    cout << "Unconstrained linear least squares problem has no solution" << endl;

  // Next, non-negative solution only (expected result [0, 2.57] (CHECK!))
  cout << endl;
  if (lsq.solve(LinearLeastSquares::Constraint::NON_NEGATIVE))
    cout << "Non-negative linear least squares solution: [" << lsq.getSolution()[0] << ", " << lsq.getSolution()[1] << "]"
         << endl;
  else
    cout << "Non-negative linear least squares problem has no solution" << endl;

  return true;
}

bool
testLogisticRegression()
{
  LogisticRegression logreg(1);
  double x, y;

  x = 1900; y = 0.075995; logreg.addObservation(&x, y);
  x = 1910; y = 0.091972; logreg.addObservation(&x, y);
  x = 1920; y = 0.105711; logreg.addObservation(&x, y);
  x = 1930; y = 0.122755; logreg.addObservation(&x, y);
  x = 1940; y = 0.131699; logreg.addObservation(&x, y);
  x = 1950; y = 0.150697; logreg.addObservation(&x, y);
  x = 1960; y = 0.179323; logreg.addObservation(&x, y);
  x = 1970; y = 0.203212; logreg.addObservation(&x, y);
  x = 1980; y = 0.226505; logreg.addObservation(&x, y);

  cout << endl;

  if (logreg.solve())
  {
    double a = logreg.getSolution()[0];
    double b = logreg.getSolution()[1];

    double r = b;
    double t_0 = 1900;
    double c = std::exp(a + r * t_0);

    cout << "Logistic regression solution: a = " << a << ", b = " << b << endl;
    cout << "             => equivalently: c = " << c << ", r = " << r << endl;
  }
  else
    cout << "Logistic regression has no solution" << endl;

  return true;
}
