#include "../Algorithms/StdLinearSolver.hpp"
#include "../Algorithms/LogisticRegression.hpp"
#include "../Array.hpp"
#include "../BinaryInputStream.hpp"
#include "../MatrixWrapper.hpp"
#include "../MatVec.hpp"
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
printMatrix(MatrixX<float> const & m)
{
  for (intx r = 0; r < m.rows(); ++r)
  {
    cout << " [";

    for (intx c = 0; c < m.cols(); ++c)
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

  MatrixXf a(5, 4); a.setZero();
  a(0, 0) = 1;
  a(1, 3) = 4;
  a(2, 2) = 3;
  a(4, 0) = 2;

  // Make sure things work with more columns than rows
  a = a.transpose();
  Eigen::JacobiSVD<MatrixXf> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
  MatrixXf u = svd.matrixU();
  MatrixXf v = svd.matrixV();
  VectorXf d = svd.singularValues();

  cout.precision(3);

  cout << "U: " << u.rows() << " x " << u.cols() << endl;
  cout << "D: " << d.size() << " elements" << endl;
  cout << "V: " << v.rows() << " x " << v.cols() << endl;

  // Print A
  cout << "\nA =" << endl;
  printMatrix(a);

  // Print U
  cout << "\nU =" << endl;
  printMatrix(u);

  // Print D
  cout << "\nD =" << endl;
  printMatrix(d);

  // Print V
  cout << "\nV =" << endl;
  printMatrix(v);

  // Print U * D * V
  cout << "\nU * D * V^T =" << endl;
  printMatrix(u * MatrixXf(d.asDiagonal()) * v.transpose());

  // Compute the pseudo-inverse of A via SVD
  MatrixXf svd_inv = v.transpose() * d.asDiagonal().inverse() * u;
  cout << "\nSVD pseudo-inverse =" << endl;
  printMatrix(svd_inv);

  // Compute the pseudo-inverse of A directly
  // cout << "\nDirect pseudo-inverse =" << endl;
  // printMatrix((a.transpose() * a).inverse() * a.transpose());

  //====================================================================
  // Linear least squares
  //====================================================================

  // First, unconstrained (expected result [-3.5, 1.4])
  StdLinearSolver lsq;
  MatrixXd coeffs; coeffs << -1, 1,
                             -1, 2,
                             -1, 3,
                             -1, 4;
  VectorXd consts; consts << 6, 5, 7, 10;

  cout << endl;
  if (lsq.solve(&asLvalue(Math::wrapMatrix(coeffs)), consts.data()))
    cout << "Unconstrained linear least squares solution: [" << lsq.getSolution()[0] << ", " << lsq.getSolution()[1] << "]"
         << endl;
  else
    cout << "Unconstrained linear least squares problem has no solution" << endl;

  // Next, non-negative solution only (expected result [0, 2.57] (CHECK!))
  cout << endl;
  StdLinearSolver lsq_nneg(StdLinearSolver::Method::DEFAULT, StdLinearSolver::Constraint::NON_NEGATIVE);
  if (lsq.solve(&asLvalue(Math::wrapMatrix(coeffs)), consts.data()))
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
