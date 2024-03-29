#include "../Common.hpp"
#include "../Application.hpp"
#include "../Array.hpp"
#include "../FilePath.hpp"
#include "../MatVec.hpp"
#include "../Options.hpp"
#include "../SparseMatrixWrapper.hpp"
#include "../SparseMatVec.hpp"
#include "../IPlugin.hpp"
#include "../Algorithms/ILinearSolver.hpp"
#include <cmath>
#include <iostream>
#include <cstdio>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testCSPARSE(int argc, char * argv[]);
int cleanup(int status);

int
main(int argc, char * argv[])
{
  bool ok;
  try
  {
    ok = testCSPARSE(argc, argv);
  }
  THEA_CATCH(return cleanup(-1);, ERROR, "%s", "An error occurred")

  if (ok)  // hooray, all tests passed
    cout << "Test passed" << endl;
  else  // test failed
    cout << "Test failed" << endl;

  return cleanup(ok ? 0 : -1);
}

bool
testCSPARSE(int argc, char * argv[])
{
  // Get the path containing the executable
  string bin_path = FilePath::parent(argv[0]);

  // Try to load the CSPARSE plugin from the same parent directory as the executable
#ifdef THEA_DEBUG_BUILD
  string plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginCSPARSEd");
#else
  string plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginCSPARSE");
#endif

  cout << "Loading plugin: " << plugin_path << endl;
  IPlugin * plugin = Application::getPluginManager().load(plugin_path);

  // Start up the plugin
  plugin->startup();

  // We should now have a CSPARSE linear solver factory
  ILinearSolverFactory * factory = Application::getLinearSolverManager().getFactory("CSPARSE");

  // Create a linear solver
  ILinearSolver * ls = factory->createLinearSolver("My CSPARSE linear solver");

  MatrixXd A(3, 3);
  A(0, 0) =   1; A(0, 1) =   3; A(0, 2) =  -2;
  A(1, 0) =   3; A(1, 1) =   5; A(1, 2) =   6;
  A(2, 0) =   2; A(2, 1) =   4; A(2, 2) =   3;
  double b[3] = { 5, 7, 8 };

  typedef SparseColumnMatrix<double> CSC;
  CSC A_sparse = A.sparseView();
  A_sparse.makeCompressed();

  Options opts;
  opts.set("method", "LU");
  ls->solve(&asLvalue(Math::wrapMatrix(A_sparse)), b, /* guess = */ nullptr, &opts);

  static double const EPSILON = 0.0001;
  double expected[3] = {-15, 8, 2};

  bool ok = true;
  if (ls->hasSolution())
  {
    cout << "Solution:" << endl;

    double const * sols = ls->getSolution();
    for (intx i = 0; i < ls->dims(); ++i)
    {
      printf("  x[%ld] = %g (expected %g)\n", (intx)i, sols[i], expected[i]);

      double allowed_error = std::max(1e-10, EPSILON * expected[i]);
      if (std::fabs(sols[i] - expected[i]) > allowed_error)
        ok = false;
    }
  }
  else
  {
    cout << "System could not be solved" << endl;
    ok = false;
  }

  // Destroy the linear solver
  factory->destroyLinearSolver(ls);

  // Cleanup and quit
  plugin->shutdown();

  return ok;
}

int
cleanup(int status)
{
  Application::getPluginManager().unloadAllPlugins();
  return status;
}
