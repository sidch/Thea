#include "../Common.hpp"
#include "../Application.hpp"
#include "../FilePath.hpp"
#include "../IPlugin.hpp"
#include "../MatrixWrapper.hpp"
#include "../MatVec.hpp"
#include "../Options.hpp"
#include "../SparseMatrixWrapper.hpp"
#include "../SparseMatVec.hpp"
#include "../Algorithms/IEigenSolver.hpp"
#include <iostream>
#include <cstdio>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testARPACK(int argc, char * argv[]);
int cleanup(int status);

int
main(int argc, char * argv[])
{
  bool ok;
  try
  {
    ok = testARPACK(argc, argv);
  }
  THEA_CATCH(return cleanup(-1);, ERROR, "%s", "An error occurred")

  if (ok)  // hooray, all tests passed
    cout << "Test passed" << endl;
  else  // test failed
    cout << "Test failed" << endl;

  return cleanup(ok ? 0 : -1);
}

bool
testARPACK(int argc, char * argv[])
{
  // Get the path containing the executable
  string bin_path = FilePath::parent(argv[0]);

  // Try to load the ARPACK plugin from the same parent directory as the executable
#ifdef THEA_DEBUG_BUILD
  string plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginARPACKd");
#else
  string plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginARPACK");
#endif

  cout << "Loading plugin: " << plugin_path << endl;
  auto plugin = Application::getPluginManager().load(plugin_path);

  // Start up the plugin
  plugin->startup();

  // We should now have a ARPACK eigensolver factory
  auto factory = Application::getEigenSolverManager().getFactory("ARPACK");

  // Create an eigensolver
  auto eigs = factory->createEigenSolver("My ARPACK eigensolver");

  MatrixXd A(4, 4);
  A(0, 0) = 3; A(0, 1) = 0; A(0, 2) = 0; A(0, 3) = 1;
  A(1, 0) = 0; A(1, 1) = 4; A(1, 2) = 0; A(1, 3) = 0;
  A(2, 0) = 0; A(2, 1) = 0; A(2, 2) = 4; A(2, 3) = 0;
  A(3, 0) = 1; A(3, 1) = 0; A(3, 2) = 0; A(3, 3) = 3;
  auto A_wrapper = Math::wrapMatrix(A);

  typedef SparseColumnMatrix<double> CSC;
  CSC A_sparse = A.sparseView();
  A_sparse.makeCompressed();
  auto A_sparse_wrapper = Math::wrapMatrix(A_sparse);

  Options opts;
  opts.set("which", "SM");
  opts.setInteger("maxit", 1000);
  auto num_eigs = eigs->solve(&A_wrapper, /* compute_eigenvectors = */ false, /* num_requested_eigenpairs = */ 2, &opts);

  bool ok = true;
  if (num_eigs >= 0)
  {
    THEA_CONSOLE << "Found " << num_eigs << " eigenpair(s)";
  }
  else
  {
    THEA_CONSOLE << "System could not be solved";
    ok = false;
  }

  // Destroy the linear solver
  factory->destroyEigenSolver(eigs);

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
