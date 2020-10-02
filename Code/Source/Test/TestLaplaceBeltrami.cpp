#include "../Common.hpp"

using namespace Thea;

#define MESH_TYPE GENERAL
// #define MESH_TYPE DCEL

#if MESH_TYPE == GENERAL
#  include "../Graphics/GeneralMesh.hpp"
   typedef Graphics::GeneralMesh<> Mesh;

#else // MESH_TYPE == DCEL
#  include "../Graphics/DcelMesh.hpp"
   typedef Graphics::DcelMesh<> Mesh;
#endif

#include "../Algorithms/LaplaceBeltrami.hpp"
#include "../Algorithms/IEigenSolver.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Application.hpp"
#include "../FilePath.hpp"
#include "../MappedMatrix.hpp"
#include "../MatVec.hpp"
#include "../Options.hpp"
#include "../IPlugin.hpp"
#include "../SparseMatrixWrapper.hpp"
#include "../SparseMatVec.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace Algorithms;
using namespace Graphics;

#ifdef THEA_DEBUG_BUILD
  static std::string const arpack_plugin = "libTheaPluginARPACKd";
#else
  static std::string const arpack_plugin = "libTheaPluginARPACK";
#endif

void testLB(int argc, char * argv[]);

int
main(int argc, char * argv[])
{
  // Do the OpenGL tests
  try
  {
    testLB(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}

Mesh::Ptr
loadMesh(std::string const & path)
{
  Array< MeshCodec<Mesh>::Ptr > codecs {
    std::make_shared< CodecObj<Mesh> >(CodecObj<Mesh>::ReadOptions().setIgnoreTexCoords(true).setFlatten(true)),
    std::make_shared< Codec3ds<Mesh> >(Codec3ds<Mesh>::ReadOptions().setIgnoreTexCoords(true).setFlatten(true))
  };

  MeshGroup<Mesh> mesh_group;
  mesh_group.load(path, codecs);

  if (mesh_group.numMeshes() <= 0)
    throw Error("Mesh group is empty");

  Mesh::Ptr mesh = *mesh_group.meshesBegin();
  mesh->setName(FilePath::baseName(path));
  return mesh;
}

void
testLB(int argc, char * argv[])
{
  if (argc < 3)
  {
    cerr << "Usage: " << argv[0] << " <mesh> <num_eigs>" << endl;
    return;
  }

  string model_path = argv[1];
  intx num_eigs = std::atoi(argv[2]);

  Mesh::Ptr mesh = loadMesh(model_path);
  intx n = mesh->numVertices();

  MappedMatrix<float64> lb;
  LaplaceBeltrami::compute(*mesh, LaplaceBeltrami::Method::XU_2006, lb);

  SparseColumnMatrix<float64> sparse_lb(n, n);
  sparse_lb.setFromTriplets(lb.tripletsBegin(), lb.tripletsEnd());

  IPlugin * plugin = Application::getPluginManager().load(arpack_plugin);
  plugin->startup();
  Algorithms::IEigenSolverFactory * factory = Application::getEigenSolverManager().getFactory("ARPACK");
  Algorithms::IEigenSolver * eig = factory->createEigenSolver("My ARPACK eigensolver");

  intx num_ret = eig->solve(&asLvalue(Math::wrapMatrix(sparse_lb)),
                            /* compute_eigenvectors = */ true,
                            /* num_requested_eigenpairs = */ num_eigs);

  float64 eigval_re, eigval_im;
  float64 const * eigvec_re, * eigvec_im;
  for (intx i = 0; i < num_ret; ++i)
  {
    eig->getEigenvalue(i, &eigval_re, &eigval_im);
    eig->getEigenvector(i, &eigvec_re, &eigvec_im);

    THEA_CONSOLE << "eig[" << i << "] = ((" << eigval_re << ", " << eigval_im << "), "
                 << "(re: " << toString(VectorXdConstMap(eigvec_re, n))
                 << ", im: " << toString(VectorXdConstMap(eigvec_im, n)) << "))";
  }
}
