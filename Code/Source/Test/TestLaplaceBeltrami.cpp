#include "../Common.hpp"

using namespace Thea;

#define MESH_TYPE GENERAL
// #define MESH_TYPE DCEL
// #define MESH_TYPE CGAL

#if MESH_TYPE == GENERAL
#  include "../Graphics/GeneralMesh.hpp"
   typedef Graphics::GeneralMesh<> Mesh;

#elif MESH_TYPE == DCEL
#  include "../Graphics/DCELMesh.hpp"
   typedef Graphics::DCELMesh<> Mesh;

#else
#  include "../Graphics/CGALMesh.hpp"
   typedef Graphics::CGALMesh<> Mesh;
#endif

#include "../Algorithms/LaplaceBeltrami.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../MappedMatrix.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace Algorithms;
using namespace Graphics;

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

void
testLB(int argc, char * argv[])
{
  if (argc < 2)
  {
    cerr << "Usage: " << argv[0] << " <mesh>" << endl;
    return;
  }

  string model_path = argv[1];
  MeshGroup<Mesh> mesh_group("Manifold");

  CodecOBJ<Mesh> const codec_obj(CodecOBJ<Mesh>::ReadOptions().setIgnoreTexCoords(true).setFlatten(true));
  Codec3DS<Mesh> const codec_3ds(Codec3DS<Mesh>::ReadOptions().setIgnoreTexCoords(true).setFlatten(true));

  if (endsWith(toLower(model_path), "obj"))
    mesh_group.load(model_path, codec_obj);
  else if (endsWith(toLower(model_path), "3ds"))
    mesh_group.load(model_path, codec_3ds);
  else
    mesh_group.load(model_path);

  if (mesh_group.numMeshes() <= 0)
    throw Error("Mesh group is empty");

  Mesh::Ptr mesh = *mesh_group.meshesBegin();

  typedef MappedMatrix<Real> Mat;
  Mat lb;
  LaplaceBeltrami::compute(*mesh, LaplaceBeltrami::Method::XU_2006, lb);

  cout << "[\n";

  for (Mat::ConstIterator mi = lb.begin(); mi != lb.end(); ++mi)
    cout << "    (" << mi->first.first << ", " << mi->first.second << ") = " << mi->second << '\n';

  cout << ']' << endl;
}
