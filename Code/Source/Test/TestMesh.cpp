#define TEST_GENERAL_MESH
#define TEST_DCEL_MESH
#define TEST_CONNECTED_COMPONENTS
#define TEST_MANIFOLD
// #define TEST_IMLS

#include "../Common.hpp"

#if defined(TEST_GENERAL_MESH) || defined(TEST_CONNECTED_COMPONENTS) || defined(TEST_MANIFOLD) || defined(TEST_IMLS)
#include "../Graphics/GeneralMesh.hpp"
typedef Thea::Graphics::GeneralMesh<> GM;
#endif

#if defined(TEST_DCEL_MESH)
#include "../Graphics/DCELMesh.hpp"
typedef Thea::Graphics::DCELMesh<> DM;
#endif

#include "../Graphics/MeshType.hpp"
#include "../Algorithms/ConnectedComponents.hpp"
#include "../Algorithms/IMLSSurface.hpp"
#include "../Algorithms/ImplicitSurfaceMesher.hpp"
#include "../Algorithms/MeshKDTree.hpp"
#include "../Application.hpp"
#include "../Array.hpp"
#include "../FilePath.hpp"
#include "../FileSystem.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

void testMesh(int argc, char * argv[]);
void testGeneralMesh(int argc, char * argv[]);
void testDCELMesh(int argc, char * argv[]);
void testManifold(int argc, char * argv[]);
void testIMLS(int argc, char * argv[]);

string data_dir;

int
main(int argc, char * argv[])
{
  try
  {
#ifdef _MSC_VER
    data_dir = FilePath::concat(FilePath::parent(FileSystem::resolve(Application::programPath())),
                                "../../../../../Data/Models");
#else
    data_dir = FilePath::concat(FilePath::parent(FileSystem::resolve(Application::programPath())),
                                "../../../../Data/Models");
#endif
    testMesh(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  // Hooray, all tests passed
  cout << "Mesh: Test completed" << endl;
  return 0;
}

void
testMesh(int argc, char * argv[])
{
  testGeneralMesh(argc, argv);
  testDCELMesh(argc, argv);
  testManifold(argc, argv);
  testIMLS(argc, argv);
}

void
testGeneralMesh(int argc, char * argv[])
{
#ifdef TEST_GENERAL_MESH

  MeshGroup<GM> mg("General Mesh Group");
  GM::Ptr mesh(new GM("General Mesh"));
  mg.addMesh(mesh);

  Algorithms::MeshKDTree<GM> kdtree;
  kdtree.add(mg);

#ifdef TEST_CONNECTED_COMPONENTS
  TheaArray< TheaArray<GM::Face *> > components;
  ConnectedComponents::findEdgeConnected(*mesh, components);
#endif

#endif
}

void
testDCELMesh(int argc, char * argv[])
{
#ifdef TEST_DCEL_MESH

  MeshGroup<DM> mg("DCEL Mesh Group");
  DM::Ptr mesh(new DM("DCEL Mesh"));
  mg.addMesh(mesh);

  Algorithms::MeshKDTree<DM> kdtree;
  kdtree.add(mg);

#ifdef TEST_CONNECTED_COMPONENTS
  TheaArray< TheaArray<DM::Face *> > components;
  ConnectedComponents::findEdgeConnected(*mesh, components);
#endif

#endif
}

#ifdef TEST_MANIFOLD
bool isManifold(GM const & mesh)
{
  bool is_closed_manifold = mesh.isManifold(true  /* require_closed */);
  bool is_open_manifold = true;
  if (!is_closed_manifold)
    is_open_manifold = mesh.isManifold(false  /* require_closed */);

  cout << "Mesh '" << mesh.getName()
       << "' (V: " << mesh.numVertices() << ", F: " << mesh.numFaces() << ", E: " << mesh.numEdges() << ") ";

  if (is_closed_manifold)
    cout << "is a closed manifold" << endl;
  else if (is_open_manifold)
    cout << "is an open manifold" << endl;
  else
    cout << "is not manifold" << endl;

  return false;
}
#endif

void
testManifold(int argc, char * argv[])
{
#ifdef TEST_MANIFOLD

  string model_path = (argc < 2 ? FilePath::concat(data_dir, "teapot.obj") : argv[1]);

  MeshGroup<GM> input_mg(FilePath::objectName(model_path));
  input_mg.load(model_path);

  input_mg.forEachMeshUntil(isManifold);

#endif
}

#ifdef TEST_IMLS
class BoundaryValuesFunctor
{
  private:
    IMLSSurface const & imls;

  public:
    BoundaryValuesFunctor(IMLSSurface const & imls_) : imls(imls_) {}

    bool operator()(GM const & mesh)
    {
      for (GM::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
        cout << "imls" << vi->getPosition().toString() << " = " << imls(vi->getPosition()) << endl;

      return false;
    }
};
#endif

inline double
log10(double x)
{
  static double const LOG_10 = log(10.0);
  return log(x) / LOG_10;
}

void
testIMLS(int argc, char * argv[])
{
#ifdef TEST_IMLS
  MeshGroup<GM> input_mg("IMLS Input");
  // string model_path = "../../Data/Models/m0.off";
  string model_path = FilePath::concat(data_dir, "teapot.obj");
  input_mg.load(model_path);
  input_mg.updateBounds();

  // Compute the function
  IMLSSurface imls(input_mg, 0, 0.001);
  cout << "IMLS function computed" << endl;

// #define EVAL_BOUNDARY
#ifdef EVAL_BOUNDARY
  BoundaryValuesFunctor bvf(imls);
  input_mg.forEachMeshUntil(&bvf);
#endif

#define VOXELIZE
#ifdef VOXELIZE

  AxisAlignedBox3 const & bounds = input_mg.getBounds().scaleCenteredCopy(1.2);
  Vector3 extent = bounds.getExtent();
  static int const NUM_STEPS = 30;
  double voxels[NUM_STEPS][NUM_STEPS][NUM_STEPS];
  Vector3 step = extent / NUM_STEPS;
  for (int i = 0; i < NUM_STEPS; ++i)
  {
    for (int j = 0; j < NUM_STEPS; ++j)
      for (int k = 0; k < NUM_STEPS; ++k)
      {
        Vector3 p = bounds.getLow() + Vector3((i + 0.5f) * step.x(), (j + 0.5f) * step.y(), (k + 0.5f) * step.z());
        voxels[i][j][k] = imls(p);
      }

    cout << "Processed slice " << i << endl;
  }

  // binvox format, see http://www.cs.princeton.edu/~min/binvox/binvox.html
  ofstream voxout((FilePath::baseName(model_path) + ".binvox").c_str(), ios::binary);
  voxout << "#binvox 1\n"
         << "dim " << NUM_STEPS << ' ' << NUM_STEPS << ' ' << NUM_STEPS << '\n'  // add a border
         << "translate 0 0 0\n"
         << "scale 1\n"
         << "data\n";

  double max_value = -1;
  for (int i = 0; i < NUM_STEPS; ++i)
    for (int j = 0; j < NUM_STEPS; ++j)
      for (int k = 0; k < NUM_STEPS; ++k)
        if (fabs(voxels[i][j][k]) > max_value)
          max_value = fabs(voxels[i][j][k]);

  static double const THRESHOLD = 0.03;
  for (int i = 0; i < NUM_STEPS; ++i)
  {
    cout << "Slice " << i << ':' << endl;

    for (int k = 0; k < NUM_STEPS; ++k)  // binvox voxel ordering has Z running slower than Y
    {
      for (int j = 0; j < NUM_STEPS; ++j)
      {
        voxels[i][j][k] /= max_value;
        // cout << "val = " << voxels[i][j][k] << endl;

        if (fabs(voxels[i][j][k]) < THRESHOLD)
        {
          int ll = (int)floor(-log10(fabs(voxels[i][j][k])));
          if (ll < 1)
            ll = 1;
          else if (ll > 9)
            ll = 9;

          cout << ' ' << ll;
          voxout << (unsigned char)1 << (unsigned char)1;
        }
        else
        {
          cout << " 0";
          voxout << (unsigned char)0 << (unsigned char)1;
        }
      }

      cout << endl;
    }
  }

  voxout.close();

#endif

// #define POLYGONIZE
#ifdef POLYGONIZE
  GM::Ptr output(new GM("IMLS Output"));

  ImplicitSurfaceMesher::Options mesh_opts = ImplicitSurfaceMesher::Options::defaults();

#  if 1
  // Get a vertex on the mesh
  Vector3 p = (*input_mg.meshesBegin())->verticesBegin()->getPosition();

  // mesh_opts.bloomenthal.cell_size = input_mg.getBounds().getExtent().max() / 15.0;
  // mesh_opts.bloomenthal.max_search_steps = 10;

  cout << "Cell size = " << mesh_opts.bloomenthal.cell_size
       << ", max search steps = " << mesh_opts.bloomenthal.max_search_steps << endl;

  ImplicitSurfaceMesher::meshBloomenthal(&imls, imls.getBoundingBallWithNegativeCenter(), p, mesh_opts.bloomenthal, *output);
#  else
  ImplicitSurfaceMesher::meshBoissonnatOudot(&imls, imls.getBoundingBallWithNegativeCenter(), mesh_opts.boissonnat_oudot,
                                             *output);
#  endif

  output->updateBounds();
  cout << "Output mesh has " << output->numVertices() << " vertices, " << output->numFaces() << " faces and bounding box "
       << output->getBounds().toString() << endl;

  MeshGroup<GM> output_mg("IMLS Output");
  output_mg.addMesh(output);
  output_mg.save((FilePath::baseName(model_path) + "_remeshed.obj").c_str(), CodecOBJ<GM>());

#endif

#endif
}
