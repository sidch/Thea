#include <Thea/Algorithms/BestFitSphereN.hpp>
#include <Thea/Algorithms/ImplicitSurfaceMesher.hpp>
#include <Thea/Algorithms/MeshKDTree.hpp>
#include <Thea/Algorithms/MetricL2.hpp>
#include <Thea/Algorithms/PointCollectorN.hpp>
#include <Thea/Graphics/GeneralMesh.hpp>
#include <Thea/Graphics/MeshGroup.hpp>
#include <Thea/Stopwatch.hpp>
#include <cstdlib>

// Namespaces and typedefs to make things easier to read below
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef GeneralMesh<> Mesh;

/** Assumes consistent outward normals, and no nested structures. */
class SignedDistance
{
  public:
    /** Constructor. */
    SignedDistance(MeshGroup<Mesh> & m) : num_calls(0)
    {
      kdtree.add(m);
      kdtree.init();
      kdtree.enableNearestNeighborAcceleration();
    }

    /** Evaluate the signed distance from a given point. */
    Real operator()(Vector3 const & p) const
    {
      num_calls++;

      double d; Vector3 cp;
      auto index = kdtree.closestElement<MetricL2>(p, -1, &d, &cp);
      Vector3 cn = kdtree.getElements()[index].getNormal();
      return cn.dot(cp - p) < 0 ? d : -d;
    }

    /** How many times was the kdtree queried? */
    intx numCalls() const { return num_calls; }

  private:
    MeshKDTree<Mesh> kdtree;
    mutable intx num_calls;
};

int
remesh(int argc, char * argv[])
{
  // Parse cmdline
  if (argc < 3)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " <mesh-in> <mesh-out> [<precision>]";
    return -1;
  }

  int precision = 20;
  if (argc > 3) precision = std::atoi(argv[3]);

  // Read the input mesh
  MeshGroup<Mesh> m("Input shape");
  m.load(argv[1]);

  // Setup the remeshing
  Mesh::Ptr remeshed(new Mesh);
  SignedDistance sd(m);

  // Compute the bounding sphere, for a remeshing range and a scale parameter
  BestFitSphereN<3> bounds;
  PointCollectorN< BestFitSphereN<3>, 3 > c(&bounds);
  c.addMeshVertices(m);

  // Do the remeshing (and time it)
  Stopwatch timer;
  timer.tick();

    ImplicitSurfaceMesher::Options::Bloomenthal opts(
      /* cell_size = */ bounds.getDiameter() / precision,
      /* max_search_steps = */ precision);

    ImplicitSurfaceMesher::meshBloomenthal(
      &sd,
      bounds.getBall(),
      (*m.meshesBegin())->verticesBegin()->getPosition(),
      opts,
      *remeshed);

  timer.tock();
  THEA_CONSOLE << sd.numCalls() << " kdtree NN queries in "
               << 1000000 * timer.elapsedTime() / sd.numCalls() << "ns each";

  // Save the output mesh
  MeshGroup<Mesh> out("Remeshed");
  out.addMesh(remeshed);
  out.save(argv[2]);

  return 0;
}

int
main(int argc, char * argv[])
{
  try
  {
    return remesh(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "Could not remesh the shape")
}
