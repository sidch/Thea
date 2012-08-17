#include "../../Common.hpp"
#include "../../Algorithms/KDTree3.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Graphics/VertexWelder.hpp"
#include "../../LineSegment3.hpp"
#include "../../Math.hpp"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Graphics;

int meshFix(int argc, char * argv[]);

int
main(int argc, char * argv[])
{
  int status = 0;
  try
  {
    status = meshFix(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return status;
}

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;

namespace Thea {
namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <>
class IsPointN<Mesh::Vertex, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <>
struct PointTraitsN<Mesh::Vertex, 3>
{
  static Vector3 const & getPosition(Mesh::Vertex const & t) { return t.getPosition(); }
};

} // namespace Algorithms
} // namespace Thea

bool verbose = false;
string infile, outfile;

bool do_del_danglers = false;

bool do_v_weld = false;
bool v_weld_boundary_only = true;
double v_weld_tolerance = -1;

bool do_zipper = false;
double zipper_tolerance = -1;

bool do_t_juncts = false;
double t_juncts_tolerance = -1;
int t_juncts_iters = -1;

int parseArgs(int argc, char * argv[]);
bool delDanglers(Mesh & mesh);
bool zipper(Mesh & mesh);
bool vWeld(Mesh & mesh);
bool tJuncts(Mesh & mesh);

int
meshFix(int argc, char * argv[])
{
  int parse_status = parseArgs(argc, argv);
  if (parse_status <= 0)
    return parse_status;

  MG mg(G3D::FilePath::baseExt(infile));
  mg.load(infile);
  mg.updateBounds();  // load() should do this, but let's play safe
  if (mg.isEmpty())
  {
    cerr << "Input mesh is empty, no output written" << endl;
    return -1;
  }

  if (do_del_danglers)
    mg.forEachMeshUntil(&delDanglers);

  if (do_t_juncts)
    mg.forEachMeshUntil(&tJuncts);

  if (do_v_weld)
    mg.forEachMeshUntil(&vWeld);

  if (do_zipper)
    mg.forEachMeshUntil(&zipper);

  string lc_out = toLower(outfile);
  if (endsWith(lc_out, ".3ds"))
    mg.save(outfile, Codec3DS<Mesh>());
  else if (endsWith(lc_out, ".obj"))
    mg.save(outfile, CodecOBJ<Mesh>());
  else if (endsWith(lc_out, ".off"))
    mg.save(outfile, CodecOFF<Mesh>());
  else
  {
    cerr << "Unrecognized mesh output format" << endl;
    return -1;
  }

  return 0;
}

int
parseArgs(int argc, char * argv[])
{
  namespace po = boost::program_options;

  static std::string const usage("\nUsage: " + string(argv[0]) + " [options] <infile> <outfile>\n"
                                   "   (all tolerances are specified as a fraction of the sub-mesh bounding box diagonal)\n");

  po::options_description visible("Allowed options");
  visible.add_options()
          ("help,h",              "Print this help message")
          ("version,v",           "Print the program version")
          ("verbose",             "Print extra status information")
          ("infile",              po::value<string>(&infile), "Path to input mesh file")
          ("outfile",             po::value<string>(&outfile), "Path to output mesh file")

          ("del-danglers",          "Delete isolated vertices and edges, and empty faces")

          ("zipper",              po::value<double>(&zipper_tolerance),
                                  "Seal pairs of boundary edges whose endpoints are equal upto the specified tolerance (implies --del-danglers)")

          ("v-weld",              po::value<double>(&v_weld_tolerance),
                                  "Merge vertices closer than the specified tolerance (implies --del-danglers)")

          ("v-weld-boundary",     po::value<double>(&v_weld_tolerance),
                                  "Merge boundary vertices closer than the specified tolerance (implies --del-danglers)")

          ("t-juncts",            po::value<double>(&t_juncts_tolerance),
                                  "Split boundary edges at t-junctions so they can be zippered properly (always executed before any welding or zippering)")

          ("tj-iters",            po::value<int>(&t_juncts_iters),
                                  "Number of iterations for fixing t-junctions")
  ;

  po::options_description desc;
  desc.add(visible) /* .add(hidden) */ ;

  po::positional_options_description pdesc;
  pdesc.add("infile", 1);
  pdesc.add("outfile", 1);

  // Read cmdline options first (overrides conflicting config file values)
  po::parsed_options cmdline_parsed = po::basic_command_line_parser<char>(argc, argv).options(desc).positional(pdesc).run();
  po::variables_map vm;
  po::store(cmdline_parsed, vm);
  po::notify(vm);

  if (vm.count("verbose") > 0)
    verbose = true;

  bool quit = false;
  if (vm.count("version") > 0)
  {
    cerr << "MeshFix version 0.1" << endl;
    cerr << "Computer Graphics Lab, Stanford University, 2011" << endl;
    quit = true;
  }

  if (argc <= 2 || vm.count("help") > 0)
  {
    if (quit) cerr << endl;
    cerr << usage << endl;
    cerr << visible << endl;
    quit = true;
  }

  if (quit)
    return 0;

  if (vm.count("infile") <= 0 || vm.count("outfile") <= 0)
  {
    cerr << "Both input and output files should be specified" << endl;
    return -1;
  }

  string lc_out = toLower(outfile);
  if (!(endsWith(lc_out, ".3ds") || endsWith(lc_out, ".obj") || endsWith(lc_out, ".off")))
  {
    cerr << "Unrecognized mesh output format" << endl;
    return -1;
  }

  bool do_something = false;

  //============================================================================================================================
  // del_danglers
  //============================================================================================================================

  if (vm.count("del-danglers") > 0)
    do_del_danglers = do_something = true;

  if (vm.count("zipper") > 0)
    do_zipper = do_something = true;

  v_weld_boundary_only = true;
  if (vm.count("v-weld") > 0)
  {
    do_v_weld = do_something = true;
    v_weld_boundary_only = false;
  }

  // Must be tested after v_weld, to ensure v_weld_boundary_only is set correctly
  if (vm.count("v-weld-boundary") > 0)
    do_v_weld = do_something = true;

  if (vm.count("t-juncts") > 0)
    do_t_juncts = do_something = true;

  if (!do_something)
  {
    cerr << "No repair operation specified, no output written" << endl;
    return 0;
  }

  return 1;
}

double
scaledTolerance(Mesh const & mesh, double relative_tolerance)
{
  return relative_tolerance * mesh.getBounds().getExtent().length();
}

bool
delDanglers(Mesh & mesh)
{
  long nv = mesh.numVertices(), ne = mesh.numEdges(), nf = mesh.numFaces();

  mesh.removeDanglers();

  if (verbose)
  {
    cout << "del-danglers('" << mesh.getName() << "'): Removed " << (nv - mesh.numVertices()) << " vertices, "
         << (ne - mesh.numEdges()) << " edges and " << (nf - mesh.numFaces()) << " faces" << endl;
  }

  return false;
}

bool
zipper(Mesh & mesh)
{
  if (zipper_tolerance < 0)
  {
    zipper_tolerance = 0.001;
    if (verbose)
      cout << "Using default zipper tolerance: " << zipper_tolerance << endl;
  }

  long nv = mesh.numVertices(), ne = mesh.numEdges(), nf = mesh.numFaces();

  mesh.sealSeams((Real)scaledTolerance(mesh, zipper_tolerance));
  mesh.removeDanglers();

  if (verbose)
  {
    cout << "zipper('" << mesh.getName() << "', " << zipper_tolerance << "): Removed " << (nv - mesh.numVertices())
         << " vertices, " << (ne - mesh.numEdges()) << " edges and " << (nf - mesh.numFaces()) << " faces" << endl;
  }

  return false;
}

bool
vWeld(Mesh & mesh)
{
  if (v_weld_tolerance < 0)
  {
    v_weld_tolerance = 0.001;
    if (verbose)
      cout << "Using default vertex welding tolerance: " << v_weld_tolerance << endl;
  }

  long nv = mesh.numVertices();

  VertexWelder welder((Real)scaledTolerance(mesh, v_weld_tolerance));
  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
  {
    Mesh::Vertex * vx = &(*vi);

    if (v_weld_boundary_only && !vx->isBoundary())
      continue;

    Mesh::Vertex * existing = (Mesh::Vertex *)welder.getVertex(vx->getPosition());
    if (existing)
      mesh.replaceVertex(vx, existing);
    else
      welder.addVertex(vx, vx->getPosition());
  }

  mesh.removeDanglers();

  if (verbose)
  {
    cout << (v_weld_boundary_only ? "v-weld-boundary('" : "v-weld('") << mesh.getName() << "', " << v_weld_tolerance
         << "): Removed " << (nv - mesh.numVertices()) << " vertices" << endl;
  }

  return false;
}

class VertexFilter : public Algorithms::Filter<Mesh::Vertex *>
{
  public:
    VertexFilter(Mesh::Edge * edge_) : edge(edge_) {}
    bool allows(Mesh::Vertex * const & vx) const { return !edge->hasEndpoint(vx); }

  private:
    Mesh::Edge * edge;

}; // class VertexFilter

bool
segFromEdge(Mesh::Edge const & edge, double tol, LineSegment3 & seg)
{
  LineSegment3 full_seg(edge.getEndpoint(0)->getPosition(), edge.getEndpoint(1)->getPosition());
  Real len = full_seg.length();
  if (len < 2 * tol)
    return false;

  double s = tol / len;
  double t = 1 - s;
  seg = LineSegment3((1 - (Real)s) * full_seg.getEndpoint(0) + (Real)s * full_seg.getEndpoint(1),
                     (1 - (Real)t) * full_seg.getEndpoint(0) + (Real)t * full_seg.getEndpoint(1));

  return true;
}

bool
tJuncts(Mesh & mesh)
{
  using namespace Algorithms;

  if (t_juncts_tolerance < 0)
  {
    t_juncts_tolerance = 0.001;
    if (verbose)
      cout << "Using default t-junction tolerance: " << t_juncts_tolerance << endl;
  }

  double tol = scaledTolerance(mesh, t_juncts_tolerance);
  double sqtol = tol * tol;

  TheaArray<Mesh::Vertex *> boundary_verts;
  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    if (vi->isBoundary())
      boundary_verts.push_back(&(*vi));

  typedef KDTree3<Mesh::Vertex *> VertexKDTree;
  VertexKDTree kdtree(boundary_verts.begin(), boundary_verts.end());

  TheaArray<Mesh::Edge *> boundary_edges;
  TheaArray<LineSegment3> boundary_segs;
  LineSegment3 seg;
  for (Mesh::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
    if (ei->isBoundary())
    {
      if (segFromEdge(*ei, tol, seg))
      {
        boundary_edges.push_back(&(*ei));
        boundary_segs.push_back(seg);
      }
    }

  if (t_juncts_iters <= 0)
    t_juncts_iters = 3;

  long num_fixed = 0;
  for (int n = 0; n < t_juncts_iters; ++n)
  {
    for (array_size_t i = 0; i < boundary_segs.size(); ++i)
    {
      VertexFilter filter(boundary_edges[i]);
      kdtree.pushFilter(&filter);
        long index = kdtree.closestElement<MetricL2>(boundary_segs[i], 2 * tol);  // add 2 for a little safety margin
      kdtree.popFilter();

      if (index < 0)
        continue;

      Mesh::Vertex * vx = kdtree.getElements()[(size_t)index];
      Vector3 cp = boundary_segs[i].closestPoint(vx->getPosition());

      if ((cp - boundary_edges[i]->getEndpoint(0)->getPosition()).squaredLength() > sqtol
       && (cp - boundary_edges[i]->getEndpoint(1)->getPosition()).squaredLength() > sqtol
       && (cp - vx->getPosition()).squaredLength() < sqtol)
      {
        Mesh::Edge * edge = boundary_edges[i];
        Mesh::Edge * new_edge = mesh.splitEdge(edge, vx);
        if (!new_edge)
        {
          THEA_WARNING << mesh.getName() << ": Could not split boundary edge";
          break;
        }

        if (segFromEdge(*edge, tol, seg))
          boundary_segs[i] = seg;

        if (segFromEdge(*new_edge, tol, seg))
        {
          boundary_edges.push_back(new_edge);
          boundary_segs.push_back(seg);
        }

        num_fixed++;
      }
    }
  }

  mesh.removeDanglers();

  if (verbose)
  {
    cout << "t-juncts('" << mesh.getName() << "', " << t_juncts_tolerance << ", " << t_juncts_iters << "): Fixed " << num_fixed
         << " t-junctions" << endl;
  }

  return false;
}
