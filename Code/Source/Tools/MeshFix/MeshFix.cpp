#include "../../Common.hpp"
#include "../../FilePath.hpp"
#include "../../Algorithms/MeshFeatures/Local/ShapeDiameter.hpp"
#include "../../Algorithms/ConnectedComponents.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Graphics/VertexWelder.hpp"
#include "../../LineSegment3.hpp"
#include "../../Math.hpp"
#include <boost/program_options.hpp>
#include <algorithm>
#include <iostream>
#include <queue>
#include <set>
#include <string>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace MeshFeatures;
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

struct FaceAttribute
{
  int flag;

  void draw(RenderSystem & render_system, RenderOptions const & options) const {}  // noop
};

typedef GeneralMesh< Graphics::NullAttribute,  // vertex attribute
                     Graphics::NullAttribute,  // edge attribute
                     FaceAttribute             // face attribute
                   > Mesh;
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

bool do_orient = false;
bool do_orient_sdf = false;
bool do_orient_majority = false;

int parseArgs(int argc, char * argv[]);
bool delDanglers(Mesh & mesh);
bool zipper(Mesh & mesh);
bool vWeld(Mesh & mesh);
bool tJuncts(Mesh & mesh);
bool orient(Mesh & mesh);
bool orientMajority(Mesh & mesh);
void orientSDF(MG & mesh_group);
bool checkProblems(Mesh & mesh);

int
meshFix(int argc, char * argv[])
{
  int parse_status = parseArgs(argc, argv);
  if (parse_status <= 0)
    return parse_status;

  MG mg(FilePath::objectName(infile));
  mg.load(infile);
  mg.updateBounds();  // load() should do this, but let's play safe
  if (mg.isEmpty())
  {
    THEA_ERROR << "Input mesh is empty, no output written";
    return -1;
  }

  // The order of operations is important here -- e.g. it's better to do orient() after welding vertices
  if (do_del_danglers)
  {
    mg.forEachMeshUntil(&delDanglers);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_t_juncts)
  {
    mg.forEachMeshUntil(&tJuncts);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_v_weld)
  {
    mg.forEachMeshUntil(&vWeld);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_zipper)
  {
    mg.forEachMeshUntil(&zipper);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_orient)
  {
    mg.forEachMeshUntil(&orient);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_orient_sdf)
  {
    orientSDF(mg);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_orient_majority)
  {
    mg.forEachMeshUntil(&orientMajority);
    mg.forEachMeshUntil(&checkProblems);
  }

  string lc_out = toLower(outfile);
  if (endsWith(lc_out, ".3ds"))
    mg.save(outfile, Codec3DS<Mesh>());
  else if (endsWith(lc_out, ".obj"))
    mg.save(outfile, CodecOBJ<Mesh>());
  else if (endsWith(lc_out, ".off"))
    mg.save(outfile, CodecOFF<Mesh>());
  else
  {
    THEA_ERROR << "Unrecognized mesh output format";
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

          ("del-danglers",        "Delete isolated vertices and edges, and empty faces")

          ("zipper",              po::value<double>(&zipper_tolerance),
                                  "Seal pairs of boundary edges whose endpoints are equal upto the specified tolerance"
                                  " (implies --del-danglers)")

          ("v-weld",              po::value<double>(&v_weld_tolerance),
                                  "Merge vertices closer than the specified tolerance (implies --del-danglers)")

          ("v-weld-boundary",     po::value<double>(&v_weld_tolerance),
                                  "Merge boundary vertices closer than the specified tolerance (implies --del-danglers)")

          ("t-juncts",            po::value<double>(&t_juncts_tolerance),
                                  "Split boundary edges at t-junctions so they can be zippered properly (always executed before"
                                  " any welding or zippering)")

          ("tj-iters",            po::value<int>(&t_juncts_iters),
                                  "Number of iterations for fixing t-junctions")

          ("orient",              "Consistently orient each edge-connected component of the mesh so that normals defined by"
                                  " counter-clockwise winding point inside-out")

          ("orient-sdf",          "Consistently orient each edge-connected component of the mesh so that normals defined by"
                                  " counter-clockwise winding point inside-out, using the direction of smaller shape diameter"
                                  " as a heuristic")

          ("orient-majority",     "Consistently orient each edge-connected component of the mesh so that normals defined by"
                                  " counter-clockwise winding point inside-out, flipping faces that disagree most with their"
                                  " neighbors first")
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
    THEA_CONSOLE << "MeshFix version 0.1";
    THEA_CONSOLE << "Computer Graphics Lab, Stanford University, 2011";
    quit = true;
  }

  if (argc <= 2 || vm.count("help") > 0)
  {
    if (quit) THEA_CONSOLE << "";
    THEA_CONSOLE << usage;
    THEA_CONSOLE << visible;
    quit = true;
  }

  if (quit)
    return 0;

  if (vm.count("infile") <= 0 || vm.count("outfile") <= 0)
  {
    THEA_ERROR << "Both input and output files should be specified";
    return -1;
  }

  string lc_out = toLower(outfile);
  if (!(endsWith(lc_out, ".3ds") || endsWith(lc_out, ".obj") || endsWith(lc_out, ".off")))
  {
    THEA_ERROR << "Unrecognized mesh output format";
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

  if (vm.count("orient") > 0)
    do_orient = do_something = true;

  if (vm.count("orient-sdf") > 0)
    do_orient_sdf = do_something = true;

  if (vm.count("orient-majority") > 0)
    do_orient_majority = do_something = true;

  if (!do_something)
  {
    THEA_WARNING << "No repair operation specified, no output written";
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
    THEA_CONSOLE << "del-danglers('" << mesh.getName() << "'): Removed " << (nv - mesh.numVertices()) << " vertices, "
                 << (ne - mesh.numEdges()) << " edges and " << (nf - mesh.numFaces()) << " faces";
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
      THEA_CONSOLE << "Using default zipper tolerance: " << zipper_tolerance;
  }

  long nv = mesh.numVertices(), ne = mesh.numEdges(), nf = mesh.numFaces();

  mesh.sealSeams((Real)scaledTolerance(mesh, zipper_tolerance));

  // checkProblems(mesh);

  mesh.removeDanglers();

  if (verbose)
  {
    THEA_CONSOLE << "zipper('" << mesh.getName() << "', " << zipper_tolerance << "): Removed " << (nv - mesh.numVertices())
                 << " vertices, " << (ne - mesh.numEdges()) << " edges and " << (nf - mesh.numFaces()) << " faces";
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
      THEA_CONSOLE << "Using default vertex welding tolerance: " << v_weld_tolerance;
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
    {
      // Don't weld vertices connected by an edge or a face
      bool are_connected = false;
      if (existing->hasEdgeTo(vx))
        are_connected = true;

      if (!are_connected)
      {
        for (Mesh::Vertex::FaceConstIterator vfi = vx->facesBegin(); vfi != vx->facesEnd(); ++vfi)
          if (existing->hasIncidentFace(*vfi))
          {
            are_connected = true;
            break;
          }
      }

      if (!are_connected)
        mesh.replaceVertex(vx, existing);
    }
    else
      welder.addVertex(vx, vx->getPosition());
  }

  // checkProblems(mesh);

  mesh.removeDanglers();

  if (verbose)
  {
    THEA_CONSOLE << (v_weld_boundary_only ? "v-weld-boundary('" : "v-weld('") << mesh.getName() << "', " << v_weld_tolerance
                 << "): Removed " << (nv - mesh.numVertices()) << " vertices";
  }

  return false;
}

class VertexFilter : public Filter<Mesh::Vertex *>
{
  public:
    VertexFilter(Mesh::Edge * edge_) : edge(edge_) {}
    bool allows(Mesh::Vertex * const & vx) const
    {
      if (edge->hasEndpoint(vx))
        return false;

      for (Mesh::Edge::FaceConstIterator efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
        if ((*efi)->hasVertex(vx))
          return false;

      return true;
    }

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
  if (t_juncts_tolerance < 0)
  {
    t_juncts_tolerance = 0.001;
    if (verbose)
      THEA_CONSOLE << "Using default t-junction tolerance: " << t_juncts_tolerance;
  }

  double tol = scaledTolerance(mesh, t_juncts_tolerance);
  double sqtol = tol * tol;

  TheaArray<Mesh::Vertex *> boundary_verts;
  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    if (vi->isBoundary())
      boundary_verts.push_back(&(*vi));

  typedef KDTreeN<Mesh::Vertex *, 3> VertexKDTree;
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
      Mesh::Edge * edge = boundary_edges[i];

      VertexFilter filter(edge);
      kdtree.pushFilter(&filter);
        long index = kdtree.closestElement<MetricL2>(boundary_segs[i], 2 * tol);  // add 2 for a little safety margin
      kdtree.popFilter();

      if (index < 0)
        continue;

      Mesh::Vertex * vx = kdtree.getElements()[(size_t)index];
      Vector3 cp = boundary_segs[i].closestPoint(vx->getPosition());

      if ((cp - edge->getEndpoint(0)->getPosition()).squaredLength() > sqtol
       && (cp - edge->getEndpoint(1)->getPosition()).squaredLength() > sqtol
       && (cp - vx->getPosition()).squaredLength() < sqtol)
      {
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

  // checkProblems(mesh);

  mesh.removeDanglers();

  if (verbose)
  {
    THEA_CONSOLE << "t-juncts('" << mesh.getName() << "', " << t_juncts_tolerance << ", " << t_juncts_iters << "): Fixed "
                 << num_fixed << " t-junctions";
  }

  return false;
}

bool
consistentWinding(Mesh::Face const * face0, Mesh::Face const * face1, Mesh::Edge const * edge)
{
  Mesh::Face const * f[2] = { face0, face1 };
  int in_order[2] = { -1, -1 };

  for (int i = 0; i < 2; ++i)
  {
    for (Mesh::Face::VertexConstIterator fvj = f[i]->verticesBegin(); fvj != f[i]->verticesEnd(); ++fvj)
    {
      Mesh::Vertex const * v0 = *fvj;
      Mesh::Vertex const * v1 = f[i]->getSuccessor(fvj);

      if (edge->getEndpoint(0) == v0 && edge->getEndpoint(1) == v1)
      {
        in_order[i] = 1;
        break;
      }
      else if (edge->getEndpoint(0) == v1 && edge->getEndpoint(1) == v0)
      {
        in_order[i] = 0;
        break;
      }
    }

    if (in_order[i] < 0)
    {
      THEA_CONSOLE << "face = " << f[i] << ", e[0] = " << edge->getEndpoint(0) << ", e[1] = " << edge->getEndpoint(1);
      for (Mesh::Face::VertexConstIterator fvj = f[i]->verticesBegin(); fvj != f[i]->verticesEnd(); ++fvj)
        THEA_CONSOLE << "v = " << *fvj;
    }

    alwaysAssertM(in_order[i] >= 0, "Edge not part of face");
  }

  return in_order[0] != in_order[1];
}

bool
orient(Mesh & mesh)
{
  TheaArray< TheaArray<Mesh::Face *> > cc;
  long num_cc = ConnectedComponents::findEdgeConnected(mesh, cc);

  enum { NOT_VISITED, VISITED };

  for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    fi->attr().flag = NOT_VISITED;

  long num_flipped = 0;
  for (array_size_t i = 0; i < cc.size(); ++i)
  {
    if (cc.empty())  // shouldn't happen, but check just in case
    {
      THEA_WARNING << mesh.getName() << ": Empty connected component";
      continue;
    }

    // Find vertex with max X
    Mesh::Vertex * xmax_vertex = NULL;
    for (array_size_t j = 0; j < cc[i].size(); ++j)
      for (Mesh::Face::VertexConstIterator fvi = cc[i][j]->verticesBegin(); fvi != cc[i][j]->verticesEnd(); ++fvi)
        if (!xmax_vertex || (*fvi)->getPosition().x() > xmax_vertex->getPosition().x())
          xmax_vertex = *fvi;

    if (!xmax_vertex)  // shouldn't happen, but check just in case
    {
      THEA_WARNING << mesh.getName() << ": No vertex with maximum X found";
      continue;
    }

    // Find the "most X-facing" face incident on xmax_vertex. This face should always have a normal with positive X-component.
    Mesh::Face * xmax_face = NULL;
    Real xmax_face_nx = -1;
    for (Mesh::Vertex::FaceConstIterator vfi = xmax_vertex->facesBegin(); vfi != xmax_vertex->facesEnd(); ++vfi)
    {
      Mesh::Face * face = *vfi;
      Real nx = fabs(face->getNormal().x());
      if (!xmax_face || nx > xmax_face_nx)
      {
        xmax_face = face;
        xmax_face_nx = nx;
      }
    }

    if (!xmax_face)  // shouldn't happen, but check just in case
    {
      THEA_WARNING << mesh.getName() << ": No face with maximum X found";
      continue;
    }

    // Fix the orientation of this face
    if (xmax_face->getNormal().x() < 0)  // assume this is never zero... that would be a degenerate case
    {
      xmax_face->reverseWinding();
      num_flipped++;
    }

    xmax_face->attr().flag = VISITED;

    queue<Mesh::Face *> q;
    q.push(xmax_face);
    while (!q.empty())
    {
      Mesh::Face * face = q.front();
      q.pop();

      if (!face)
      {
        THEA_WARNING << mesh.getName() << ": Null face";
        continue;
      }

      alwaysAssertM(face->attr().flag == VISITED, string(mesh.getName()) + ": Unvisited face in queue");

      for (Mesh::Face::EdgeConstIterator fei = face->edgesBegin(); fei != face->edgesEnd(); ++fei)
      {
        Mesh::Edge * edge = *fei;
        for (Mesh::Edge::FaceConstIterator efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
        {
           Mesh::Face * adj_face = *efi;
           if (adj_face->attr().flag == VISITED)
             continue;

           if (!consistentWinding(face, adj_face, edge))
           {
             adj_face->reverseWinding();
             num_flipped++;
           }

           adj_face->attr().flag = VISITED;
           q.push(adj_face);
        }
      }
    }
  }

  for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    fi->updateNormal();

  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    if (!vi->hasPrecomputedNormal())
      vi->updateNormal();

  // checkProblems(mesh);

  if (verbose)
  {
    THEA_CONSOLE << "orient('" << mesh.getName() << "): Flipped " << num_flipped << '/' << mesh.numFaces() << " faces in "
                 << num_cc << " connected components";
  }

  return false;
}

struct SDFOrienter
{
  SDFOrienter(Local::ShapeDiameter<Mesh> * sdf_) : sdf(sdf_) {}

  bool operator()(Mesh & mesh)
  {
    long num_flipped = 0;
    for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      if (fi->numVertices() < 3)
        continue;

      Vector3 centroid = Vector3::zero();
      for (Mesh::Face::VertexConstIterator fvi = fi->verticesBegin(); fvi != fi->verticesEnd(); ++fvi)
        centroid += (*fvi)->getPosition();

      centroid /= fi->numVertices();

      Vector3 n = fi->getNormal();
      double sdf_pos = sdf->compute(centroid,  n, false);
      double sdf_neg = sdf->compute(centroid, -n, false);
#if 0
      if ((sdf_neg < 0 || sdf_neg > 0.4) && (sdf_pos < 0 || sdf_pos > 0.4))
      {
        if (fi->getNormal().z() < 0)
        {
          fi->reverseWinding();
          num_flipped++;
        }
      }
      else
#endif
      if (sdf_neg >= 0 && (sdf_pos < 0 || sdf_neg < sdf_pos))
      {
        fi->reverseWinding();
        num_flipped++;
      }
    }

    if (verbose)
    {
      THEA_CONSOLE << "orient-sdf('" << mesh.getName() << "): Flipped " << num_flipped << '/' << mesh.numFaces() << " faces";
    }

    return false;
  }

  Local::ShapeDiameter<Mesh> * sdf;

}; // struct SDFOrienter

void
orientSDF(MG & mesh_group)
{
  Local::ShapeDiameter<Mesh> sdf(mesh_group);
  SDFOrienter func(&sdf);
  mesh_group.forEachMeshUntil(&func);
}

struct CountComparator
{
  bool operator()(Mesh::Face const * f0, Mesh::Face const * f1)
  {
    return f0->attr().flag > f1->attr().flag || (f0->attr().flag == f1->attr().flag && f0 < f1);
  }
};

bool
orientMajority(Mesh & mesh)
{
  for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    fi->attr().flag = 0;

  for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
  {
    Mesh::Face * face = &(*fi);
    for (Mesh::Face::EdgeIterator fei = face->edgesBegin(); fei != face->edgesEnd(); ++fei)
    {
      Mesh::Edge * edge = *fei;
      for (Mesh::Edge::FaceIterator efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
      {
        Mesh::Face * nbr = *efi;
        if (face == nbr)
          continue;

        if (!consistentWinding(face, nbr, edge))
        {
          face->attr().flag++;
          nbr->attr().flag++;
        }
      }
    }
  }

  typedef set<Mesh::Face *, CountComparator> FaceMaxHeap;
  FaceMaxHeap heap;
  for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    heap.insert(&(*fi));

  long num_flipped = false;
  while (!heap.empty())
  {
    Mesh::Face * face = *heap.begin();
    heap.erase(face);

    if (face->attr().flag <= 0)
    {
      face->attr().flag = -1;  // visited
      continue;
    }

    face->attr().flag = -1;  // visited

    for (Mesh::Face::EdgeIterator fei = face->edgesBegin(); fei != face->edgesEnd(); ++fei)
    {
      Mesh::Edge * edge = *fei;
      for (Mesh::Edge::FaceIterator efi = edge->facesBegin(); efi != edge->facesEnd(); ++efi)
      {
        Mesh::Face * nbr = *efi;
        if (face == nbr || nbr->attr().flag < 0)
          continue;

        if (!consistentWinding(face, nbr, edge))
        {
          alwaysAssertM(nbr->attr().flag > 0,
                        "orient-majority('" + string(mesh.getName()) + "): Expected > 0 neighbor consistency count");

          heap.erase(nbr);
          nbr->attr().flag--;
          heap.insert(nbr);
        }
      }
    }

    face->reverseWinding();
    num_flipped++;
  }

  if (verbose)
  {
    THEA_CONSOLE << "orient-majority('" << mesh.getName() << "): Flipped " << num_flipped << '/' << mesh.numFaces() << " faces";
  }

  return false;
}

bool
checkProblems(Mesh & mesh)
{
  long num_faces_with_repeated_vertices = 0;
  long num_faces_with_repeated_edges = 0;

  long num_edges_with_repeated_vertices = 0;
  long num_edges_with_repeated_faces = 0;

  long num_vertices_with_repeated_edges = 0;
  long num_vertices_with_repeated_faces = 0;

  for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
  {
    Mesh::Face const & face = *fi;

    for (Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
      if (std::count(face.verticesBegin(), face.verticesEnd(), *vi) > 1)
      {
        num_faces_with_repeated_vertices++;
        break;
      }

    for (Mesh::Face::EdgeConstIterator ei = face.edgesBegin(); ei != face.edgesEnd(); ++ei)
      if (std::count(face.edgesBegin(), face.edgesEnd(), *ei) > 1)
      {
        num_faces_with_repeated_edges++;
        break;
      }
  }

  for (Mesh::EdgeConstIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
  {
    Mesh::Edge const & edge = *ei;

    if (edge.getEndpoint(0) == edge.getEndpoint(1))
      num_edges_with_repeated_vertices++;

    for (Mesh::Edge::FaceConstIterator fi = edge.facesBegin(); fi != edge.facesEnd(); ++fi)
      if (std::count(edge.facesBegin(), edge.facesEnd(), *fi) > 1)
      {
        num_edges_with_repeated_faces++;
        break;
      }
  }

  for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
  {
    Mesh::Vertex const & vertex = *vi;

    for (Mesh::Vertex::EdgeConstIterator ei = vertex.edgesBegin(); ei != vertex.edgesEnd(); ++ei)
      if (std::count(vertex.edgesBegin(), vertex.edgesEnd(), *ei) > 1)
      {
        num_vertices_with_repeated_edges++;
        break;
      }

    for (Mesh::Vertex::FaceConstIterator fi = vertex.facesBegin(); fi != vertex.facesEnd(); ++fi)
      if (std::count(vertex.facesBegin(), vertex.facesEnd(), *fi) > 1)
      {
        num_vertices_with_repeated_faces++;
        break;
      }
  }

  if (num_faces_with_repeated_vertices > 0)
    THEA_CONSOLE << mesh.getName() << ": " << num_faces_with_repeated_vertices << " faces with repeated vertices";
  if (num_faces_with_repeated_edges > 0)
    THEA_CONSOLE << mesh.getName() << ": " << num_faces_with_repeated_edges << " faces with repeated edges";
  if (num_edges_with_repeated_vertices > 0)
    THEA_CONSOLE << mesh.getName() << ": " << num_edges_with_repeated_vertices << " self-loop edges";
  if (num_edges_with_repeated_faces > 0)
    THEA_CONSOLE << mesh.getName() << ": " << num_edges_with_repeated_faces << " edges with repeated faces";
  if (num_vertices_with_repeated_edges > 0)
    THEA_CONSOLE << mesh.getName() << ": " << num_vertices_with_repeated_edges << " vertices with repeated edges";
  if (num_vertices_with_repeated_faces > 0)
    THEA_CONSOLE << mesh.getName() << ": " << num_vertices_with_repeated_faces << " vertices with repeated faces";

  return false;
}
