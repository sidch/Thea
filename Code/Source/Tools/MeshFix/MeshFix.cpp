#include "../../Common.hpp"
#include "../../Algorithms/SurfaceFeatures/Local/ShapeDiameter.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/ConnectedComponents.hpp"
#include "../../Algorithms/BvhN.hpp"
#include "../../Algorithms/MeshBvh.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Algorithms/RayIntersectionTester.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Graphics/VertexWelder.hpp"
#include "../../FilePath.hpp"
#include "../../Hash.hpp"
#include "../../LineSegment3.hpp"
#include "../../Math.hpp"
#include "../../Random.hpp"
#include "../../UnorderedSet.hpp"
#include "../../UnorderedMap.hpp"
#include "../../ThirdParty/CLI11/CLI11.hpp"
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <queue>
#include <set>
#include <string>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace SurfaceFeatures;
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
  THEA_CATCH(return -1;, ERROR, "%s", "An error occurred")

  return status;
}

typedef Array<string> LabelArray;
typedef UnorderedMap<string, int> LabelIndexMap;

struct FaceAttribute
{
  int label_index;
  int flag;

  FaceAttribute() : label_index(-1), flag(0) {}

  void draw(IRenderSystem & render_system, IRenderOptions const & options) const {}  // noop
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

bool no_normals = false;
bool no_texcoords = false;
bool no_empty = false;

bool do_flatten = false;
bool do_del_danglers = false;

bool do_del_dup_faces = false;
bool del_dup_faces_sorted = true;

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
bool do_orient_visibility = false;
bool orient_visibility_hi_qual = false;

bool do_triangulate = false;
double triangulate_tolerance = -1;

bool export_face_labels = false;

int parseArgs(int argc, char * argv[]);
void flatten(MG & mesh_group);
bool delDanglers(Mesh & mesh);
void delDuplicateFaces(MG & mesh_group, bool sorted);
bool zipper(Mesh & mesh);
bool vWeld(Mesh & mesh);
bool tJuncts(Mesh & mesh);
bool orient(Mesh & mesh);
bool orientMajority(Mesh & mesh);
void orientSDF(MG & mesh_group);
void orientVisibility(MG & mesh_group);
bool triangulate(Mesh & mesh);
bool checkProblems(Mesh & mesh);

struct ReadCallback : public MG::ReadCallback
{
  ReadCallback(LabelArray & labels_) : labels(labels_) {}

  void faceRead(Mesh * mesh, intx index, Mesh::FaceHandle face)
  {
    string mesh_name = mesh->getName();
    if (mesh_name.empty())
      face->attr().label_index = -1;
    else
    {
      LabelIndexMap::const_iterator existing = label_indices.find(mesh_name);
      if (existing == label_indices.end())
      {
        face->attr().label_index = (int)labels.size();
        label_indices[mesh_name] = face->attr().label_index;
        labels.push_back(mesh_name);
      }
      else
        face->attr().label_index = existing->second;
    }
  }

  LabelArray & labels;
  LabelIndexMap label_indices;

}; // struct ReadCallback

struct WriteCallback : public MG::WriteCallback
{
  WriteCallback(LabelArray const & labels_) : labels(labels_), faces_per_label(labels_.size()) {}

  void faceWritten(Mesh const * mesh, intx index, Mesh::FaceConstHandle face)
  {
    intx lab = face->attr().label_index;
    alwaysAssertM(lab < (intx)labels.size(), "Face label index out of bounds");

    if (lab < 0)
      unlabeled_faces.push_back(index);
    else
      faces_per_label[(size_t)lab].push_back(index);
  }

  LabelArray const & labels;
  Array< Array<intx> > faces_per_label;
  Array<intx> unlabeled_faces;

}; // struct WriteCallback

int
meshFix(int argc, char * argv[])
{
  int parse_status = parseArgs(argc, argv);
  if (parse_status <= 0)
    return parse_status;

  Codec3ds<Mesh>::ReadOptions read_opts_3ds = Codec3ds<Mesh>::ReadOptions::defaults();
  read_opts_3ds.setReadTexCoords(!no_texcoords)
               .setSkipEmptyMeshes(no_empty);
  Codec3ds<Mesh>::Ptr codec_3ds(new Codec3ds<Mesh>(read_opts_3ds, Codec3ds<Mesh>::WriteOptions::defaults()));

  CodecObj<Mesh>::ReadOptions read_opts_obj = CodecObj<Mesh>::ReadOptions::defaults();
  read_opts_obj.setReadTexCoords(!no_texcoords)
               .setReadNormals(!no_normals)
               .setSkipEmptyMeshes(no_empty);
  CodecObj<Mesh>::WriteOptions write_opts_obj = CodecObj<Mesh>::WriteOptions::defaults();
  write_opts_obj.setWriteTexCoords(!no_texcoords)
                .setWriteNormals(!no_normals)
                .setSkipEmptyMeshes(no_empty);
  CodecObj<Mesh>::Ptr codec_obj(new CodecObj<Mesh>(read_opts_obj, write_opts_obj));

  CodecOff<Mesh>::ReadOptions opts_off = CodecOff<Mesh>::ReadOptions::defaults();
  opts_off.setSkipEmptyMeshes(no_empty);
  CodecOff<Mesh>::Ptr codec_off(new CodecOff<Mesh>(opts_off));

  CodecPly<Mesh>::ReadOptions opts_ply = CodecPly<Mesh>::ReadOptions::defaults();
  opts_ply.setSkipEmptyMeshes(no_empty);
  CodecPly<Mesh>::Ptr codec_ply(new CodecPly<Mesh>(opts_ply));

  Array<MeshCodec<Mesh>::Ptr> codecs;
  codecs.push_back(codec_3ds);
  codecs.push_back(codec_obj);
  codecs.push_back(codec_off);
  codecs.push_back(codec_3ds);

  LabelArray labels;
  ReadCallback read_callback(labels);

  MG mg(FilePath::objectName(infile));
  mg.load(infile, codecs, &read_callback);

  mg.updateBounds();  // load() should do this, but let's play safe
  if (mg.empty())
  {
    THEA_ERROR << "Input mesh is empty, no output written";
    return -1;
  }

  // The order of operations is important here -- e.g. it's better to do orient() after welding vertices. flatten() must happen
  // right at the beginning.

  if (do_flatten)
  {
    flatten(mg);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_del_danglers)
  {
    mg.forEachMeshUntil(delDanglers);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_del_dup_faces)
  {
    delDuplicateFaces(mg, del_dup_faces_sorted);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_t_juncts)
  {
    mg.forEachMeshUntil(tJuncts);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_v_weld)
  {
    mg.forEachMeshUntil(vWeld);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_zipper)
  {
    mg.forEachMeshUntil(zipper);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_orient)
  {
    mg.forEachMeshUntil(orient);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_orient_sdf)
  {
    orientSDF(mg);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_orient_majority)
  {
    mg.forEachMeshUntil(orientMajority);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_orient_visibility)
  {
    orientVisibility(mg);
    mg.forEachMeshUntil(checkProblems);
  }

  if (do_triangulate)
  {
    if (mg.forEachMeshUntil(triangulate)) return -1;
    mg.forEachMeshUntil(checkProblems);
  }

  string lc_out = toLower(outfile);
  if (endsWith(lc_out, ".off.bin"))
  {
    CodecOff<Mesh>::WriteOptions write_opts_off = CodecOff<Mesh>::WriteOptions().setBinary(true);
    *codec_off = CodecOff<Mesh>(opts_off, write_opts_off);
  }

  WriteCallback write_callback(labels);
  mg.save(outfile, codecs, (export_face_labels ? &write_callback : nullptr));

  if (export_face_labels)
  {
    string lab_path = FilePath::changeCompleteExtension(outfile, "lab");
    ofstream out(lab_path.c_str(), ios::binary);
    if (!out)
    {
      THEA_ERROR << "Could not open labels output file '" << lab_path << '\'';
      return -1;
    }

    bool first = true;
    for (size_t i = 0; i < labels.size(); ++i)
    {
      if (write_callback.faces_per_label[i].empty())
        continue;

      if (first) first = false;
      else out << '\n';

      out << labels[i] << '\n';
      for (size_t j = 0; j < write_callback.faces_per_label[i].size(); ++j)
      {
        if (j > 0) out << ' ';
        out << write_callback.faces_per_label[i][j] + 1;
      }
      out << '\n';
    }

    if (!write_callback.unlabeled_faces.empty())
    {
      if (first) first = false;
      else out << '\n';

      out << "Unlabeled\n";
      for (size_t j = 0; j < write_callback.unlabeled_faces.size(); ++j)
      {
        if (j > 0) out << ' ';
        out << write_callback.unlabeled_faces[j] + 1;
      }
      out << '\n';
    }
  }

  return 0;
}

int
parseArgs(int argc, char * argv[])
{
  CLI::App app{"MeshFix"};

  string s_orient_visibility_qual;

  app.add_option("infile", infile, "Path to input mesh file")->required()->check(CLI::ExistingFile);
  app.add_option("outfile", outfile, "Path to output mesh file")->required();
  app.add_flag("--verbose", verbose, "Print extra debugging output");
  app.add_flag("--no-normals", no_normals, "Ignore vertex and face normals in the input, and don't write any to output");
  app.add_flag("--no-texcoords", no_texcoords,
      "Ignore vertex and face texture coordinates in the input, and don't write any to output");
  app.add_flag("--no-empty", no_empty, "Ignore empty sub-meshes");
  app.add_flag("--flatten", do_flatten, "Flatten the object to a single sub-mesh when reading");
  app.add_flag("--del-danglers", do_del_danglers, "Delete isolated vertices and edges, and empty faces");
  app.add_flag("--del-dup-faces", do_del_dup_faces, "Delete duplicate faces, regardless of orientation");
  app.add_option("--zipper", zipper_tolerance,
      "Seal pairs of boundary edges whose endpoints are equal upto the given tolerance, specified here and for subsequent "
      "options as a fraction of the sub-mesh bounding box diagonal (implies --del-danglers)");
  app.add_option("--v-weld", v_weld_tolerance, "Merge boundary vertices closer than the specified tolerance (implies --del-danglers)");
  app.add_option("--v-weld-boundary", v_weld_tolerance,
      "Merge vertices closer than the specified tolerance (implies --del-danglers)");
  app.add_option("--t-juncts", t_juncts_tolerance,
      "Split boundary edges at t-junctions so they can be zippered properly (always executed before any welding or zippering)");
  app.add_option("--tj-iters", t_juncts_iters, "Number of iterations for fixing t-junctions");
  app.add_flag("--orient",
      "Consistently orient each edge-connected component of the mesh so that normals defined by counter-clockwise winding "
      "point inside-out");
  app.add_flag("--orient-sdf",
      "Consistently orient each face of the mesh so that the normal defined by counter-clockwise winding points inside-out, "
      "using the direction of smaller shape diameter as a heuristic");
  app.add_flag("--orient-majority",
      "Consistently orient each edge-connected component of the mesh so that normals defined by counter-clockwise winding "
      "point inside-out, first flipping faces that disagree most with their neighbors");
  app.add_option("--orient-visibility", s_orient_visibility_qual,
      "Consistently orient each face of the mesh so that the normal defined by counter-clockwise winding points towards the "
      "direction of an external camera from which the face is visible");
  app.add_option("--triangulate", triangulate_tolerance, "Triangulate every face with more than 3 vertices");
  app.add_flag("--export-face-labels", export_face_labels,
      "Along with the output mesh, export a .lab file that associates every face with the name of the submesh it belonged to "
      "in the input. Works even if --flatten is specified.");

  try
  {
    app.parse(argc, argv);
  }
  catch (CLI::ParseError const & e)
  {
    app.exit(e);
    return -1;
  }

  string lc_out = toLower(outfile);
  if (!(endsWith(lc_out, ".3ds") || endsWith(lc_out, ".obj") || endsWith(lc_out, ".off") || endsWith(lc_out, ".off.bin")))
  {
    THEA_ERROR << "Unrecognized mesh output format";
    return -1;
  }

  do_zipper = (app.count("--zipper") > 0);
  do_v_weld = (app.count("--v-weld") > 0 || app.count("--v-weld-boundary") > 0);
  v_weld_boundary_only = (do_v_weld && app.count("--v-weld") <= 0);
  do_t_juncts = (app.count("--t-juncts") > 0);
  do_orient = (app.count("--orient") > 0);
  do_orient_sdf = (app.count("--orient-sdf") > 0);
  do_orient_majority = (app.count("--orient-majority") > 0);
  do_orient_visibility = (app.count("--orient-visibility") > 0);
  do_triangulate = (app.count("--triangulate") > 0);
  bool do_something = no_normals || no_texcoords || no_empty || do_flatten || do_del_danglers || do_del_dup_faces || do_zipper
                   || do_v_weld || do_t_juncts || do_orient || do_orient_sdf || do_orient_majority || do_orient_visibility
                   || do_triangulate || export_face_labels;

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
  return relative_tolerance * mesh.getBounds().getExtent().norm();
}

struct Flattener
{
  Flattener(string const & mesh_name = "Flattened") { flattened = Mesh::Ptr(new Mesh(mesh_name)); }

  bool operator()(Mesh const & mesh)
  {
    typedef UnorderedMap<Mesh::Vertex const *, Mesh::Vertex *> VertexMap;
    VertexMap vmap;

    for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    {
      Mesh::Vertex * new_vertex = flattened->addVertex(vi->getPosition(), -1, &vi->getNormal());
      if (!new_vertex)
        throw Error("Could not add copy of vertex to flattened mesh");

      new_vertex->attr() = vi->attr();
      vmap[&(*vi)] = new_vertex;
    }

    Array<Mesh::Vertex *> face_vertices;
    for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      face_vertices.clear();
      for (Mesh::Face::VertexConstIterator fvi = fi->verticesBegin(); fvi != fi->verticesEnd(); ++fvi)
      {
        VertexMap::const_iterator existing = vmap.find(*fvi);
        if (existing == vmap.end())
          throw Error("Could not find copy of vertex in flattened mesh");

        face_vertices.push_back(existing->second);
      }

      Mesh::Face * new_face = flattened->addFace(face_vertices.begin(), face_vertices.end());
      if (new_face)
        new_face->attr() = fi->attr();
      else
        THEA_WARNING << "Could not add copy of face to flattened mesh";  // face might be malformed, this is not a fatal error
    }

    return false;
  }

  Mesh::Ptr flattened;
};

void
flatten(MG & mesh_group)
{
  Flattener flattener;
  mesh_group.forEachMeshUntil(ref(flattener));
  mesh_group.clear();
  mesh_group.addMesh(flattener.flattened);
}

bool
delDanglers(Mesh & mesh)
{
  intx nv = mesh.numVertices(), ne = mesh.numEdges(), nf = mesh.numFaces();

  mesh.removeDanglers();

  if (verbose)
  {
    THEA_CONSOLE << "del-danglers('" << mesh.getName() << "'): Removed " << (nv - mesh.numVertices()) << " vertices, "
                 << (ne - mesh.numEdges()) << " edges and " << (nf - mesh.numFaces()) << " faces";
  }

  return false;
}

struct FaceSeq
{
  Array<Mesh::Vertex const *> seq;

  FaceSeq(Mesh::Face const & face, bool sorted)
  {
    for (Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
      seq.push_back(*vi);

    if (sorted)
      sort(seq.begin(), seq.end());
  }
};

bool
operator==(FaceSeq const & lhs, FaceSeq const & rhs)
{
  return lhs.seq.size() == rhs.seq.size()
      && equal(lhs.seq.begin(), lhs.seq.end(), rhs.seq.begin());
}

template <>
struct Hasher<FaceSeq>
{
  size_t operator()(FaceSeq const & f) const
  {
    return hashRange(f.seq.begin(), f.seq.end());
  }
};

typedef UnorderedSet< FaceSeq, Hasher<FaceSeq> > FaceSet;

struct DupFaceDeleter
{
  bool sorted;

  DupFaceDeleter(bool sorted_) : sorted(sorted_) {}

  bool operator()(Mesh & mesh)
  {
    FaceSet seqs;
    intx num_deleted = 0;

    for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); )
    {
      FaceSeq s(*fi, sorted);
      if (seqs.find(s) == seqs.end())
      {
        seqs.insert(s);
        ++fi;
      }
      else  // delete the duplicated face
      {
        Mesh::FaceIterator curr = fi; ++fi;
        mesh.removeFace(curr);
        num_deleted++;
      }
    }

    if (verbose)
    {
      THEA_CONSOLE << "del-dup-faces('" << mesh.getName() << "'): Removed " << num_deleted << " faces";
    }

    return false;
  }
};

void
delDuplicateFaces(MG & mesh_group, bool sorted)
{
  mesh_group.forEachMeshUntil(DupFaceDeleter(sorted));
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

  intx nv = mesh.numVertices(), ne = mesh.numEdges(), nf = mesh.numFaces();

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

  intx nv = mesh.numVertices();

  VertexWelder welder((Real)scaledTolerance(mesh, v_weld_tolerance));
  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
  {
    Mesh::Vertex * vx = &(*vi);

    if (v_weld_boundary_only && !vx->isBoundaryVertex())
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

  Array<Mesh::Vertex *> boundary_verts;
  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    if (vi->isBoundaryVertex())
      boundary_verts.push_back(&(*vi));

  typedef BvhN<Mesh::Vertex *, 3> VertexBvh;
  VertexBvh bvh(boundary_verts.begin(), boundary_verts.end());

  Array<Mesh::Edge *> boundary_edges;
  Array<LineSegment3> boundary_segs;
  LineSegment3 seg;
  for (Mesh::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
    if (ei->isBoundaryEdge())
    {
      if (segFromEdge(*ei, tol, seg))
      {
        boundary_edges.push_back(&(*ei));
        boundary_segs.push_back(seg);
      }
    }

  if (t_juncts_iters <= 0)
    t_juncts_iters = 3;

  intx num_fixed = 0;
  for (int n = 0; n < t_juncts_iters; ++n)
  {
    for (size_t i = 0; i < boundary_segs.size(); ++i)
    {
      Mesh::Edge * edge = boundary_edges[i];

      VertexFilter filter(edge);
      bvh.pushFilter(&filter);
        intx index = bvh.closestElement<MetricL2>(boundary_segs[i], /* dist_bound = */ 2 * tol);  // a little safety margin
      bvh.popFilter();

      if (index < 0)
        continue;

      Mesh::Vertex * vx = bvh.getElements()[(size_t)index];
      Vector3 cp = boundary_segs[i].closestPoint(vx->getPosition());

      if ((cp - edge->getEndpoint(0)->getPosition()).squaredNorm() > sqtol
       && (cp - edge->getEndpoint(1)->getPosition()).squaredNorm() > sqtol
       && (cp - vx->getPosition()).squaredNorm() < sqtol)
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
  Array< Array<Mesh::Face *> > cc;
  intx num_cc = ConnectedComponents::findEdgeConnected(mesh, cc);

  enum { NOT_VISITED, VISITED };

  for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    fi->attr().flag = NOT_VISITED;

  intx num_flipped = 0;
  for (size_t i = 0; i < cc.size(); ++i)
  {
    if (cc.empty())  // shouldn't happen, but check just in case
    {
      THEA_WARNING << mesh.getName() << ": Empty connected component";
      continue;
    }

    // Find vertex with max X
    Mesh::Vertex * xmax_vertex = nullptr;
    for (size_t j = 0; j < cc[i].size(); ++j)
      for (Mesh::Face::VertexConstIterator fvi = cc[i][j]->verticesBegin(); fvi != cc[i][j]->verticesEnd(); ++fvi)
        if (!xmax_vertex || (*fvi)->getPosition().x() > xmax_vertex->getPosition().x())
          xmax_vertex = *fvi;

    if (!xmax_vertex)  // shouldn't happen, but check just in case
    {
      THEA_WARNING << mesh.getName() << ": No vertex with maximum X found";
      continue;
    }

    // Find the "most X-facing" face incident on xmax_vertex. This face should always have a normal with positive X-component.
    Mesh::Face * xmax_face = nullptr;
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
    THEA_CONSOLE << "orient('" << mesh.getName() << "'): Flipped " << num_flipped << '/' << mesh.numFaces() << " faces in "
                 << num_cc << " connected components";
  }

  return false;
}

struct SDFOrienter
{
  SDFOrienter(Local::ShapeDiameter<Mesh> * sdf_) : sdf(sdf_) {}

  bool operator()(Mesh & mesh)
  {
    intx num_flipped = 0;
    for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      if (fi->numVertices() < 3)
        continue;

      Vector3 centroid = Vector3::Zero();
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
      THEA_CONSOLE << "orient-sdf('" << mesh.getName() << "'): Flipped " << num_flipped << '/' << mesh.numFaces() << " faces";
    }

    return false;
  }

  Local::ShapeDiameter<Mesh> * sdf;

}; // struct SDFOrienter

void
orientSDF(MG & mesh_group)
{
  Local::ShapeDiameter<Mesh> sdf(mesh_group);
  mesh_group.forEachMeshUntil(SDFOrienter(&sdf));
}

struct CountComparator
{
  bool operator()(Mesh::Face const * f0, Mesh::Face const * f1) const
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

  intx num_flipped = false;
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
                        "orient-majority('" + string(mesh.getName()) + "'): Expected > 0 neighbor consistency count");

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
    THEA_CONSOLE << "orient-majority('" << mesh.getName() << "'): Flipped " << num_flipped << '/' << mesh.numFaces()
                 << " faces";
  }

  return false;
}

struct VisibilityOrienter
{
  typedef MeshBvh<Mesh> Bvh;

  VisibilityOrienter(MG & mg, bool hi_qual_ = false)
  : hi_qual(hi_qual_)
  {
    bvh.add(mg);
    bvh.init();

    mesh_center = mg.getBounds().getCenter();
    camera_distance = 2 * mg.getBounds().getExtent().norm();
  }

  bool operator()(Mesh & mesh)
  {
    static Real PHI = (Real)((1.0 + sqrt(5.0)) / 2.0);
    static Vector3 const CAMERAS_LO[] = {
      Vector3( 1,  1,  1).normalized(),
      Vector3( 1,  1, -1).normalized(),
      Vector3( 1, -1,  1).normalized(),
      Vector3( 1, -1, -1).normalized(),
      Vector3(-1,  1,  1).normalized(),
      Vector3(-1,  1, -1).normalized(),
      Vector3(-1, -1,  1).normalized(),
      Vector3(-1, -1, -1).normalized(),

      Vector3(0,  1 / PHI,  PHI).normalized(),
      Vector3(0,  1 / PHI, -PHI).normalized(),
      Vector3(0, -1 / PHI,  PHI).normalized(),
      Vector3(0, -1 / PHI, -PHI).normalized(),

      Vector3( PHI, 0,  1 / PHI).normalized(),
      Vector3( PHI, 0, -1 / PHI).normalized(),
      Vector3(-PHI, 0,  1 / PHI).normalized(),
      Vector3(-PHI, 0, -1 / PHI).normalized(),

      Vector3( 1 / PHI,  PHI, 0).normalized(),
      Vector3(-1 / PHI,  PHI, 0).normalized(),
      Vector3( 1 / PHI, -PHI, 0).normalized(),
      Vector3(-1 / PHI, -PHI, 0).normalized(),
    };

    static Vector3 const CAMERAS_HI[] = {
      Vector3(      0,    0.5257,    0.8507),
      Vector3(      0,    0.5257,   -0.8507),
      Vector3(      0,   -0.5257,   -0.8507),
      Vector3(      0,   -0.5257,    0.8507),
      Vector3( 0.5257,    0.8507,         0),
      Vector3( 0.5257,   -0.8507,         0),
      Vector3(-0.5257,   -0.8507,         0),
      Vector3(-0.5257,    0.8507,         0),
      Vector3( 0.8507,         0,    0.5257),
      Vector3(-0.8507,         0,    0.5257),
      Vector3(-0.8507,         0,   -0.5257),
      Vector3( 0.8507,         0,   -0.5257),
      Vector3( 0.8090,    0.5000,    0.3090),
      Vector3( 0.5000,   -0.3090,    0.8090),
      Vector3(-0.5000,   -0.3090,    0.8090),
      Vector3( 0.5000,    0.3090,    0.8090),
      Vector3( 0.3090,    0.8090,   -0.5000),
      Vector3(-0.5000,    0.3090,    0.8090),
      Vector3(      0,    1.0000,         0),
      Vector3( 0.8090,   -0.5000,    0.3090),
      Vector3( 1.0000,         0,         0),
      Vector3( 0.8090,    0.5000,   -0.3090),
      Vector3( 0.3090,   -0.8090,   -0.5000),
      Vector3( 0.5000,    0.3090,   -0.8090),
      Vector3(      0,   -1.0000,         0),
      Vector3(-0.8090,   -0.5000,    0.3090),
      Vector3( 0.3090,   -0.8090,    0.5000),
      Vector3(-0.5000,    0.3090,   -0.8090),
      Vector3(      0,         0,   -1.0000),
      Vector3(-0.8090,    0.5000,    0.3090),
      Vector3(-1.0000,         0,         0),
      Vector3(-0.3090,   -0.8090,   -0.5000),
      Vector3( 0.3090,    0.8090,    0.5000),
      Vector3(      0,         0,    1.0000),
      Vector3(-0.3090,    0.8090,   -0.5000),
      Vector3(-0.3090,    0.8090,    0.5000),
      Vector3( 0.8090,   -0.5000,   -0.3090),
      Vector3( 0.5000,   -0.3090,   -0.8090),
      Vector3(-0.3090,   -0.8090,    0.5000),
      Vector3(-0.8090,    0.5000,   -0.3090),
      Vector3(-0.5000,   -0.3090,   -0.8090),
      Vector3(-0.8090,   -0.5000,   -0.3090),
    };

    static size_t const NUM_CAMERAS_LO = sizeof(CAMERAS_LO) / sizeof(Vector3);
    static size_t const NUM_CAMERAS_HI = sizeof(CAMERAS_HI) / sizeof(Vector3);

    static int const NUM_JITTERED = 5;
    static Real const JITTER_SCALE = 0.25;

    Vector3 const * cameras = hi_qual ? CAMERAS_HI : CAMERAS_LO;
    size_t num_cameras = hi_qual ? NUM_CAMERAS_HI : NUM_CAMERAS_LO;

    Array<Vector3> face_pts;
    face_pts.reserve(32);

    for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      face_pts.clear();
      {
        Vector3 c = CentroidN<Mesh::Vertex, 3>::compute(fi->verticesBegin(), fi->verticesEnd());
        face_pts.push_back(c);
        for (Mesh::Face::VertexConstIterator fvi = fi->verticesBegin(); fvi != fi->verticesEnd(); ++fvi)
          face_pts.push_back((*fvi)->getPosition());
      }

      int best_camera = -1;
      Real best_exposure = -1;
      Vector3 best_dir = Vector3::Zero();

      for (size_t j = 0; j < face_pts.size(); ++j)
      {
        Vector3 const & p = face_pts[j];

        for (size_t k = 0; k < num_cameras; ++k)
        {
#define ORIENT_VISIBILITY_RELATIVE_CAMERAS
#ifdef ORIENT_VISIBILITY_RELATIVE_CAMERAS
          Vector3 camera_pos = p + camera_distance * cameras[k];
#else
          Vector3 camera_pos = mesh_center + camera_distance * cameras[k];
#endif
          Ray3 ray(p, camera_pos - p);
          ray.setOrigin(ray.getPoint(0.0001));
          if (bvh.rayIntersects<RayIntersectionTester>(ray))
            continue;

          // Shoot a few jittered rays to test "openness" of the face w.r.t. this camera
          int num_open = 0;
          for (int h = 0; h < NUM_JITTERED; ++h)
          {
            Real jx = Random::common().uniform(-JITTER_SCALE, JITTER_SCALE);
            Real jy = Random::common().uniform(-JITTER_SCALE, JITTER_SCALE);
            Real jz = Random::common().uniform(-JITTER_SCALE, JITTER_SCALE);
            Ray3 jray(p, camera_pos + Vector3(jx, jy, jz) - p);
            jray.setOrigin(ray.getPoint(0.0001));

            if (!bvh.rayIntersects<RayIntersectionTester>(jray))
              num_open++;
          }

          Vector3 dir = ray.getDirection().normalized();
          Real exposure = num_open / (Real)NUM_JITTERED + fabs(fi->getNormal().dot(dir));
          if (best_camera < 0 || exposure > best_exposure)
          {
            best_camera = (int)k;
            best_exposure = exposure;
            best_dir = dir;
          }
        }

        if (best_camera >= 0)
          break;
      }

      if (best_camera >= 0 && fi->getNormal().dot(best_dir) < 0)
        fi->reverseWinding();
    }

    return false;
  }

  Bvh bvh;
  bool hi_qual;
  Vector3 mesh_center;
  Real camera_distance;
};

void
orientVisibility(MG & mesh_group)
{
  VisibilityOrienter func(mesh_group, orient_visibility_hi_qual);
  mesh_group.forEachMeshUntil(ref(func));
}

bool
triangulate(Mesh & mesh)
{
  intx before_num_faces = mesh.numFaces();

  Real tol = (triangulate_tolerance < 0 ? 1.0e-6f : (Real)triangulate_tolerance);
  Real epsilon = max(tol * mesh.getBounds().getExtent().norm(), tol);
  intx num_triangulated_faces = mesh.triangulate(epsilon);
  if (num_triangulated_faces < 0)
    return true;  // error, stop recursion through mesh group

  if (num_triangulated_faces == 0)  // already triangulated
    return false;

  intx after_num_faces = mesh.numFaces();

  THEA_CONSOLE << "Mesh " << mesh.getName() << ": " << num_triangulated_faces << " face(s) triangulated with "
               << (after_num_faces - before_num_faces + num_triangulated_faces) << " triangles";

  for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
  {
    if (fi->numVertices() > 3)
    {
      THEA_ERROR << "Mesh was triangulated but still has a face with " << fi->numVertices() << " vertices";
      return true;
    }
  }

  return false;
}

bool
checkProblems(Mesh & mesh)
{
  intx num_faces_with_repeated_vertices = 0;
  intx num_faces_with_repeated_edges = 0;

  intx num_edges_with_repeated_vertices = 0;
  intx num_edges_with_repeated_faces = 0;

  intx num_vertices_with_repeated_edges = 0;
  intx num_vertices_with_repeated_faces = 0;

  for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
  {
    Mesh::Face const & face = *fi;

    for (Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
      if (count(face.verticesBegin(), face.verticesEnd(), *vi) > 1)
      {
        num_faces_with_repeated_vertices++;
        break;
      }

    for (Mesh::Face::EdgeConstIterator ei = face.edgesBegin(); ei != face.edgesEnd(); ++ei)
      if (count(face.edgesBegin(), face.edgesEnd(), *ei) > 1)
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
      if (count(edge.facesBegin(), edge.facesEnd(), *fi) > 1)
      {
        num_edges_with_repeated_faces++;
        break;
      }
  }

  for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
  {
    Mesh::Vertex const & vertex = *vi;

    for (Mesh::Vertex::EdgeConstIterator ei = vertex.edgesBegin(); ei != vertex.edgesEnd(); ++ei)
      if (count(vertex.edgesBegin(), vertex.edgesEnd(), *ei) > 1)
      {
        num_vertices_with_repeated_edges++;
        break;
      }

    for (Mesh::Vertex::FaceConstIterator fi = vertex.facesBegin(); fi != vertex.facesEnd(); ++fi)
      if (count(vertex.facesBegin(), vertex.facesEnd(), *fi) > 1)
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
