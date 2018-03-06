#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/Local/ShapeDiameter.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/ConnectedComponents.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Algorithms/RayIntersectionTester.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Graphics/VertexWelder.hpp"
#include "../../FilePath.hpp"
#include "../../LineSegment3.hpp"
#include "../../Math.hpp"
#include "../../Random.hpp"
#include "../../UnorderedSet.hpp"
#include "../../UnorderedMap.hpp"
#include <boost/program_options.hpp>
#include <algorithm>
#include <fstream>
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

typedef TheaArray<string> LabelArray;
typedef TheaUnorderedMap<string, int> LabelIndexMap;

struct FaceAttribute
{
  int label_index;
  int flag;

  FaceAttribute() : label_index(-1), flag(0) {}

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

  void faceRead(Mesh * mesh, long index, Mesh::FaceHandle face)
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

  void faceWritten(Mesh const * mesh, long index, Mesh::FaceConstHandle face)
  {
    long lab = face->attr().label_index;
    alwaysAssertM(lab < (long)labels.size(), "Face label index out of bounds");

    if (lab < 0)
      unlabeled_faces.push_back(index);
    else
      faces_per_label[(size_t)lab].push_back(index);
  }

  LabelArray const & labels;
  TheaArray< TheaArray<long> > faces_per_label;
  TheaArray<long> unlabeled_faces;

}; // struct WriteCallback

int
meshFix(int argc, char * argv[])
{
  int parse_status = parseArgs(argc, argv);
  if (parse_status <= 0)
    return parse_status;

  Codec3DS<Mesh>::ReadOptions read_opts_3ds = Codec3DS<Mesh>::ReadOptions::defaults();
  read_opts_3ds.setIgnoreTexCoords(no_texcoords)
               .setSkipEmptyMeshes(no_empty);
  Codec3DS<Mesh>::Ptr codec_3ds(new Codec3DS<Mesh>(read_opts_3ds, Codec3DS<Mesh>::WriteOptions::defaults()));

  CodecOBJ<Mesh>::ReadOptions read_opts_obj = CodecOBJ<Mesh>::ReadOptions::defaults();
  read_opts_obj.setIgnoreTexCoords(no_texcoords)
               .setIgnoreNormals(no_normals)
               .setSkipEmptyMeshes(no_empty);
  CodecOBJ<Mesh>::WriteOptions write_opts_obj = CodecOBJ<Mesh>::WriteOptions::defaults();
  write_opts_obj.setIgnoreTexCoords(no_texcoords)
                .setIgnoreNormals(no_normals)
                .setSkipEmptyMeshes(no_empty);
  CodecOBJ<Mesh>::Ptr codec_obj(new CodecOBJ<Mesh>(read_opts_obj, write_opts_obj));

  CodecOFF<Mesh>::ReadOptions opts_off = CodecOFF<Mesh>::ReadOptions::defaults();
  opts_off.setSkipEmptyMeshes(no_empty);
  CodecOFF<Mesh>::Ptr codec_off(new CodecOFF<Mesh>(opts_off));

  CodecPLY<Mesh>::ReadOptions opts_ply = CodecPLY<Mesh>::ReadOptions::defaults();
  opts_ply.setSkipEmptyMeshes(no_empty);
  CodecPLY<Mesh>::Ptr codec_ply(new CodecPLY<Mesh>(opts_ply));

  TheaArray<MeshCodec<Mesh>::Ptr> codecs;
  codecs.push_back(codec_3ds);
  codecs.push_back(codec_obj);
  codecs.push_back(codec_off);
  codecs.push_back(codec_3ds);

  LabelArray labels;
  ReadCallback read_callback(labels);

  MG mg(FilePath::objectName(infile));
  mg.load(infile, codecs, &read_callback);

  mg.updateBounds();  // load() should do this, but let's play safe
  if (mg.isEmpty())
  {
    THEA_ERROR << "Input mesh is empty, no output written";
    return -1;
  }

  // The order of operations is important here -- e.g. it's better to do orient() after welding vertices. flatten() must happen
  // right at the beginning.

  if (do_flatten)
  {
    flatten(mg);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_del_danglers)
  {
    mg.forEachMeshUntil(&delDanglers);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_del_dup_faces)
  {
    delDuplicateFaces(mg, del_dup_faces_sorted);
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

  if (do_orient_visibility)
  {
    orientVisibility(mg);
    mg.forEachMeshUntil(&checkProblems);
  }

  if (do_triangulate)
  {
    if (mg.forEachMeshUntil(&triangulate)) return -1;
    mg.forEachMeshUntil(&checkProblems);
  }

  string lc_out = toLower(outfile);
  if (endsWith(lc_out, ".off.bin"))
  {
    CodecOFF<Mesh>::WriteOptions write_opts_off = CodecOFF<Mesh>::WriteOptions().setBinary(true);
    *codec_off = CodecOFF<Mesh>(opts_off, write_opts_off);
  }

  WriteCallback write_callback(labels);
  mg.save(outfile, codecs, (export_face_labels ? &write_callback : NULL));

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
  namespace po = boost::program_options;

  static std::string const usage("\nUsage: " + string(argv[0]) + " [options] <infile> <outfile>\n"
                                   "   (all tolerances are specified as a fraction of the sub-mesh bounding box diagonal)\n");

  string s_orient_visibility_qual;

  po::options_description visible("Allowed options");
  visible.add_options()
          ("help,h",              "Print this help message")
          ("version,v",           "Print the program version")
          ("verbose",             "Print extra status information")
          ("infile",              po::value<string>(&infile), "Path to input mesh file")
          ("outfile",             po::value<string>(&outfile), "Path to output mesh file")

          ("no-normals",          "Ignore vertex and face normals in the input, and don't write any to output")

          ("no-texcoords",        "Ignore vertex and face texture coordinates in the input, and don't write any to output")

          ("no-empty",            "Ignore empty sub-meshes")

          ("flatten",             "Flatten the object to a single sub-mesh when reading")

          ("del-danglers",        "Delete isolated vertices and edges, and empty faces")

          ("del-dup-faces",       "Delete duplicate faces, regardless of orientation")

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

          ("orient-sdf",          "Consistently orient each face of the mesh so that the normal defined by counter-clockwise"
                                  " winding points inside-out, using the direction of smaller shape diameter as a heuristic")

          ("orient-majority",     "Consistently orient each edge-connected component of the mesh so that normals defined by"
                                  " counter-clockwise winding point inside-out, flipping faces that disagree most with their"
                                  " neighbors first")

          ("orient-visibility",   po::value<string>(&s_orient_visibility_qual)->implicit_value(""),
                                  "Consistently orient each face of the mesh so that the normal defined by counter-clockwise"
                                  " winding points towards the direction of an external camera from which the face is visible")

          ("triangulate",         po::value<double>(&triangulate_tolerance),
                                  "Triangulate every face with more than 3 vertices")

          ("export-face-labels",  "Along with the output mesh, export a .lab file that associates every face with the name of"
                                  " the submesh it belonged to in the input. Works even if --flatten is specified.")
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
  if (!(endsWith(lc_out, ".3ds") || endsWith(lc_out, ".obj") || endsWith(lc_out, ".off") || endsWith(lc_out, ".off.bin")))
  {
    THEA_ERROR << "Unrecognized mesh output format";
    return -1;
  }

  bool do_something = false;

  if (vm.count("no-normals") > 0)
    no_normals = do_something = true;

  if (vm.count("no-texcoords") > 0)
    no_texcoords = do_something = true;

  if (vm.count("no-empty") > 0)
    no_empty = do_something = true;

  if (vm.count("flatten") > 0)
    do_flatten = do_something = true;

  if (vm.count("del-danglers") > 0)
    do_del_danglers = do_something = true;

  if (vm.count("del-dup-faces") > 0)
    do_del_dup_faces = do_something = true;

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

  if (vm.count("orient-visibility") > 0)
  {
    do_orient_visibility = do_something = true;
    string s = toLower(s_orient_visibility_qual);
    if (s == "hi" || s == "high")
      orient_visibility_hi_qual = true;
    else if (s.empty() || s == "lo" || s == "low")
      orient_visibility_hi_qual = false;
    else
    {
      THEA_ERROR << "Orient by visibility quality setting not recognized: " << s_orient_visibility_qual;
      return -1;
    }
  }

  if (vm.count("triangulate") > 0)
    do_triangulate = do_something = true;

  if (vm.count("export-face-labels") > 0)
    export_face_labels = do_something = true;

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

struct Flattener
{
  Flattener(string const & mesh_name = "Flattened") { flattened = Mesh::Ptr(new Mesh(mesh_name)); }

  bool operator()(Mesh const & mesh)
  {
    typedef TheaUnorderedMap<Mesh::Vertex const *, Mesh::Vertex *> VertexMap;
    VertexMap vmap;

    for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    {
      Mesh::Vertex * new_vertex = flattened->addVertex(vi->getPosition(), -1, &vi->getNormal());
      if (!new_vertex)
        throw Error("Could not add copy of vertex to flattened mesh");

      new_vertex->attr() = vi->attr();
      vmap[&(*vi)] = new_vertex;
    }

    TheaArray<Mesh::Vertex *> face_vertices;
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
  mesh_group.forEachMeshUntil(&flattener);
  mesh_group.clear();
  mesh_group.addMesh(flattener.flattened);
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

struct FaceSeq
{
  TheaArray<Mesh::Vertex const *> seq;

  FaceSeq(Mesh::Face const & face, bool sorted)
  {
    for (Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
      seq.push_back(*vi);

    if (sorted)
      std::sort(seq.begin(), seq.end());
  }
};

bool
operator==(FaceSeq const & lhs, FaceSeq const & rhs)
{
  return lhs.seq.size() == rhs.seq.size()
      && std::equal(lhs.seq.begin(), lhs.seq.end(), rhs.seq.begin());
}

std::size_t
hash_value(FaceSeq const & f)
{
  return boost::hash_range(f.seq.begin(), f.seq.end());
}

typedef TheaUnorderedSet<FaceSeq> FaceSet;

struct DupFaceDeleter
{
  bool sorted;

  DupFaceDeleter(bool sorted_) : sorted(sorted_) {}

  bool operator()(Mesh & mesh)
  {
    FaceSet seqs;
    long num_deleted = 0;

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
  DupFaceDeleter func(sorted);
  mesh_group.forEachMeshUntil(&func);
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
    for (size_t i = 0; i < boundary_segs.size(); ++i)
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
  for (size_t i = 0; i < cc.size(); ++i)
  {
    if (cc.empty())  // shouldn't happen, but check just in case
    {
      THEA_WARNING << mesh.getName() << ": Empty connected component";
      continue;
    }

    // Find vertex with max X
    Mesh::Vertex * xmax_vertex = NULL;
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
  SDFOrienter func(&sdf);
  mesh_group.forEachMeshUntil(&func);
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
    THEA_CONSOLE << "orient-majority('" << mesh.getName() << "'): Flipped " << num_flipped << '/' << mesh.numFaces() << " faces";
  }

  return false;
}

struct VisibilityOrienter
{
  typedef MeshKDTree<Mesh> KDTree;

  VisibilityOrienter(MG & mg, bool hi_qual_ = false)
  : hi_qual(hi_qual_)
  {
    kdtree.add(mg);
    kdtree.init();

    mesh_center = mg.getBounds().getCenter();
    camera_distance = 2 * mg.getBounds().getExtent().length();
  }

  bool operator()(Mesh & mesh)
  {
    static Real PHI = (Real)((1.0 + std::sqrt(5.0)) / 2.0);
    static Vector3 const CAMERAS_LO[] = {
      Vector3( 1,  1,  1).unit(),
      Vector3( 1,  1, -1).unit(),
      Vector3( 1, -1,  1).unit(),
      Vector3( 1, -1, -1).unit(),
      Vector3(-1,  1,  1).unit(),
      Vector3(-1,  1, -1).unit(),
      Vector3(-1, -1,  1).unit(),
      Vector3(-1, -1, -1).unit(),

      Vector3(0,  1 / PHI,  PHI).unit(),
      Vector3(0,  1 / PHI, -PHI).unit(),
      Vector3(0, -1 / PHI,  PHI).unit(),
      Vector3(0, -1 / PHI, -PHI).unit(),

      Vector3( PHI, 0,  1 / PHI).unit(),
      Vector3( PHI, 0, -1 / PHI).unit(),
      Vector3(-PHI, 0,  1 / PHI).unit(),
      Vector3(-PHI, 0, -1 / PHI).unit(),

      Vector3( 1 / PHI,  PHI, 0).unit(),
      Vector3(-1 / PHI,  PHI, 0).unit(),
      Vector3( 1 / PHI, -PHI, 0).unit(),
      Vector3(-1 / PHI, -PHI, 0).unit(),
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

    TheaArray<Vector3> face_pts;
    face_pts.reserve(32);

    for (Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      face_pts.clear();
      {
        Vector3 c = CentroidN<Mesh::Vertex const *, 3>::compute(fi->verticesBegin(), fi->verticesEnd());
        face_pts.push_back(c);
        for (Mesh::Face::VertexConstIterator fvi = fi->verticesBegin(); fvi != fi->verticesEnd(); ++fvi)
          face_pts.push_back((*fvi)->getPosition());
      }

      int best_camera = -1;
      Real best_exposure = -1;
      Vector3 best_dir = Vector3::zero();

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
          if (kdtree.rayIntersects<RayIntersectionTester>(ray))
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

            if (!kdtree.rayIntersects<RayIntersectionTester>(jray))
              num_open++;
          }

          Vector3 dir = ray.getDirection().unit();
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

  KDTree kdtree;
  bool hi_qual;
  Vector3 mesh_center;
  Real camera_distance;
};

void
orientVisibility(MG & mesh_group)
{
  VisibilityOrienter func(mesh_group, orient_visibility_hi_qual);
  mesh_group.forEachMeshUntil(&func);
}

bool
triangulate(Mesh & mesh)
{
  long before_num_faces = mesh.numFaces();

  Real tol = (triangulate_tolerance < 0 ? 1.0e-6f : (Real)triangulate_tolerance);
  Real epsilon = max(tol * mesh.getBounds().getExtent().length(), tol);
  long num_triangulated_faces = mesh.triangulate(epsilon);
  if (num_triangulated_faces < 0)
    return true;  // error, stop recursion through mesh group

  if (num_triangulated_faces == 0)  // already triangulated
    return false;

  long after_num_faces = mesh.numFaces();

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
