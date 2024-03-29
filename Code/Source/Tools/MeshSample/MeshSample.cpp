#include "../../Common.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/FurthestPointSampling.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../FileSystem.hpp"
#include "../../FilePath.hpp"
#include "../../UnorderedMap.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

struct IndexAttribute
{
  intx index;

  IndexAttribute() : index(-1) {}
  void draw(IRenderSystem & render_system, IRenderOptions const & options) const {}
};

typedef GeneralMesh<IndexAttribute, Graphics::NullAttribute, IndexAttribute> Mesh;
typedef MeshGroup<Mesh> MG;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [options] <mesh> <out-pts>";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << " -nN       : Generate N samples [=5000]";
  THEA_CONSOLE << " -s[F]     : Generate approximately uniformly separated samples";
  THEA_CONSOLE << "             with an initial oversampling factor of F (omit for";
  THEA_CONSOLE << "             default). F is ignored if -i is also present";
  THEA_CONSOLE << " -v        : Generate samples at mesh vertices (ignores -n)";
  THEA_CONSOLE << " -f        : Generate samples at face centers (ignores -n)";
  THEA_CONSOLE << " -id       : Write the index of the face (or if -v, the vertex)";
  THEA_CONSOLE << "             from which each sample was drawn";
  THEA_CONSOLE << " -l <file> : Load face labels from file and write sample labels";
  THEA_CONSOLE << " -i <file> : Load a set of initial point samples for the -s option";

  return 0;
}

struct ReadCallback : public MeshCodec<Mesh>::ReadCallback
{
  void vertexRead(Mesh * mesh, intx index, Mesh::VertexHandle vertex)
  {
    vertex->attr().index = index;
  }

  void faceRead(Mesh * mesh, intx index, Mesh::FaceHandle face)
  {
    face->attr().index = index;
  }
};

struct VertexCollector
{
  VertexCollector(Array<Vector3> * positions_, Array<Vector3> * normals_, Array<intx> * indices_)
  : positions(positions_), normals(normals_), indices(indices_) {}

  bool operator()(Mesh const & mesh)
  {
    for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    {
      positions->push_back(vi->getPosition());
      normals->push_back(vi->getNormal());
      indices->push_back(vi->attr().index);
    }

    return false;
  }

  Array<Vector3> * positions;
  Array<Vector3> * normals;
  Array<intx> * indices;
};

struct FaceCenterCollector
{
  FaceCenterCollector(Array<Vector3> * positions_, Array<Vector3> * normals_, Array<intx> * indices_)
  : positions(positions_), normals(normals_), indices(indices_) {}

  bool operator()(Mesh const & mesh)
  {
    for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      positions->push_back(CentroidN<Mesh::Vertex, 3>::compute(fi->verticesBegin(), fi->verticesEnd()));
      normals->push_back(fi->getNormal());
      indices->push_back(fi->attr().index);
    }

    return false;
  }

  Array<Vector3> * positions;
  Array<Vector3> * normals;
  Array<intx> * indices;
};

bool
loadLabels_Lab(string const & path, Array<string> & labels, Array<intx> & face_labels)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not open labels file: " << path;
    return false;
  }

  labels.clear();
  face_labels.clear();

  typedef map<string, intx> LabelIndexMap;
  LabelIndexMap labmap;

  string line;
  intx label_index;
  while (getline(in, line))
  {
    string label = trimWhitespace(line);
    if (label.empty())
      continue;

    LabelIndexMap::const_iterator existing = labmap.find(label);
    if (existing == labmap.end())
    {
      label_index = (intx)labmap.size();
      labmap[label] = label_index;
    }
    else
      label_index = existing->second;

    if (!getline(in, line))
    {
      THEA_ERROR << "Could not read list of faces for label '" << label << '\'';
      return false;
    }

    istringstream line_in(line);
    intx face_index;
    while (line_in >> face_index)
    {
      if (face_index < 1)
      {
        THEA_ERROR << "Face index out of bounds: " << face_index;
        return false;
      }

      if (face_index > (intx)face_labels.size())
      {
        face_labels.reserve((size_t)(2 * face_index));
        face_labels.resize((size_t)face_index, -1);
      }

      face_labels[(size_t)face_index - 1] = label_index;
    }
  }

  labels.resize((size_t)labmap.size());
  for (LabelIndexMap::const_iterator li = labmap.begin(); li != labmap.end(); ++li)
    labels[(size_t)li->second] = li->first;

  return true;
}

struct FaceToMeshMapper
{
  bool operator()(Mesh const & mesh)
  {
    for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      intx face_index = fi->attr().index;
      if (face_index >= (intx)mapping.size())
      {
        mapping.reserve((size_t)(2 * (face_index + 1)));
        mapping.resize((size_t)(face_index + 1), nullptr);
      }

      mapping[(size_t)face_index] = &mesh;
    }

    return false;
  }

  Array<Mesh const *> mapping;
};

bool
loadLabels_Labels(string const & path, MG const & mg, Array<string> & labels, Array<intx> & face_labels)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not open labels file: " << path;
    return false;
  }

  labels.clear();
  face_labels.clear();

  typedef map<string, intx> LabelIndexMap;
  LabelIndexMap labmap;

  FaceToMeshMapper f2m;
  mg.forEachMeshUntil(std::ref(f2m));
  face_labels.resize(f2m.mapping.size(), -1);

  string line;
  intx label_index;
  while (getline(in, line))
  {
    string label = trimWhitespace(line);
    if (label.empty())
      continue;

    LabelIndexMap::const_iterator existing = labmap.find(label);
    if (existing == labmap.end())
    {
      label_index = (intx)labmap.size();
      labmap[label] = label_index;
    }
    else
      label_index = existing->second;

    if (!getline(in, line))
    {
      THEA_ERROR << "Could not read representative faces for label '" << label << '\'';
      return false;
    }

    istringstream line_in(line);
    intx face_index;
    while (line_in >> face_index)
    {
      if (face_index < 0 || face_index >= (intx)face_labels.size())
      {
        THEA_ERROR << "Face index out of bounds: " << face_index;
        return false;
      }

      Mesh const * selected_mesh = f2m.mapping[(size_t)face_index];
      for (Mesh::FaceConstIterator fi = selected_mesh->facesBegin(); fi != selected_mesh->facesEnd(); ++fi)
        face_labels[(size_t)fi->attr().index] = label_index;
    }
  }

  labels.resize((size_t)labmap.size());
  for (LabelIndexMap::const_iterator li = labmap.begin(); li != labmap.end(); ++li)
    labels[(size_t)li->second] = li->first;

  return true;
}

bool
loadLabels_FaceLabels(string const & path, Array<string> & labels, Array<intx> & face_labels)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not open labels file: " << path;
    return false;
  }

  labels.clear();
  face_labels.clear();

  typedef map<string, intx> LabelIndexMap;
  LabelIndexMap labmap;

  string line;
  intx label_index;
  while (getline(in, line))
  {
    string label = trimWhitespace(line);
    LabelIndexMap::const_iterator existing = labmap.find(label);
    if (existing == labmap.end())
    {
      label_index = (intx)labmap.size();
      labmap[label] = label_index;
    }
    else
      label_index = existing->second;

    face_labels.push_back(label_index);
  }

  labels.resize((size_t)labmap.size());
  for (LabelIndexMap::const_iterator li = labmap.begin(); li != labmap.end(); ++li)
    labels[(size_t)li->second] = li->first;

  return true;
}

bool
loadLabels(string const & path, MG const & mg, Array<string> & labels, Array<intx> & face_labels)
{
  string ext = toLower(FilePath::extension(path));
  if (ext == "lab")
    return loadLabels_Lab(path, labels, face_labels);
  else if (ext == "labels")
    return loadLabels_Labels(path, mg, labels, face_labels);
  else
    return loadLabels_FaceLabels(path, labels, face_labels);
}

bool
loadSamples(string const & path, Array<Vector3> & positions, Array<string> & lines)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Couldn't open pre-sampled points file '" << path << '\'';
    return false;
  }

  positions.clear();
  lines.clear();

  string line;
  Vector3 p, n;
  while (getline(in, line))
  {
    line = trimWhitespace(line);
    if (line.empty() || line[0] == '#')
      continue;

    istringstream line_in(line);
    if (!(line_in >> p[0] >> p[1] >> p[2]))
    {
      THEA_ERROR << "Couldn't read pre-sampled point from line '" << line << '\'';
      return false;
    }

    positions.push_back(p);
    lines.push_back(line);
  }

  THEA_CONSOLE << "Read " << positions.size() << " pre-sampled point(s) from '" << path << '\'';

  return true;
}

int
main(int argc, char * argv[])
{
  if (argc < 3)
    return usage(argc, argv);

  string mesh_path;
  string out_path;
  long num_samples = 5000;
  bool uniformly_separated = false;
  float oversampling_factor = -1;
  bool vertex_samples = false;
  bool face_samples = false;
  bool output_ids = false;
  string labels_path;
  string presampled_path;

  int curr_pos_arg = 0;
  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];
    if (arg.length() >= 2 && arg[0] == '-')
    {
      if (arg == "-v")
        vertex_samples = true;
      else if (arg == "-f")
        face_samples = true;
      else if (arg == "-id")
        output_ids = true;
      else if (beginsWith(arg, "-s"))
      {
        uniformly_separated = true;

        if (arg.length() > 2)
        {
          if (sscanf(arg.substr(2).c_str(), "%f", &oversampling_factor) != 1)
          {
            THEA_ERROR << "Invalid oversampling factor";
            return -1;
          }
        }
      }
      else if (beginsWith(arg, "-n"))
      {
        if (arg.length() <= 2)
          return usage(argc, argv);

        if (sscanf(arg.substr(2).c_str(), "%ld", &num_samples) != 1 || num_samples < 0)
        {
          THEA_ERROR << "Invalid number of samples";
          return -1;
        }
      }
      else if (arg == "-l")
      {
        if (i >= argc - 1)
          return usage(argc, argv);

        labels_path = argv[++i];
        if (!FileSystem::fileExists(labels_path))
        {
          THEA_ERROR << "Labels file '" << labels_path << "' does not exist";
          return -1;
        }
      }
      else if (arg == "-i")
      {
        if (i >= argc - 1)
          return usage(argc, argv);

        presampled_path = argv[++i];
        if (!FileSystem::fileExists(presampled_path))
        {
          THEA_ERROR << "Pre-sampled points file '" << presampled_path << "' does not exist";
          return -1;
        }
      }
      else
        return usage(argc, argv);
    }
    else
    {
      switch (curr_pos_arg)
      {
        case 0: mesh_path = arg; break;
        case 1: out_path = arg; break;
        default: return usage(argc, argv);
      }

      curr_pos_arg++;
    }
  }

  if (curr_pos_arg < 2)
    return usage(argc, argv);

  try
  {
    if (!presampled_path.empty())
    {
      //=======================================================================================================================
      // Subsampling a previously sampled set of points
      //=======================================================================================================================

      if (!uniformly_separated)
      {
        THEA_ERROR << "Pre-sampled points require the -s option as well";
        return -1;
      }

      Array<Vector3> orig_pos;
      Array<string> orig_lines;
      if (!loadSamples(presampled_path, orig_pos, orig_lines))
        return -1;

      Array<intx> selected;
      if (!orig_pos.empty())
      {
        selected.resize((size_t)num_samples);
        if (FurthestPointSampling::subsample((intx)orig_pos.size(), &orig_pos[0], num_samples, &selected[0],
                                             DistanceType::GEODESIC, /* verbose = */ true) < num_samples)
          return -1;
      }

      ofstream out(out_path.c_str());
      if (!out)
      {
        THEA_ERROR << "Could not open output file '" << out_path << "' for writing";
        return -1;
      }

      for (size_t i = 0; i < selected.size(); ++i)
        out << orig_lines[(size_t)selected[i]] << '\n';
    }
    else
    {
      //=======================================================================================================================
      // Sampling points from scratch
      //=======================================================================================================================

      ReadCallback read_callback;
      Codec3ds<Mesh>::Ptr codec_3ds(new Codec3ds<Mesh>(Codec3ds<Mesh>::ReadOptions().setReadTexCoords(false)));
      CodecObj<Mesh>::Ptr codec_obj(new CodecObj<Mesh>(CodecObj<Mesh>::ReadOptions().setReadNormals(false)
                                                                                    .setReadTexCoords(false)));
      Array<MeshCodec<Mesh>::Ptr> read_codecs;
      read_codecs.push_back(codec_3ds);
      read_codecs.push_back(codec_obj);

      MG mg("MeshGroup");
      mg.load(mesh_path, read_codecs, &read_callback);

      Array<string> labels;
      Array<intx> face_labels;
      bool output_labels = (!vertex_samples && !labels_path.empty());
      if (output_labels)
      {
        if (!loadLabels(labels_path, mg, labels, face_labels))
          return -1;
      }

      Array<Vector3> positions;
      Array<Vector3> normals;
      Array<intx> indices;

      bool need_face_ids = (output_ids || output_labels);
      bool do_random_sampling = true;

      if (vertex_samples)
      {
        VertexCollector collector(&positions, &normals, &indices);
        mg.forEachMeshUntil(std::ref(collector));
        do_random_sampling = false;
      }

      if (face_samples)
      {
        FaceCenterCollector collector(&positions, &normals, &indices);
        mg.forEachMeshUntil(std::ref(collector));
        do_random_sampling = false;
      }

      if (do_random_sampling)
      {
        MeshSampler<Mesh> sampler(mg);
        Array< MeshSampler<Mesh>::Triangle const * > tris;

        if (uniformly_separated)
        {
          sampler.sampleEvenlyBySeparation(num_samples, positions, &normals, (need_face_ids ? &tris : nullptr),
                                           MeshSampler<Mesh>::CountMode::EXACT, oversampling_factor, true);
        }
        else
        {
          sampler.sampleEvenlyByArea(num_samples, positions, &normals, (need_face_ids ? &tris : nullptr),
                                     MeshSampler<Mesh>::CountMode::EXACT, true);
        }

        if (need_face_ids)
        {
          alwaysAssertM(tris.size() >= positions.size(), "Triangle IDs not initialized");

          indices.resize(positions.size());
          for (size_t i = 0; i < positions.size(); ++i)
          {
            Mesh::Face const * face = tris[i]->getVertices().getMeshFace();
            indices[i] = face->attr().index;
          }
        }
      }

      // Sanitize outputs, sometimes this is a problem with values smaller than 32-bit float precision like 4.01752e-42 getting
      // written to the output stream. This seems to be a GCC bug (feature?) where the optimizer is allowed to use higher
      // precision registers for floating point numbers unless -ffloat-store is enabled
      //
      // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=10644
      // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=323
      //
      for (size_t i = 0; i < positions.size(); ++i)
      {
        if (std::fabs(positions[i].x()) < 1e-35) positions[i].x() = 0;
        if (std::fabs(positions[i].y()) < 1e-35) positions[i].y() = 0;
        if (std::fabs(positions[i].z()) < 1e-35) positions[i].z() = 0;

        if (std::fabs(normals[i].x()) < 1e-35) normals[i].x() = 0;
        if (std::fabs(normals[i].y()) < 1e-35) normals[i].y() = 0;
        if (std::fabs(normals[i].z()) < 1e-35) normals[i].z() = 0;
      }

      ofstream out(out_path.c_str());
      if (!out)
      {
        THEA_ERROR << "Could not open output file '" << out_path << "' for writing";
        return -1;
      }

      for (size_t i = 0; i < positions.size(); ++i)
      {
        out << positions[i].x() << ' ' << positions[i].y() << ' ' << positions[i].z() << ' '
            << normals[i].x() << ' ' << normals[i].y() << ' ' << normals[i].z();

        if (output_ids)
          out << ' ' << indices[i];

        if (output_labels)
        {
          intx label_index = (indices[i] < 0 || indices[i] >= (intx)face_labels.size()
                            ? -1 : face_labels[(size_t)indices[i]]);
          out << " \"" << (label_index < 0 ? "" : labels[(size_t)label_index]) << '"';
        }

        out << '\n';
      }
      out.flush();

      THEA_CONSOLE << "Computed " << positions.size() << " samples from '" << mesh_path << '\'';
    }
  }
  THEA_CATCH(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
