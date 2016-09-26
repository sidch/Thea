#include "../../Common.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

struct IndexAttribute
{
  long index;

  IndexAttribute() : index(-1) {}
  void draw(RenderSystem &render_system, RenderOptions const &options) const {}
};

typedef GeneralMesh<IndexAttribute, Graphics::NullAttribute, IndexAttribute> Mesh;
typedef MeshGroup<Mesh> MG;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [options] <mesh> <out-pts>";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << " -nN   : Generate N samples [=5000]";
  THEA_CONSOLE << " -s[F] : Generate approximately uniformly separated samples";
  THEA_CONSOLE << "         with an initial oversampling factor of F";
  THEA_CONSOLE << " -v    : Generate samples at mesh vertices (ignores -n)";
  THEA_CONSOLE << " -f    : Generate samples at face centers (ignores -n)";
  THEA_CONSOLE << " -id   : Write the index of the face (or if -v, the vertex)";
  THEA_CONSOLE << "         from which each sample was drawn";

  return 0;
}

struct ReadCallback : public MeshCodec<Mesh>::ReadCallback
{
  void vertexAdded(Mesh * mesh, long index, IncrementalMeshBuilder<Mesh>::VertexHandle vertex)
  {
    vertex->attr().index = index;
  }

  void faceAdded(Mesh * mesh, long index, IncrementalMeshBuilder<Mesh>::FaceHandle face)
  {
    face->attr().index = index;
  }
};

struct VertexCollector
{
  VertexCollector(TheaArray<Vector3> * positions_, TheaArray<Vector3> * normals_, TheaArray<long> * indices_)
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

  TheaArray<Vector3> * positions;
  TheaArray<Vector3> * normals;
  TheaArray<long> * indices;
};

struct FaceCenterCollector
{
  FaceCenterCollector(TheaArray<Vector3> * positions_, TheaArray<Vector3> * normals_, TheaArray<long> * indices_)
  : positions(positions_), normals(normals_), indices(indices_) {}

  bool operator()(Mesh const & mesh)
  {
    for (Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      positions->push_back(CentroidN<Mesh::Vertex const *, 3>::compute(fi->verticesBegin(), fi->verticesEnd()));
      normals->push_back(fi->getNormal());
      indices->push_back(fi->attr().index);
    }

    return false;
  }

  TheaArray<Vector3> * positions;
  TheaArray<Vector3> * normals;
  TheaArray<long> * indices;
};

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
    ReadCallback read_callback;
    Codec3DS<Mesh>::Ptr codec_3ds(new Codec3DS<Mesh>(Codec3DS<Mesh>::ReadOptions().setIgnoreTexCoords(true)));
    CodecOBJ<Mesh>::Ptr codec_obj(new CodecOBJ<Mesh>(CodecOBJ<Mesh>::ReadOptions().setIgnoreNormals(true)
                                                                                  .setIgnoreTexCoords(true)));
    TheaArray<MeshCodec<Mesh>::Ptr> read_codecs;
    read_codecs.push_back(codec_3ds);
    read_codecs.push_back(codec_obj);

    MG mg("MeshGroup");
    mg.load(mesh_path, read_codecs, &read_callback);

    TheaArray<Vector3> positions;
    TheaArray<Vector3> normals;
    TheaArray<long> indices;

    bool do_random_sampling = true;

    if (vertex_samples)
    {
      VertexCollector collector(&positions, &normals, &indices);
      mg.forEachMeshUntil(&collector);
      do_random_sampling = false;
    }

    if (face_samples)
    {
      FaceCenterCollector collector(&positions, &normals, &indices);
      mg.forEachMeshUntil(&collector);
      do_random_sampling = false;
    }

    if (do_random_sampling)
    {
      MeshSampler<Mesh> sampler(mg);
      TheaArray< MeshSampler<Mesh>::Triangle const * > tris;

      if (uniformly_separated)
      {
        sampler.sampleEvenlyBySeparation(num_samples, positions, &normals, (output_ids ? &tris : NULL),
                                         MeshSampler<Mesh>::CountMode::EXACT, oversampling_factor, true);
      }
      else
      {
        sampler.sampleEvenlyByArea(num_samples, positions, &normals, (output_ids ? &tris : NULL),
                                   MeshSampler<Mesh>::CountMode::EXACT, true);
      }

      if (output_ids)
      {
        alwaysAssertM(tris.size() >= positions.size(), "Triangle ID's not initialized");

        indices.resize(positions.size());
        for (array_size_t i = 0; i < positions.size(); ++i)
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
    for (array_size_t i = 0; i < positions.size(); ++i)
    {
      if (std::fabs(positions[i].x()) < 1e-35) positions[i].x() = 0;
      if (std::fabs(positions[i].y()) < 1e-35) positions[i].y() = 0;
      if (std::fabs(positions[i].z()) < 1e-35) positions[i].z() = 0;

      if (std::fabs(normals[i].x()) < 1e-35) normals[i].x() = 0;
      if (std::fabs(normals[i].y()) < 1e-35) normals[i].y() = 0;
      if (std::fabs(normals[i].z()) < 1e-35) normals[i].z() = 0;
    }

    ofstream out(out_path.c_str());
    for (array_size_t i = 0; i < positions.size(); ++i)
    {
      out << positions[i].x() << ' ' << positions[i].y() << ' ' << positions[i].z() << ' '
          << normals[i].x() << ' ' << normals[i].y() << ' ' << normals[i].z();

      if (output_ids)
        out << ' ' << indices[i];

      out << '\n';
    }
    out.flush();

    THEA_CONSOLE << "Computed " << positions.size() << " samples from '" << mesh_path << '\'';
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
