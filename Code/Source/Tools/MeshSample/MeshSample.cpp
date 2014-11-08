#include "../../Common.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include <cstdio>
#include <fstream>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " [options] <mesh> <out-pts>";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << " -nN   : Generate N samples [=5000]";
  THEA_CONSOLE << " -s[F] : Generate approximately uniformly separated samples";
  THEA_CONSOLE << "         with an initial oversampling factor of F";
  THEA_CONSOLE << " -v    : Generate samples only at mesh vertices (ignores -n)";
  return 0;
}

struct VertexCollector
{
  VertexCollector(TheaArray<Vector3> * positions_, TheaArray<Vector3> * normals_) : positions(positions_), normals(normals_) {}

  bool operator()(Mesh const & mesh)
  {
    for (Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    {
      positions->push_back(vi->getPosition());
      normals->push_back(vi->getNormal());
    }

    return false;
  }

  TheaArray<Vector3> * positions;
  TheaArray<Vector3> * normals;
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

  int curr_pos_arg = 0;
  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];
    if (arg.length() >= 2 && arg[0] == '-')
    {
      if (arg == "-v")
        vertex_samples = true;
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
    MG mg("MeshGroup");
    mg.load(mesh_path);

    TheaArray<Vector3> positions;
    TheaArray<Vector3> normals;

    if (vertex_samples)
    {
      VertexCollector collector(&positions, &normals);
      mg.forEachMeshUntil(&collector);
    }
    else
    {
      MeshSampler<Mesh> sampler(mg);
      if (uniformly_separated)
      {
        sampler.sampleEvenlyBySeparation(num_samples, positions, &normals, NULL, MeshSampler<Mesh>::CountMode::EXACT,
                                         oversampling_factor, true);
      }
      else
      {
        sampler.sampleEvenlyByArea(num_samples, positions, &normals, NULL, MeshSampler<Mesh>::CountMode::EXACT, true);
      }
    }

    ofstream out(out_path.c_str());
    for (size_t i = 0; i < positions.size(); ++i)
    {
      out << positions[i].x() << ' ' << positions[i].y() << ' ' << positions[i].z() << ' '
          << normals[i].x() << ' ' << normals[i].y() << ' ' << normals[i].z() << '\n';
    }
    out.flush();

    THEA_CONSOLE << "Computed " << positions.size() << " samples from '" << mesh_path << '\'';
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
