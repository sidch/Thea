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
main(int argc, char * argv[])
{
  if (argc < 4)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " <mesh> <approx-num-samples> <out-pts>";
    return 0;
  }

  string mesh_path = argv[1];
  string nsamples_str = argv[2];
  string out_path = argv[3];

  long approx_num_samples = 0;
  if (sscanf(nsamples_str.c_str(), "%ld", &approx_num_samples) != 1 || approx_num_samples < 0)
  {
    THEA_ERROR << "Invalid number of samples";
    return -1;
  }

  try
  {
    MG mg("MeshGroup");
    mg.load(mesh_path);

    TheaArray<Vector3> positions;
    TheaArray<Vector3> face_normals;

    MeshSampler<Mesh> sampler(mg);
    sampler.sampleEvenlyByArea(approx_num_samples, positions, &face_normals);

    ofstream out(out_path.c_str());
    for (size_t i = 0; i < positions.size(); ++i)
    {
      out << positions[i].x() << ' ' << positions[i].y() << ' ' << positions[i].z() << ' '
          << face_normals[i].x() << ' ' << face_normals[i].y() << ' ' << face_normals[i].z() << '\n';
    }
    out.flush();

    THEA_CONSOLE << "Computed " << positions.size() << " samples from '" << mesh_path << '\'';
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
