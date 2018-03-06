#include "../../Common.hpp"
#include "../../Algorithms/GeodesicSphere3.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include <cstdio>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef DisplayMesh Mesh;
typedef MeshGroup<Mesh> MG;

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "Usage: " << argv[0] << " <#subdivs> <outfile>";
  return 0;
}

int
main(int argc, char * argv[])
{
  if (argc < 3)
    return usage(argc, argv);

  long num_subdivs;
  if (sscanf(argv[1], " %ld", &num_subdivs) != 1)
  {
    THEA_ERROR << "Could not parse number of subdivisions";
    return -1;
  }

  string out_path = argv[2];

  TheaArray<Vector3> vertices;
  TheaArray<long> triangles;
  if (!GeodesicSphere3::compute(num_subdivs, vertices, &triangles))
    return -1;

  Mesh::Ptr mesh(new Mesh);
  for (size_t i = 0; i < vertices.size(); ++i)
    mesh->addVertex(vertices[i]);

  for (size_t i = 0; i < triangles.size(); i += 3)
    mesh->addTriangle(triangles[i], triangles[i + 1], triangles[i + 2]);

  MG mg;
  mg.addMesh(mesh);

  try
  {
    mg.save(out_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not save geodesic sphere to '%s'", out_path.c_str())

  return 0;
}
