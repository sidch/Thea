#include "../Algorithms/SurfaceParametrization.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../FileSystem.hpp"
#include "../FilePath.hpp"
#include "../Math.hpp"

// #define DCEL_MESH
#ifdef DCEL_MESH
#  include "../Graphics/DcelMesh.hpp"
#else
#  include "../Graphics/GeneralMesh.hpp"
#endif

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

#ifdef DCEL_MESH
  typedef DcelMesh<> Mesh;
#else
  typedef GeneralMesh<> Mesh;
#endif

int
main(int argc, char * argv[])
{
  try
  {
    string in_path;
    if (argc > 1)
      in_path = argv[1];
    else
    {
      string app_dir = FilePath::parent(FileSystem::resolve(argv[0]));
      in_path = FilePath::concat(app_dir, "../../../../Data/Models/cube_open.obj");
    }

    SurfaceParametrization::WeightType weight_type = SurfaceParametrization::WeightType::MEAN_VALUE;
    if (argc > 2)
    {
      string w_str = toLower(argv[2]);
      if (w_str == "uniform")
        weight_type = SurfaceParametrization::WeightType::UNIFORM;
      else if (w_str == "cotangent")
        weight_type = SurfaceParametrization::WeightType::COTANGENT;
      else if (w_str == "mean_value")
        weight_type = SurfaceParametrization::WeightType::MEAN_VALUE;
      else
      {
        THEA_ERROR << "Unsupported weight type: " << argv[2];
        return -1;
      }
    }

    MeshGroup<Mesh> mg;
    mg.load(in_path);
    auto mesh = *mg.meshesBegin();

#ifndef DCEL_MESH
    mesh->triangulate();
#endif

    Mesh::Vertex const * vx = nullptr;
    while (true)
    {
      for (auto vi = mesh->verticesBegin(); vi != mesh->verticesEnd(); ++vi)
        if (vi->isBoundaryVertex())
        {
          vx = &(*vi);
          break;
        }

      if (vx || mesh->numFaces() <= 1)
            break;
      else
      {
#ifdef DCEL_MESH
        THEA_ERROR << "Mesh has no open boundary";
        return -1;
#else
        mesh->removeFace(&(*mesh->facesBegin()));
#endif
      }
    }

    Array<Mesh::Vertex const *> boundary_loop;
    UnorderedSet<Mesh::Vertex const *> visited;
    while (vx)
    {
      boundary_loop.push_back(vx);
      visited.insert(vx);

      // THEA_CONSOLE << "Vertex = " << toString(vx->getPosition());

      Mesh::Vertex const * next = nullptr;
      for (auto vei = vx->edgesBegin(); vei != vx->edgesEnd(); ++vei)
      {
        // THEA_CONSOLE << "Edge " << toString((*vei)->getEndpoint(0)->getPosition()) << " - "
        //                         << toString((*vei)->getEndpoint(1)->getPosition())
        //                         << ": boundary = " << (*vei)->isBoundaryEdge();

        if (!(*vei)->isBoundaryEdge())
          continue;

        auto other = (*vei)->getOtherEndpoint(vx);
        if (visited.find(other) == visited.end())
        {
          next = other;
          break;
        }
      }

      vx = next;
    }

    THEA_CONSOLE << "Detected boundary loop with " << boundary_loop.size() << " vertices";

    Array<Real> edge_lengths(boundary_loop.size());
    Real total_len = 0;
    for (size_t i = 0; i < boundary_loop.size(); ++i)
    {
      size_t j = (i + 1) % boundary_loop.size();
      edge_lengths[i] = (boundary_loop[j]->getPosition() - boundary_loop[i]->getPosition()).norm();
      total_len += edge_lengths[i];
    }

    typedef UnorderedMap<Mesh::Vertex const *, Vector2> VertexParameterMap;
    VertexParameterMap vertex_params;
    Real angle = 0;
    for (size_t i = 0; i < boundary_loop.size(); ++i)
    {
      vertex_params[boundary_loop[i]] = Vector2(std::cos(angle), std::sin(angle));
      angle += Math::twoPi() * edge_lengths[i] / total_len;
    }

    THEA_CONSOLE << "Mapped boundary to circle";

    if (!SurfaceParametrization::parametrize(mesh->facesBegin(), mesh->facesEnd(), weight_type, vertex_params, vertex_params))
      return -1;

    for (auto vi = mesh->verticesBegin(); vi != mesh->verticesEnd(); ++vi)
    {
      auto v = &(*vi);
      auto params = vertex_params.find(v)->second;
      v->setPosition(Vector3(params[0], params[1], 0));
    }

    mg.save("flattened.obj");
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
