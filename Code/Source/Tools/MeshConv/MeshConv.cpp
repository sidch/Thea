#include "../../Common.hpp"
#include "../../Algorithms/ConnectedComponents.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../FilePath.hpp"
#include "../../UnorderedMap.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;

bool
splitMesh(MG::Ptr mg)
{
  typedef TheaUnorderedMap<Mesh::Vertex const *, Mesh::Vertex *> VertexMap;

  TheaArray<Mesh::Ptr> new_meshes;
  bool has_new = false;
  for (MG::MeshIterator mi = mg->meshesBegin(); mi != mg->meshesEnd(); ++mi)
  {
    TheaArray< TheaArray<Mesh::Face *> > cc;
    ConnectedComponents::findEdgeConnected(**mi, cc);

    if (cc.size() <= 1)
    {
      new_meshes.push_back(*mi);
      continue;
    }

    THEA_CONSOLE << "Splitting submesh " << (*mi)->getName() << " into " << cc.size() << " connected components";

    has_new = true;
    for (array_size_t j = 0; j < cc.size(); ++j)
    {
      Mesh::Ptr m(new Mesh(format("%s/%ld", (*mi)->getName(), (long)j)));
      VertexMap vmap;
      TheaArray<Mesh::Vertex *> new_face_vertices;
      Mesh::Vertex * new_vertex = NULL;

      for (array_size_t k = 0; k < cc[j].size(); ++k)
      {
        Mesh::Face const & face = *cc[j][k];
        new_face_vertices.clear();
        for (Mesh::Face::VertexConstIterator fvh = face.verticesBegin(); fvh != face.verticesEnd(); ++fvh)
        {
          VertexMap::const_iterator existing = vmap.find(*fvh);
          if (existing == vmap.end())
          {
            new_vertex = m->addVertex((*fvh)->getPosition(), &(*fvh)->getNormal());
            vmap[*fvh] = new_vertex;
          }
          else
            new_vertex = existing->second;

          new_face_vertices.push_back(new_vertex);
        }

        m->addFace(new_face_vertices.begin(), new_face_vertices.end());
      }

      new_meshes.push_back(m);
    }
  }

  if (has_new)
  {
    mg->clearMeshes();
    for (array_size_t i = 0; i < new_meshes.size(); ++i)
      mg->addMesh(new_meshes[i]);
  }

  for (MG::GroupIterator ci = mg->childrenBegin(); ci != mg->childrenEnd(); ++ci)
    if (!splitMesh(*ci))
      return false;

  return true;
}

int
main(int argc, char * argv[])
{
  if (argc < 3)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " [--binary] [--split] <infile> [<infile> ...] <outfile>";
    THEA_CONSOLE << "";
    THEA_CONSOLE << "Options:";
    THEA_CONSOLE << "  --binary  :  Force a binary output encoding wherever possible";
    THEA_CONSOLE << "  --split   :  Make each connected component a separate submesh";
    return -1;
  }

  CodecOBJ<Mesh>::Ptr codec_obj(new CodecOBJ<Mesh>(CodecOBJ<Mesh>::ReadOptions().setIgnoreNormals(true)
                                                                                .setIgnoreTexCoords(true)));
  Codec3DS<Mesh>::Ptr codec_3ds(new Codec3DS<Mesh>(Codec3DS<Mesh>::ReadOptions().setIgnoreTexCoords(true)));
  TheaArray<MeshCodec<Mesh>::Ptr> read_codecs;
  read_codecs.push_back(codec_obj);
  read_codecs.push_back(codec_3ds);

  CodecOFF<Mesh>::Ptr codec_off_bin(new CodecOFF<Mesh>(CodecOFF<Mesh>::ReadOptions(),
                                                       CodecOFF<Mesh>::WriteOptions().setBinary(true)));
  CodecPLY<Mesh>::Ptr codec_ply_bin(new CodecPLY<Mesh>(CodecPLY<Mesh>::ReadOptions(),
                                                       CodecPLY<Mesh>::WriteOptions().setBinary(true)));

  try
  {
    MG::Ptr main_group;
    bool force_binary = false;
    bool split = false;

    long num_mg = 0;
    for (int i = 1; i < argc - 1; ++i)
    {
      std::string arg = argv[i];
      if (arg == "--binary")
      {
        force_binary = true;
        continue;
      }
      else if (arg == "--split")
      {
        split = true;
        continue;
      }

      MG::Ptr mg(new MG(FilePath::baseName(arg)));

      std::string arg_lc = toLower(arg);
      mg->load(arg, read_codecs);

      if (!main_group)
        main_group = mg;
      else
      {
        if (num_mg == 1)  // create an aggregating node at the root level and reparent the original mesh group
        {
          MG::Ptr first_mg = main_group;
          main_group = MG::Ptr(new MG("MeshGroup"));
          main_group->addChild(first_mg);
        }

        main_group->addChild(mg);
      }

      num_mg++;
    }

    if (split)
    {
      if (!splitMesh(main_group))
        return -1;
    }

    std::string outfile = toLower(argv[argc - 1]);

    if (endsWith(outfile, ".off.bin") || (force_binary && endsWith(outfile, ".off")))
      main_group->save(argv[argc - 1], *codec_off_bin);
    else if (force_binary && endsWith(outfile, ".ply"))
      main_group->save(argv[argc - 1], *codec_ply_bin);
    else
      main_group->save(argv[argc - 1]);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
