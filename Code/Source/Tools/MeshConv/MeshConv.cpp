#include "../../Common.hpp"
#include "../../FilePath.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Graphics;

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;

int
main(int argc, char * argv[])
{
  if (argc < 3)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " [--binary] <infile> [<infile> ...] <outfile>";
    THEA_CONSOLE << "";
    THEA_CONSOLE << "Options:";
    THEA_CONSOLE << "  --binary : Force a binary output encoding wherever possible";
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

    long num_mg = 0;
    for (int i = 1; i < argc - 1; ++i)
    {
      std::string arg = argv[i];
      if (arg == "--binary")
      {
        force_binary = true;
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
