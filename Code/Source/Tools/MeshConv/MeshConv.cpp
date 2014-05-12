#include "../../Common.hpp"
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
    std::cerr << "Usage: " << argv[0] << " <infile> [<infile> ...] <outfile>" << std::endl;
    return -1;
  }

  CodecOBJ<Mesh> codec_obj(NULL, CodecOBJ<Mesh>::ReadOptions().setIgnoreNormals(true).setIgnoreTexCoords(true));
  CodecOFF<Mesh> codec_off;
  CodecOFF<Mesh> codec_off_bin(NULL, CodecOFF<Mesh>::ReadOptions(), CodecOFF<Mesh>::WriteOptions().setBinary(true));
  Codec3DS<Mesh> codec_3ds(NULL, Codec3DS<Mesh>::ReadOptions().setIgnoreTexCoords(true));

  try
  {
    MG::Ptr main_group;

    for (int i = 1; i < argc - 1; ++i)
    {
      MG::Ptr mg(new MG("MeshGroup"));

      std::string infile = toLower(argv[i]);
      if (endsWith(infile, ".obj"))
        mg->load(argv[i], codec_obj);
      else if (endsWith(infile, ".3ds"))
        mg->load(argv[i], codec_3ds);
      else
        mg->load(argv[i]);  // use whatever other fallback we may have

      if (argc == 3)
        main_group = mg;
      else
      {
        if (!main_group)
          main_group = MG::Ptr(new MG("MeshGroup"));

        main_group->addChild(mg);
      }
    }

    std::string outfile = toLower(argv[argc - 1]);
    if (endsWith(outfile, ".obj"))
      main_group->save(argv[argc - 1], codec_obj);
    else if (endsWith(outfile, ".off.bin"))
      main_group->save(argv[argc - 1], codec_off_bin);
    else if (endsWith(outfile, ".off"))
      main_group->save(argv[argc - 1], codec_off);
    else if (endsWith(outfile, ".3ds"))
      main_group->save(argv[argc - 1], codec_3ds);
    else
      throw Error("Output format not recognized");
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
