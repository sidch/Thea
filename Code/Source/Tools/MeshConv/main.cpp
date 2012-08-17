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
    std::cerr << "Usage: " << argv[0] << " <infile> <outfile>" << std::endl;
    return -1;
  }

  CodecOBJ<Mesh> codec_obj(NULL, CodecOBJ<Mesh>::ReadOptions().setIgnoreNormals(true).setIgnoreTexCoords(true));
  CodecOFF<Mesh> codec_off;
  CodecOFF<Mesh> codec_off_bin(NULL, CodecOFF<Mesh>::ReadOptions(), CodecOFF<Mesh>::WriteOptions().setBinary(true));
  Codec3DS<Mesh> codec_3ds(NULL, Codec3DS<Mesh>::ReadOptions().setIgnoreTexCoords(true));

  try
  {
    MG mg("MeshGroup");

    std::string infile = toLower(argv[1]);
    if (endsWith(infile, ".obj"))
      mg.load(argv[1], codec_obj);
    else if (endsWith(infile, ".3ds"))
      mg.load(argv[1], codec_3ds);
    else
      mg.load(argv[1]);  // use whatever other fallback we may have

    std::string outfile = toLower(argv[2]);
    if (endsWith(outfile, ".obj"))
      mg.save(argv[2], codec_obj);
    else if (endsWith(outfile, ".off.bin"))
      mg.save(argv[2], codec_off_bin);
    else if (endsWith(outfile, ".off"))
      mg.save(argv[2], codec_off);
    else if (endsWith(outfile, ".3ds"))
      mg.save(argv[2], codec_3ds);
    else
      throw Error("Output format not recognized");
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
