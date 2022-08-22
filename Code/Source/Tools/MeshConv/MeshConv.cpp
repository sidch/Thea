#include "../../Common.hpp"
#include "../../Algorithms/ConnectedComponents.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../AffineTransform3.hpp"
#include "../../Array.hpp"
#include "../../FilePath.hpp"
#include "../../UnorderedMap.hpp"
#include <cstdlib>
#include <fstream>
#include <functional>
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
  typedef UnorderedMap<Mesh::Vertex const *, Mesh::Vertex *> VertexMap;

  Array<Mesh::Ptr> new_meshes;
  bool has_new = false;
  for (MG::MeshIterator mi = mg->meshesBegin(); mi != mg->meshesEnd(); ++mi)
  {
    Array< Array<Mesh::Face *> > cc;
    ConnectedComponents::findEdgeConnected(**mi, cc);

    if (cc.size() <= 1)
    {
      new_meshes.push_back(*mi);
      continue;
    }

    THEA_CONSOLE << "Splitting submesh " << (*mi)->getName() << " into " << cc.size() << " connected components";

    has_new = true;
    for (size_t j = 0; j < cc.size(); ++j)
    {
      Mesh::Ptr m(new Mesh(format("%s/%ld", (*mi)->getName(), (intx)j)));
      VertexMap vmap;
      Array<Mesh::Vertex *> new_face_vertices;
      Mesh::Vertex * new_vertex = nullptr;

      for (size_t k = 0; k < cc[j].size(); ++k)
      {
        Mesh::Face const & face = *cc[j][k];
        new_face_vertices.clear();
        for (Mesh::Face::VertexConstIterator fvh = face.verticesBegin(); fvh != face.verticesEnd(); ++fvh)
        {
          VertexMap::const_iterator existing = vmap.find(*fvh);
          if (existing == vmap.end())
          {
            new_vertex = m->addVertex((*fvh)->getPosition(), -1, &(*fvh)->getNormal());
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
    for (size_t i = 0; i < new_meshes.size(); ++i)
      mg->addMesh(new_meshes[i]);
  }

  for (MG::GroupIterator ci = mg->childrenBegin(); ci != mg->childrenEnd(); ++ci)
    if (!splitMesh(*ci))
      return false;

  return true;
}

MG::Ptr
extractFaces(MG::ConstPtr mg, Array<bool> const & selection)
{
  auto extracted = make_shared<MG>(mg->getName());
  Array<Mesh::Face const *> selected_faces;
  for (auto mi = mg->meshesBegin(); mi != mg->meshesEnd(); ++mi)
  {
    selected_faces.clear();

    auto const & mesh = **mi;
    for (auto fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
    {
      intx index = fi->getIndex();
      if (index >= 0 && (size_t)index < selection.size() && selection[(size_t)index])
        selected_faces.push_back(&(*fi));
    }

    if (!selected_faces.empty())
    {
      auto extracted_mesh = make_shared<Mesh>(mesh.getName());
      mesh.extractFaces(selected_faces.begin(), selected_faces.end(), *extracted_mesh);
      extracted->addMesh(extracted_mesh);
    }
  }

  for (auto ci = mg->childrenBegin(); ci != mg->childrenEnd(); ++ci)
  {
    auto ec = extractFaces(*ci, selection);
    if (ec && !ec->empty())
      extracted->addChild(ec);
  }

  return extracted;
}

MG::Ptr
extractFaces(MG::ConstPtr mg, string labels_path, string label)
{
  Array<bool> selection;
  ifstream in(labels_path);
  if (!in) { throw Error("Could not open labels file: " + labels_path); }

  string line;
  while (getline(in, line))
    selection.push_back(trimWhitespace(line) == label);

  return extractFaces(mg, selection);
}

struct Transformer
{
  Transformer(AffineTransform3 const & t) : transform(t) {}

  bool operator()(Mesh & mesh) const
  {
    for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      vi->setPosition(transform * vi->getPosition());

    return false;
  }

  AffineTransform3 transform;
};

bool
centerMesh(MG & mg)
{
  mg.updateBounds();
  Vector3 c = mg.getBounds().getCenter();
  Transformer tr(AffineTransform3::translation(-c));
  mg.forEachMeshUntil(cref(tr));

  return true;
}

bool
rescaleMesh(MG & mg, CoordinateAxis axis, Real len)
{
  mg.updateBounds();
  Real ext = mg.getBounds().getExtent()[axis];
  if (ext > 0)
  {
    Real s = len / ext;
    Transformer tr(AffineTransform3::scaling(s));
    mg.forEachMeshUntil(cref(tr));
  }

  return true;
}

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "\nUsage: " << argv[0] << " [OPTIONS] <infile> [<infile> ...] <out_path>\n"
               << '\n'
               << "Options:\n"
               << "  --binary                         :  Force a binary output encoding wherever possible\n"
               << "  --split                          :  Make each connected component a separate submesh\n"
               << "  --center                         :  Center the mesh bounding box at the origin\n"
                  "                                      (always precedes rescale)\n"
               << "  --rescale <x|y|z> <len>          :  Rescale the mesh to a given length along an axis\n"
               << "  --extract <labels-file> <label>  :  Extract faces with a given label from the input\n"
               << '\n';

  return -1;
}

int
main(int argc, char * argv[])
{
  CodecObj<Mesh>::Ptr codec_obj(new CodecObj<Mesh>(CodecObj<Mesh>::ReadOptions().setReadNormals(false)
                                                                                .setReadTexCoords(false)));
  Codec3ds<Mesh>::Ptr codec_3ds(new Codec3ds<Mesh>(Codec3ds<Mesh>::ReadOptions().setReadTexCoords(false)));
  Array<MeshCodec<Mesh>::Ptr> read_codecs;
  read_codecs.push_back(codec_obj);
  read_codecs.push_back(codec_3ds);

  CodecOff<Mesh>::Ptr codec_off_bin(new CodecOff<Mesh>(CodecOff<Mesh>::ReadOptions(),
                                                       CodecOff<Mesh>::WriteOptions().setBinary(true)));
  CodecPly<Mesh>::Ptr codec_ply_bin(new CodecPly<Mesh>(CodecPly<Mesh>::ReadOptions(),
                                                       CodecPly<Mesh>::WriteOptions().setBinary(true)));

  try
  {
    Array<string> paths;
    bool force_binary = false;
    bool do_split = false;
    bool do_center = false;
    bool do_rescale = false;
    CoordinateAxis rescale_axis = CoordinateAxis::X;
    Real rescale_len = 1;
    bool do_extract = false;
    string labels_path;
    string label_to_extract;

    for (int i = 1; i < argc; ++i)
    {
      string arg = argv[i];
      if (arg.empty())
        continue;
      else if (arg[0] == '-')
      {
        if (arg == "--binary")
          force_binary = true;
        else if (arg == "--split")
          do_split = true;
        else if (arg == "--center")
          do_center = true;
        else if (arg == "--rescale")
        {
          if (i + 3 > argc)
            return usage(argc, argv);

          arg = toLower(argv[++i]);
          if      (arg == "x") rescale_axis = CoordinateAxis::X;
          else if (arg == "y") rescale_axis = CoordinateAxis::Y;
          else if (arg == "z") rescale_axis = CoordinateAxis::Z;
          else
          {
            THEA_ERROR << "Invalid rescale axis: " << arg;
            return -1;
          }

          arg = argv[++i];
          rescale_len = atof(arg.c_str());
          if (rescale_len <= 0)
          {
            THEA_ERROR << "Invalid rescale length: " << arg;
            return -1;
          }

          do_rescale = true;
        }
        else if (arg == "--extract")
        {
          if (i + 3 > argc)
            return usage(argc, argv);

          labels_path = argv[++i];
          label_to_extract = argv[++i];
          do_extract = true;
        }
      }
      else
        paths.push_back(arg);
    }

    if (paths.size() < 2)
      return usage(argc, argv);

    if (do_extract && paths.size() != 2)
    {
      THEA_ERROR << "Face extraction requires exactly one input mesh";
      return -1;
    }

    string out_path = paths.back();
    paths.pop_back();

    MG::Ptr main_group;
    for (size_t i = 0; i < paths.size(); ++i)
    {
      auto mg = make_shared<MG>(FilePath::baseName(paths[i]));
      mg->load(paths[i], read_codecs);

      if (paths.size() == 1)
        main_group = mg;
      else
      {
        if (!main_group)
          main_group = make_shared<MG>("MeshGroup");

        main_group->addChild(mg);
      }
    }

    if (do_extract)
      main_group = extractFaces(main_group, labels_path, label_to_extract);  // done first, affects all subsequent operations

    if (do_split)
    {
      if (!splitMesh(main_group))
        return -1;
    }

    if (do_center)
    {
      if (!centerMesh(*main_group))
        return -1;
    }

    if (do_rescale)
    {
      if (!rescaleMesh(*main_group, rescale_axis, rescale_len))
        return -1;
    }

    string out_path_lc = toLower(out_path);
    if (endsWith(out_path_lc, ".off.bin") || (force_binary && endsWith(out_path_lc, ".off")))
      main_group->save(out_path, *codec_off_bin);
    else if (force_binary && endsWith(out_path_lc, ".ply"))
      main_group->save(out_path, *codec_ply_bin);
    else
      main_group->save(out_path);
  }
  THEA_CATCH(return -1;, ERROR, "%s", "An error occurred")

  return 0;
}
