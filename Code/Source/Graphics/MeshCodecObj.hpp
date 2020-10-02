//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2009
//
//============================================================================

#ifndef __Thea_Graphics_MeshCodecObj_hpp__
#define __Thea_Graphics_MeshCodecObj_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../UnorderedMap.hpp"
#include "MeshGroup.hpp"
#include "MeshCodec.hpp"
#include "MeshType.hpp"
#include <functional>
#include <sstream>
#include <type_traits>
#include <utility>

namespace Thea {

namespace CodecObjInternal {

class VTN
{
  public:
    intx operator[](size_t i) const { return elems[i]; }
    intx & operator[](size_t i) { return elems[i]; }

    bool operator==(VTN const & rhs) const
    {
      return elems[0] == rhs.elems[0] && elems[1] == rhs.elems[1] && elems[2] == rhs.elems[2];
    }

    size_t hash() const { return boost::hash_range(elems, elems + 3); }

  private:
    intx elems[3];
};

template <typename MeshT, typename Enable = void>
struct VertexIndexMap
{
  typedef UnorderedMap<typename MeshT::Vertex const *, intx> type;
};

template <typename MeshT>
struct VertexIndexMap<MeshT, typename std::enable_if< Graphics::IsDisplayMesh<MeshT>::value >::type>
{
  typedef UnorderedMap<std::pair<MeshT const *, intx>, intx> type;
};

} // namespace CodecObjInternal

} // namespace Thea

namespace std {

template <> struct hash<Thea::CodecObjInternal::VTN>
{
  size_t operator()(Thea::CodecObjInternal::VTN const & v) const { return v.hash(); }
};

} // namespace std

namespace Thea {

/** %Codec for reading and writing OBJ files. */
template <typename MeshT, typename BuilderT>
class CodecObj : public CodecObjBase<MeshT>
{
  private:
    typedef CodecObjBase<MeshT> BaseT;
    typedef typename CodecObjInternal::VertexIndexMap<MeshT>::type VertexIndexMap;

  public:
    typedef MeshT Mesh;                                   ///< The type of mesh processed by the codec.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;          ///< A group of meshes.
    typedef typename MeshGroup::MeshPtr MeshPtr;          ///< A shared pointer to a mesh.
    typedef BuilderT Builder;                             ///< The mesh builder class used by the codec.
    typedef typename BaseT::ReadCallback ReadCallback;    ///< Called when a mesh element is read.
    typedef typename BaseT::WriteCallback WriteCallback;  ///< Called when a mesh element is written.
    using BaseT::getName;

    /** %Options for deserializing meshes. */
    class ReadOptions
    {
      private:
        bool ignore_texcoords;
        bool ignore_normals;
        bool skip_empty_meshes;
        bool flatten;
        bool store_vertex_indices;
        bool store_face_indices;
        bool strict;
        bool verbose;

        friend class CodecObj;

      public:
        /* Constructor. Sets default values. */
        ReadOptions()
        : ignore_texcoords(false), ignore_normals(false), skip_empty_meshes(true), flatten(false), store_vertex_indices(true),
          store_face_indices(true), strict(false), verbose(false) {}

        /**
         * Ignore texture coordinates when reading from/writing to the OBJ file? If false, each unique vertex/texcoord pair
         * referenced by faces is treated as a distinct vertex. Note that not all mesh types can store texture coordinates by
         * default.
         */
        ReadOptions & setIgnoreTexCoords(bool value) { ignore_texcoords = value; return *this; }

        /**
         * Ignore normals when reading from/writing to the OBJ file? If false, each unique vertex/normal pair referenced by
         * faces is treated as a distinct vertex. Note that not all mesh types can store normals by default.
         */
        ReadOptions & setIgnoreNormals(bool value) { ignore_normals = value; return *this; }

        /** Skip meshes with no vertices or faces? */
        ReadOptions & setSkipEmptyMeshes(bool value) { skip_empty_meshes = value; return *this; }

        /** Flatten mesh hierarchy into a single mesh? */
        ReadOptions & setFlatten(bool value) { flatten = value; return *this; }

        /** Store vertex indices in mesh? */
        ReadOptions & setStoreVertexIndices(bool value) { store_vertex_indices = value; return *this; }

        /** Store face indices in mesh? */
        ReadOptions & setStoreFaceIndices(bool value) { store_face_indices = value; return *this; }

        /** Treat warnings as errors */
        ReadOptions & setStrict(bool value) { strict = value; return *this; }

        /** Print debugging information? */
        ReadOptions & setVerbose(bool value) { verbose = value; return *this; }

        /**
         * The set of default options. The default options correspond to
         * ReadOptions().setIgnoreTexCoords(false).setIgnoreNormals(false).setSkipEmptyMeshes(true).setFlatten(false)
         *              .setStoreVertexIndices(true).setStoreFaceIndices(true).setVerbose(false).
         */
        static ReadOptions const & defaults() { static ReadOptions const def; return def; }

    }; // class ReadOptions

    /** %Options for serializing meshes. */
    class WriteOptions
    {
      private:
        bool ignore_texcoords;
        bool ignore_normals;
        bool skip_empty_meshes;
        bool flatten;
        bool verbose;

        friend class CodecObj;

      public:
        /* Constructor. Sets default values. */
        WriteOptions() : ignore_texcoords(false), ignore_normals(false), skip_empty_meshes(true), flatten(false), verbose(false)
        {}

        /**
         * Ignore texture coordinates when writing to the OBJ file? If false, each unique vertex/texcoord pair referenced by
         * faces is treated as a distinct vertex. Note that not all mesh types can store texture coordinates by default.
         */
        WriteOptions & setIgnoreTexCoords(bool value) { ignore_texcoords = value; return *this; }

        /**
         * Ignore normals when writing to the OBJ file? If false, each unique vertex/normal pair referenced by faces is treated
         * as a distinct vertex. Note that not all mesh types can store normals by default.
         */
        WriteOptions & setIgnoreNormals(bool value) { ignore_normals = value; return *this; }

        /** Skip meshes with no vertices or faces? */
        WriteOptions & setSkipEmptyMeshes(bool value) { skip_empty_meshes = value; return *this; }

        /** Flatten mesh hierarchy into a single mesh? */
        WriteOptions & setFlatten(bool value) { flatten = value; return *this; }

        /** Print debugging information? */
        WriteOptions & setVerbose(bool value) { verbose = value; return *this; }

        /**
         * The set of default options. The default options correspond to
         * WriteOptions().setIgnoreTexCoords(false).setIgnoreNormals(false).setSkipEmptyMeshes(true).setFlatten(false).setVerbose(false).
         */
        static WriteOptions const & defaults() { static WriteOptions const def; return def; }

    }; // class WriteOptions

    /** Constructor. */
    CodecObj(ReadOptions const & read_opts_ = ReadOptions::defaults(),
             WriteOptions const & write_opts_ = WriteOptions::defaults())
    : read_opts(read_opts_), write_opts(write_opts_) {}

    void readMeshGroup(MeshGroup & mesh_group, BinaryInputStream & input, Codec::BlockHeader const * block_header,
                       ReadCallback * callback) const
    {
      mesh_group.clear();

      BinaryInputStream * in = &input;
      BinaryInputStream::Ptr tmp_in;
      if (block_header)
      {
        if (block_header->data_size <= 0) { return; }
        tmp_in = std::make_shared<BinaryInputStream>(input, (int64)block_header->data_size);
        in = tmp_in.get();
      }

      using CodecObjInternal::VTN;
      typedef UnorderedMap<VTN, typename Builder::VertexHandle> VTNVertexMap;
      typedef UnorderedMap<intx, typename Builder::VertexHandle> IndexVertexMap;
      VTNVertexMap vtn_refs;
      IndexVertexMap vrefs;
      Array<typename Builder::VertexHandle> face;

      // OBJ is not neatly divided into separate meshes (e.g. *all* the vertices can be put at the beginning), so we need to
      // cache the vertices and add them to meshes on-demand.
      Array<Vector3> vertices;
      Array<Vector2> texcoords;
      Array<Vector3> normals;

      std::string line;
      double x, y, z;
      intx index;

      std::string group_name = std::string(mesh_group.getName()) + (read_opts.flatten ? "/FlattenedMesh" : "/AnonymousMesh0");
      int anon_index = 0;

      MeshPtr mesh;
      std::shared_ptr<Builder> bp;
      Builder * builder = nullptr;

      intx num_faces = 0;

      while (in->hasMore())
      {
        std::string line = trimWhitespace(in->readLine());
        bool done = false;
        while (line.empty() || !(line[0] == 'v' || line[0] == 'f' || line[0] == 'p' || line[0] == 'g' || line[0] == 'o'))
        {
          if (in->hasMore())
            line = trimWhitespace(in->readLine());
          else
          {
            done = true;
            break;
          }
        }

        if (done) break;

        if (line[0] == 'v' && line.length() >= 2)
        {
          std::istringstream vstr(line); vstr.get();

          if (line[1] == 't')
          {
            if (!read_opts.ignore_texcoords)  // texcoord
            {
              vstr.get();
              if (!(vstr >> x >> y))
                throw Error(std::string(getName()) + ": Could not read texture coordinate on line '" + line + '\'');

              texcoords.push_back(Vector2((Real)x, (Real)y));
            }
          }
          else if (line[1] == 'n')
          {
            if (!read_opts.ignore_normals)  // normal
            {
              vstr.get();
              if (!(vstr >> x >> y >> z))
                throw Error(std::string(getName()) + ": Could not read normal on line '" + line + '\'');

              normals.push_back(Vector3((Real)x, (Real)y, (Real)z));
            }
          }
          else if (line[1] == ' ' || line[1] == '\t')  // vertex
          {
            if (!(vstr >> x >> y >> z))
              throw Error(std::string(getName()) + ": Could not read vertex position on line '" + line + '\'');

            vertices.push_back(Vector3((Real)x, (Real)y, (Real)z));
          }
        }
        else if ((line[0] == 'f' || line[0] == 'p') && line.length() >= 2 && (line[1] == ' ' || line[1] == '\t'))  // face
        {
          // If no mesh+builder have been created yet, create them
          if (!builder)
          {
            mesh = MeshPtr(new Mesh(group_name));
            bp = std::shared_ptr<Builder>(new Builder(mesh));
            builder = bp.get();
            builder->begin();
          }

          face.clear();
          size_t field_begin = line.find_first_not_of("fp \t"), field_end = 0;
          bool bad_face = false;
          while (field_end != std::string::npos)
          {
            field_end = line.find_first_of(" \t", field_begin);

            std::istringstream fstr(line.substr(field_begin, (field_end == std::string::npos ? std::string::npos
                                                                                             : field_end - field_begin)));
            fstr.setf(std::ios::skipws);

            if (!read_opts.ignore_texcoords || !read_opts.ignore_normals)  // use the VTN map
            {
              // OBJ stores a vertex reference as VertexIndex[/[TexCoordIndex][/NormalIndex]]
              VTN vtn; vtn[0] = 0, vtn[1] = 0; vtn[2] = 0;

              bool bad_index = false;
              if (!(fstr >> vtn[0]))
                bad_index = true;

              if (!bad_index && fstr.get() == '/')
              {
                if (fstr.peek() == '/')
                {
                  if (!read_opts.ignore_normals)
                  {
                    fstr.get();
                    if (!(fstr >> vtn[2]))
                      bad_index = true;
                  }
                }
                else
                {
                  if (!(fstr >> vtn[1]))
                    bad_index = true;

                  if (!bad_index)
                  {
                    if (read_opts.ignore_texcoords)  // reset field
                      vtn[1] = 0;

                    if (!read_opts.ignore_normals && fstr.get() == '/')
                    {
                      if (!(fstr >> vtn[2]))
                        bad_index = true;
                    }
                  }
                }
              }

              if (bad_index)
              {
                THEA_WARNING << getName() << ": Could not read index";
                bad_face = true; break;
              }

              if (std::abs(vtn[0]) < 1 || std::abs(vtn[0]) > (intx)vertices.size())
              {
                THEA_WARNING << getName() << ": Vertex index " << vtn[0] << " out of bounds (#vertices = "
                             << vertices.size() << ')';
                bad_face = true; break;
              }

              if (!read_opts.ignore_texcoords && std::abs(vtn[1]) > (intx)texcoords.size())
              {
                THEA_WARNING << getName() << ": ITexture coordinate index " << vtn[1] << " out of bounds (#texcoords = "
                             << texcoords.size() << ')';
                vtn[1] = 0;
              }

              if (!read_opts.ignore_normals && std::abs(vtn[2]) > (intx)normals.size())
              {
                THEA_WARNING << getName() << ": Normal index " << vtn[2] << " out of bounds (#normals = "
                             << normals.size() << ')';
                vtn[2] = 0;
              }

              // Negative indices indicate counting from the last element. We'll subtract one eventually, so -1 should currently
              // map to vertices.size().
              if (vtn[0] < 0) vtn[0] = (intx)vertices.size()  + 1 + vtn[0];
              if (vtn[1] < 0) vtn[1] = (intx)texcoords.size() + 1 + vtn[1];
              if (vtn[2] < 0) vtn[2] = (intx)normals.size()   + 1 + vtn[2];

              // Add the vertex referenced by the triple to the mesh builder if it has not already been added
              typename VTNVertexMap::const_iterator existing = vtn_refs.find(vtn);
              if (existing == vtn_refs.end())
              {
                typename Builder::VertexHandle vref = builder->addVertex(vertices[(size_t)vtn[0] - 1],
                                                                         read_opts.store_vertex_indices ? vtn[0] - 1 : -1,
                                                                         vtn[2] > 0 ? &normals[(size_t)vtn[2] - 1] : nullptr,
                                                                         nullptr,  // color
                                                                         vtn[1] > 0 ? &texcoords[(size_t)vtn[1] - 1] : nullptr);
                if (callback)
                  callback->vertexRead(mesh.get(), vtn[0] - 1, vref);

                vtn_refs[vtn] = vref;
                face.push_back(vref);
              }
              else
                face.push_back(existing->second);
            }
            else
            {
              if (!(fstr >> index))
              {
                THEA_WARNING << getName() << ": Could not read index";
                bad_face = true; break;
              }

              if (std::abs(index) < 1 || std::abs(index) > (intx)vertices.size())
              {
                THEA_WARNING << getName() << ": Vertex index " << index << " out of bounds (#vertices = "
                             << vertices.size() << ')';
                bad_face = true; break;
              }

              // OBJ indices start from 1. Negative indices indicate counting from the last element.
              index = (index < 0 ? (intx)vertices.size() + index : index - 1);

              // Add the referenced vertex to the mesh builder if it has not already been added
              typename IndexVertexMap::const_iterator existing = vrefs.find(index);
              if (existing == vrefs.end())
              {
                typename Builder::VertexHandle vref = builder->addVertex(vertices[(size_t)index],
                                                                         (read_opts.store_vertex_indices ? index : -1));
                if (callback)
                  callback->vertexRead(mesh.get(), index, vref);

                vrefs[index] = vref;
                face.push_back(vref);
              }
              else
                face.push_back(existing->second);
            }

            if (field_end != std::string::npos)
              field_begin = line.find_first_not_of(" \t", field_end);
          }

          if (!bad_face)
          {
            typename Builder::FaceHandle fref = builder->addFace(face.begin(), face.end(),
                                                                 (read_opts.store_face_indices ? num_faces : -1));
            if (callback)
              callback->faceRead(mesh.get(), num_faces, fref);

            num_faces++;
          }
          else
          {
            if (read_opts.strict)
              throw Error(std::string(getName()) + ": Malformed face: '" + line + '\'');
            else
              THEA_WARNING << getName() << ": Skipping malformed face: '" << line << '\'';
          }
        }
        else if (!read_opts.flatten
              && ((line[0] == 'g' || line[0] == 'o') && (line.length() < 2 || line[1] == ' ' || line[1] == '\t')))  // group
        {
          // Add the previous mesh to the mesh group
          if (builder)
          {
            builder->end();
            if (builder->numFaces() > 0 || !read_opts.skip_empty_meshes)
            {
              if (read_opts.verbose)
              {
                THEA_CONSOLE << getName() << ": Mesh " << mesh->getName() << " has " << builder->numVertices()
                             << " vertices and " << builder->numFaces() << " faces";
              }

              mesh_group.addMesh(mesh);
            }
          }

          // Read the new group name
          group_name = trimWhitespace(line.substr(1));
          if (group_name.empty())
            group_name = format("%s/AnonymousMesh%d", mesh_group.getName(), ++anon_index);

          // Create a new mesh and a builder for it
          mesh = MeshPtr(new Mesh(group_name));
          vrefs.clear();  // start a new set of vertex handles for the new mesh
          vtn_refs.clear();
          bp = std::shared_ptr<Builder>(new Builder(mesh));  // old builder gets destroyed here
          builder = bp.get();
          builder->begin();
        }
        // Else ignore the line
      }

      // Add the final mesh to the mesh group
      if (builder)
      {
        builder->end();
        if (builder->numVertices() > 0 || !read_opts.skip_empty_meshes)
        {
          if (read_opts.verbose)
          {
            THEA_CONSOLE << getName() << ": Mesh " << mesh->getName() << " has " << builder->numVertices() << " vertices and "
                         << builder->numFaces() << " faces";
          }

          mesh_group.addMesh(mesh);
        }
      }

      THEA_CONSOLE << getName() << ": Read " << mesh_group.numMeshes() << " submesh(es) with a total of " << vertices.size()
                   << " vertices and " << num_faces << " faces";
    }

    void writeMeshGroup(MeshGroup const & mesh_group, BinaryOutputStream & output, bool write_block_header,
                        WriteCallback * callback) const
    {
      Codec::BlockHeader bh(this->getMagic());
      if (write_block_header)
        bh.markAndSkip(output);

      // No need to set endianness, OBJ is a purely text-based format

      VertexIndexMap vertex_indices;
      writeVertices(mesh_group, output, vertex_indices, callback);

      intx next_index = 0;
      writeFaces(mesh_group, vertex_indices, output, callback, next_index);

      if (write_block_header)
        bh.calcAndWrite(output);
    }

  private:
    /** Write out all the vertices from a mesh group and map them to indices. */
    void writeVertices(MeshGroup const & mesh_group, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
      {
        writeVertices(**mi, output, vertex_indices, callback);
      }

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
      {
        writeVertices(**ci, output, vertex_indices, callback);
      }
    }

    /** Write out all the vertices from a general mesh or DCEL mesh and map them to indices. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                      || Graphics::IsDcelMesh<_MeshT>::value, int >::type = 0 >
    void writeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      intx vertex_index = (intx)vertex_indices.size() + 1;  // OBJ numbers vertices starting from 1
      for (typename Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++vertex_index)
      {
        output.printf("v %f %f %f\n", vi->getPosition().x(), vi->getPosition().y(), vi->getPosition().z());
        vertex_indices[&(*vi)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index - 1, &(*vi));
      }

      if (!write_opts.ignore_normals)
      {
        for (typename Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
          output.printf("vn %f %f %f\n", vi->getNormal().x(), vi->getNormal().y(), vi->getNormal().z());
      }
    }

    /** Write out all the vertices from a display mesh and map them to indices. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value, int >::type = 0 >
    void writeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      typedef std::pair<_MeshT const *, intx> DisplayMeshVRef;
      typename Mesh::VertexArray const & vertices = mesh.getVertices();
      intx vertex_index = (intx)vertex_indices.size() + 1;  // OBJ numbers vertices starting from 1

      for (size_t i = 0; i < vertices.size(); ++i, ++vertex_index)
      {
        Vector3 const & v = vertices[i];
        output.printf("v %f %f %f\n", v.x(), v.y(), v.z());
        vertex_indices[DisplayMeshVRef(&mesh, (intx)i)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index - 1, (intx)i);
      }

      if (!write_opts.ignore_texcoords && mesh.hasTexCoords())
      {
        typename Mesh::TexCoordArray const & texcoords = mesh.getTexCoords();
        alwaysAssertM(texcoords.size() == vertices.size(),
                      std::string(getName()) + ": Mesh has unequal numbers of vertices and texture coordinates");

        for (size_t i = 0; i < texcoords.size(); ++i)
        {
          Vector2 const & t = texcoords[i];
          output.printf("vt %f %f\n", t.x(), t.y());
        }
      }

      if (!write_opts.ignore_normals && mesh.hasNormals())
      {
        typename Mesh::NormalArray const & normals = mesh.getNormals();
        alwaysAssertM(normals.size() == vertices.size(),
                      std::string(getName()) + ": Mesh has unequal numbers of vertices and normals");

        for (size_t i = 0; i < normals.size(); ++i)
        {
          Vector3 const & n = normals[i];
          output.printf("vn %f %f %f\n", n.x(), n.y(), n.z());
        }
      }
    }

    /** Write out all the faces from a mesh group. Returns the number of faces written. */
    void writeFaces(MeshGroup const & mesh_group, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
      {
        writeFaces(**mi, vertex_indices, output, callback, next_index);
      }

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
      {
        writeFaces(**ci, vertex_indices, output, callback, next_index);
      }
    }

    /** Write out all the faces from a general or DCEL mesh. Returns the number of faces written. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                      || Graphics::IsDcelMesh<_MeshT>::value, int >::type = 0 >
    void writeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      if (write_opts.skip_empty_meshes && mesh.numFaces() <= 0)
        return;

      if (!write_opts.flatten)
        output.printf("\ng %s\n", mesh.getName());

      for (typename Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
      {
        typename Mesh::Face const & face = *fi;
        if (face.numVertices() < 3) continue;

        std::ostringstream os; os << 'f';
        for (typename Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
        {
          typename VertexIndexMap::const_iterator ii = vertex_indices.find(*vi);
          alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

          if (!write_opts.ignore_normals)
            os << ' ' << ii->second << "//" << ii->second;
          else
            os << ' ' << ii->second;
        }

        os << '\n';
        output.writeBytes((int64)os.str().length(), os.str().data());

        if (callback) callback->faceWritten(&mesh, next_index++, &face);
      }
    }

    /** Write out all the faces from a display mesh. Returns the number of faces written. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value, int >::type = 0 >
    void writeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      if (write_opts.skip_empty_meshes && mesh.numFaces() <= 0)
        return;

      if (!write_opts.flatten)
        output.printf("\ng %s\n", mesh.getName());

      typedef std::pair<_MeshT const *, intx> DisplayMeshVRef;

      for (int type = 0; type < 2; ++type)  // 0: triangles, 1: quads
      {
        typename Mesh::IndexArray const & indices = (type == 0 ? mesh.getTriangleIndices() : mesh.getQuadIndices());
        size_t degree = (type == 0 ? 3 : 4);

        for (size_t i = 0; i < indices.size(); i += degree)
        {
          std::ostringstream os; os << 'f';

          for (size_t j = 0; j < degree; ++j)
          {
            typename VertexIndexMap::const_iterator ii = vertex_indices.find(DisplayMeshVRef(&mesh, (intx)indices[i + j]));
            alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

            if (!write_opts.ignore_texcoords && mesh.hasTexCoords())
            {
              if (!write_opts.ignore_normals && mesh.hasNormals())
                os << ' ' << ii->second << '/' << ii->second << '/' << ii->second;
              else
                os << ' ' << ii->second << '/' << ii->second;
            }
            else
            {
              if (!write_opts.ignore_normals && mesh.hasNormals())
                os << ' ' << ii->second << "//" << ii->second;
              else
                os << ' ' << ii->second;
            }
          }

          os << '\n';
          output.writeBytes((int64)os.str().length(), os.str().data());

          if (callback)
          {
            typename Mesh::Face face(const_cast<Mesh *>(&mesh), (int)degree, (type == 0), (intx)i, 1);
            callback->faceWritten(&mesh, next_index++, face);
          }
        }
      }
    }

    ReadOptions read_opts;
    WriteOptions write_opts;

}; // class CodecObj

} // namespace Thea

#endif
