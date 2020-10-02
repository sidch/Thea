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

#ifndef __Thea_Graphics_MeshCodecOff_hpp__
#define __Thea_Graphics_MeshCodecOff_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "MeshGroup.hpp"
#include "MeshCodec.hpp"
#include <algorithm>

namespace Thea {

namespace CodecOffInternal {

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

} // namespace CodecOffInternal

/** %Codec for reading and writing OFF files. */
template <typename MeshT, typename BuilderT>
class CodecOff : public CodecOffBase<MeshT>
{
  private:
    typedef CodecOffBase<MeshT> BaseT;
    typedef typename CodecOffInternal::VertexIndexMap<MeshT>::type VertexIndexMap;

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
        bool skip_empty_meshes;
        bool store_vertex_indices;
        bool store_face_indices;
        bool strict;
        bool verbose;

        friend class CodecOff;

      public:
        /** Constructor. Sets default values. */
        ReadOptions()
        : skip_empty_meshes(true), store_vertex_indices(true), store_face_indices(true), strict(false), verbose(false)
        {}

        /** Skip meshes with no faces? */
        ReadOptions & setSkipEmptyMeshes(bool value) { skip_empty_meshes = value; return *this; }

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
         * ReadOptions().setSkipEmptyMeshes(true).setStoreVertexIndices(true).setStoreFaceIndices(true).setVerbose(false).
         */
        static ReadOptions const & defaults() { static ReadOptions const def; return def; }

    }; // class ReadOptions

    /** %Options for serializing meshes. */
    class WriteOptions
    {
      private:
        bool binary;
        bool verbose;

        friend class CodecOff;

      public:
        /** Constructor. Sets default values. */
        WriteOptions() : binary(false), verbose(false) {}

        /** Write in the binary format? */
        WriteOptions & setBinary(bool value) { binary = value; return *this; }

        /** Print debugging information? */
        WriteOptions & setVerbose(bool value) { verbose = value; return *this; }

        /** The set of default options. The default options correspond to WriteOptions().setBinary(false).setVerbose(false). */
        static WriteOptions const & defaults() { static WriteOptions const def; return def; }

    }; // class WriteOptions

    /** Constructor. */
    CodecOff(ReadOptions const & read_opts_ = ReadOptions::defaults(),
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

      std::string header = trimWhitespace(in->readLine());
      if (header != "OFF" && !beginsWith(header, "OFF "))
        throw Error(std::string(getName()) + ": Invalid OFF stream (does not start with 'OFF')");

      bool binary = (header == "OFF BINARY" || beginsWith(header, "OFF BINARY "));
      if (binary)
        readBinary(mesh_group, *in, callback);
      else
        readAscii(mesh_group, *in, callback);
    }

    void writeMeshGroup(MeshGroup const & mesh_group, BinaryOutputStream & output, bool write_block_header,
                        WriteCallback * callback) const
    {
      Codec::BlockHeader bh(this->getMagic());
      if (write_block_header)
        bh.markAndSkip(output);

      { BinaryOutputStream::EndiannessScope scope(output, Endianness::BIG);  // binary OFF uses big-endian storage

        intx num_vertices = 0, num_faces = 0;
        getStats(mesh_group, num_vertices, num_faces);

        VertexIndexMap vertex_indices;
        if (write_opts.binary)
        {
          output.printf("OFF BINARY\n");
          output.writeInt32((int32)num_vertices);
          output.writeInt32((int32)num_faces);
          output.writeInt32(0);  // num_edges
        }
        else
          output.printf("OFF\n%ld %ld 0\n", num_vertices, num_faces);

        writeVertices(mesh_group, output, vertex_indices, callback);

        intx next_index = 0;
        writeFaces(mesh_group, vertex_indices, output, callback, next_index);
      }

      if (write_block_header)
        bh.calcAndWrite(output);
    }

  private:
    /** Read a mesh group in ASCII format. */
    void readAscii(MeshGroup & mesh_group, BinaryInputStream & in, ReadCallback * callback) const
    {
      std::string line = trimWhitespace(in.readLine());
      while (line.empty() || line[0] == '#')
      {
        if (in.hasMore())
          line = trimWhitespace(in.readLine());
        else
          throw Error(std::string(getName()) + ": Unexpected end of input");
      }

      std::istringstream counts(line);
      intx num_vertices, num_faces, num_edges;
      if (!(counts >> num_vertices >> num_faces >> num_edges))
        throw Error(std::string(getName()) + ": Could not read mesh statistics on line '" + line + '\'');

      THEA_CONSOLE << getName() << ": Mesh has " << num_vertices << " vertices, " << num_faces << " faces and " << num_edges
                   << " edges";

      if (read_opts.skip_empty_meshes && num_vertices <= 0)
        return;

      // Create new mesh
      MeshPtr mesh(new Mesh(std::string(mesh_group.getName()) + "/Mesh"));

      // Create a builder for the mesh
      Builder builder(mesh);
      builder.begin();

      Array<typename Builder::VertexHandle> vrefs;
      Array<typename Builder::VertexHandle> face;

      // Read list of vertices
      double x, y, z;
      for (intx v = 0; v < num_vertices; ++v)
      {
        std::string line = trimWhitespace(in.readLine());
        while (line.empty() || line[0] == '#')
        {
          if (in.hasMore())
            line = trimWhitespace(in.readLine());
          else
            throw Error(std::string(getName()) + ": Unexpected end of input");
        }

        std::istringstream vstr(line);
        if (!(vstr >> x >> y >> z))
          throw Error(std::string(getName()) + ": Could not read vertex on line '" + line + '\'');

        typename Builder::VertexHandle vref = builder.addVertex(Vector3((Real)x, (Real)y, (Real)z),
                                                                (read_opts.store_vertex_indices ? v : -1));
        if (callback)
          callback->vertexRead(mesh.get(), v, vref);

        vrefs.push_back(vref);
      }

      // Read list of faces
      intx index;
      intx num_face_vertices;
      for (intx f = 0; f < num_faces; ++f)
      {
        std::string line = trimWhitespace(in.readLine());
        while (line.empty() || line[0] == '#')
        {
          if (in.hasMore())
            line = trimWhitespace(in.readLine());
          else
            throw Error(std::string(getName()) + ": Unexpected end of input");
        }

        std::istringstream vstr(line);
        if (!(vstr >> num_face_vertices))
          throw Error(std::string(getName()) + ": Could not read number of vertices in face on line '" + line + '\'');

        if (num_face_vertices > 0)
        {
          face.resize(num_face_vertices);

          bool skip = false;
          for (intx v = 0; v < num_face_vertices && !skip; ++v)
          {
            if (!(vstr >> index))
            {
              if (read_opts.strict)
                throw Error(std::string(getName()) + ": Could not read vertex index on line '" + line + '\'');
              else
              {
                THEA_WARNING << getName() << ": Skipping face, could not read vertex index on line '" << line << '\'';
                skip = true; break;
              }
            }

            if (index < 0 || index >= (intx)vrefs.size())
            {
              if (read_opts.strict)
                throw Error(getName() + format(": Vertex index %ld out of bounds (#vertices = %ld) on line '%s'",
                                               index, (intx)vrefs.size(), line.c_str()));
              else
              {
                THEA_WARNING << getName() << ": Skipping face, vertex index " << index << " out of bounds (#vertices = "
                             << vrefs.size() << ") on line '" << line << '\'';
                skip = true; break;
              }
            }

            face[(size_t)v] = vrefs[(size_t)index];

            for (int w = 0; w < v; ++w)
              if (face[w] == face[v])  // face has repeated vertices
              {
                if (read_opts.strict)
                  throw Error(std::string(getName()) + ": Face has repeated vertices on line '" + line + '\'');
                else
                {
                  THEA_WARNING << getName() << ": Skipping face with repeated vertices on line '" << line << '\'';
                  skip = true; break;
                }
              }
          }

          if (!skip)
          {
            typename Builder::FaceHandle fref = builder.addFace(face.begin(), face.end(),
                                                                (read_opts.store_face_indices ? f : -1));
            if (callback)
              callback->faceRead(mesh.get(), f, fref);
          }
        }
      }

      builder.end();
      mesh_group.addMesh(mesh);
    }

    /** Read a mesh group in binary format. */
    void readBinary(MeshGroup & mesh_group, BinaryInputStream & in, ReadCallback * callback) const
    {
      BinaryInputStream::EndiannessScope scope(in, Endianness::BIG);  // binary OFF uses big-endian

      intx num_vertices = (intx)in.readInt32();
      intx num_faces = (intx)in.readInt32();
      intx num_edges = (intx)in.readInt32();

      THEA_CONSOLE << getName() << ": Mesh has " << num_vertices << " vertices, " << num_faces << " faces and " << num_edges
                   << " edges";

      if (read_opts.skip_empty_meshes && num_vertices <= 0)
        return;

      // Create new mesh
      MeshPtr mesh(new Mesh(std::string(mesh_group.getName()) + "/Mesh"));

      // Create a builder for the mesh
      Builder builder(mesh);
      builder.begin();

      Array<typename Builder::VertexHandle> vrefs;
      Array<typename Builder::VertexHandle> face;

      // Read list of vertices
      Vector3 vertex;
      for (intx v = 0; v < num_vertices; ++v)
      {
        vertex[0] = (Real)in.readFloat32();
        vertex[1] = (Real)in.readFloat32();
        vertex[2] = (Real)in.readFloat32();

        typename Builder::VertexHandle vref = builder.addVertex(vertex, (read_opts.store_vertex_indices ? v : -1));
        if (callback)
          callback->vertexRead(mesh.get(), v, vref);

        vrefs.push_back(vref);
      }

      // Read list of faces
      for (intx f = 0; f < num_faces; ++f)
      {
        int num_face_vertices = (int)in.readInt32();
        if (num_face_vertices > 0)
        {
          face.resize(num_face_vertices);

          bool skip = false;
          for (int v = 0; v < num_face_vertices; ++v)
          {
            if (skip)
            {
              in.skip(4);
              continue;
            }

            intx index = (intx)in.readInt32();
            if (index < 0 || index >= (intx)vrefs.size())
            {
              if (read_opts.strict)
                throw Error(getName() + format(": Vertex index %ld out of bounds (#vertices = %ld) in face %ld",
                                               index, (intx)vrefs.size(), f));
              else
              {
                THEA_WARNING << getName() << ": Skipping face " << f << ", vertex index " << index
                                          << " out of bounds (#vertices = " << vrefs.size() << ')';
                skip = true; break;
              }
            }

            face[v] = vrefs[(size_t)index];

            for (int w = 0; w < v; ++w)
              if (face[w] == face[v])  // face has repeated vertices
              {
                if (read_opts.strict)
                  throw Error(getName() + format(": Face %ld has repeated vertices", f));
                else
                {
                  THEA_WARNING << getName() << ": Skipping face " << f << " with repeated vertices";
                  skip = true; break;
                }
              }
          }

          int num_color_components = (int)in.readInt32();
          if (num_color_components > 0)
            in.skip(4 * num_color_components);

          if (!skip)
          {
            typename Builder::FaceHandle fref = builder.addFace(face.begin(), face.end(),
                                                                (read_opts.store_face_indices ? f : -1));
            if (callback)
              callback->faceRead(mesh.get(), f, fref);
          }
        }
      }

      builder.end();
      mesh_group.addMesh(mesh);
    }

    /** Count the number of vertices and faces in a mesh (increments parameters). */
    template <typename _MeshT>
    static void getStats(_MeshT const & mesh, intx & num_vertices, intx & num_faces)
    {
      num_vertices += mesh.numVertices();
      num_faces += mesh.numFaces();
    }

    /** Count the number of vertices and faces in a mesh group (increments parameters). */
    static void getStats(MeshGroup const & mesh_group, intx & num_vertices, intx & num_faces)
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
        getStats(**mi, num_vertices, num_faces);

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
        getStats(**ci, num_vertices, num_faces);
    }

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

    /** Write out all the vertices from a general or DCEL mesh and map them to indices. */
    template <typename _MeshT, typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                     || Graphics::IsDcelMesh<_MeshT>::value, int >::type = 0>
    void writeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      intx vertex_index = (intx)vertex_indices.size();
      for (typename Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++vertex_index)
      {
        if (write_opts.binary)
        {
          output.writeFloat32((float32)vi->getPosition().x());
          output.writeFloat32((float32)vi->getPosition().y());
          output.writeFloat32((float32)vi->getPosition().z());
        }
        else
          output.printf("%f %f %f\n", vi->getPosition().x(), vi->getPosition().y(), vi->getPosition().z());

        vertex_indices[&(*vi)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index, &(*vi));
      }
    }

    /** Write out all the vertices from a display mesh and map them to indices. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value, int >::type = 0 >
    void writeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      typedef std::pair<_MeshT const *, intx> DisplayMeshVRef;
      typename Mesh::VertexArray const & vertices = mesh.getVertices();
      intx vertex_index = (intx)vertex_indices.size();

      for (size_t i = 0; i < vertices.size(); ++i, ++vertex_index)
      {
        Vector3 const & v = vertices[i];

        if (write_opts.binary)
        {
          output.writeFloat32((float32)v.x());
          output.writeFloat32((float32)v.y());
          output.writeFloat32((float32)v.z());
        }
        else
          output.printf("%f %f %f\n", v.x(), v.y(), v.z());

        vertex_indices[DisplayMeshVRef(&mesh, (intx)i)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index, (intx)i);
      }
    }

    /** Write out all the faces from a mesh group. */
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

    /** Write out all the faces from a general mesh. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                      || Graphics::IsDcelMesh<_MeshT>::value, int >::type = 0 >
    void writeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      for (typename Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
      {
        typename Mesh::Face const & face = *fi;
        if (face.numVertices() < 3) continue;

        if (write_opts.binary)
        {
          output.writeInt32((int32)face.numVertices());
          for (typename Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
          {
            typename VertexIndexMap::const_iterator ii = vertex_indices.find(*vi);
            alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

            output.writeInt32((int32)ii->second);
          }

          output.writeInt32(0);  // no color components
        }
        else
        {
          std::ostringstream os; os << face.numVertices();
          for (typename Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
          {
            typename VertexIndexMap::const_iterator ii = vertex_indices.find(*vi);
            alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

            os << ' ' << ii->second;
          }

          os << '\n';
          output.writeBytes((int64)os.str().length(), os.str().data());
        }

        if (callback) callback->faceWritten(&mesh, next_index++, &face);
      }
    }

    /** Write out all the faces from a display mesh. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value, int >::type = 0 >
    void writeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      typedef std::pair<_MeshT const *, intx> DisplayMeshVRef;

      for (int type = 0; type < 2; ++type)  // 0: triangles, 1: quads
      {
        typename Mesh::IndexArray indices = (type == 0 ? mesh.getTriangleIndices() : mesh.getQuadIndices());
        size_t degree = (type == 0 ? 3 : 4);

        for (size_t i = 0; i < indices.size(); i += degree)
        {
          if (write_opts.binary)
          {
            output.writeInt32((int32)degree);
            for (size_t j = 0; j < degree; ++j)
            {
              typename VertexIndexMap::const_iterator ii = vertex_indices.find(DisplayMeshVRef(&mesh, (intx)indices[i + j]));
              alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

              output.writeInt32((int32)ii->second);
            }

            output.writeInt32(0);  // no color components
          }
          else
          {
            std::ostringstream os; os << degree;

            for (size_t j = 0; j < degree; ++j)
            {
              typename VertexIndexMap::const_iterator ii = vertex_indices.find(DisplayMeshVRef(&mesh, (intx)indices[i + j]));
              alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

              os << ' ' << ii->second;
            }

            os << '\n';
            output.writeBytes((int64)os.str().length(), os.str().data());
          }

          if (callback)
          {
            typename Mesh::Face face(const_cast<Mesh *>(&mesh), degree, (type == 0), (intx)i, 1);
            callback->faceWritten(&mesh, next_index++, face);
          }
        }
      }
    }

    ReadOptions read_opts;
    WriteOptions write_opts;

}; // class CodecOff

} // namespace Thea

#endif
