//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#ifndef __Thea_Graphics_MeshCodecOFF_hpp__
#define __Thea_Graphics_MeshCodecOFF_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "MeshGroup.hpp"
#include "MeshCodec.hpp"
#include <algorithm>

namespace Thea {

namespace CodecOFFInternal {

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

} // namespace CodecOFFInternal

/** %Codec for reading and writing OFF files. */
template <typename MeshT, typename BuilderT>
class CodecOFF : public CodecOFFBase<MeshT>
{
  private:
    typedef CodecOFFBase<MeshT> BaseT;
    typedef typename CodecOFFInternal::VertexIndexMap<MeshT>::type VertexIndexMap;

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

        friend class CodecOFF;

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

        friend class CodecOFF;

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
    CodecOFF(ReadOptions const & read_opts_ = ReadOptions::defaults(),
             WriteOptions const & write_opts_ = WriteOptions::defaults())
    : read_opts(read_opts_), write_opts(write_opts_) {}

    intx serializeMeshGroup(MeshGroup const & mesh_group, BinaryOutputStream & output, bool prefix_info,
                            WriteCallback * callback) const
    {
      output.setEndianness(Endianness::LITTLE);
      int64 initial_pos = output.getPosition();

      int64 size_pos = 0;
      if (prefix_info)
      {
        output.writeBytes(BaseT::MAGIC_LENGTH, BaseT::getMagic());

        // Placeholder for the size field
        size_pos = output.getPosition();
        output.writeUInt32(0);
      }

      int64 enc_start = output.getPosition();

        intx num_vertices = 0, num_faces = 0;
        getStats(mesh_group, num_vertices, num_faces);

        VertexIndexMap vertex_indices;
        if (write_opts.binary)
        {
          output.setEndianness(Endianness::BIG);  // binary OFF uses big-endian storage
          writeString("OFF BINARY\n", output);
          output.writeInt32((int32)num_vertices);
          output.writeInt32((int32)num_faces);
          output.writeInt32(0);  // num_edges
        }
        else
        {
          writeString("OFF\n", output);
          writeString(format("%ld %ld 0\n", num_vertices, num_faces), output);
        }

        serializeVertices(mesh_group, output, vertex_indices, callback);

        intx next_index = 0;
        serializeFaces(mesh_group, vertex_indices, output, callback, next_index);

      int64 enc_end = output.getPosition();

      if (prefix_info)
      {
        output.setEndianness(Endianness::LITTLE);
        output.setPosition(size_pos);
        output.writeUInt32((uint32)(enc_end - enc_start));
      }

      return (intx)(enc_end - initial_pos);
    }

    void deserializeMeshGroup(MeshGroup & mesh_group, BinaryInputStream & input, bool read_prefixed_info,
                              ReadCallback * callback) const
    {
      mesh_group.clear();

      BinaryInputStream * in = &input;
      Array<uint8> enc_block;
      BinaryInputStream::Ptr tmp_in;

      if (read_prefixed_info)
      {
        input.setEndianness(Endianness::LITTLE);
        input.skip(BaseT::MAGIC_LENGTH);
        uint32 encoding_size = input.readUInt32();

        if (encoding_size <= 0)
          return;

        enc_block.resize((size_t)encoding_size);
        input.readBytes((int64)encoding_size, &enc_block[0]);

        tmp_in = BinaryInputStream::Ptr(new BinaryInputStream(&enc_block[0], (int64)encoding_size, Endianness::BIG, false));
                                                              // shared pointer ensures deallocation on return
        in = tmp_in.get();
      }

      std::string header = trimWhitespace(in->readLine());
      if (header != "OFF" && !beginsWith(header, "OFF "))
        throw Error(std::string(getName()) + ": Invalid OFF stream (does not start with 'OFF')");

      bool binary = (header == "OFF BINARY" || beginsWith(header, "OFF BINARY "));
      if (binary)
        deserializeBinary(mesh_group, *in, callback);
      else
        deserializeAscii(mesh_group, *in, callback);
    }

  private:
    /** Deserialize a mesh group in ASCII format. */
    void deserializeAscii(MeshGroup & mesh_group, BinaryInputStream & in, ReadCallback * callback) const
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

    /** Deserialize a mesh group in binary format. */
    void deserializeBinary(MeshGroup & mesh_group, BinaryInputStream & in, ReadCallback * callback) const
    {
      in.setEndianness(Endianness::BIG);  // binary OFF uses big-endian

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

    /** Write the bytes of a string (without any trailing zero) to a binary output stream. */
    static void writeString(std::string const & str, BinaryOutputStream & output)
    {
      output.writeBytes((int64)str.length(), str.data());
    }

    /** Write out all the vertices from a mesh group and map them to indices. */
    void serializeVertices(MeshGroup const & mesh_group, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                           WriteCallback * callback) const
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
      {
        serializeVertices(**mi, output, vertex_indices, callback);
      }

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
      {
        serializeVertices(**ci, output, vertex_indices, callback);
      }
    }

    /** Write out all the vertices from a general or DCEL mesh and map them to indices. */
    template <typename _MeshT>
    void serializeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                           WriteCallback * callback,
                           typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                 || Graphics::IsDCELMesh<_MeshT>::value >::type * dummy = NULL) const
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
          writeString(format("%f %f %f\n", vi->getPosition().x(), vi->getPosition().y(), vi->getPosition().z()), output);

        vertex_indices[&(*vi)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index, &(*vi));
      }
    }

    /** Write out all the vertices from a display mesh and map them to indices. */
    template <typename _MeshT>
    void serializeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                           WriteCallback * callback,
                           typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value >::type * dummy = NULL) const
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
          writeString(format("%f %f %f\n", v.x(), v.y(), v.z()), output);

        vertex_indices[DisplayMeshVRef(&mesh, (intx)i)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index, (intx)i);
      }
    }

    /** Write out all the faces from a mesh group. */
    void serializeFaces(MeshGroup const & mesh_group, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                        WriteCallback * callback, intx & next_index) const
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
      {
        serializeFaces(**mi, vertex_indices, output, callback, next_index);
      }

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
      {
        serializeFaces(**ci, vertex_indices, output, callback, next_index);
      }
    }

    /** Write out all the faces from a general mesh. */
    template <typename _MeshT>
    void serializeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                        WriteCallback * callback, intx & next_index,
                        typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                              || Graphics::IsDCELMesh<_MeshT>::value>::type * dummy = NULL) const
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
          writeString(os.str(), output);
        }

        if (callback) callback->faceWritten(&mesh, next_index++, &face);
      }
    }

    /** Write out all the faces from a display mesh. */
    template <typename _MeshT>
    void serializeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                        WriteCallback * callback, intx & next_index,
                        typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value >::type * dummy = NULL) const
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
            writeString(os.str(), output);
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

}; // class CodecOFF

} // namespace Thea

#endif
