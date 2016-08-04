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

#ifndef __Thea_Graphics_MeshCodec3DS_hpp__
#define __Thea_Graphics_MeshCodec3DS_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "MeshGroup.hpp"
#include "MeshCodec.hpp"
#include <algorithm>

#if THEA_LIB3DS_VERSION_MAJOR >= 2
#  include <lib3ds.h>
#else
#  include <lib3ds/file.h>
#  include <lib3ds/io.h>
#  include <lib3ds/mesh.h>
#  include <lib3ds/node.h>
#endif

namespace Thea {

/** %Codec for reading and writing 3DS Max files. */
template <typename MeshT, typename BuilderT>
class Codec3DS : public Codec3DSBase<MeshT>
{
  private:
    typedef Codec3DSBase<MeshT> BaseT;

  public:
    typedef MeshT Mesh;  ///< The type of mesh processed by the codec.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;  ///< A group of meshes.
    typedef typename MeshGroup::MeshPtr MeshPtr;  ///< A shared pointer to a mesh.
    typedef BuilderT Builder;  ///< The mesh builder class used by the codec.
    typedef typename BaseT::ReadCallback ReadCallback;  ///< Called when a mesh undergoes an incremental update.
    using BaseT::getName;
    using BaseT::_getName;

  private:
    typedef shared_ptr<Builder> BuilderPtr;

  public:
    /** Read/write options for 3DS codec. */
    class Options
    {
      private:
        bool use_transforms;
        bool ignore_texcoords;
        bool skip_empty_meshes;
        bool flatten;
        bool verbose;

        friend class Codec3DS;

      public:
        /** Constructor. Sets default values. */
        Options() : use_transforms(false), ignore_texcoords(false), skip_empty_meshes(true), flatten(false), verbose(false) {}

        /** Apply node transforms embedded in the 3DS file? */
        Options & setUseTransforms(bool value) { use_transforms = value; return *this; }

        /** Ignore texture coordinates when reading from/writing to the 3DS file? */
        Options & setIgnoreTexCoords(bool value) { ignore_texcoords = value; return *this; }

        /** Skip meshes with no vertices or faces? */
        Options & setSkipEmptyMeshes(bool value) { skip_empty_meshes = value; return *this; }

        /** Flatten mesh hierarchy into a single mesh? */
        Options & setFlatten(bool value) { flatten = value; return *this; }

        /** Print debugging information? */
        Options & setVerbose(bool value) { verbose = value; return *this; }

        /**
         * The set of default options. The default options correspond to
         * Options().setUseTransforms(false).setIgnoreTexCoords(false).setSkipEmptyMeshes(true).setFlatten(false).setVerbose(false).
         */
        static Options const & defaults() { static Options const def; return def; }

    }; // class Options

    typedef Options ReadOptions;   ///< %Options for deserializing meshes.
    typedef Options WriteOptions;  ///< %Options for serializing meshes.

    /** Constructor. */
    Codec3DS(ReadOptions const & read_opts_ = ReadOptions::defaults(),
             WriteOptions const & write_opts_ = WriteOptions::defaults())
    : read_opts(read_opts_), write_opts(write_opts_) {}

    long serializeMeshGroup(MeshGroup const & mesh_group, BinaryOutputStream & output, bool prefix_info) const
    {
      throw Error(std::string(getName()) + ": Not implemented");
    }

    void deserializeMeshGroup(MeshGroup & mesh_group, BinaryInputStream & input, bool read_prefixed_info,
                              ReadCallback * callback) const
    {
      mesh_group.clear();

      BinaryInputStream * in = &input;
      TheaArray<uint8> enc_block;
      BinaryInputStream::Ptr tmp_in;

      if (read_prefixed_info)
      {
        input.setEndianness(Endianness::LITTLE);
        input.skip(BaseT::MAGIC_LENGTH);
        uint32 encoding_size = input.readUInt32();

        if (encoding_size <= 0)
          return;

        enc_block.resize((array_size_t)encoding_size);
        input.readBytes((int64)encoding_size, &enc_block[0]);

        tmp_in = BinaryInputStream::Ptr(new BinaryInputStream(&enc_block[0], (int64)encoding_size, Endianness::LITTLE, false));
                                                              // shared pointer ensures deallocation on return
        in = tmp_in.get();
      }

      Lib3dsFile * file3ds = lib3ds_file_new();
      if (!file3ds)
        throw Error(std::string(getName()) + ": Couldn't create lib3ds file object");

#if THEA_LIB3DS_VERSION_MAJOR >= 2

      Lib3dsIo io;

      io.self        =  static_cast<void *>(in);
      io.impl        =  NULL;
      io.seek_func   =  seekInput;
      io.tell_func   =  tellInput;
      io.read_func   =  read;
      io.write_func  =  NULL;
      io.log_func    =  log;

      if (!lib3ds_file_read(file3ds, &io))
        throw Error(std::string(getName()) + ": Error reading 3DS file");

#else

      Lib3dsIo * io = lib3ds_io_new(static_cast<void *>(in), errorInput, seekInput, tellInput, read, NULL);
      if (!io)
        throw Error(std::string(getName()) + ": Couldn't create lib3ds IO object");

      if (lib3ds_file_read(file3ds, io) != LIB3DS_TRUE)
      {
        lib3ds_io_free(io);
        throw Error(std::string(getName()) + ": Error reading 3DS file");
      }

      lib3ds_io_free(io);

#endif

      vertex_count = face_count = 0;

      if (read_opts.flatten)
      {
        MeshPtr mesh(new Mesh("FlattenedMesh"));
        Builder flat_builder(mesh);
        flat_builder.begin();

          convert3DSSubtree(file3ds->meshes, file3ds->nodes, NULL, Matrix4::identity(), &flat_builder, callback);

        flat_builder.end();
        mesh_group.addMesh(mesh);
      }
      else
        convert3DSSubtree(file3ds->meshes, file3ds->nodes, &mesh_group, Matrix4::identity(), NULL, callback);

      THEA_CONSOLE << getName() << ": Read a total of " << vertex_count << " vertices and " << face_count
                   << " faces (including malformed faces which may have been dropped from the final mesh)";

      lib3ds_file_free(file3ds);
    }

  private:
    /**
     * Convert a subtree of the 3DS scene graph to a hierarchical mesh group. Exactly one of \a mesh_group and \a flat_builder
     * must be null.
     */
    void convert3DSSubtree(Lib3dsMesh * meshes, Lib3dsNode * nodes, MeshGroup * mesh_group, Matrix4 const & accum_transform,
                           Builder * flat_builder, ReadCallback * callback) const
    {
      for (Lib3dsMesh * m = meshes; m; m = m->next)
      {
        long num_vertices = (long)m->points;
        if (num_vertices <= 0 && read_opts.skip_empty_meshes)
          continue;

        long num_faces = (long)m->faces;
        if (read_opts.verbose)
        {
          THEA_CONSOLE << getName() << ": Mesh " << m->name << " has " << num_vertices << " vertices and " << num_faces
                       << " faces";
        }

        MeshPtr mesh;
        BuilderPtr local_builder;
        Builder * builder = flat_builder;
        if (!builder)
        {
          std::string name = m->name;
          mesh = MeshPtr(new Mesh(std::string(mesh_group->getName()) + '/' + name));

          local_builder = BuilderPtr(new Builder(mesh));
          builder = local_builder.get();
          builder->begin();
        }

        Matrix4 transform = read_opts.use_transforms
                          ? accum_transform * Matrix4(m->matrix[0][0], m->matrix[1][0], m->matrix[2][0], m->matrix[3][0],
                                                      m->matrix[0][1], m->matrix[1][1], m->matrix[2][1], m->matrix[3][1],
                                                      m->matrix[0][2], m->matrix[1][2], m->matrix[2][2], m->matrix[3][2],
                                                      m->matrix[0][3], m->matrix[1][3], m->matrix[2][3], m->matrix[3][3])
                          : Matrix4::identity();

        // Read list of vertices
        TheaArray<typename Builder::VertexHandle> vrefs;
        typename Builder::VertexHandle vref;
        for (array_size_t i = 0; i < (array_size_t)num_vertices; ++i)
        {
          Vector3 vertex(m->pointL[i].pos[0], m->pointL[i].pos[1], m->pointL[i].pos[2]);

          if (!read_opts.ignore_texcoords && (long)i < (long)m->texels)
          {
            Vector2 texcoord(m->texelL[i][0], m->texelL[i][1]);
            vref = builder->addVertex(read_opts.use_transforms ? transform * vertex : vertex, NULL, NULL, &texcoord);
          }
          else
            vref = builder->addVertex(read_opts.use_transforms ? transform * vertex : vertex);

          if (callback)
            callback->vertexAdded(mesh.get(), vertex_count, vref);

          vrefs.push_back(vref);
          vertex_count++;
        }

        // Read list of faces
        typename Builder::VertexHandle face[3];
        long indices[3];
        for (long i = 0; i < num_faces; ++i)
        {
          indices[0] = (long)m->faceL[i].points[0];
          indices[1] = (long)m->faceL[i].points[1];
          indices[2] = (long)m->faceL[i].points[2];

          debugAssertM(indices[0] >= 0 && indices[0] < num_vertices
                    && indices[1] >= 0 && indices[1] < num_vertices
                    && indices[2] >= 0 && indices[2] < num_vertices, std::string(getName()) + ": Vertex index out of bounds");

          if (indices[0] == indices[1] || indices[1] == indices[2] || indices[2] == indices[0])
          {
            if (read_opts.verbose)
            {
              THEA_WARNING << getName() << ": Skipping malformed face " << i << " == " << indices[0] << ' ' << indices[1]
                           << ' ' << indices[2];
            }

            face_count++;  // increment in any case to keep face indexing correct
            continue;
          }

          face[0] = vrefs[(array_size_t)indices[0]];
          face[1] = vrefs[(array_size_t)indices[1]];
          face[2] = vrefs[(array_size_t)indices[2]];

          typename Builder::FaceHandle fref = builder->addFace(face, face + 3);
          if (callback)
            callback->faceAdded(mesh.get(), face_count, fref);

          face_count++;
        }

        if (!flat_builder)
        {
          builder->end();
          mesh_group->addMesh(mesh);
        }
      }

      // Recurse over child nodes, if any
      int i = 0;
      for (Lib3dsNode * n = nodes; n && i < 10; n = n->next, ++i)
      {
        if (n->type == LIB3DS_OBJECT_NODE)
        {
          Matrix4 transform = read_opts.use_transforms
                            ? accum_transform * Matrix4(n->matrix[0][0], n->matrix[1][0], n->matrix[2][0], n->matrix[3][0],
                                                        n->matrix[0][1], n->matrix[1][1], n->matrix[2][1], n->matrix[3][1],
                                                        n->matrix[0][2], n->matrix[1][2], n->matrix[2][2], n->matrix[3][2],
                                                        n->matrix[0][3], n->matrix[1][3], n->matrix[2][3], n->matrix[3][3])
                            : Matrix4::identity();

          if (flat_builder)
            convert3DSSubtree(n->user.mesh, n->childs, NULL, transform, flat_builder, callback);
          else
          {
            typename MeshGroup::Ptr child_group(new MeshGroup(mesh_group->getName() + format("/Child%d", i)));
            convert3DSSubtree(n->user.mesh, n->childs, child_group.get(), transform, NULL, callback);
            mesh_group->addChild(child_group);
          }
        }
      }
    }

    /** Seek in an input stream. */
    static long seekInput(void * self, long offset, Lib3dsIoSeek origin)
    {
      if (!self) return 0;
      BinaryInputStream & in = *static_cast<BinaryInputStream *>(self);

      int64 abs_pos;
      switch (origin)
      {
        case LIB3DS_SEEK_SET: abs_pos = (int64)offset; break;
        case LIB3DS_SEEK_CUR: abs_pos = (int64)(in.getPosition() + offset); break;
        case LIB3DS_SEEK_END: abs_pos = (int64)(in.size() + offset); break;
        default: throw Error(std::string(_getName()) + ": Invalid enum value specified to seekInput");
      }

      if (abs_pos >= 0 && abs_pos <= in.size())
      {
        in.setPosition(abs_pos);
        return 0;
      }
      else
        return -1;
    }

    /** Seek in an output stream. */
    static long seekOutput(void * self, long offset, Lib3dsIoSeek origin)
    {
      if (!self) return 0;
      BinaryOutputStream & out = *static_cast<BinaryOutputStream *>(self);

      int64 abs_pos;
      switch (origin)
      {
        case LIB3DS_SEEK_SET: abs_pos = (int64)offset; break;
        case LIB3DS_SEEK_CUR: abs_pos = (int64)(out.getPosition() + offset); break;
        case LIB3DS_SEEK_END: abs_pos = (int64)(out.size() + offset); break;
        default: throw Error(std::string(_getName()) + ": Invalid enum value specified to seekOutput");
      }

      if (abs_pos >= 0 && abs_pos <= out.size())
      {
        out.setPosition(abs_pos);
        return 0;
      }
      else
        return -1;
    }

    /** Get current position in an input stream. */
    static long tellInput(void * self)
    {
      return self ? (long)static_cast<BinaryInputStream *>(self)->getPosition() : -1;
    }

    /** Get current position in an output stream. */
    static long tellOutput(void * self)
    {
      return self ? (long)static_cast<BinaryOutputStream *>(self)->getPosition() : -1;
    }

    /** Read bytes from an input stream. */
    static size_t read(void * self, void * buffer, size_t size)
    {
      if (!self) return 0;
      BinaryInputStream & in = *static_cast<BinaryInputStream *>(self);

      int64 max_bytes = in.size() - in.getPosition();
      int64 bytes_to_read = std::max(std::min(max_bytes, (int64)size), (int64)0);
      if (bytes_to_read > 0) in.readBytes(bytes_to_read, buffer);

      return static_cast<size_t>(bytes_to_read);
    }

    /** Write bytes to an output stream. */
    static size_t write(void * self, void const * buffer, size_t size)
    {
      if (!self) return 0;
      BinaryOutputStream & out = *static_cast<BinaryOutputStream *>(self);

      out.writeBytes((int64)size, buffer);
      return size;
    }

#if THEA_LIB3DS_VERSION_MAJOR >= 2
    /** Print logging data. */
    static void log(void * self, Lib3dsLogLevel level, int indent, const char * msg)
    {
      std::string indent_str((size_t)indent, ' ');
      switch (level)
      {
        case LIB3DS_LOG_ERROR: THEA_ERROR   << indent_str << msg; break;
        case LIB3DS_LOG_WARN:  THEA_WARNING << indent_str << msg; break;
        case LIB3DS_LOG_DEBUG: THEA_DEBUG   << indent_str << msg; break;
        case LIB3DS_LOG_INFO:  THEA_CONSOLE << indent_str << msg; break;
        default: throw Error(std::string(_getName()) + ": Invalid enum value specified to log");
      }
    }
#else
    /** Check if an input stream has any errors. */
    static Lib3dsBool errorInput(void * self) { return LIB3DS_FALSE; }

    /** Check if an output stream has any errors. */
    static Lib3dsBool errorOutput(void * self)
    {
      if (self)
        return static_cast<BinaryOutputStream *>(self)->ok() ? LIB3DS_FALSE : LIB3DS_TRUE;
      else
        return LIB3DS_FALSE;
    }
#endif

    ReadOptions read_opts;
    WriteOptions write_opts;

    mutable long vertex_count;
    mutable long face_count;

    friend struct ::Lib3dsIo;

}; // class Codec3DS

} // namespace Thea

#endif
