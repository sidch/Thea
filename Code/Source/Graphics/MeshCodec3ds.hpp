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

#ifndef __Thea_Graphics_MeshCodec3ds_hpp__
#define __Thea_Graphics_MeshCodec3ds_hpp__

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
class Codec3ds : public Codec3dsBase<MeshT>
{
  private:
    typedef Codec3dsBase<MeshT> BaseT;

  public:
    typedef MeshT Mesh;                                   ///< The type of mesh processed by the codec.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;          ///< A group of meshes.
    typedef typename MeshGroup::MeshPtr MeshPtr;          ///< A shared pointer to a mesh.
    typedef BuilderT Builder;                             ///< The mesh builder class used by the codec.
    typedef typename BaseT::ReadCallback ReadCallback;    ///< Called when a mesh element is read.
    typedef typename BaseT::WriteCallback WriteCallback;  ///< Called when a mesh element is written.
    using BaseT::getName;
    using BaseT::_getName;

  private:
    typedef std::shared_ptr<Builder> BuilderPtr;

  public:
    /** %Options for deserializing meshes. */
    class ReadOptions
    {
      private:
        bool use_transforms;
        bool ignore_texcoords;
        bool skip_empty_meshes;
        bool flatten;
        bool store_vertex_indices;
        bool store_face_indices;
        bool strict;
        bool verbose;

        friend class Codec3ds;

      public:
        /** Constructor. Sets default values. */
        ReadOptions()
        : use_transforms(false), ignore_texcoords(false), skip_empty_meshes(true), flatten(false), store_vertex_indices(true),
          store_face_indices(true), strict(false), verbose(false)
        {}

        /** Apply node transforms embedded in the 3DS file? */
        ReadOptions & setUseTransforms(bool value) { use_transforms = value; return *this; }

        /** Ignore texture coordinates when reading from/writing to the 3DS file? */
        ReadOptions & setIgnoreTexCoords(bool value) { ignore_texcoords = value; return *this; }

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
         * ReadOptions().setUseTransforms(false).setIgnoreTexCoords(false).setSkipEmptyMeshes(true).setFlatten(false).
         *              .setStoreVertexIndices(true).setStoreFaceIndices(true).setVerbose(false).
         */
        static ReadOptions const & defaults() { static ReadOptions const def; return def; }

    }; // class ReadOptions

    /** %Options for serializing meshes. */
    class WriteOptions
    {
      private:
        bool verbose;

        friend class Codec3ds;

      public:
        /** Constructor. Sets default values. */
        WriteOptions() : verbose(false)
        {}

        /** Print debugging information? */
        WriteOptions & setVerbose(bool value) { verbose = value; return *this; }

        /**
         * The set of default options. The default options correspond to
         * WriteOptions().setVerbose(false).
         */
        static WriteOptions const & defaults() { static WriteOptions const def; return def; }

    }; // class WriteOptions

    /** Constructor. */
    Codec3ds(ReadOptions const & read_opts_ = ReadOptions::defaults(),
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

      Lib3dsFile * file3ds = lib3ds_file_new();
      if (!file3ds)
        throw Error(std::string(getName()) + ": Couldn't create lib3ds file object");

#if THEA_LIB3DS_VERSION_MAJOR >= 2

      Lib3dsIo io;

      io.self        =  static_cast<void *>(in);
      io.impl        =  nullptr;
      io.seek_func   =  seekInput;
      io.tell_func   =  tellInput;
      io.read_func   =  read;
      io.write_func  =  nullptr;
      io.log_func    =  log;

      if (!lib3ds_file_read(file3ds, &io))
        throw Error(std::string(getName()) + ": Error reading 3DS file");

#else

      Lib3dsIo * io = lib3ds_io_new(static_cast<void *>(in), errorInput, seekInput, tellInput, read, nullptr);
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
        MeshPtr mesh(new Mesh(std::string(mesh_group.getName()) + "/FlattenedMesh"));
        Builder flat_builder(mesh);
        flat_builder.begin();

          convert3DSSubtree(file3ds->meshes, file3ds->nodes, nullptr, Matrix4::Identity(), &flat_builder, callback);

        flat_builder.end();
        mesh_group.addMesh(mesh);
      }
      else
        convert3DSSubtree(file3ds->meshes, file3ds->nodes, &mesh_group, Matrix4::Identity(), nullptr, callback);

      THEA_CONSOLE << getName() << ": Read a total of " << vertex_count << " vertices and " << face_count
                   << " faces (including malformed faces which may have been dropped from the final mesh)";

      lib3ds_file_free(file3ds);
    }

    void writeMeshGroup(MeshGroup const & mesh_group, BinaryOutputStream & output, bool write_block_header,
                        WriteCallback * callback) const
    {
      throw Error(std::string(getName()) + ": Not implemented");
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
        intx num_vertices = (intx)m->points;
        if (num_vertices <= 0 && read_opts.skip_empty_meshes)
          continue;

        intx num_faces = (intx)m->faces;
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

        Matrix4 transform;
        if (read_opts.use_transforms)
          transform = accum_transform
                    * (Matrix4() << m->matrix[0][0], m->matrix[1][0], m->matrix[2][0], m->matrix[3][0],
                                    m->matrix[0][1], m->matrix[1][1], m->matrix[2][1], m->matrix[3][1],
                                    m->matrix[0][2], m->matrix[1][2], m->matrix[2][2], m->matrix[3][2],
                                    m->matrix[0][3], m->matrix[1][3], m->matrix[2][3], m->matrix[3][3]).finished();
        else
          transform = Matrix4::Identity();

        // Read list of vertices
        Array<typename Builder::VertexHandle> vrefs;
        typename Builder::VertexHandle vref;
        for (size_t i = 0; i < (size_t)num_vertices; ++i)
        {
          Vector3 vertex(m->pointL[i].pos[0], m->pointL[i].pos[1], m->pointL[i].pos[2]);

          intx vindex = (read_opts.store_vertex_indices ? vertex_count : -1);
          if (!read_opts.ignore_texcoords && (intx)i < (intx)m->texels)
          {
            Vector2 texcoord(m->texelL[i][0], m->texelL[i][1]);
            vref = builder->addVertex(read_opts.use_transforms ? Math::hmul(transform, vertex)
                                                               : vertex, vindex, nullptr, nullptr, &texcoord);
          }
          else
            vref = builder->addVertex(read_opts.use_transforms ? Math::hmul(transform, vertex) : vertex, vindex);

          if (callback)
            callback->vertexRead(mesh.get(), vertex_count, vref);

          vrefs.push_back(vref);
          vertex_count++;
        }

        // Read list of faces
        typename Builder::VertexHandle face[3];
        intx indices[3];
        for (intx i = 0; i < num_faces; ++i)
        {
          indices[0] = (intx)m->faceL[i].points[0];
          indices[1] = (intx)m->faceL[i].points[1];
          indices[2] = (intx)m->faceL[i].points[2];

          if (indices[0] < 0 || indices[0] >= num_vertices
           || indices[1] < 0 || indices[1] >= num_vertices
           || indices[2] < 0 || indices[2] >= num_vertices)
          {
            if (read_opts.strict)
              throw Error(getName() + format(": Vertex index out of bounds in face %ld: (%ld, %ld, %ld)",
                                             i, indices[0], indices[1], indices[2]));
            else
            {
              THEA_WARNING << getName() << ": Skipping face " << i << ", vertex index out of bounds: ("
                                        << indices[0] << ", " << indices[1] << ", " << indices[2] << ')';
              face_count++;  // increment in any case to keep face indexing correct
              continue;
            }
          }

          if (indices[0] == indices[1] || indices[1] == indices[2] || indices[2] == indices[0])
          {
            if (read_opts.strict)
              throw Error(getName() + format(": Repeated vertices in face %ld: (%ld, %ld, %ld)",
                                             i, indices[0], indices[1], indices[2]));
            else
            {
              THEA_WARNING << getName() << ": Skipping face " << i << " with repeated vertices: ("
                                        << indices[0] << ", " << indices[1] << ", " << indices[2] << ')';
              face_count++;  // increment in any case to keep face indexing correct
              continue;
            }
          }

          face[0] = vrefs[(size_t)indices[0]];
          face[1] = vrefs[(size_t)indices[1]];
          face[2] = vrefs[(size_t)indices[2]];

          typename Builder::FaceHandle fref = builder->addFace(face, face + 3,
                                                               (read_opts.store_face_indices ? face_count : -1));
          if (callback)
            callback->faceRead(mesh.get(), face_count, fref);

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
          Matrix4 transform;
          if (read_opts.use_transforms)
            transform = accum_transform
                      * (Matrix4() << n->matrix[0][0], n->matrix[1][0], n->matrix[2][0], n->matrix[3][0],
                                      n->matrix[0][1], n->matrix[1][1], n->matrix[2][1], n->matrix[3][1],
                                      n->matrix[0][2], n->matrix[1][2], n->matrix[2][2], n->matrix[3][2],
                                      n->matrix[0][3], n->matrix[1][3], n->matrix[2][3], n->matrix[3][3]).finished();
          else
            transform = Matrix4::Identity();

          if (flat_builder)
            convert3DSSubtree(n->user.mesh, n->childs, nullptr, transform, flat_builder, callback);
          else
          {
            typename MeshGroup::Ptr child_group(new MeshGroup(mesh_group->getName() + format("/Child%d", i)));
            convert3DSSubtree(n->user.mesh, n->childs, child_group.get(), transform, nullptr, callback);
            mesh_group->addChild(child_group);
          }
        }
      }
    }

    /** Seek in an input stream. */
    static intx seekInput(void * self, intx offset, Lib3dsIoSeek origin)
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
    static intx seekOutput(void * self, intx offset, Lib3dsIoSeek origin)
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
    static intx tellInput(void * self)
    {
      return self ? (intx)static_cast<BinaryInputStream *>(self)->getPosition() : -1;
    }

    /** Get current position in an output stream. */
    static intx tellOutput(void * self)
    {
      return self ? (intx)static_cast<BinaryOutputStream *>(self)->getPosition() : -1;
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

    mutable intx vertex_count;
    mutable intx face_count;

    friend struct ::Lib3dsIo;

}; // class Codec3ds

} // namespace Thea

#endif
