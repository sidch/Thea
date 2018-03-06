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

#ifndef __Thea_Graphics_MeshGroup_hpp__
#define __Thea_Graphics_MeshGroup_hpp__

#include "../Common.hpp"
#include "../FilePath.hpp"
#include "../NamedObject.hpp"
#include "../Serializable.hpp"
#include "../Set.hpp"
#include "DrawableObject.hpp"
#include "MeshCodec.hpp"

namespace Thea {
namespace Graphics {

/** A collection of meshes and subgroups. */
template <typename MeshT>
class MeshGroup : public virtual NamedObject, public DrawableObject, public Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(MeshGroup, shared_ptr, weak_ptr)

    typedef MeshT Mesh;  ///< The type of mesh in the group.
    typedef shared_ptr<Mesh> MeshPtr;  ///< A shared pointer to a mesh.

    /** Interface for callback functions that are called when a vertex or face is deserialized. */
    class ReadCallback
    {
      public:
        /** Destructor. */
        virtual ~ReadCallback() = 0;

        /**
         * Called after a vertex has been read and added to the mesh. \a index is the sequential index of the vertex in the mesh
         * file, and need not correspond to the sequence in which vertices are added to the mesh.
         */
        virtual void vertexRead(Mesh * mesh, long index, typename Mesh::VertexHandle vertex) {}

        /**
         * Called after a face has been read and added to the mesh. \a index is the sequential index of the face in the mesh
         * file, and need not correspond to the sequence in which faces are added to the mesh.
         */
        virtual void faceRead(Mesh * mesh, long index, typename Mesh::FaceHandle face) {}

    }; // class ReadCallback

    /** Interface for callback functions that are called when a vertex or face is serialized. */
    class WriteCallback
    {
      public:
        /** Destructor. */
        virtual ~WriteCallback() = 0;

        /**
         * Called after a vertex has been written to an output stream. \a index is the sequential index of the vertex in the
         * output order, or -1 if the codec does not provide this information.
         */
        virtual void vertexWritten(Mesh const * mesh, long index, typename Mesh::VertexConstHandle vertex) {}

        /**
         * Called after a face has been written to an output stream. \a index is the sequential index of the face in the output
         * order, or -1 if the codec does not provide this information.
         */
        virtual void faceWritten(Mesh const * mesh, long index, typename Mesh::FaceConstHandle face) {}

    }; // class WriteCallback

  private:
    typedef TheaSet<MeshPtr>  MeshSet;
    typedef TheaSet<Ptr>      GroupSet;

  public:
    typedef typename MeshSet::iterator         MeshIterator;
    typedef typename MeshSet::const_iterator   MeshConstIterator;
    typedef typename GroupSet::iterator        GroupIterator;
    typedef typename GroupSet::const_iterator  GroupConstIterator;

    /** Default constructor. */
    MeshGroup() : parent(NULL) {}

    /** Constructor. */
    MeshGroup(std::string const & name_) : NamedObject(name_), parent(NULL) {}

    /**
     * Construct a flat group (no children) from a set of meshes. The input iterator must dereference to a pointer to a mesh.
     * The range used is [first, last).
     */
    template <typename MeshInputIterator> MeshGroup(MeshInputIterator first, MeshInputIterator last) : meshes(first, last) {}

    /** Check if the mesh group is empty. */
    bool isEmpty() const { return meshes.empty() && children.empty(); }

    /** Get the parent group of this mesh group, or null if this is the root of the hierarchy. */
    MeshGroup * getParent() const { return parent; }

    /**
     * Get the depth of the mesh group in the hierarchy (the root is at depth 0, its children at depth 1, and so on). The
     * function takes time proportional to the depth.
     */
    long getDepth() const
    {
      long depth = 0;
      MeshGroup const * p = parent;
      while (p)
      {
        p = p->parent;
        depth++;
      }

      return depth;
    }

    /** Number of meshes in the group. */
    long numMeshes() const { return (long)meshes.size(); }

    /** Get the first mesh in the group. */
    MeshIterator meshesBegin() { return meshes.begin(); }

    /** Get the first mesh in the group. */
    MeshConstIterator meshesBegin() const { return meshes.begin(); }

    /** Get the position beyond the last mesh in the group. */
    MeshIterator meshesEnd() { return meshes.end(); }

    /** Get the position beyond the last mesh in the group. */
    MeshConstIterator meshesEnd() const { return meshes.end(); }

    /** Add a mesh to the group. No-op if the pointer is null. */
    void addMesh(MeshPtr mesh) { if (mesh) meshes.insert(mesh); }

    /** Remove a mesh from the group. */
    void removeMesh(MeshPtr mesh) { meshes.erase(mesh); }

    /** Remove all meshes from the group. */
    void clearMeshes() { meshes.clear(); }

    /** Number of children of the group. */
    long numChildren() const { return (long)children.size(); }

    /** Get the first child of the group. */
    GroupIterator childrenBegin() { return children.begin(); }

    /** Get the first child of the group. */
    GroupConstIterator childrenBegin() const { return children.begin(); }

    /** Get the position beyond the last child of the group. */
    GroupIterator childrenEnd() { return children.end(); }

    /** Get the position beyond the last child in the group. */
    GroupConstIterator childrenEnd() const { return children.end(); }

    /** Add a child to the group. No-op if the pointer is null. */
    void addChild(Ptr child)
    {
      if (child)
      {
        children.insert(child);
        child->parent = this;
      }
    }

    /** Remove a mesh from the group. */
    void removeChild(Ptr child)
    {
      if (child)
      {
        children.erase(child);
        child->parent = NULL;
      }
    }

    /** Remove all children of the group. */
    void clearChildren()
    {
      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
        (*ci)->parent = NULL;

      children.clear();
    }

    /** Remove all meshes and children of the group. */
    void clear() { clearMeshes(); clearChildren(); }

    /**
     * Apply a functor to each mesh in the group, at any level, until the functor returns true. The functor may not alter a
     * mesh. The order of processing the meshes is not specified.
     *
     * The functor should overload the () operator as follows (or be a function pointer with the equivalent signature):
     *
     * \code
     * bool operator()(Mesh const & mesh)
     * {
     *   // Do something with the mesh
     *   // Return true if we should stop after seeing this mesh, else false
     * }
     * \endcode
     *
     * @return True if the functor evaluated to true on any mesh (and hence stopped after processing this mesh), else false.
     */
    template <typename MeshFunctorT> bool forEachMeshUntil(MeshFunctorT * functor) const
    {
      // Need to explicitly const_cast to ensure the functor can't change the meshes
      for (MeshConstIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if ((*functor)(const_cast<Mesh const &>(**mi)))
          return true;

      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
        if (const_cast<MeshGroup const &>(**ci).forEachMeshUntil(functor))
          return true;

      return false;
    }

    /**
     * Apply a functor to each mesh in the group, at any level, until the functor returns true. The order of processing the
     * meshes is not specified.
     *
     * The functor should overload the () operator as follows (or be a function pointer with the equivalent signature):
     *
     * \code
     * bool operator()(Mesh & mesh)
     * {
     *   // Do something with the mesh
     *   // Return true if we should stop after seeing this mesh, else false
     * }
     * \endcode
     *
     * @return True if the functor evaluated to true on any mesh (and hence stopped after processing this mesh), else false.
     */
    template <typename MeshFunctorT> bool forEachMeshUntil(MeshFunctorT * functor)
    {
      for (MeshIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if ((*functor)(**mi))
          return true;

      for (GroupIterator ci = children.begin(); ci != children.end(); ++ci)
        if ((*ci)->forEachMeshUntil(functor))
          return true;

      return false;
    }

    void uploadToGraphicsSystem(RenderSystem & render_system)
    {
      for (MeshConstIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if (*mi) (*mi)->uploadToGraphicsSystem(render_system);

      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
        if (*ci) (*ci)->uploadToGraphicsSystem(render_system);
    }

    void draw(RenderSystem & render_system, RenderOptions const & options = RenderOptions::defaults()) const
    {
      for (MeshConstIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if (*mi) (*mi)->draw(render_system, options);

      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
        if (*ci) (*ci)->draw(render_system, options);
    }

    void updateBounds()
    {
      bounds = AxisAlignedBox3();

      for (MeshIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if (*mi)
        {
          Mesh & mesh = **mi;
          mesh.updateBounds();
          bounds.merge(mesh.getBounds());
        }

      for (GroupIterator ci = children.begin(); ci != children.end(); ++ci)
        if (*ci)
        {
          MeshGroup & child = **ci;
          child.updateBounds();
          bounds.merge(child.getBounds());
        }
    }

    AxisAlignedBox3 const & getBounds() const { return bounds; }

    using Serializable::serialize;    // Added to suppress a Clang warning about hidden overloaded virtual function, but is this
    using Serializable::deserialize;  // really necessary? Cannot be reproduced in toy example -- maybe a compiler bug? In any
                                      // case the behavior appears to be as expected with this line.

    /**
     * {@inheritDoc}
     *
     * The serialized mesh group is prefixed with a header indicating the type and size of the encoding.
     */
    void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const
    {
      serialize(output, codec, NULL);
    }

    /**
     * {@inheritDoc}
     *
     * The serialized mesh group is prefixed with a header indicating the type and size of the encoding.
     */
    void serialize(BinaryOutputStream & output, Codec const & codec, WriteCallback * callback) const
    {
      if (codec == Codec_AUTO())
        throw Error(getNameStr() + ": You must explicitly choose a codec for serializing mesh groups");

      try
      {
        MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
        mesh_codec.serializeMeshGroup(*this, output, true, callback);
      }
      catch (std::bad_cast &)
      {
        // Serious programming error
        throw FatalError(getNameStr() + ": Codec specified for mesh group serialization is not a mesh codec.");
      }
    }

    /**
     * {@inheritDoc}
     *
     * The mesh group <b>must</b> have been serialized using the layout specified in serialize().
     */
    void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO())
    {
      deserialize(input, codec, NULL);
    }

    /**
     * Read the mesh group from an input stream.
     *
     * The mesh group <b>must</b> have been serialized using the layout specified in serialize().
     */
    void deserialize(BinaryInputStream & input, Codec const & codec, ReadCallback * callback)
    {
      if (codec == Codec_AUTO())
        deserialize_AUTO(input, true, callback);
      else
      {
        try
        {
          MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
          mesh_codec.deserializeMeshGroup(*this, input, true, callback);
        }
        catch (std::bad_cast &)
        {
          // Serious programming error
          throw FatalError(getNameStr() + ": Codec specified for mesh group deserialization is not a mesh codec.");
        }
      }

      updateBounds();
    }

    /**
     * Save the mesh group to a file. Unlike serialize(), the file will <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be saved.
     */
    void save(std::string const & path, Codec const & codec = Codec_AUTO(), WriteCallback * callback = NULL) const
    {
      save(path, codec, NULL, callback);
    }

    /**
     * Save the mesh group to a file, choosing one of a given set of codecs that best matches the output path, else falling back
     * on a default option. Unlike serialize(), the file will <b>not</b> have a prefixed header. An exception will be thrown if
     * the mesh group cannot be saved.
     */
    void save(std::string const & path, TheaArray< typename MeshCodec<Mesh>::Ptr > const & codecs,
              WriteCallback * callback = NULL) const
    {
      save(path, Codec_AUTO(), &codecs, callback);
    }

    /**
     * Load the mesh from a file. Unlike deserialize(), the file should <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be loaded.
     */
    void load(std::string const & path, Codec const & codec = Codec_AUTO(), ReadCallback * callback = NULL)
    {
      load(path, codec, NULL, callback);
    }

    /**
     * Load the mesh from a file, choosing one of a given set of codecs that best matches the input path, else falling back on a
     * default option. Unlike deserialize(), the file should <b>not</b> have a prefixed header. An exception will be thrown if
     * the mesh group cannot be loaded.
     */
    void load(std::string const & path, TheaArray< typename MeshCodec<Mesh>::Ptr > const & codecs,
              ReadCallback * callback = NULL)
    {
      load(path, Codec_AUTO(), &codecs, callback);
    }

  private:
    /**
     * Save the mesh group to a file. Unlike serialize(), the file will <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be saved.
     */
    void save(std::string const & path, Codec const & codec, TheaArray< typename MeshCodec<Mesh>::Ptr > const * codecs,
              WriteCallback * callback) const
    {
      MeshCodec<Mesh> const * mesh_codec = NULL;

      if (codec == Codec_AUTO())
      {
        if (codecs)
          mesh_codec = codecFromPath(path, *codecs);
        if (!mesh_codec)  // fallback
          mesh_codec = codecFromPath(path);
        if (!mesh_codec)
          throw Error(getNameStr() + ": Could not autodetect codec for saving mesh group");
      }
      else
      {
        mesh_codec = dynamic_cast<MeshCodec<Mesh> const *>(&codec);
        if (!mesh_codec)
          throw Error(getNameStr() + ": Codec specified for saving mesh group is not a mesh codec");
      }

      BinaryOutputStream out(path, Endianness::LITTLE);
      mesh_codec->serializeMeshGroup(*this, out, false, callback);

      out.commit();
      if (!out.ok())
        throw Error(getNameStr() + ": Could not save mesh file");
    }

    /**
     * Load the mesh from a file. Unlike deserialize(), the file should <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be loaded.
     */
    void load(std::string const & path, Codec const & codec, TheaArray< typename MeshCodec<Mesh>::Ptr > const * codecs,
              ReadCallback * callback)
    {
      MeshCodec<Mesh> const * mesh_codec = NULL;

      if (codec == Codec_AUTO())
      {
        if (codecs)
          mesh_codec = codecFromPath(path, *codecs);
        if (!mesh_codec)  // fallback
          mesh_codec = codecFromPath(path);
        if (!mesh_codec)
          throw Error(getNameStr() + ": Could not autodetect codec for loading mesh group");
      }
      else
      {
        mesh_codec = dynamic_cast<MeshCodec<Mesh> const *>(&codec);
        if (!mesh_codec)
          throw Error(getNameStr() + ": Codec specified for loading mesh group is not a mesh codec");
      }

      BinaryInputStream in(path, Endianness::LITTLE);
      mesh_codec->deserializeMeshGroup(*this, in, false, callback);

      setName(FilePath::objectName(path));

      updateBounds();
    }

    /**
     * Automatically detect the type of the encoded mesh group and deserialize it appropriately. Limited to a hard-coded set
     * of standard mesh codecs in the default specialization.
     */
    void deserialize_AUTO(BinaryInputStream & input, bool read_prefixed_info, ReadCallback * callback)
    {
      if (read_prefixed_info)
      {
        // Try to identify the codec by the magic string
        int64 pos = input.getPosition();
        std::string magic = input.readString(MeshCodec<Mesh>::MAGIC_LENGTH);
        input.setPosition(pos);

        MeshCodec<Mesh> const * codec = NULL;
        long codec_index = 0;
        while ((codec = getDefaultCodec(codec_index++)))
          if (codec->getMagic() == magic)
          {
            codec->deserializeMeshGroup(*this, input, true, callback);
            return;
          }
      }

      // Try to identify by filename extension
      MeshCodec<Mesh> const * codec = codecFromPath(input.getPath());
      if (codec)
      {
        codec->deserializeMeshGroup(*this, input, false, callback);
        return;
      }

      throw Error(getNameStr() + ": Could not detect mesh encoding from input. Please specify the encoding explicitly.");
    }

    /** Try to get the appropriate codec for a mesh, given the path to the mesh. */
    static MeshCodec<Mesh> const * codecFromPath(std::string path)
    {
      path = toLower(path);
      MeshCodec<Mesh> const * codec = NULL;
      long codec_index = 0;
      while ((codec = getDefaultCodec(codec_index++)))
      {
        for (size_t j = 0; j < codec->getExtensions().size(); ++j)
          if (endsWith(path, '.' + codec->getExtensions()[j]))
            return codec;
      }

      return NULL;
    }

    /** Get one of the default mesh codecs. */
    static MeshCodec<Mesh> const * getDefaultCodec(long index)
    {
      // A set of default codecs that should be implemented for each mesh type
      static Codec3DS<Mesh> const codec_3DS;
      static CodecOBJ<Mesh> const codec_OBJ;
      static CodecOFF<Mesh> const codec_OFF;
      static CodecPLY<Mesh> const codec_PLY;
      static MeshCodec<Mesh> const * codecs[] = { &codec_3DS, &codec_OBJ, &codec_OFF, &codec_PLY };
      static long NUM_CODECS = (long)(sizeof(codecs) / sizeof(MeshCodec<Mesh> const *));

      if (index >= 0 && index < NUM_CODECS)
        return codecs[index];
      else
        return NULL;
    }

    /** Try to get the appropriate codec for a mesh, given a collection of codecs. */
    static MeshCodec<Mesh> const * codecFromPath(std::string path, TheaArray< typename MeshCodec<Mesh>::Ptr > const & codecs)
    {
      path = toLower(path);
      for (size_t i = 0; i < codecs.size(); ++i)
      {
        MeshCodec<Mesh> const * codec = codecs[i].get();
        if (!codec)
          continue;

        for (size_t j = 0; j < codec->getExtensions().size(); ++j)
          if (endsWith(path, '.' + codec->getExtensions()[j]))
            return codec;
      }

      return NULL;
    }

    MeshGroup * parent;  // not a smart pointer, to avoid circular references
    MeshSet meshes;
    GroupSet children;
    AxisAlignedBox3 bounds;

}; // class MeshGroup

// Pure virtual destructor should have a body
// http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
template <typename MeshT> inline MeshGroup<MeshT>::ReadCallback::~ReadCallback() {}
template <typename MeshT> inline MeshGroup<MeshT>::WriteCallback::~WriteCallback() {}

} // namespace Graphics
} // namespace Thea

#endif
