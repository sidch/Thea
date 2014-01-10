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

  private:
    typedef TheaSet<MeshPtr>  MeshSet;
    typedef TheaSet<Ptr>      GroupSet;

  public:
    typedef typename MeshSet::iterator         MeshIterator;
    typedef typename MeshSet::const_iterator   MeshConstIterator;
    typedef typename GroupSet::iterator        GroupIterator;
    typedef typename GroupSet::const_iterator  GroupConstIterator;

    /** Default constructor. */
    MeshGroup() {}

    /** Constructor. */
    MeshGroup(std::string const & name_) : NamedObject(name_) {}

    /**
     * Construct a flat group (no children) from a set of meshes. The input iterator must dereference to a pointer to a mesh.
     * The range used is [first, last).
     */
    template <typename MeshInputIterator> MeshGroup(MeshInputIterator first, MeshInputIterator last) : meshes(first, last) {}

    /** Check if the mesh group is empty. */
    bool isEmpty() const { return meshes.empty() && children.empty(); }

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
    void addChild(Ptr child) { if (child) children.insert(child); }

    /** Remove a mesh from the group. */
    void removeChild(Ptr child) { children.erase(child); }

    /** Remove all children of the group. */
    void clearChildren() { children.clear(); }

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

    /**
     * {@inheritDoc}
     *
     * The serialized mesh group is prefixed with a header indicating the type and size of the encoding.
     */
    void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const
    {
      if (codec == Codec_AUTO())
        throw Error(getNameStr() + ": You must explicitly choose a codec for serializing mesh groups");

      try
      {
        MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
        mesh_codec.serializeMeshGroup(*this, output, true);
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
      if (codec == Codec_AUTO())
        deserialize_AUTO(input, true);
      else
      {
        try
        {
          MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
          mesh_codec.deserializeMeshGroup(*this, input, true);
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
    void save(std::string const & path, Codec const & codec = Codec_AUTO()) const
    {
      MeshCodec<Mesh> const * mesh_codec = NULL;

      if (codec == Codec_AUTO())
      {
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
      mesh_codec->serializeMeshGroup(*this, out, false);

      out.commit();
      if (!out.ok())
        throw Error(getNameStr() + ": Could not save mesh file");
    }

    /**
     * Load the mesh from a file. Unlike deserialize(), the file should <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be loaded.
     */
    void load(std::string const & path, Codec const & codec = Codec_AUTO())
    {
      MeshCodec<Mesh> const * mesh_codec = NULL;

      if (codec == Codec_AUTO())
      {
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
      mesh_codec->deserializeMeshGroup(*this, in, false);

      updateBounds();
    }

  private:
    /**
     * Automatically detect the type of the encoded mesh group and deserialize it appropriately. Limited to a hard-coded set
     * of standard mesh codecs in the default specialization.
     */
    void deserialize_AUTO(BinaryInputStream & input, bool read_prefixed_info)
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
            codec->deserializeMeshGroup(*this, input, true);
            return;
          }
      }

      // Try to identify by filename extension
      MeshCodec<Mesh> const * codec = codecFromPath(input.getPath());
      if (codec)
      {
        codec->deserializeMeshGroup(*this, input, false);
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
        for (array_size_t j = 0; j < codec->getExtensions().size(); ++j)
          if (endsWith(path, '.' + codec->getExtensions()[j]))
            return codec;
      }

      return NULL;
    }

    /** Get one of the default mesh codecs. */
    static MeshCodec<Mesh> const * getDefaultCodec(long index)
    {
      // A set of default codecs that should be implemented for each mesh type
      static CodecOBJ<Mesh> const codec_OBJ;
      static CodecOFF<Mesh> const codec_OFF;
      static Codec3DS<Mesh> const codec_3DS;
      static long const NUM_CODECS = 3;
      static MeshCodec<Mesh> const * codecs[NUM_CODECS] = { &codec_OBJ, &codec_OFF, &codec_3DS };

      if (index >= 0 && index < NUM_CODECS)
        return codecs[index];
      else
        return NULL;
    }

    MeshSet meshes;
    GroupSet children;
    AxisAlignedBox3 bounds;

}; // class MeshGroup

} // namespace Graphics
} // namespace Thea

#endif
