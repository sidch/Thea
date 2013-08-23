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
        throw Error(getName() + ": You must explicitly choose a codec for serializing mesh groups");

      try
      {
        MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
        mesh_codec.serializeMeshGroup(*this, output, true);
      }
      catch (std::bad_cast &)
      {
        // Serious programming error
        throw FatalError(getName() + ": Codec specified for mesh group serialization is not a mesh codec.");
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
          throw FatalError(getName() + ": Codec specified for mesh group deserialization is not a mesh codec.");
        }
      }

      updateBounds();
    }

    /**
     * Save the mesh group to a file. Unlike serialize(), the file will <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be saved.
     */
    void save(std::string const & filename, Codec const & codec = Codec_AUTO()) const
    {
      if (codec == Codec_AUTO())
        throw Error(getName() + ": You must explicitly choose a codec for saving mesh groups");

      BinaryOutputStream out(filename, Endianness::LITTLE);
      if (!out.ok())
        throw Error(getName() + ": Could not open mesh file for writing");

      try
      {
        MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
        mesh_codec.serializeMeshGroup(*this, out, false);
      }
      catch (std::bad_cast &)
      {
        // Serious programming error
        throw FatalError(getName() + ": Codec specified for saving mesh group is not a mesh codec.");
      }

      out.commit();
      if (!out.ok())
        throw Error(getName() + ": Could not save mesh file");
    }

    /**
     * Load the mesh from a file. Unlike deserialize(), the file should <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be loaded.
     */
    void load(std::string const & filename, Codec const & codec = Codec_AUTO())
    {
      BinaryInputStream in(filename, Endianness::LITTLE);
      int64 file_size = in.size();
      if (file_size <= 0)
        throw Error(getName() + ": Mesh file does not exist or is empty");

      if (codec == Codec_AUTO())
        deserialize_AUTO(in, false);
      else
      {
        try
        {
          MeshCodec<Mesh> const & mesh_codec = dynamic_cast< MeshCodec<Mesh> const & >(codec);
          mesh_codec.deserializeMeshGroup(*this, in, false);
        }
        catch (std::bad_cast &)
        {
          // Serious programming error
          throw FatalError(getName() + ": Codec specified for loading mesh group is not a mesh codec.");
        }
      }

      updateBounds();
    }

  private:
    /**
     * Automatically detect the type of the encoded mesh group and deserialize it appropriately. Limited to a hard-coded set
     * of standard mesh codecs in the default specialization.
     */
    void deserialize_AUTO(BinaryInputStream & input, bool read_prefixed_info)
    {
      // A set of default codecs that should be implemented for each mesh type
      static CodecOBJ<Mesh> const codec_OBJ;
      static CodecOFF<Mesh> const codec_OFF;
      static Codec3DS<Mesh> const codec_3DS;
      static int const NUM_CODECS = 3;
      static MeshCodec<Mesh> const * codecs[NUM_CODECS] = { &codec_OBJ, &codec_OFF, &codec_3DS };

      if (read_prefixed_info)
      {
        // Try to identify the codec by the magic string
        int64 pos = input.getPosition();
        std::string magic = input.readString(MeshCodec<Mesh>::MAGIC_LENGTH);
        input.setPosition(pos);

        for (int i = 0; i < NUM_CODECS; ++i)
          if (codecs[i]->getMagic() == magic)
          {
            codecs[i]->deserializeMeshGroup(*this, input, true);
            return;
          }
      }

      // Try to identify by filename extension
      std::string path = toLower(input.getPath());
      for (int i = 0; i < NUM_CODECS; ++i)
        for (array_size_t j = 0; j < codecs[i]->getExtensions().size(); ++j)
          if (endsWith(path, '.' + codecs[i]->getExtensions()[j]))
          {
            codecs[i]->deserializeMeshGroup(*this, input, false);
            return;
          }

      throw Error(getName() + ": Could not detect mesh encoding from input. Please specify the encoding explicitly.");
    }

    MeshSet meshes;
    GroupSet children;
    AxisAlignedBox3 bounds;

}; // class MeshGroup

} // namespace Graphics
} // namespace Thea

#endif
