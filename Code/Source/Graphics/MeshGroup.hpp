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

#ifndef __Thea_Graphics_MeshGroup_hpp__
#define __Thea_Graphics_MeshGroup_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../FilePath.hpp"
#include "../NamedObject.hpp"
#include "../Serializable.hpp"
#include "../Set.hpp"
#include "IDrawable.hpp"
#include "MeshCodec.hpp"

namespace Thea {
namespace Graphics {

/** A collection of meshes and subgroups. */
template <typename MeshT>
class MeshGroup : public virtual NamedObject, public IDrawable, public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(MeshGroup)

    typedef MeshT Mesh;                                ///< The type of mesh in the group.
    typedef std::shared_ptr<Mesh> MeshPtr;             ///< A shared pointer to a mesh.
    typedef std::shared_ptr<Mesh const> MeshConstPtr;  ///< A shared pointer to a const mesh.

    /** Interface for callback functions that are called when a vertex or face is read. */
    class ReadCallback
    {
      public:
        /** Destructor. */
        virtual ~ReadCallback() = 0;

        /**
         * Called after a vertex has been read and added to the mesh. \a index is the sequential index of the vertex in the mesh
         * file, and need not correspond to the sequence in which vertices are added to the mesh.
         */
        virtual void vertexRead(Mesh * mesh, intx index, typename Mesh::VertexHandle vertex) {}

        /**
         * Called after a face has been read and added to the mesh. \a index is the sequential index of the face in the mesh
         * file, and need not correspond to the sequence in which faces are added to the mesh.
         */
        virtual void faceRead(Mesh * mesh, intx index, typename Mesh::FaceHandle face) {}

    }; // class ReadCallback

    /** Interface for callback functions that are called when a vertex or face is written. */
    class WriteCallback
    {
      public:
        /** Destructor. */
        virtual ~WriteCallback() = 0;

        /**
         * Called after a vertex has been written to an output stream. \a index is the sequential index of the vertex in the
         * output order, or -1 if the codec does not provide this information.
         */
        virtual void vertexWritten(Mesh const * mesh, intx index, typename Mesh::VertexConstHandle vertex) {}

        /**
         * Called after a face has been written to an output stream. \a index is the sequential index of the face in the output
         * order, or -1 if the codec does not provide this information.
         */
        virtual void faceWritten(Mesh const * mesh, intx index, typename Mesh::FaceConstHandle face) {}

    }; // class WriteCallback

  private:
    typedef Set<MeshPtr>  MeshSet;
    typedef Set<Ptr>      GroupSet;

  public:
    typedef typename MeshSet::iterator         MeshIterator;
    typedef typename MeshSet::const_iterator   MeshConstIterator;
    typedef typename GroupSet::iterator        GroupIterator;
    typedef typename GroupSet::const_iterator  GroupConstIterator;

    /** Default constructor. */
    MeshGroup() : parent(nullptr) {}

    /** Constructor. */
    MeshGroup(std::string const & name_) : NamedObject(name_), parent(nullptr) {}

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
    intx getDepth() const
    {
      intx depth = 0;
      MeshGroup const * p = parent;
      while (p)
      {
        p = p->parent;
        depth++;
      }

      return depth;
    }

    /** Number of meshes in the group. */
    intx numMeshes() const { return (intx)meshes.size(); }

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
    intx numChildren() const { return (intx)children.size(); }

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
        child->parent = nullptr;
      }
    }

    /** Remove all children of the group. */
    void clearChildren()
    {
      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
        (*ci)->parent = nullptr;

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
     * To pass a functor by reference (e.g. if the functor has modifiable state), wrap it in <tt>std::ref</tt>.
     *
     * @return The mesh, if any, for which the functor evaluated to true (and hence processing stopped after this mesh), else a
     *   null pointer.
     */
    template <typename MeshFunctorT> MeshConstPtr forEachMeshUntil(MeshFunctorT functor) const
    {
      // Need to explicitly const_cast to ensure the functor can't change the meshes
      for (MeshConstIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if (functor(const_cast<Mesh const &>(**mi)))
          return *mi;

      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
      {
        auto m = const_cast<MeshGroup const &>(**ci).forEachMeshUntil(functor);
        if (m) return m;
      }

      return MeshConstPtr();
    }

    /**
     * Apply a functor to each mesh in the group, at any level, until the functor returns true. The order of processing the
     * meshes is not specified.
     *
     * The functor should overload the () operator as follows (or be a function pointer with the equivalent signature):
     *
     * \code
     * bool operator()(Mesh [const] & mesh)
     * {
     *   // Do something with the mesh
     *   // Return true if we should stop after seeing this mesh, else false
     * }
     * \endcode
     *
     * To pass a functor by reference (e.g. if the functor has modifiable state), wrap it in <tt>std::ref</tt>.
     *
     * @return The mesh, if any, for which the functor evaluated to true (and hence processing stopped after this mesh), else a
     *   null pointer.
     */
    template <typename MeshFunctorT> MeshPtr forEachMeshUntil(MeshFunctorT functor)
    {
      for (MeshIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if (functor(**mi))
          return *mi;

      for (GroupIterator ci = children.begin(); ci != children.end(); ++ci)
      {
        auto m = (*ci)->forEachMeshUntil(functor);
        if (m) return m;
      }

      return MeshPtr();
    }

    /** Recompute and cache the bounding box for the mesh group. Make sure this has been called before calling getBounds(). */
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

    /**
     * Get the cached bounding box of the mesh group. Will be out-of-date unless updateBounds() has been called after all
     * modifications.
     */
    AxisAlignedBox3 const & getBounds() const { return bounds; }

    void THEA_ICALL draw(IRenderSystem * render_system, IRenderOptions const * options = nullptr) const
    {
      for (MeshConstIterator mi = meshes.begin(); mi != meshes.end(); ++mi)
        if (*mi) (*mi)->draw(render_system, options);

      for (GroupConstIterator ci = children.begin(); ci != children.end(); ++ci)
        if (*ci) (*ci)->draw(render_system, options);
    }

    using Serializable::read;   // Added to suppress a Clang warning about hidden overloaded virtual function, but is this
    using Serializable::write;  // really necessary? Cannot be reproduced in toy example -- maybe a compiler bug? In any case
                                // the behavior appears to be as expected with this line.

    void read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false)
    {
      read(input, codec, read_block_header, nullptr);
    }

    /** Read the mesh group from an input stream, with a callback function called after every element is read. */
    void read(BinaryInputStream & input, Codec const & codec, bool read_block_header, ReadCallback * callback)
    {
      if (codec == CodecAuto())
        readAuto(input, read_block_header, callback);
      else
      {
        Codec::BlockHeader bh_obj, * bh = nullptr;
        if (read_block_header)
        {
          bh_obj.read(input);
          bh = &bh_obj;
        }

        MeshCodec<Mesh> const * mesh_codec = dynamic_cast< MeshCodec<Mesh> const * >(&codec);
        if (!mesh_codec)
          throw Error(getNameStr() + ": Codec specified for reading mesh group is not a mesh codec");

        if (bh && bh->magic != mesh_codec->getMagic())
          throw Error(getNameStr() + ": Magic string mismatch");

        mesh_codec->readMeshGroup(*this, input, bh, callback);
      }

      updateBounds();
    }

    void write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false) const
    {
      write(output, codec, write_block_header, nullptr);
    }

    /** Write the mesh group to an output stream, with a callback function called after every element is written. */
    void write(BinaryOutputStream & output, Codec const & codec, bool write_block_header, WriteCallback * callback) const
    {
      if (codec == CodecAuto())
        throw Error(getNameStr() + ": You must explicitly choose a codec for writing mesh groups");

      MeshCodec<Mesh> const * mesh_codec = dynamic_cast< MeshCodec<Mesh> const * >(&codec);
      if (!mesh_codec)
        throw Error(getNameStr() + ": Codec specified for writing mesh group is not a mesh codec");

      mesh_codec->writeMeshGroup(*this, output, write_block_header, callback);
    }

    /**
     * Save the mesh group to a file. The file will <b>not</b> have a prefixed Codec::BlockHeader. An exception will be thrown
     * if the mesh group cannot be saved.
     */
    void save(std::string const & path, Codec const & codec = CodecAuto(), WriteCallback * callback = nullptr) const
    {
      save(path, codec, nullptr, callback);
    }

    /**
     * Save the mesh group to a file, choosing one of a given set of codecs that best matches the output path, else falling back
     * on a default option. The file will <b>not</b> have a prefixed header. An exception will be thrown if
     * the mesh group cannot be saved.
     */
    void save(std::string const & path, Array< typename MeshCodec<Mesh>::Ptr > const & codecs,
              WriteCallback * callback = nullptr) const
    {
      save(path, CodecAuto(), &codecs, callback);
    }

    /**
     * Load the mesh from a file. Unlike read(), the file should <b>not</b> have a prefixed header. An exception will be thrown
     * if the mesh group cannot be loaded.
     */
    void load(std::string const & path, Codec const & codec = CodecAuto(), ReadCallback * callback = nullptr)
    {
      load(path, codec, nullptr, callback);
    }

    /**
     * Load the mesh from a file, choosing one of a given set of codecs that best matches the input path, else falling back on a
     * default option. Unlike read(), the file should <b>not</b> have a prefixed header. An exception will be thrown if
     * the mesh group cannot be loaded.
     */
    void load(std::string const & path, Array< typename MeshCodec<Mesh>::Ptr > const & codecs,
              ReadCallback * callback = nullptr)
    {
      load(path, CodecAuto(), &codecs, callback);
    }

  private:
    /**
     * Save the mesh group to a file. Unlike write(), the file will <b>not</b> have a prefixed header. An exception will be
     * thrown if the mesh group cannot be saved.
     */
    void save(std::string const & path, Codec const & codec, Array< typename MeshCodec<Mesh>::Ptr > const * codecs,
              WriteCallback * callback) const
    {
      MeshCodec<Mesh> const * mesh_codec = nullptr;

      if (codec == CodecAuto())
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
      mesh_codec->writeMeshGroup(*this, out, false, callback);

      out.commit();
      if (!out.ok())
        throw Error(getNameStr() + ": Could not save mesh file");
    }

    /**
     * Load the mesh from a file. Unlike read(), the file should <b>not</b> have a prefixed header. An exception will be thrown
     * if the mesh group cannot be loaded.
     */
    void load(std::string const & path, Codec const & codec, Array< typename MeshCodec<Mesh>::Ptr > const * codecs,
              ReadCallback * callback)
    {
      MeshCodec<Mesh> const * mesh_codec = nullptr;

      if (codec == CodecAuto())
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
      mesh_codec->readMeshGroup(*this, in, nullptr, callback);

      setName(FilePath::objectName(path));

      updateBounds();
    }

    /**
     * Automatically detect the type of the encoded mesh group and read it appropriately. Limited to a hard-coded set of
     * standard mesh codecs in the default specialization.
     */
    void readAuto(BinaryInputStream & input, bool read_block_header, ReadCallback * callback)
    {
      if (read_block_header)
      {
        // Try to identify the codec by the magic string
        Codec::BlockHeader header; header.read(input);
        MeshCodec<Mesh> const * codec = nullptr;
        intx codec_index = 0;
        while ((codec = getDefaultCodec(codec_index++)))
          if (codec->getMagic() == header.magic)
          {
            codec->readMeshGroup(*this, input, &header, callback);
            return;
          }
      }
      else
      {
        // Try to identify by filename extension
        MeshCodec<Mesh> const * codec = codecFromPath(input.getPath());
        if (codec)
        {
          codec->readMeshGroup(*this, input, nullptr, callback);
          return;
        }
      }

      throw Error(getNameStr() + ": Could not detect mesh encoding from input. Please specify the encoding explicitly.");
    }

    /** Try to get the appropriate codec for a mesh, given the path to the mesh. */
    static MeshCodec<Mesh> const * codecFromPath(std::string path)
    {
      path = toLower(path);
      MeshCodec<Mesh> const * codec = nullptr;
      intx codec_index = 0;
      while ((codec = getDefaultCodec(codec_index++)))
      {
        for (size_t j = 0; j < codec->getExtensions().size(); ++j)
          if (endsWith(path, '.' + codec->getExtensions()[j]))
            return codec;
      }

      return nullptr;
    }

    /** Get one of the default mesh codecs. */
    static MeshCodec<Mesh> const * getDefaultCodec(intx index)
    {
      // A set of default codecs that should be implemented for each mesh type
      static Codec3ds<Mesh> const codec_3DS;
      static CodecObj<Mesh> const codec_OBJ;
      static CodecOff<Mesh> const codec_OFF;
      static CodecPly<Mesh> const codec_PLY;
      static MeshCodec<Mesh> const * codecs[] = { &codec_3DS, &codec_OBJ, &codec_OFF, &codec_PLY };
      static int NUM_CODECS = (intx)(sizeof(codecs) / sizeof(MeshCodec<Mesh> const *));

      if (index >= 0 && index < NUM_CODECS)
        return codecs[index];
      else
        return nullptr;
    }

    /** Try to get the appropriate codec for a mesh, given a collection of codecs. */
    static MeshCodec<Mesh> const * codecFromPath(std::string path, Array< typename MeshCodec<Mesh>::Ptr > const & codecs)
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

      return nullptr;
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
