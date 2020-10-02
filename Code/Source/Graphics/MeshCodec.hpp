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

#ifndef __Thea_Graphics_MeshCodec_hpp__
#define __Thea_Graphics_MeshCodec_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Codec.hpp"
#include "IncrementalMeshBuilder.hpp"

namespace Thea {

namespace Graphics {

// Forward declaration
template <typename MeshT> class MeshGroup;

/** Abstract base class for codecs that read/write meshes to/from streams. */
template <typename MeshT>
class MeshCodec : public Codec
{
  public:
    THEA_DECL_SMART_POINTERS(MeshCodec)

    typedef MeshT Mesh;  ///< The mesh type.
    typedef typename MeshGroup<MeshT>::ReadCallback ReadCallback;    ///< Called when a mesh element is read.
    typedef typename MeshGroup<MeshT>::WriteCallback WriteCallback;  ///< Called when a mesh element is written.

    /** Destructor. */
    virtual ~MeshCodec() {}

    /**
     * Read a mesh group from a binary output stream. Any block header, if present, is assumed to have already been extracted
     * from the stream and stored in the object pointed to by \a block_header. Else, \a block_header should be a null pointer.
     */
    virtual void readMeshGroup(MeshGroup<Mesh> & mesh_group, BinaryInputStream & input, Codec::BlockHeader const * block_header,
                               ReadCallback * callback) const = 0;

    /**
     * Write a mesh group to a binary output stream. Optionally (if \a write_block_header is true) prefixes extra information
     * about the mesh block such as its size and type (which may have not been specified in the encoding format itself).
     */
    virtual void writeMeshGroup(MeshGroup<Mesh> const & mesh_group, BinaryOutputStream & output, bool write_block_header,
                                WriteCallback * callback) const = 0;

    /** Get the filename extensions for the codec. */
    virtual Array<std::string> const & getExtensions() const = 0;

}; // class MeshCodec

} // namespace Graphics

// Define some convenience base classes that hard-code properties of codecs that are invariant of the mesh type used, and
// forward declare the codec templates.
#define THEA_DEF_MESH_CODEC(name, basename, desc, magic, ...)                                                                 \
  template <typename MeshT>                                                                                                   \
  class basename : public Graphics::MeshCodec<MeshT>                                                                          \
  {                                                                                                                           \
    public:                                                                                                                   \
      typedef MeshT Mesh;                                                                                                     \
                                                                                                                              \
      char const * getName() const { return _getName(); }                                                                     \
      Codec::MagicString const & getMagic() const { return _getMagic(); }                                                     \
      Array<std::string> const & getExtensions() const { return _getExtensions(); }                                           \
                                                                                                                              \
    protected:                                                                                                                \
      static char const * _getName() { static std::string const name_ = desc + std::string(" codec"); return name_.c_str(); } \
      static Codec::MagicString const &  _getMagic()                                                                          \
      { static Codec::MagicString const magic_ = Codec::toMagic(magic); return magic_; }                                      \
      static Array<std::string> const & _getExtensions()                                                                      \
      { static Array<std::string> const exts_ = initExts(); return exts_; }                                                   \
                                                                                                                              \
    private:                                                                                                                  \
      static Array<std::string> initExts()                                                                                    \
      {                                                                                                                       \
        Array<std::string> e;                                                                                                 \
        std::string const s[] = { __VA_ARGS__ };                                                                              \
        for (size_t i = 0; i < sizeof(s) / sizeof(std::string); ++i)                                                          \
          e.push_back(s[i]);                                                                                                  \
                                                                                                                              \
        return e;                                                                                                             \
      }                                                                                                                       \
  };                                                                                                                          \
                                                                                                                              \
  template < typename MeshT, typename BuilderT = Graphics::IncrementalMeshBuilder<MeshT> > class name;

THEA_DEF_MESH_CODEC(Codec3ds, Codec3dsBase, "3D Studio Max",             "3DS",  "3ds")
THEA_DEF_MESH_CODEC(CodecObj, CodecObjBase, "Wavefront OBJ",             "OBJ",  "obj")
THEA_DEF_MESH_CODEC(CodecOff, CodecOffBase, "Object File Format (OFF)",  "OFF",  "off", "off.bin")
THEA_DEF_MESH_CODEC(CodecPly, CodecPlyBase, "Polygon File Format (PLY)", "PLY",  "ply")

#undef THEA_DEF_MESH_CODEC

} // namespace Thea

#endif
