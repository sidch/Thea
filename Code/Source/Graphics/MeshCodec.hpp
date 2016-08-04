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

#ifndef __Thea_Graphics_MeshCodec_hpp__
#define __Thea_Graphics_MeshCodec_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Serializable.hpp"
#include "IncrementalMeshBuilder.hpp"

namespace Thea {

namespace Graphics {

// Forward declaration
template <typename MeshT> class MeshGroup;

/** Abstract base class for all mesh codecs. */
template <typename MeshT>
class MeshCodec : public Codec
{
  public:
    THEA_DEF_POINTER_TYPES(MeshCodec, shared_ptr, weak_ptr)

    typedef MeshT Mesh;  ///< The mesh type.
    typedef typename MeshGroup<MeshT>::ReadCallback ReadCallback;  ///< Called when a mesh undergoes an incremental update.

    /** Destructor. */
    virtual ~MeshCodec() {}

    /**
     * Serialize a mesh group to a binary output stream. Optionally prefixes extra information about the mesh block such as its
     * size and type (which may have not been specified in the encoding format itself).
     */
    virtual long serializeMeshGroup(MeshGroup<Mesh> const & mesh_group, BinaryOutputStream & output, bool prefix_info)
                 const = 0;

    /**
     * Deserialize a mesh group from a binary output stream. If the <code>read_prefixed_info</code> parameter is true, extra
     * information about the mesh block (such as its size and type) will be read first from the input stream. Else, the entire
     * input will be treated as the mesh block (the size() function of the stream must return the correct value in this case).
     *
     * @see serializeMeshGroup
     */
    virtual void deserializeMeshGroup(MeshGroup<Mesh> & mesh_group, BinaryInputStream & input, bool read_prefixed_info,
                                      ReadCallback * callback) const = 0;

    /**
     * Get the <b>4-byte</b> magic string for the codec.
     *
     * @see MAGIC_LENGTH
     */
    virtual char const * getMagic() const = 0;

    /** Get the filename extensions for the codec. */
    virtual TheaArray<std::string> const & getExtensions() const = 0;

    /** The standard length of the magic string. getMagic() must conform to this value. */
    static int const MAGIC_LENGTH = 4;

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
      char const * getMagic() const { return _getMagic(); }                                                                   \
      TheaArray<std::string> const & getExtensions() const { return _getExtensions(); }                                       \
                                                                                                                              \
    protected:                                                                                                                \
      static char const * _getName() { static std::string const name_ = desc + std::string(" codec"); return name_.c_str(); } \
      static char const * _getMagic() { static char const * magic_ = magic; return magic_; }                                  \
      static TheaArray<std::string> const & _getExtensions()                                                                  \
      { static TheaArray<std::string> const exts_ = initExts(); return exts_; }                                               \
                                                                                                                              \
    private:                                                                                                                  \
      static TheaArray<std::string> initExts()                                                                                \
      {                                                                                                                       \
        TheaArray<std::string> e;                                                                                             \
        std::string const s[] = { __VA_ARGS__ };                                                                              \
        for (size_t i = 0; i < sizeof(s) / sizeof(std::string); ++i)                                                          \
          e.push_back(s[i]);                                                                                                  \
                                                                                                                              \
        return e;                                                                                                             \
      }                                                                                                                       \
  };                                                                                                                          \
                                                                                                                              \
  template < typename MeshT, typename BuilderT = Graphics::IncrementalMeshBuilder<MeshT> > class name;

THEA_DEF_MESH_CODEC(Codec3DS, Codec3DSBase, "3D Studio Max",               "MM  ", "3ds")
THEA_DEF_MESH_CODEC(CodecOBJ, CodecOBJBase, "Wavefront OBJ",               "OBJ ", "obj")
THEA_DEF_MESH_CODEC(CodecOFF, CodecOFFBase, "Object File Format (OFF)",    "OFF ", "off", "off.bin")
THEA_DEF_MESH_CODEC(CodecPLY, CodecPLYBase, "Polygon File Format (PLY)",   "PLY ", "ply")

#undef THEA_DEF_MESH_CODEC

} // namespace Thea

#endif
