//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Browse3D_Mesh_hpp__
#define __Browse3D_Mesh_hpp__

#include "Common.hpp"
#include "MeshFwd.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"

namespace Browse3D {

class Mesh : public Graphics::DisplayMesh
{
  private:
    typedef Graphics::DisplayMesh BaseType;

  public:
    THEA_DEF_POINTER_TYPES(Mesh, shared_ptr, weak_ptr)

    Mesh(std::string const & name = "AnonymousMesh") : NamedObject(name), BaseType(name), parent(NULL), valid_features(false) {}

    typedef BaseType::Vertex Vertex;
    typedef BaseType::Face Face;

    template <typename IndexIterator> Face addFace(IndexIterator vi_begin, IndexIterator vi_end)
    {
      Face face = BaseType::addFace(vi_begin, vi_end);

      // Increment regardless of whether the face was successfully added or not, since we want the indices to correspond exactly
      // to the input file
      long next_index = nextFaceIndex(face);

      if (face)
      {
        if (face.hasTriangles())
        {
          long num_tris = face.numTriangles();
          for (long i = 0; i < num_tris; ++i)
            tri_face_indices.push_back(next_index);

          alwaysAssertM((long)tri_face_indices.size() == numTriangles(),
                        std::string(getName()) + ": Face indices and triangle list out of sync");
        }

        if (face.hasQuads())
        {
          long num_quads = face.numQuads();
          for (long i = 0; i < num_quads; ++i)
            quad_face_indices.push_back(next_index);

          alwaysAssertM((long)quad_face_indices.size() == numQuads(),
                        std::string(getName()) + ": Face indices and quad list out of sync");
        }
      }

      return face;
    }

    long getTriFaceIndex(long tri) const { return tri_face_indices[(array_size_t)tri]; }
    long getQuadFaceIndex(long quad) const { return quad_face_indices[(array_size_t)quad]; }

    static long nextFaceIndex(Face const & face, bool reset = false)
    {
      if (reset)
      {
        indexToFace().clear();
        return 0;
      }
      else
      {
        indexToFace().push_back(face);
        return (long)indexToFace().size() - 1;
      }
    }

    static Face const & mapIndexToFace(long face)
    {
      static Face const INVALID;

      if (face < 0 || face >= (long)indexToFace().size())
        return INVALID;

      return indexToFace()[(array_size_t)face];
    }

    void setParent(MeshGroup * p) { parent = p; }
    MeshGroup * getParent() const { return parent; }

    // Ancestor at generations = 1 is the parent. If the hierarchy is not deep enough, returns the root.
    MeshGroup * getAncestor(long generations) const;
    bool hasAncestor(MeshGroup const * anc) const;

    TheaArray<double> const & getFeatures() const { updateFeatures(); return features; }
    void invalidateFeatures() { valid_features = false; }

  private:
    static TheaArray<Face> & indexToFace()
    {
      static TheaArray<Face> index_to_face;
      return index_to_face;
    }

    void updateFeatures() const;

    TheaArray<long> tri_face_indices;
    TheaArray<long> quad_face_indices;
    MeshGroup * parent;

    mutable bool valid_features;
    mutable TheaArray<double> features;

}; // class Mesh

bool isSimilarTo(Mesh const & lhs, Mesh const & rhs);
bool isSimilarTo(Mesh const & lhs, MeshGroup const & rhs);
bool isSimilarTo(MeshGroup const & lhs, Mesh const & rhs);
bool isSimilarTo(MeshGroup const & lhs, MeshGroup const & rhs);

} // namespace Browse3D

namespace Thea {
namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <>
class IsPointN<Browse3D::MeshVertex, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <>
struct PointTraitsN<Browse3D::MeshVertex, 3>
{
  static Vector3 const & getPosition(Browse3D::MeshVertex const & t) { return t.getPosition(); }
};

} // namespace Algorithms
} // namespace Thea

#endif
