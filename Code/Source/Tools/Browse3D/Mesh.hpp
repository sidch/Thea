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
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/GraphicsAttributes.hpp"
#include "../../Graphics/MeshGroup.hpp"

namespace Browse3D {

class VertexAttribute : public Graphics::ColorAttribute<ColorRGBA>
{
  public:
    VertexAttribute() : Graphics::ColorAttribute<ColorRGBA>(), index(-1), parent(NULL) {}

    void draw(Graphics::RenderSystem & render_system, Graphics::RenderOptions const & options) const
    {
      if (options.sendColors() && options.useVertexData())
        Graphics::ColorAttribute<ColorRGBA>::draw(render_system, options);
    }

    void setIndex(long i) { index = i; }
    long getIndex() const { return index; }

    void setParent(Mesh * p) { parent = p; }
    Mesh * getParent() const { return parent; }

  private:
    long index;
    Mesh * parent;
};

class FaceAttribute : public Graphics::ColorAttribute<ColorRGBA>
{
  public:
    FaceAttribute() : Graphics::ColorAttribute<ColorRGBA>(), index(-1), parent(NULL) {}

    void draw(Graphics::RenderSystem & render_system, Graphics::RenderOptions const & options) const
    {
      if (options.sendColors() && !options.useVertexData())
        Graphics::ColorAttribute<ColorRGBA>::draw(render_system, options);
    }

    void setIndex(long i) { index = i; }
    long getIndex() const { return index; }

    void setParent(Mesh * p) { parent = p; }
    Mesh * getParent() const { return parent; }

  private:
    long index;
    Mesh * parent;
};

// FIXME: In a huge hack, we make the assumption there is only one set of vertex indices and one set of face indices in the
// entire program, and declare these as static variables. This avoids having to pass separate structures around.
class Mesh : public Graphics::GeneralMesh<VertexAttribute, Graphics::NullAttribute, FaceAttribute>
{
  private:
    typedef Graphics::GeneralMesh<VertexAttribute, Graphics::NullAttribute, FaceAttribute> BaseType;

  public:
    THEA_DEF_POINTER_TYPES(Mesh, shared_ptr, weak_ptr)

    Mesh(std::string const & name = "AnonymousMesh") : NamedObject(name), BaseType(name), parent(NULL), valid_features(false) {}

    typedef BaseType::Vertex Vertex;
    typedef BaseType::Face Face;

    Vertex * addVertex(Vector3 const & point, Vector3 const * normal = NULL, ColorRGBA const * color = NULL,
                       Vector2 const * texcoord = NULL)
    {
      Vertex * vertex = BaseType::addVertex(point, normal, color, texcoord);

      // Increment regardless of whether the vertex was successfully added or not, since we want the indices to correspond
      // exactly to the input file
      long next_index = nextVertexIndex(vertex);
      if (vertex)
      {
        vertex->attr().setIndex(next_index);
        vertex->attr().setParent(this);
      }

      return vertex;
    }

    template <typename VertexInputIterator> Face * addFace(VertexInputIterator vi_begin, VertexInputIterator vi_end)
    {
      Face * face = BaseType::addFace(vi_begin, vi_end);

      // Increment regardless of whether the face was successfully added or not, since we want the indices to correspond exactly
      // to the input file
      long next_index = nextFaceIndex(face);
      if (face)
      {
        face->attr().setIndex(next_index);
        face->attr().setParent(this);
      }

      return face;
    }

    static void resetVertexIndices()
    {
      index_to_vertex.clear();
    }

    static Vertex * mapIndexToVertex(long index)
    {
      if (index < 0 || index >= (long)index_to_vertex.size())
        return NULL;

      return index_to_vertex[(array_size_t)index];
    }

    static void resetFaceIndices()
    {
      index_to_face.clear();
    }

    static Face * mapIndexToFace(long index)
    {
      if (index < 0 || index >= (long)index_to_face.size())
        return NULL;

      return index_to_face[(array_size_t)index];
    }

    void setParent(MeshGroup * p) { parent = p; }
    MeshGroup * getParent() const { return parent; }

    // Ancestor at generations = 1 is the parent. If the hierarchy is not deep enough, returns the root.
    MeshGroup * getAncestor(long generations) const;
    bool hasAncestor(MeshGroup const * anc) const;

    TheaArray<double> const & getFeatures() const { updateFeatures(); return features; }
    void invalidateFeatures() { valid_features = false; }

  private:
    static long nextVertexIndex(Vertex * vertex)
    {
      index_to_vertex.push_back(vertex);
      return (long)index_to_vertex.size() - 1;
    }

    static long nextFaceIndex(Face * face)
    {
      index_to_face.push_back(face);
      return (long)index_to_face.size() - 1;
    }

    void updateFeatures() const;

    MeshGroup * parent;

    mutable bool valid_features;
    mutable TheaArray<double> features;

    static TheaArray<Vertex *>  index_to_vertex;  // horrible hack
    static TheaArray<Face *>    index_to_face;  // horrible hack

}; // class Mesh

bool isSimilarTo(Mesh const & lhs, Mesh const & rhs);
bool isSimilarTo(Mesh const & lhs, MeshGroup const & rhs);
bool isSimilarTo(MeshGroup const & lhs, Mesh const & rhs);
bool isSimilarTo(MeshGroup const & lhs, MeshGroup const & rhs);

} // namespace Browse3D

#endif
