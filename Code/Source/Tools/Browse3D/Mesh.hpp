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
#include "../../UnorderedMap.hpp"

namespace Browse3D {

class VertexAttribute : public Graphics::ColorAttribute<ColorRGBA>
{
  public:
    VertexAttribute() : Graphics::ColorAttribute<ColorRGBA>(), parent(NULL) {}

    void draw(Graphics::RenderSystem & render_system, Graphics::AbstractRenderOptions const & options) const
    {
      if (options.sendColors() && options.useVertexData())
        Graphics::ColorAttribute<ColorRGBA>::draw(render_system, options);
    }

    void setParent(Mesh * p) { parent = p; }
    Mesh * getParent() const { return parent; }

  private:
    Mesh * parent;
};

class FaceAttribute : public Graphics::ColorAttribute<ColorRGBA>
{
  public:
    FaceAttribute() : Graphics::ColorAttribute<ColorRGBA>(), parent(NULL) {}

    void draw(Graphics::RenderSystem & render_system, Graphics::AbstractRenderOptions const & options) const
    {
      if (options.sendColors() && !options.useVertexData())
        Graphics::ColorAttribute<ColorRGBA>::draw(render_system, options);
    }

    void setParent(Mesh * p) { parent = p; }
    Mesh * getParent() const { return parent; }

  private:
    Mesh * parent;
};

// FIXME: In a huge hack, we make the assumption there is only one set of vertex indices and one set of face indices in the
// entire program, and declare these as static variables. This avoids having to pass separate structures around.
class Mesh : public Graphics::GeneralMesh<VertexAttribute, Graphics::NullAttribute, FaceAttribute>
{
  private:
    typedef Graphics::GeneralMesh<VertexAttribute, Graphics::NullAttribute, FaceAttribute> BaseType;
    typedef UnorderedMap<intx, Vertex *> IndexVertexMap;
    typedef UnorderedMap<intx, Face *> IndexFaceMap;

  public:
    THEA_DECL_SMART_POINTERS(Mesh)

    Mesh(std::string const & name = "AnonymousMesh") : NamedObject(name), BaseType(name), parent(NULL), valid_features(false) {}

    typedef BaseType::Vertex Vertex;
    typedef BaseType::Face Face;

    Vertex * addVertex(Vector3 const & point, intx index = -1, Vector3 const * normal = NULL, ColorRGBA const * color = NULL,
                       Vector2 const * texcoord = NULL)
    {
      Vertex * vertex = BaseType::addVertex(point, index, normal, color, texcoord);

      if (vertex)
      {
        index_to_vertex[index] = vertex;
        vertex->attr().setParent(this);
      }

      return vertex;
    }

    template <typename VertexInputIterator>
    Face * addFace(VertexInputIterator vi_begin, VertexInputIterator vi_end, intx index = -1)
    {
      Face * face = BaseType::addFace(vi_begin, vi_end, index);

      if (face)
      {
        index_to_face[index] = face;
        face->attr().setParent(this);
      }

      return face;
    }

    static void resetVertexIndices()
    {
      index_to_vertex.clear();
    }

    static Vertex * mapIndexToVertex(intx index)
    {
      IndexVertexMap::const_iterator existing = index_to_vertex.find(index);
      return existing == index_to_vertex.end() ? NULL : existing->second;
    }

    static void resetFaceIndices()
    {
      index_to_face.clear();
    }

    static Face * mapIndexToFace(intx index)
    {
      IndexFaceMap::const_iterator existing = index_to_face.find(index);
      return existing == index_to_face.end() ? NULL : existing->second;
    }

    void setParent(MeshGroup * p) { parent = p; }
    MeshGroup * getParent() const { return parent; }

    // Ancestor at generations = 1 is the parent. If the hierarchy is not deep enough, returns the root.
    MeshGroup * getAncestor(intx generations) const;
    bool hasAncestor(MeshGroup const * anc) const;

    Array<double> const & getFeatures() const { updateFeatures(); return features; }
    void invalidateFeatures() { valid_features = false; }

  private:
    void updateFeatures() const;

    MeshGroup * parent;

    mutable bool valid_features;
    mutable Array<double> features;

    static IndexVertexMap index_to_vertex;  // horrible hack
    static IndexFaceMap index_to_face;  // horrible hack

}; // class Mesh

bool isSimilarTo(Mesh const & lhs, Mesh const & rhs);
bool isSimilarTo(Mesh const & lhs, MeshGroup const & rhs);
bool isSimilarTo(MeshGroup const & lhs, Mesh const & rhs);
bool isSimilarTo(MeshGroup const & lhs, MeshGroup const & rhs);

} // namespace Browse3D

#endif
