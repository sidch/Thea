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
// First version: 2011
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

class VertexAttribute : public Graphics::ColorAttribute<ColorRgba>
{
  public:
    VertexAttribute() : Graphics::ColorAttribute<ColorRgba>(), parent(nullptr) {}

    void draw(Graphics::IRenderSystem & render_system, Graphics::IRenderOptions const & options) const
    {
      if (options.sendColors() && options.useVertexData())
        Graphics::ColorAttribute<ColorRgba>::draw(render_system, options);
    }

    void setParent(Mesh * p) { parent = p; }
    Mesh * getParent() const { return parent; }

  private:
    Mesh * parent;
};

class FaceAttribute : public Graphics::ColorAttribute<ColorRgba>
{
  public:
    FaceAttribute() : Graphics::ColorAttribute<ColorRgba>(), parent(nullptr) {}

    void draw(Graphics::IRenderSystem & render_system, Graphics::IRenderOptions const & options) const
    {
      if (options.sendColors() && !options.useVertexData())
        Graphics::ColorAttribute<ColorRgba>::draw(render_system, options);
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

    Mesh(std::string const & name = "AnonymousMesh")
    : BaseType(name), parent(nullptr), valid_features(false) {}

    typedef BaseType::Vertex Vertex;
    typedef BaseType::Face Face;

    Vertex * addVertex(Vector3 const & point, intx index = -1, Vector3 const * normal = nullptr,
                       ColorRgba const * color = nullptr, Vector2 const * texcoord = nullptr)
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
      return existing == index_to_vertex.end() ? nullptr : existing->second;
    }

    static void resetFaceIndices()
    {
      index_to_face.clear();
    }

    static Face * mapIndexToFace(intx index)
    {
      IndexFaceMap::const_iterator existing = index_to_face.find(index);
      return existing == index_to_face.end() ? nullptr : existing->second;
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
