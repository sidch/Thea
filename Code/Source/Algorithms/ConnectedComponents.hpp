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

#ifndef __Thea_Algorithms_ConnectedComponents_hpp__
#define __Thea_Algorithms_ConnectedComponents_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshType.hpp"
#include "../Array.hpp"
#include "../UnionFind.hpp"
#include "../UnorderedMap.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

/** Compute connected components of graphs and meshes. */
class THEA_API ConnectedComponents
{
  public:
    /**
     * Compute all maximal edge-connected components of a general mesh. Two faces are edge-connected if they share an edge. An
     * edge-connected path is a sequence of faces such that every consecutive pair is edge-connected. An edge-connected
     * component is a set of faces such that every pair is connected by an edge-connected path consisting only of faces from the
     * set.
     *
     * @return The number of components found (== components.size())
     *
     * @see GeneralMesh
     */
    template < typename MeshT, typename FaceT, typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value, int >::type = 0 >
    static intx findEdgeConnected(MeshT & mesh, Array< Array<FaceT *> > & components)
    {
      return findEdgeConnectedDefault(mesh, components);
    }

    /**
     * Compute all maximal edge-connected components of a DCEL mesh. Two faces are edge-connected if they share an edge. An
     * edge-connected path is a sequence of faces such that every consecutive pair is edge-connected. An edge-connected
     * component is a set of faces such that every pair is connected by an edge-connected path consisting only of faces from the
     * set.
     *
     * @return The number of components found (== components.size())
     *
     * @see DcelMesh
     */
    template < typename MeshT, typename FaceT, typename std::enable_if< Graphics::IsDcelMesh<MeshT>::value, int >::type = 0 >
    static intx findEdgeConnected(MeshT & mesh, Array< Array<FaceT *> > & components)
    {
      return findEdgeConnectedDefault(mesh, components);
    }

  private:
    /**
     * Compute all maximal edge-connected components of a general mesh or DCEL mesh. Two faces are edge-connected if they share
     * an edge. An edge-connected path is a sequence of faces such that every consecutive pair is edge-connected.
     * An edge-connected component is a set of faces such that every pair is connected by an edge-connected path consisting only
     * of faces from the set.
     *
     * @return The number of components found (== components.size())
     */
    template <typename MeshT, typename FaceT>
    static intx findEdgeConnectedDefault(MeshT & mesh, Array< Array<FaceT *> > & components)
    {
      // Begin with all faces as separate components
      Array<FaceT *> faces;
      collectFaces(mesh, faces);

      typedef UnionFind<FaceT *> FaceUnionFind;
      FaceUnionFind uf(faces.begin(), faces.end());

      // Now go over all edges, connecting the faces on either side
      unifyAdjacentFaces(mesh, uf);

      // Reserve space for the result
      components.resize(uf.numSets());

      typedef UnorderedMap<intx, size_t> ComponentIndexMap;
      ComponentIndexMap component_indices;

      // Loop over faces, adding each to the appropriate result subarray
      size_t component_index;
      for (size_t i = 0; i < faces.size(); ++i)
      {
        FaceT * face = faces[i];
        intx rep = uf.find(uf.getObjectId(face));
        typename ComponentIndexMap::iterator existing = component_indices.find(rep);
        if (existing == component_indices.end())
        {
          component_index = component_indices.size();
          component_indices[rep] = component_index;
          components[component_index].clear();
          components[component_index].reserve(uf.sizeOfSet(rep));
        }
        else
          component_index = existing->second;

        components[component_index].push_back(face);
      }

      return (intx)components.size();
    }

    /**
     * Unify sets containing each pair of adjacent edges of a general mesh.
     *
     * @see GeneralMesh
     */
    template < typename MeshT, typename FaceUnionFind,
               typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value, int >::type = 0 >
    static void unifyAdjacentFaces(MeshT & mesh, FaceUnionFind & uf)
    {
      for (typename MeshT::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
      {
        for (typename MeshT::Edge::FaceIterator fi = ei->facesBegin(); fi != ei->facesEnd(); ++fi)
          for (typename MeshT::Edge::FaceIterator fj = ei->facesBegin(); fj != fi; ++fj)
          {
            intx handle1 = uf.getObjectId(*fi);
            intx handle2 = uf.getObjectId(*fj);
            uf.merge(handle1, handle2);
          }
      }
    }

    /**
     * Unify sets containing each pair of adjacent edges of a DCEL mesh.
     *
     * @see DcelMesh
     */
    template < typename MeshT, typename FaceUnionFind,
               typename std::enable_if< Graphics::IsDcelMesh<MeshT>::value, int >::type = 0 >
    static void unifyAdjacentFaces(MeshT & mesh, FaceUnionFind & uf)
    {
      for (typename MeshT::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
      {
        typename MeshT::Halfedge * edge = &(*ei);
        if (!edge->isBoundaryEdge())
        {
          intx handle1 = uf.getObjectId(edge->getFace());
          intx handle2 = uf.getObjectId(edge->twin()->getFace());
          uf.merge(handle1, handle2);
        }
      }
    }

    /** Collect all the faces of a GeneralMesh or DcelMesh. */
    template < typename MeshT, typename FaceT,
               typename std::enable_if< Graphics::IsGeneralMesh<MeshT>::value
                                     || Graphics::IsDcelMesh<MeshT>::value, int >::type = 0 >
    static void collectFaces(MeshT & mesh, Array<FaceT *> & faces)
    {
      faces.clear();
      for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); fi++)
        faces.push_back(&(*fi));
    }

}; // class ConnectedComponents

} // namespace Algorithms
} // namespace Thea

#endif
