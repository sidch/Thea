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

#ifndef __Thea_Algorithms_ConnectedComponents_hpp__
#define __Thea_Algorithms_ConnectedComponents_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshType.hpp"
#include "../Array.hpp"
#include "../UnionFind.hpp"
#include "../UnorderedMap.hpp"
#include <boost/utility/enable_if.hpp>

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
    template <typename MeshT, typename FaceT>
    static long findEdgeConnected(MeshT & mesh, TheaArray< TheaArray<FaceT *> > & components,
                                  typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type * dummy = NULL)
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
     * @see DCELMesh
     */
    template <typename MeshT, typename FaceT>
    static long findEdgeConnected(MeshT & mesh, TheaArray< TheaArray<FaceT *> > & components,
                                  typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type * dummy = NULL)
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
    static long findEdgeConnectedDefault(MeshT & mesh, TheaArray< TheaArray<FaceT *> > & components)
    {
      // Begin with all faces as separate components
      TheaArray<FaceT *> faces;
      collectFaces(mesh, faces);

      typedef UnionFind<FaceT *> FaceUnionFind;
      FaceUnionFind uf(faces.begin(), faces.end());

      // Now go over all edges, connecting the faces on either side
      unifyAdjacentFaces(mesh, uf);

      // Reserve space for the result
      components.resize(uf.numSets());

      typedef TheaUnorderedMap<long, size_t> ComponentIndexMap;
      ComponentIndexMap component_indices;

      // Loop over faces, adding each to the appropriate result subarray
      size_t component_index;
      for (size_t i = 0; i < faces.size(); ++i)
      {
        FaceT * face = faces[i];
        long rep = uf.find(uf.getObjectID(face));
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

      return (long)components.size();
    }

    /**
     * Unify sets containing each pair of adjacent edges of a general mesh.
     *
     * @see GeneralMesh
     */
    template <typename MeshT, typename FaceUnionFind>
    static void unifyAdjacentFaces(MeshT & mesh, FaceUnionFind & uf,
                                   typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type * dummy = NULL)
    {
      for (typename MeshT::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
      {
        for (typename MeshT::Edge::FaceIterator fi = ei->facesBegin(); fi != ei->facesEnd(); ++fi)
          for (typename MeshT::Edge::FaceIterator fj = ei->facesBegin(); fj != fi; ++fj)
          {
            long handle1 = uf.getObjectID(*fi);
            long handle2 = uf.getObjectID(*fj);
            uf.merge(handle1, handle2);
          }
      }
    }

    /**
     * Unify sets containing each pair of adjacent edges of a DCEL mesh.
     *
     * @see DCELMesh
     */
    template <typename MeshT, typename FaceUnionFind>
    static void unifyAdjacentFaces(MeshT & mesh, FaceUnionFind & uf,
                                   typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type * dummy = NULL)
    {
      for (typename MeshT::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
      {
        typename MeshT::Halfedge * edge = *ei;
        if (!edge->isBoundaryEdge())
        {
          long handle1 = uf.getObjectID(edge->getFace());
          long handle2 = uf.getObjectID(edge->twin()->getFace());
          uf.merge(handle1, handle2);
        }
      }
    }

    /** Collect all the faces of a GeneralMesh. */
    template <typename MeshT, typename FaceT>
    static void collectFaces(MeshT & mesh, TheaArray<FaceT *> & faces,
                             typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type * dummy = NULL)
    {
      faces.clear();
      for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); fi++)
        faces.push_back(&(*fi));
    }

    /** Collect all the faces of a DCELMesh. */
    template <typename MeshT, typename FaceT>
    static void collectFaces(MeshT & mesh, TheaArray<FaceT *> & faces,
                             typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type * dummy = NULL)
    {
      for (typename MeshT::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); fi++)
        faces.push_back(*fi);
    }

};  // class ConnectedComponents

} // namespace Algorithms
} // namespace Thea

#endif
