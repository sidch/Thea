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
#include "../Array.hpp"
#include "../UnorderedMap.hpp"
#include "../Graphics/MeshType.hpp"

// TODO: Remove dependence on CGAL union-find by rolling one of our own
#include <CGAL/Union_find.h>
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

    /**
     * Compute all maximal edge-connected components of a CGAL mesh. Two faces are edge-connected if they share an edge. An
     * edge-connected path is a sequence of faces such that every consecutive pair is edge-connected. An edge-connected
     * component is a set of faces such that every pair is connected by an edge-connected path consisting only of faces from the
     * set.
     *
     * The face handle type of the mesh should be identical (or at least implicitly convertible) to FaceHandleT.
     *
     * @return The number of components found (== components.size())
     *
     * @see CGALMesh, CGAL::Polyhedron_3
     */
    template <typename MeshT, typename FaceHandleT>
    static long findEdgeConnected(MeshT & mesh, TheaArray< TheaArray<FaceHandleT> > & components,
                                  typename boost::enable_if< Graphics::IsCGALMesh<MeshT> >::type * dummy = NULL)
    {
      typedef CGAL::Union_find<typename MeshT::Facet *> FacetUnionFind;
      FacetUnionFind uf;

      // Begin with all faces as separate components
      TheaUnorderedMap<typename MeshT::Facet *, typename FacetUnionFind::handle> facet_handles;
      for (typename MeshT::Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); fi++)
      {
        typename MeshT::Facet * fp = &(*fi);
        typename FacetUnionFind::handle handle = uf.push_back(fp);
        facet_handles[fp] = handle;
      }

      // Now go over all edges, connecting the faces on either side
      for (typename MeshT::Edge_iterator ei = mesh.edges_begin(); ei != mesh.edges_end(); ei++)
      {
        if (!ei->is_border_edge())
        {
          typename FacetUnionFind::handle handle1 = facet_handles[&(*ei->facet())];
          typename FacetUnionFind::handle handle2 = facet_handles[&(*ei->opposite()->facet())];
          uf.unify_sets(handle1, handle2);
        }
      }

      // Reserve space for the result
      components.resize(uf.number_of_sets());

      typedef TheaUnorderedMap<typename MeshT::Facet *, array_size_t> ComponentIndexMap;
      ComponentIndexMap component_indices;

      // Loop over facets, adding each to the appropriate result subarray
      array_size_t component_index;
      for (typename MeshT::Facet_iterator fi = mesh.facets_begin(); fi != mesh.facets_end(); fi++)
      {
        typename FacetUnionFind::handle rep = uf.find(facet_handles[&(*fi)]);

        typename ComponentIndexMap::iterator existing = component_indices.find(*rep);
        if (existing == component_indices.end())
        {
          component_index = component_indices.size();
          component_indices[*rep] = component_index;
          components[component_index].clear();
          components[component_index].reserve(uf.size(rep));
        }
        else
          component_index = existing->second;

        components[component_index].push_back(fi);
      }

      return (long)components.size();
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
      typedef MeshT Mesh;
      typedef FaceT Face;

      typedef CGAL::Union_find<Face *> FaceUnionFind;
      FaceUnionFind uf;

      // Begin with all faces as separate components
      TheaUnorderedMap<Face *, typename FaceUnionFind::handle> face_handles;
      for (typename Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); fi++)
      {
        Face * face = &(*fi);
        typename FaceUnionFind::handle handle = uf.push_back(face);
        face_handles[face] = handle;
      }

      // Now go over all edges, connecting the faces on either side
      unifyAdjacentFaces(mesh, face_handles, uf);

      // Reserve space for the result
      components.resize(uf.number_of_sets());

      typedef TheaUnorderedMap<Face *, array_size_t> ComponentIndexMap;
      ComponentIndexMap component_indices;

      // Loop over faces, adding each to the appropriate result subarray
      array_size_t component_index;
      for (typename Mesh::FaceIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); fi++)
      {
        Face * face = &(*fi);
        typename FaceUnionFind::handle rep = uf.find(face_handles[face]);

        typename ComponentIndexMap::iterator existing = component_indices.find(*rep);
        if (existing == component_indices.end())
        {
          component_index = component_indices.size();
          component_indices[*rep] = component_index;
          components[component_index].clear();
          components[component_index].reserve(uf.size(rep));
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
    template <typename MeshT, typename FaceHandleMap, typename FaceUnionFind>
    static void unifyAdjacentFaces(MeshT & mesh, FaceHandleMap & face_handles, FaceUnionFind & uf,
                                   typename boost::enable_if< Graphics::IsGeneralMesh<MeshT> >::type * dummy = NULL)
    {
      for (typename MeshT::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
      {
        for (typename MeshT::Edge::FaceIterator fi = ei->facesBegin(); fi != ei->facesEnd(); ++fi)
          for (typename MeshT::Edge::FaceIterator fj = ei->facesBegin(); fj != fi; ++fj)
          {
            typename FaceUnionFind::handle handle1 = face_handles[*fi];
            typename FaceUnionFind::handle handle2 = face_handles[*fj];
            uf.unify_sets(handle1, handle2);
          }
      }
    }

    /**
     * Unify sets containing each pair of adjacent edges of a DCEL mesh.
     *
     * @see DCELMesh
     */
    template <typename MeshT, typename FaceHandleMap, typename FaceUnionFind>
    static void unifyAdjacentFaces(MeshT & mesh, FaceHandleMap & face_handles, FaceUnionFind & uf,
                                   typename boost::enable_if< Graphics::IsDCELMesh<MeshT> >::type * dummy = NULL)
    {
      for (typename MeshT::EdgeIterator ei = mesh.edgesBegin(); ei != mesh.edgesEnd(); ++ei)
      {
        if (!ei->isBoundaryEdge())
        {
          typename FaceUnionFind::handle handle1 = face_handles[ei->getFace()];
          typename FaceUnionFind::handle handle2 = face_handles[ei->twin()->getFace()];
          uf.unify_sets(handle1, handle2);
        }
      }
    }

};  // class ConnectedComponents

} // namespace Algorithms
} // namespace Thea

#endif
