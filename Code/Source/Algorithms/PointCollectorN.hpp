//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_PointCollectorN_hpp__
#define __Thea_Algorithms_PointCollectorN_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Graphics/MeshType.hpp"
#include "PointTraitsN.hpp"
#include "../VectorN.hpp"
#include <boost/utility/enable_if.hpp>
#include <iterator>

namespace Thea {
namespace Algorithms {

namespace PointCollectorNInternal {

// Add vertices of general or DCEL mesh to a collection
template <typename MeshT, typename CollectionT, typename ScalarT, typename Enable = void>
struct MeshVertexCollector
{
  void addVertices(MeshT const & mesh, CollectionT & collection)
  {
    for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      collection.addPoint(vi->getPosition());
  }
};

// Add vertices of display mesh to a collection
template <typename MeshT, typename CollectionT, typename ScalarT>
struct MeshVertexCollector< MeshT, CollectionT, ScalarT, typename boost::enable_if< Graphics::IsDisplayMesh<MeshT> >::type >
{
  void addVertices(MeshT const & mesh, CollectionT & collection)
  {
    typename MeshT::VertexArray const & vertices = mesh.getVertices();
    for (array_size_t i = 0; i < vertices.size(); ++i)
      collection.addPoint(vertices[i]);
  }
};

} // namespace PointCollectorNInternal

/**
 * Add points from different objects to a collection. The function <tt>CollectionT::addPoint(VectorN<N, ScalarT> const &)</tt>
 * (or equivalent) must exist.
 */
template <typename CollectionT, long N, typename ScalarT = Real>
class /* THEA_API */ PointCollectorN
{
  public:
    THEA_DEF_POINTER_TYPES(PointCollectorN, shared_ptr, weak_ptr)

    typedef VectorN<N, ScalarT> Point;  ///< The point type.
    typedef CollectionT Collection;  ///< The collection type.

    /** Default constructor. Call setCollection() subsequently. */
    PointCollectorN() : collection(NULL) {}

    /** Constructor. */
    PointCollectorN(CollectionT * collection_)
    {
      setCollection(collection_);
    }

    /** Set the collection to which points will be added. */
    void setCollection(CollectionT * collection_)
    {
      alwaysAssertM(collection_, "PointCollectorN: Collection must be non-null");
      collection = collection_;
    }

    /** Get the collection to which points will be added. */
    CollectionT * getCollection() const { return collection; }

    /** Add a group of points to the collection. */
    template <typename PointInputIterator> void addPoints(PointInputIterator points_begin, PointInputIterator points_end)
    {
      alwaysAssertM(collection, "PointCollectorN: Collection must be non-null");

      for (PointInputIterator pi = points_begin; pi != points_end; ++pi)
      {
        collection->addPoint(PointTraitsN< typename std::iterator_traits<PointInputIterator>::value_type,
                                           N, ScalarT >::getPosition(*pi));
      }
    }

    /** Add all vertices of a mesh to the collection. */
    template <typename MeshT> void addMeshVertices(MeshT const & mesh)
    {
      alwaysAssertM(collection, "PointCollectorN: Collection must be non-null");

      PointCollectorNInternal::MeshVertexCollector<MeshT, CollectionT, ScalarT> mesh_vertex_collector;
      mesh_vertex_collector.addVertices(mesh, *collection);
    }

    /** Add all vertices of a mesh group to the collection. */
    template <typename MeshT> void addMeshVertices(Graphics::MeshGroup<MeshT> const & mesh_group)
    {
      for (typename Graphics::MeshGroup<MeshT>::MeshConstIterator mi = mesh_group.meshesBegin();
           mi != mesh_group.meshesEnd(); ++mi)
        addMeshVertices(**mi);

      for (typename Graphics::MeshGroup<MeshT>::GroupConstIterator ci = mesh_group.childrenBegin();
           ci != mesh_group.childrenEnd(); ++ci)
        addMeshVertices(**ci);
    }

  private:
    CollectionT * collection;  ///< The collection to which points are added.

}; // class PointCollectorN

} // namespace Algorithms
} // namespace Thea

#endif
