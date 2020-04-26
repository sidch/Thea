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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_PointCollectorN_hpp__
#define __Thea_Algorithms_PointCollectorN_hpp__

#include "../Common.hpp"
#include "PointTraitsN.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "../Graphics/MeshType.hpp"
#include "../MatVec.hpp"
#include <iterator>
#include <type_traits>

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
struct MeshVertexCollector< MeshT, CollectionT, ScalarT,
                            typename std::enable_if< Graphics::IsDisplayMesh<MeshT>::value >::type >
{
  void addVertices(MeshT const & mesh, CollectionT & collection)
  {
    typename MeshT::VertexArray const & vertices = mesh.getVertices();
    for (size_t i = 0; i < vertices.size(); ++i)
      collection.addPoint(vertices[i]);
  }
};

} // namespace PointCollectorNInternal

/**
 * Add points from different objects to a collection. The function <tt>CollectionT::addPoint(Vector<N, ScalarT> const &)</tt>
 * (or equivalent) must exist.
 */
template <typename CollectionT, int N, typename ScalarT = Real>
class /* THEA_API */ PointCollectorN
{
  public:
    THEA_DECL_SMART_POINTERS(PointCollectorN)

    typedef Vector<N, ScalarT> Point;  ///< The point type.
    typedef CollectionT Collection;  ///< The collection type.

    /** Default constructor. Call setCollection() subsequently. */
    PointCollectorN() : collection(nullptr) {}

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
