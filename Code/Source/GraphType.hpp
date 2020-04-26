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
// First version: 2012
//
//============================================================================

#ifndef __Thea_GraphType_hpp__
#define __Thea_GraphType_hpp__

#include "Common.hpp"
#include "Concept.hpp"

namespace Thea {

/**
 * Checks if a class is a graph that allows iterating over all vertices. The class T must define the following types:
 *
 * \code
 *   VertexHandle         // Handle to a vertex of the graph.
 *   VertexConstHandle    // Const handle to a vertex of the graph.
 *   VertexIterator       // Iterator over vertices.
 *   VertexConstIterator  // Const iterator over vertices.
 * \endcode
 *
 * and implement the following functions:
 *
 * \code
 *   intx numVertices() const;
 *   Vertex[Const]Iterator verticesBegin() [const];
 *   Vertex[Const]Iterator verticesEnd() [const];
 *   Vertex[Const]Handle getVertex(Vertex[Const]Iterator) [const];
 * \endcode
 */
template <typename T>
class IsGraph
{
  private:
    THEA_HAS_TYPE(HasVertexHandle,         VertexHandle)
    THEA_HAS_TYPE(HasVertexConstHandle,    VertexConstHandle)
    THEA_HAS_TYPE(HasVertexIterator,       VertexIterator)
    THEA_HAS_TYPE(HasVertexConstIterator,  VertexConstIterator)

    THEA_HAS_MEMBER(HasNumVertices,    numVertices)
    THEA_HAS_MEMBER(HasVerticesBegin,  verticesBegin)
    THEA_HAS_MEMBER(HasVerticesEnd,    verticesEnd)
    THEA_HAS_MEMBER(HasGetVertex,      getVertex)

  public:
    static bool const value = HasVertexHandle        <T>::value
                           && HasVertexConstHandle   <T>::value
                           && HasVertexIterator      <T>::value
                           && HasVertexConstIterator <T>::value
                           && HasNumVertices         <T>::value
                           && HasVerticesBegin       <T>::value
                           && HasVerticesEnd         <T>::value
                           && HasGetVertex           <T>::value;

}; // class IsGraph

/**
 * Checks if a class is a graph that, in addition to satisfying IsGraph, allows iterating over the neighbors of a vertex. The
 * class T must define the following types:
 *
 * \code
 *   NeighborIterator       // Iterator over neighbors of a vertex.
 *   NeighborConstIterator  // Const iterator over neighbors of a vertex.
 * \endcode
 *
 * and implement the following functions:
 *
 * \code
 *   intx numNeighbors(VertexConstHandle vertex) const;
 *   Neighbor[Const]Iterator neighborsBegin(Vertex[Const]Handle vertex) [const];
 *   Neighbor[Const]Iterator neighborsEnd(Vertex[Const]Handle vertex) [const];
 *   Vertex[Const]Handle getVertex(Neighbor[Const]Iterator) [const];
 *   [numerical-type] distance(VertexConstHandle, NeighborConstIterator) const;
 * \endcode
 */
template <typename T>
class IsAdjacencyGraph
{
  private:
    THEA_HAS_TYPE(HasNeighborIterator,       NeighborIterator)
    THEA_HAS_TYPE(HasNeighborConstIterator,  NeighborConstIterator)

    THEA_HAS_MEMBER(HasNumNeighbors,    numNeighbors)
    THEA_HAS_MEMBER(HasNeighborsBegin,  neighborsBegin)
    THEA_HAS_MEMBER(HasNeighborsEnd,    neighborsEnd)
    THEA_HAS_MEMBER(HasGetVertex,       getVertex)
    THEA_HAS_MEMBER(HasDistance,        distance)

  public:
    static bool const value = IsGraph<T>::value
                           && HasNeighborIterator      <T>::value
                           && HasNeighborConstIterator <T>::value
                           && HasNumNeighbors          <T>::value
                           && HasNeighborsBegin        <T>::value
                           && HasNeighborsEnd          <T>::value
                           && HasGetVertex             <T>::value
                           && HasDistance              <T>::value;

}; // class IsAdjacencyGraph

} // namespace Thea

#endif
