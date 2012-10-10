//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_GraphType_hpp__
#define __Thea_GraphType_hpp__

#include "Common.hpp"
#include "Concept.hpp"

namespace Thea {

/**
 * Checks if a class is a graph that allows iterating over all vertices. The class T must define the following types:
 *
 * \code
 *   Vertex               // Vertex of the graph.
 *   VertexIterator       // Iterator over vertices.
 *   VertexConstIterator  // Const iterator over vertices.
 * \endcode
 * and implement the following functions:
 * type.
 *
 * \code
 *   long numVertices() const;
 *   Vertex[Const]Iterator verticesBegin() [const];
 *   Vertex[Const]Iterator verticesEnd() [const];
 *   Vertex [const] & deref(Vertex[Const]Iterator) [const];
 * \endcode
 */
template <typename T>
class IsGraph
{
  private:
    THEA_HAS_TYPE(HasVertex,               Vertex)
    THEA_HAS_TYPE(HasVertexIterator,       VertexIterator)
    THEA_HAS_TYPE(HasVertexConstIterator,  VertexConstIterator)

    THEA_HAS_MEMBER(HasDeref,          deref)
    THEA_HAS_MEMBER(HasNumVertices,    numVertices)
    THEA_HAS_MEMBER(HasVerticesBegin,  verticesBegin)
    THEA_HAS_MEMBER(HasVerticesEnd,    verticesEnd)

  public:
    static bool const value = HasVertex              <T>::value
                           && HasVertexIterator      <T>::value
                           && HasVertexConstIterator <T>::value
                           && HasDeref               <T>::value
                           && HasNumVertices         <T>::value
                           && HasVerticesBegin       <T>::value
                           && HasVerticesEnd         <T>::value;

}; // class IsGraph

/**
 * Checks if a class is a graph that, in addition to satisfying IsGraph<T>, allows iterating over the neighbors of a vertex. The
 * class T must define the following types:
 *
 * \code
 *   NeighborIterator       // Iterator over neighbors of a vertex.
 *   NeighborConstIterator  // Const iterator over neighbors of a vertex.
 * \endcode
 * and implement the following functions:
 * type.
 *
 * \code
 *   long numNeighbors(VertexConstIterator vertex) [const];
 *   Neighbor[Const]Iterator neighborsBegin(VertexConstIterator vertex) [const];
 *   Neighbor[Const]Iterator neighborsEnd(VertexConstIterator vertex) [const];
 *   [numerical-type] distance(VertexConstIterator, NeighborConstIterator) const;
 *   Vertex [const] & deref(Neighbor[Const]Iterator) [const];
 * \endcode
 */
template <typename T>
class IsAdjacencyGraph
{
  private:
    THEA_HAS_TYPE(HasNeighborIterator,       NeighborIterator)
    THEA_HAS_TYPE(HasNeighborConstIterator,  NeighborConstIterator)

    THEA_HAS_MEMBER(HasDeref,           deref)
    THEA_HAS_MEMBER(HasNumNeighbors,    numNeighbors)
    THEA_HAS_MEMBER(HasNeighborsBegin,  neighborsBegin)
    THEA_HAS_MEMBER(HasNeighborsEnd,    neighborsEnd)
    THEA_HAS_MEMBER(HasDistance,        distance)

  public:
    static bool const value = IsGraph<T>::value
                           && HasNeighborIterator      <T>::value
                           && HasNeighborConstIterator <T>::value
                           && HasDeref                 <T>::value
                           && HasNumNeighbors          <T>::value
                           && HasNeighborsBegin        <T>::value
                           && HasNeighborsEnd          <T>::value;
                           && HasDistance              <T>::value;

}; // class IsAdjacencyGraph

} // namespace Thea

#endif
