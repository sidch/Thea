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

#ifndef __Thea_Graphics_EdgeWelder_hpp__
#define __Thea_Graphics_EdgeWelder_hpp__

#include "../Common.hpp"
#include "../MatVec.hpp"
#include "../Noncopyable.hpp"

namespace Thea {
namespace Graphics {

// Forward declaration
class EdgeWelderImpl;

/**
 * Maintains a set of edges, with associated positions, without duplication. Two edges are considered identical if both pairs of
 * corresponding endpoints are approximately co-located.
 */
class THEA_API EdgeWelder : private Noncopyable
{
  public:
    /**
     * Constructor. Two edges are considered identical if both pairs of corresponding endpoints are separated by at most the
     * welding radius.
     */
    EdgeWelder(Real weld_radius);

    /** Destructor. */
    ~EdgeWelder();

    /** Add an edge to the set, unless a coincident edge already exists. */
    void addEdge(void * edge, Vector3 const & e0, Vector3 const & e1);

    /** If a directed edge with the specified endpoints exists, return it, else return null. */
    void * getDirectedEdge(Vector3 const & e0, Vector3 const & e1) const;

    /** If an undirected edge with the specified endpoints exists, return it, else return null. */
    void * getUndirectedEdge(Vector3 const & e0, Vector3 const & e1) const;

  private:
    EdgeWelderImpl * impl;

}; // class EdgeWelder

} // namespace Graphics
} // namespace Thea

#endif
