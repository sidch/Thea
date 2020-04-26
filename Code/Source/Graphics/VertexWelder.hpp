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

#ifndef __Thea_Graphics_VertexWelder_hpp__
#define __Thea_Graphics_VertexWelder_hpp__

#include "../Common.hpp"
#include "../Noncopyable.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Graphics {

// Forward declaration
class VertexWelderImpl;

/**
 * Maintains a set of vertices, with associated positions, without duplication. Two vertices are considered identical if their
 * positions are approximately the same.
 */
class THEA_API VertexWelder : private Noncopyable
{
  public:
    /** Constructor. Two vertices are considered identical if their positions are separated by at most the welding radius. */
    VertexWelder(Real weld_radius);

    /** Destructor. */
    ~VertexWelder();

    /** Add a vertex to the set, unless a coincident vertex already exists. */
    void addVertex(void * vertex, Vector3 const & position);

    /** If a vertex exists at the specified position, return it, else return null. */
    void * getVertex(Vector3 const & position) const;

  private:
    VertexWelderImpl * impl;

}; // class VertexWelder

} // namespace Graphics
} // namespace Thea

#endif
