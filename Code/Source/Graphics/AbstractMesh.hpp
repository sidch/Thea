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
// First version: 2019
//
//============================================================================

#ifndef __Thea_Graphics_AbstractMesh_hpp__
#define __Thea_Graphics_AbstractMesh_hpp__

#include "../Common.hpp"
#include "../AbstractDenseMatrix.hpp"
#include "../NamedObject.hpp"
#include "Drawable.hpp"

namespace Thea {
namespace Graphics {

/** A minimalistic abstract base class for passing polygon mesh data across shared library boundaries. */
class AbstractMesh : public virtual AbstractNamedObject, public virtual Drawable
{
  public:
    THEA_DECL_SMART_POINTERS(AbstractMesh)

    /** Destructor. */
    virtual ~AbstractMesh() = 0;

    /**
     * Get the vertices of the mesh, as a dense 3xN column-major matrix. Each column of the matrix is the (x, y, z) coordinate
     * triplet of a vertex. The function is guaranteed to return a non-null matrix, even if it has 0 columns.
     *
     * @note This function is safe to call across a shared library boundary only if the Real type has a consistent size on both
     *   sides of the boundary.
     */
    virtual AbstractDenseMatrix<Real> const * getVertexMatrix() const = 0;

    /**
     * Get the triangles of the mesh, as a dense 3xN column-major matrix. Each column of the matrix is the set of indices of a
     * triangle's vertices. The function is guaranteed to return a non-null matrix, even if it has 0 columns.
     */
    virtual AbstractDenseMatrix<uint32> const * getTriangleMatrix() const = 0;

    /**
     * Get the quads of the mesh, as a dense 4xN column-major matrix. Each column of the matrix is the set of indices of a
     * quad's vertices. The function is guaranteed to return a non-null matrix, even if it has 0 columns.
     */
    virtual AbstractDenseMatrix<uint32> const * getQuadMatrix() const = 0;

}; // class AbstractMesh

inline AbstractMesh::~AbstractMesh() {}

} // namespace Graphics
} // namespace Thea

#endif
