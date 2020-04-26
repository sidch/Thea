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

#ifndef __Thea_Graphics_VAR_hpp__
#define __Thea_Graphics_VAR_hpp__

#include "../Common.hpp"
#include "../Colors.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Graphics {

/**
 * An interface for a Vertex Area Range object, which may be in main or GPU memory. The VAR should be constructed via the
 * functions in VARArea.
 *
 * @todo Make this safe for passing across shared library boundaries.
 */
class THEA_API VAR
{
  public:
    /** Destructor. */
    virtual ~VAR() {}

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(int64 start_elem, int64 num_elems_to_update, float32 const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(int64 start_elem, int64 num_elems_to_update, float64 const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(int64 start_elem, int64 num_elems_to_update, Vector2 const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(int64 start_elem, int64 num_elems_to_update, Vector3 const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(int64 start_elem, int64 num_elems_to_update, Vector4 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorL const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorL8 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorL16 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGB const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGB8 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGBA const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(int64 start_elem, int64 num_elems_to_update, ColorRGBA8 const * array) = 0;

    /**
     * Update a section of an index array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateIndices(int64 start_elem, int64 num_elems_to_update, uint8 const * array) = 0;

    /**
     * Update a section of an index array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateIndices(int64 start_elem, int64 num_elems_to_update, uint16 const * array) = 0;

    /**
     * Update a section of an index array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as int64 as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateIndices(int64 start_elem, int64 num_elems_to_update, uint32 const * array) = 0;

    /**
     * Clear all data in the buffer, without deallocating it. You must call this function before updating the VAR with data of
     * a different type.
     */
    virtual void clear() = 0;

    /** Get the number of elements in the VAR. */
    virtual int64 numElements() const = 0;

    /** Get the capacity of the VAR in bytes. */
    virtual int64 getCapacityInBytes() const = 0;

    /** Check if the VAR is still valid. */
    virtual int8 isValid() const = 0;

}; // class VAR

} // namespace Graphics
} // namespace Thea

#endif
