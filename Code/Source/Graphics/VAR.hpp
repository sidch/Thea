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

#ifndef __Thea_Graphics_VAR_hpp__
#define __Thea_Graphics_VAR_hpp__

#include "../Common.hpp"
#include "../Colors.hpp"
#include "../VectorN.hpp"

namespace Thea {
namespace Graphics {

/**
 * An interface for a Vertex Area Range object, which may be in main or GPU memory. The VAR should be constructed via the
 * functions in VARArea.
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
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(long start_elem, long num_elems_to_update, float const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(long start_elem, long num_elems_to_update, Vector2 const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(long start_elem, long num_elems_to_update, Vector3 const * array) = 0;

    /**
     * Update a section of a vector array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateVectors(long start_elem, long num_elems_to_update, Vector4 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorL const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorL8 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorL16 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorRGB const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorRGB8 const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorRGBA const * array) = 0;

    /**
     * Update a section of a color array. If the VAR is currently non-empty, the new and old data types must match. Call clear()
     * before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateColors(long start_elem, long num_elems_to_update, ColorRGBA8 const * array) = 0;

    /**
     * Update a section of an index array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateIndices(long start_elem, long num_elems_to_update, uint8 const * array) = 0;

    /**
     * Update a section of an index array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateIndices(long start_elem, long num_elems_to_update, uint16 const * array) = 0;

    /**
     * Update a section of an index array. If the VAR is currently non-empty, the new and old data types must match. Call
     * clear() before changing types.
     *
     * The update may increase the number of elements in the VAR, as long as it does not exceed the buffer capacity. The number
     * of elements (numElements()) is measured from the start of the VAR to the last initialized element: any uninitialized
     * portions in the middle will be included in the tally!
     *
     * @note Do <b>not</b> call this function within a RenderSystem::beginIndexedPrimitives() /
     *   RenderSystem::endIndexedPrimitives() block.
     */
    virtual void updateIndices(long start_elem, long num_elems_to_update, uint32 const * array) = 0;

    /**
     * Clear all data in the buffer, without deallocating it. You must call this function before updating the VAR with data of
     * a different type.
     */
    virtual void clear() = 0;

    /** Get the number of elements in the VAR. */
    virtual long numElements() const = 0;

    /** Get the capacity of the VAR in bytes. */
    virtual long getCapacityInBytes() const = 0;

    /** Check if the VAR is still valid. */
    virtual bool isValid() const = 0;

}; // class VAR

} // namespace Graphics
} // namespace Thea

#endif
