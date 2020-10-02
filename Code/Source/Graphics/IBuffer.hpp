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

#ifndef __Thea_Graphics_IBuffer_hpp__
#define __Thea_Graphics_IBuffer_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Graphics {

/**
 * Interface for a graphics buffer, which is in CPU or GPU memory as part of an IBufferPool. The buffer should be constructed
 * via the functions in IBufferPool.
 */
class THEA_API IBuffer
{
  public:
    /** Destructor. */
    virtual ~IBuffer() = 0;

    /** Get the dimensionality of each value in the buffer. */
    virtual int32 THEA_ICALL numComponents() const = 0;

    /** Get the type of each component of a value, as an entry from the NumericType enum. */
    virtual int32 THEA_ICALL getComponentType() const = 0;

    /** Get the number of values in the buffer. */
    virtual int64 THEA_ICALL numValues() const = 0;

    /** Get the capacity of the buffer in bytes. */
    virtual int64 THEA_ICALL getCapacityInBytes() const = 0;

    /** Check if the buffer is still valid. */
    virtual int8 THEA_ICALL isValid() const = 0;

    /**
     * Clear all data in the buffer, without deallocating it. You must call this function before updating the buffer with data
     * of a different type.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL clear() = 0;

    /**
     * Update a section of the buffer containing an array of attributes such vertex or normal coordinates. For an array of
     * indices, use updateIndices() instead. If the buffer is currently non-empty, the new and old data types and component
     * counts must match. Call clear() before changing types.
     *
     * The update may increase the number of values in the buffer, as long as it does not exceed the buffer capacity. The
     * number of values (numValues()) is measured from the start of the buffer to the last initialized value: any
     * uninitialized portions in the middle will be included in the tally!
     *
     * @param start_value The index of the starting value of the buffer section to update.
     * @param num_values_to_update The number of values to copy.
     * @param num_components The dimensionality of each value.
     * @param type The type of each component of a value, as an enum value from NumericType.
     * @param src The source buffer to copy from.
     *
     * @note Do <b>not</b> call this function within a IRenderSystem::beginIndexedPrimitives() /
     *   IRenderSystem::endIndexedPrimitives() block.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL updateAttributes(int64 start_value, int64 num_values_to_update, int32 num_components, int32 type,
                                  void const * src) = 0;

    /**
     * Update a section of the buffer containing an array of value indices. For an array of attributes such as vertex or
     * normal coordinates, use updateAttributes() instead. Usage is otherwise identical to updateAttributes().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL updateIndices(int64 start_value, int64 num_values_to_update, int32 type, void const * src) = 0;

}; // class IBuffer

inline IBuffer::~IBuffer() {}

} // namespace Graphics
} // namespace Thea

#endif
