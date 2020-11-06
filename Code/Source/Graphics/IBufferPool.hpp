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

#ifndef __Thea_Graphics_IBufferPool_hpp__
#define __Thea_Graphics_IBufferPool_hpp__

#include "../Common.hpp"
#include "../NamedObject.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
class IBuffer;

/**
 * Interface for a storage area in CPU or GPU memory that can store one or more graphics buffers.
 *
 * @see IBuffer
 */
class THEA_API IBufferPool : public virtual INamedObject
{
  public:
    /** How the storage will be used (enum class). */
    struct THEA_API Usage
    {
      /** Supported values. */
      enum Value
      {
        WRITE_EVERY_FRAME,   ///< The pool will be written to frequently.
        WRITE_OCCASIONALLY,  ///< The pool will be written to occasionally.
        WRITE_ONCE           ///< The pool will be written to at most once.
      };

      THEA_ENUM_CLASS_BODY(Usage)
    };

    /** Destructor. Automatically destroys all associated buffers. */
    virtual ~IBufferPool() = 0;

    /**
     * Destroy all buffers allocated in this pool. Don't try to dereference pointers to previously allocated buffers after
     * calling this function!
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL reset() = 0;

    /** Get the total capacity of the area in bytes. */
    virtual int64 THEA_ICALL getCapacity() const = 0;

    /** Get the number of bytes already allocated. */
    virtual int64 THEA_ICALL getAllocatedSize() const = 0;

    /** Get the number of available bytes. */
    virtual int64 THEA_ICALL getAvailableSize() const = 0;

    /** Check if the area is stored in GPU memory or not. */
    virtual int8 THEA_ICALL inGpuMemory() const = 0;

    /**
     * Create a buffer with a specified capacity in bytes. The buffer must be destroyed with destroyBuffer(). It will be
     * automatically destroyed when the IBufferPool is destroyed or reset().
     */
    virtual IBuffer * THEA_ICALL createBuffer(int64 num_bytes) = 0;

    /**
     * Destroy a buffer. The buffer must have been created with createBuffer().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL destroyBuffer(IBuffer * buf) = 0;

}; // class IBufferPool

inline IBufferPool::~IBufferPool() {}

} // namespace Graphics
} // namespace Thea

#endif
