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

#ifndef __Thea_Graphics_VARArea_hpp__
#define __Thea_Graphics_VARArea_hpp__

#include "../Common.hpp"
#include "../NamedObject.hpp"
#include "VAR.hpp"

namespace Thea {
namespace Graphics {

/**
 * An interface for a Vertex Area Range storage area.
 *
 * @see VAR
 *
 * @todo Make this safe for passing across shared library boundaries.
 */
class THEA_API VARArea : public AbstractNamedObject
{
  public:
    /** How the storage will be used (enum class). */
    struct THEA_API Usage
    {
      /** Supported values. */
      enum Value
      {
        WRITE_EVERY_FRAME,   ///< The buffer will be written to frequently.
        WRITE_OCCASIONALLY,  ///< The buffer will be written to occasionally.
        WRITE_ONCE           ///< The buffer will be written to at most once.
      };

      THEA_ENUM_CLASS_BODY(Usage)
    };

    /** Destructor. Automatically destroys all associated vertex arrays. */
    virtual ~VARArea() {}

    /**
     * Destroy all vertex arrays allocated in this block. Don't try to dereference pointers to previously allocated vertex
     * arrays (VARs) after calling this function!
     */
    virtual void reset() = 0;

    /** Get the total capacity of the area in bytes. */
    virtual int64 getCapacity() const = 0;

    /** Get the number of bytes already allocated. */
    virtual int64 getAllocatedSize() const = 0;

    /** Get the number of available bytes. */
    virtual int64 getAvailableSize() const { return getCapacity() - getAllocatedSize(); }

    /** Check if the area is stored in GPU memory or not. */
    virtual int8 inGPUMemory() const = 0;

    /**
     * Create a VAR with a specified capacity in bytes. The VAR must be destroyed with destroyArray(). It will be
     * automatically destroyed when the VARArea is destroyed or reset().
     */
    virtual VAR * createArray(int64 num_bytes) = 0;

    /** Destroy a VAR. The VAR must have been created with createArray(). */
    virtual void destroyArray(VAR * array) = 0;

}; // class VARArea

} // namespace Graphics
} // namespace Thea

#endif
