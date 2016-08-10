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
    virtual long getCapacity() const = 0;

    /** Get the number of bytes already allocated. */
    virtual long getAllocatedSize() const = 0;

    /** Get the number of available bytes. */
    virtual long getAvailableSize() const { return getCapacity() - getAllocatedSize(); }

    /** Check if the area is stored in GPU memory or not. */
    virtual bool inGPUMemory() const = 0;

    /**
     * Create a VAR with a specified capacity in bytes. The VAR must be destroyed with destroyArray(). It will be
     * automatically destroyed when the VARArea is destroyed or reset().
     */
    virtual VAR * createArray(long num_bytes) = 0;

    /** Destroy a VAR. The VAR must have been created with createArray(). */
    virtual void destroyArray(VAR * array) = 0;

}; // class VARArea

} // namespace Graphics
} // namespace Thea

#endif
