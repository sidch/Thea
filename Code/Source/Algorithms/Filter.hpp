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

#ifndef __Thea_Algorithms_Filter_hpp__
#define __Thea_Algorithms_Filter_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for all object filters. A filter has an allows() function that takes an object argument and returns a boolean value
 * indicating whether the filter allows the object to pass through or not.
 */
template <typename T>
class Filter
{
  public:
    THEA_DEF_POINTER_TYPES(Filter, shared_ptr, weak_ptr)

    /** Destructor. */
    virtual ~Filter() {}

    /**
     * Check if the filter allows an object to pass through or not.
     *
     * @return True if the object \a t is allowed through, false if it is blocked.
     */
    virtual bool allows(T const & t) const = 0;

}; // class Filter

/** A filter that allows everything to pass through. */
template <typename T>
class AlwaysPassFilter : public Filter<T>
{
  public:
    bool allows(T const & t) const { return true; }

}; // class AlwaysPassFilter

/** A filter that allows nothing to pass through. */
template <typename T>
class AlwaysBlockFilter : public Filter<T>
{
  public:
    bool allows(T const & t) const { return false; }

}; // class AlwaysPassFilter

} // namespace Algorithms
} // namespace Thea

#endif
