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

#ifndef __Thea_Algorithms_RangeQueryStructure_hpp__
#define __Thea_Algorithms_RangeQueryStructure_hpp__

#include "../Common.hpp"
#include "../Array.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for a structure that supports range queries. None of the functions are virtual, this just defines a concept
 * subclasses must implement.
 */
template <typename T>
class /* THEA_API */ RangeQueryStructure
{
  public:
    THEA_DEF_POINTER_TYPES(RangeQueryStructure, shared_ptr, weak_ptr)

    /**
     * Get all objects intersecting a range.
     *
     * @param range The range to search in.
     * @param result The objects intersecting the range are stored here.
     * @param discard_prior_results If true, the contents of \a results are cleared before the range query proceeds. If false,
     *   the previous results are retained and new objects are appended to the array (this is useful for range queries over a
     *   union of simpler ranges).
     */
    template <typename IntersectionTesterT, typename RangeT>
    void rangeQuery(RangeT const & range, TheaArray<T> & result, bool discard_prior_results = true) const;

    /**
     * Get the indices of all objects intersecting a range.
     *
     * @param range The range to search in.
     * @param result The indices of objects intersecting the range are stored here.
     * @param discard_prior_results If true, the contents of \a results are cleared before the range query proceeds. If false,
     *   the previous results are retained and indices of new objects are appended to the array (this is useful for range
     *   queries over a union of simpler ranges).
     */
    template <typename IntersectionTesterT, typename RangeT>
    void rangeQueryIndices(RangeT const & range, TheaArray<long> & result, bool discard_prior_results = true) const;

    /**
     * Apply a functor to all objects in a range, until the functor returns true. The functor should provide the member function
     * (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(long index, T & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object).
     *
     * @return True if the functor evaluated to true on any object in the range (and hence stopped immediately after processing
     *   this object), else false.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    bool processRangeUntil(RangeT const & range, FunctorT * functor);

}; // class RangeQueryStructure

} // namespace Algorithms
} // namespace Thea

#endif
