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
    THEA_DECL_SMART_POINTERS(RangeQueryStructure)

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
    void rangeQuery(RangeT const & range, Array<T> & result, bool discard_prior_results = true) const;

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
    void rangeQueryIndices(RangeT const & range, Array<intx> & result, bool discard_prior_results = true) const;

    /**
     * Apply a functor to all objects in a range, until the functor returns true. The functor should provide the member function
     * (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T const & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * @return The index of the first object in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this object), else a negative value.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    intx processRangeUntil(RangeT const & range, FunctorT functor) const;

    /**
     * Apply a functor to all objects in a range, until the functor returns true. The functor should provide the member function
     * (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T [const] & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * @return The index of the first object in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this object), else a negative value.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    intx processRangeUntil(RangeT const & range, FunctorT functor);

}; // class RangeQueryStructure

} // namespace Algorithms
} // namespace Thea

#endif
