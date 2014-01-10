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

#ifndef __Thea_Algorithms_ProximityQueryStructureN_hpp__
#define __Thea_Algorithms_ProximityQueryStructureN_hpp__

#include "../Common.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../VectorN.hpp"
#include <cmath>
#include <sstream>
#include <utility>

namespace Thea {
namespace Algorithms {

/**
 * Interface for a structure that supports proximity queries in N-space. None of the functions are virtual, this just defines a
 * concept subclasses must implement.
 */
template <long N, typename T = Real>
class /* THEA_API */ ProximityQueryStructureN
{
  public:
    THEA_DEF_POINTER_TYPES(ProximityQueryStructureN, shared_ptr, weak_ptr)

    typedef VectorN<N, T> VectorT;  ///< N-dimensional vector.

    /**
     * A return value of a k-nearest-neighbors query, specified by a monotone approximation to (for L2, square of) the distance
     * between neighboring elements, the indices of the query and target elements, and the positions of the closest pair of
     * points on them.
     */
    class /* THEA_API */ NeighborPair
    {
      public:
        /** Default constructor. */
        NeighborPair() {}

        /**
         * Construct from a pair of query and target indices, and a monotone approximation to the distance between the
         * elements.
         */
        NeighborPair(long query_index_, long target_index_ = -1, double mon_approx_dist_ = 0)
        : query_index(query_index_), target_index(target_index_), mon_approx_dist(mon_approx_dist_),
          query_point(VectorT::zero()), target_point(VectorT::zero())
        {}

        /**
         * Construct from a pair of query and target indices, a monotone approximation to the distance between the elements, and
         * the positions of the closest pair of points on them.
         */
        NeighborPair(long query_index_, long target_index_, double mon_approx_dist_, VectorT const & query_point_,
                     VectorT const & target_point_)
        : query_index(query_index_), target_index(target_index_), mon_approx_dist(mon_approx_dist_), query_point(query_point_),
          target_point(target_point_)
        {}

        /** Check if the pair has valid indices. */
        bool isValid() const { return query_index >= 0 && target_index >= 0; }

        /** Get the index of the query element. */
        long getQueryIndex() const { return query_index; }

        /** Set the index of the query element. */
        void setQueryIndex(long query_index_) { query_index = query_index_; }

        /** Get the index of the target element. */
        long getTargetIndex() const { return target_index; }

        /** Set the index of the target element. */
        void setTargetIndex(long target_index_) { target_index = target_index_; }

        /** Get a monotone approximation of the distance between the neighbors. */
        double getMonotoneApproxDistance() const { return mon_approx_dist; }

        /** Set the monotone approximation of the distance between the neighbors. */
        void setMonotoneApproxDistance(double mon_approx_dist_) { mon_approx_dist = mon_approx_dist_; }

        /**
         * Get the distance between the neighbors, assuming the stored monotone approximation was computed using the supplied
         * metric.
         */
        template <typename MetricT> double getDistance() const { return MetricT::invertMonotoneApprox(mon_approx_dist); }

        /** Get the point of the query element. */
        VectorT const & getQueryPoint() const { return query_point; }

        /** Set the point of the query element. */
        void setQueryPoint(VectorT const & query_point_) { query_point = query_point_; }

        /** Get the point of the target element. */
        VectorT const & getTargetPoint() const { return target_point; }

        /** Set the point of the target element. */
        void setTargetPoint(VectorT const & target_point_) { target_point = target_point_; }

        /** Returns a pair with the query and target indices swapped. */
        NeighborPair swapped() const
        { return NeighborPair(target_index, query_index, mon_approx_dist, target_point, query_point); }

        /** Less-than comparator, sorts neighbors by increasing distance. Uses the indices to break ties. */
        bool operator<(NeighborPair const & other) const
        {
          return mon_approx_dist < other.mon_approx_dist
              || (mon_approx_dist == other.mon_approx_dist
               && (query_index < other.query_index
                || (query_index == other.query_index && target_index == other.target_index)));
        }

        /** Equality comparator. Only compares indices. */
        bool operator==(NeighborPair const & other) const
        { return query_index == other.query_index && target_index == other.target_index; }

        /** Get a string representation of the neighboring pair. */
        std::string toString() const
        {
          std::ostringstream oss;
          oss << "NeighborPair(q: " << query_index << ", t: " << target_index << ", mad: " << mon_approx_dist << ')';
          return oss.str();
        }

      private:
        long query_index;
        long target_index;
        double mon_approx_dist;
        VectorT query_point;
        VectorT target_point;

    }; // class NeighborPair

    /** Get a bounding box for the structure. */
    AxisAlignedBox3 const & getBounds() const;

    /** Get the minimum distance between this structure and a query object. */
    template <typename MetricT, typename QueryT> double distance(QueryT const & query, double dist_bound = -1) const;

    /**
     * Get the closest element in this structure to a query object, within a specified distance bound.
     *
     * @param query Query object.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param dist The distance to the query object is placed here. Ignored if null.
     * @param closest_point The coordinates of the closest point are placed here. Ignored if null.
     *
     * @return A non-negative handle to the closest element, if one was found, else a negative number.
     */
    template <typename MetricT, typename QueryT>
    long closestElement(QueryT const & query, double dist_bound = -1, double * dist = NULL, VectorT * closest_point = NULL)
         const;

    /**
     * Get the closest pair of elements between this structure and another structure, whose separation is less than a specified
     * upper bound.
     *
     * @param query Query object.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param get_closest_points If true, the coordinates of the closest pair of points on the respective elements is computed
     *   and stored in the returned structure.
     *
     * @return Non-negative handles to the closest pair of elements in their respective objects, if such a pair was found. Else
     *   returns a pair of negative numbers.
     */
    template <typename MetricT, typename QueryT>
    NeighborPair closestPair(QueryT const & query, double dist_bound = -1, bool get_closest_points = false) const;

    /**
     * Get the k elements closest to a query object. The returned elements are placed in a set of bounded size (k). The template
     * type BoundedNeighborPairSetT should typically be BoundedSortedArray<NeighborPair> or BoundedSortedArrayN<k, NeighborPair>
     * if only a few neighbors are requested.
     *
     * @param query Query object.
     * @param k_closest_pairs The k (or fewer) nearest neighbors are placed here.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param get_closest_points If true, the coordinates of the closest pair of points on each pair of neighboring elements is
     *   computed and stored in the returned pairs.
     * @param clear_set If true (default), this function discards prior data in \a k_closest_pairs. This is chiefly for internal
     *   use and the default value of true should normally be left as is.
     * @param use_as_query_index_and_swap If non-negative, the supplied index is used as the index of the query object (instead
     *   of the default 0), following which query and target indices/points are swapped in the returned pairs of neighbors. This
     *   is chiefly for internal use and the default value of -1 should normally be left as is.
     *
     * @return The number of neighbors found (i.e. the size of \a k_closest_pairs).
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSetT>
    long kClosestPairs(QueryT const & query, BoundedNeighborPairSetT & k_closest_pairs, double dist_bound = -1,
                       bool get_closest_points = false, bool clear_set = true, long use_as_query_index_and_swap = -1) const;

}; // class ProximityQueryStructureN

/** Interface for a structure that supports proximity queries in 2-space. */
typedef ProximityQueryStructureN<2, Real> ProximityQueryStructure2;

/** Interface for a structure that supports proximity queries in 3-space. */
typedef ProximityQueryStructureN<3, Real> ProximityQueryStructure3;

} // namespace Algorithms
} // namespace Thea

#endif
