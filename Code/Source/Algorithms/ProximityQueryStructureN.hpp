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

#ifndef __Thea_Algorithms_ProximityQueryStructureN_hpp__
#define __Thea_Algorithms_ProximityQueryStructureN_hpp__

#include "../Common.hpp"
#include "Filter.hpp"
#include "../AxisAlignedBoxN.hpp"
#include "../MatVec.hpp"
#include <cmath>
#include <sstream>
#include <utility>

namespace Thea {
namespace Algorithms {

/**
 * Interface for a structure that supports proximity queries in N-space. None of the functions are virtual, this just defines a
 * concept subclasses must implement.
 */
template <int N, typename T = Real>
class /* THEA_API */ ProximityQueryStructureN
{
  private:

  public:
    THEA_DECL_SMART_POINTERS(ProximityQueryStructureN)

    typedef Vector<N, T> VectorT;  ///< N-dimensional vector.

    /**
     * A return value of a k-nearest-neighbors query, specified by a monotone approximation to (for L2, square of) the distance
     * between neighboring elements, the indices of the query and target elements, and the positions of the closest pair of
     * points on them.
     */
    class /* THEA_API */ NeighborPair
    {
      public:
        /**
         * Construct from a pair of query and target indices, and a monotone approximation to the distance between the
         * elements. A negative value for any field implies it has not been initialized.
         */
        NeighborPair(intx query_index_ = -1, intx target_index_ = -1, double mon_approx_dist_ = -1)
        : query_index(query_index_), target_index(target_index_), mon_approx_dist(mon_approx_dist_)
        {}

        /**
         * Construct from a pair of query and target indices, a monotone approximation to the distance between the elements, and
         * the positions of the closest pair of points on them. A negative value for any of the scalar fields implies it has not
         * been initialized.
         */
        NeighborPair(intx query_index_, intx target_index_, double mon_approx_dist_, VectorT const & query_point_,
                     VectorT const & target_point_)
        : query_index(query_index_), target_index(target_index_), mon_approx_dist(mon_approx_dist_), query_point(query_point_),
          target_point(target_point_)
        {}

        /** Check if the pair has valid indices. */
        bool isValid() const { return query_index >= 0 && target_index >= 0; }

        /** Get the index of the query element. */
        intx getQueryIndex() const { return query_index; }

        /** Set the index of the query element. */
        void setQueryIndex(intx query_index_) { query_index = query_index_; }

        /** Get the index of the target element. */
        intx getTargetIndex() const { return target_index; }

        /** Set the index of the target element. */
        void setTargetIndex(intx target_index_) { target_index = target_index_; }

        /**
         * Get a monotone approximation of the distance between the neighbors. A negative value indicates the field has not been
         * initialized.
         */
        double getMonotoneApproxDistance() const { return mon_approx_dist; }

        /** Set the monotone approximation of the distance between the neighbors. */
        void setMonotoneApproxDistance(double mon_approx_dist_) { mon_approx_dist = mon_approx_dist_; }

        /**
         * Get the distance between the neighbors, assuming the stored monotone approximation was computed using the supplied
         * metric. The return value is undefined if the monotone approximation has an invalid value (e.g. negative).
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
        intx query_index;
        intx target_index;
        double mon_approx_dist;
        VectorT query_point;
        VectorT target_point;

    }; // class NeighborPair

    /** Get a bounding box for the structure. */
    AxisAlignedBox3 const & getBounds() const;

    /**
     * Get the minimum distance between this structure and a query object.
     *
     * @param query Query object.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     */
    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    double distance(QueryT const & query, double dist_bound = -1, CompatibilityFunctorT compatibility = CompatibilityFunctorT())
           const;

    /**
     * Get the closest element in this structure to a query object, within a specified distance bound.
     *
     * @param query Query object.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param dist The distance to the query object is placed here. Ignored if null.
     * @param closest_point The coordinates of the closest point are placed here. Ignored if null.
     *
     * @return A non-negative handle to the closest element, if one was found, else a negative number.
     */
    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    intx closestElement(QueryT const & query, double dist_bound = -1,
                        CompatibilityFunctorT compatibility = CompatibilityFunctorT(),
                        double * dist = nullptr, VectorT * closest_point = nullptr) const;

    /**
     * Get the closest pair of elements between this structure and another structure, whose separation is less than a specified
     * upper bound.
     *
     * @param query Query object.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param get_closest_points If true, the coordinates of the closest pair of points on the respective elements is computed
     *   and stored in the returned structure.
     *
     * @return Non-negative handles to the closest pair of elements in their respective objects, if such a pair was found. Else
     *   returns a pair of negative numbers.
     */
    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    NeighborPair closestPair(QueryT const & query, double dist_bound = -1,
                             CompatibilityFunctorT compatibility = CompatibilityFunctorT(),
                             bool get_closest_points = false) const;

    /**
     * Get the k elements closest to a query object. The returned elements are placed in a set of bounded size (k). The template
     * type BoundedNeighborPairSetT should typically be BoundedSortedArray<NeighborPair> or BoundedSortedArrayN<k, NeighborPair>
     * if only a few neighbors are requested.
     *
     * @param query Query object.
     * @param k_closest_pairs The k (or fewer) nearest neighbors are placed here.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param get_closest_points If true, the coordinates of the closest pair of points on each pair of neighboring elements is
     *   computed and stored in the returned pairs.
     *
     * @return The number of neighbors found (i.e. the size of \a k_closest_pairs).
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSetT,
              typename CompatibilityFunctorT = UniversalCompatibility>
    intx kClosestPairs(QueryT const & query, BoundedNeighborPairSetT & k_closest_pairs, double dist_bound = -1,
                       CompatibilityFunctorT compatibility = CompatibilityFunctorT(), bool get_closest_points = false) const;

    /**
     * Apply a functor to all elements neighboring a query object, until the functor returns true. The functor should provide
     * the following two member functions:
     * \code
     * bool allows(NeighborPair const &) const;
     * bool operator()(NeighborPair const &)
     * \endcode
     * The first function will be used to check if a potential pair of neighbors will be accepted by the functor or not. Note
     * that this pair may be <i>incomplete</i>, e.g. one or both source/target element indices, or the element separation, could
     * be negative. Any such invalid fields will be ignored. This feature is used for early pruning of elements.
     *
     * The second function will be passed information about every valid pair of elements that passes the <tt>allows()</tt>
     * function. If the query is itself a proximity query structure, the corresponding field of the NeighborPair will point to
     * the relevant element in the structure. Otherwise, the corresponding field will always contain the index 0 and refer to
     * the entire query object. This function will always receive <i>complete</i> NeighborPair objects (no negative indices or
     * distances), though the closest point positions will be uninitialized if \a get_closest_points is false.
     *
     * This is a generic function that can be used to mimic the behavior of closestPair() or kClosestPairs(), though those
     * functions can be individually optimized in specific implementations.
     *
     * @note The same pair of elements may be found twice and passed to <tt>FunctorT::operator()</tt>.
     *   <tt>FunctorT::allows()</tt> should check for repetitions if it wants to avoid this.
     *
     * @param query Query object.
     * @param functor The functor that will be called for each pair of neighbors.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param get_closest_points If true, the coordinates of the closest pair of points on each pair of neighboring elements is
     *   computed and passed to the functor.
     * @param query_index The supplied index is passed to the functor as the index of the query object. This is chiefly used for
     *   internal processing and the default value of 0 should normally be left as is.
     * @param swap_query_and_target If true, the indices and other properties of neighboring query and target objects are
     *   swapped when they are passed to the functor. This is chiefly used for internal processing and the default value of
     *   false should normally be left as is.
     *
     * If the functor returns true on any object, the search will terminate immediately (this is useful for searching for a
     * particular pair). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     */
    template <typename MetricT, typename QueryT, typename FunctorT, typename CompatibilityFunctorT = UniversalCompatibility>
    void processNeighbors(QueryT const & query, FunctorT functor, CompatibilityFunctorT compatibility = CompatibilityFunctorT(),
                          bool get_closest_points = false, intx query_index = 0, bool swap_query_and_target = false) const;

}; // class ProximityQueryStructureN

/** Interface for a structure that supports proximity queries in 2-space. */
typedef ProximityQueryStructureN<2, Real> ProximityQueryStructure2;

/** Interface for a structure that supports proximity queries in 3-space. */
typedef ProximityQueryStructureN<3, Real> ProximityQueryStructure3;

} // namespace Algorithms
} // namespace Thea

#endif
