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

#ifndef __Thea_Algorithms_RayQueryStructureN_hpp__
#define __Thea_Algorithms_RayQueryStructureN_hpp__

#include "../Common.hpp"
#include "../RayIntersectableN.hpp"

namespace Thea {
namespace Algorithms {

/**
 * A description of the intersection point of a ray with a structure. Specifies the hit time, the normal at the intersection
 * point, and the index of the intersected element.
 */
template <int N, typename T = Real>
class /* THEA_API */ RayStructureIntersectionN : public RayIntersectionN<N, T>
{
  private:
    typedef RayIntersectionN<N, T> BaseT;

    intx element_index;

  public:
    /** Constructor. */
    RayStructureIntersectionN(Real time_ = -1, typename BaseT::VectorT const * normal_ = nullptr, intx element_index_ = -1)
    : BaseT(time_, normal_), element_index(element_index_)
    {}

    /** Constructor. */
    RayStructureIntersectionN(BaseT const & isec, intx element_index_ = -1)
    : BaseT(isec), element_index(element_index_)
    {}

    /** Get the index of the intersected element. */
    intx getElementIndex() const { return element_index; }

    /** Set the index of the intersected element. */
    void setElementIndex(intx element_index_) { element_index = element_index_; }

}; // class RayStructureIntersectionN

/**
 * Interface for a structure that supports ray intersection queries in N-space. None of the functions are virtual, this just
 * defines a concept subclasses must implement.
 */
template <int N, typename T = Real>
class /* THEA_API */ RayQueryStructureN
{
  public:
    THEA_DECL_SMART_POINTERS(RayQueryStructureN)

    typedef RayN<N, T> RayT;
    typedef RayStructureIntersectionN<N, T> RayStructureIntersectionT;

    /** Check if a ray intersects the structure in the forward direction. */
    template <typename RayIntersectionTesterT> bool rayIntersects(RayT const & ray, Real max_time = -1) const;

    /**
     * Get the time taken for a ray to intersect the structure, or a negative value if there was no intersection in the forward
     * direction.
     */
    template <typename RayIntersectionTesterT> Real rayIntersectionTime(RayT const & ray, Real max_time = -1) const;

    /**
     * Get the intersection of a ray with the structure, including the hit time, the normal at the intersection point, and the
     * index of the intersected element. A negative time is returned if there was no intersection in the forward direction. A
     * zero normal and a negative index are returned if those quantities are not known.
     */
    template <typename RayIntersectionTesterT>
    RayStructureIntersectionT rayStructureIntersection(RayT const & ray, Real max_time = -1) const;

}; // class RayQueryStructureN

/**
 * A description of the intersection point of a ray with a structure in 2-space. Specifies the hit time, the normal at the
 * intersection point, and the index of the intersected element.
 */
typedef RayStructureIntersectionN<2, Real> RayStructureIntersection2;

/** Interface for a structure that supports ray intersection queries in 2-space. */
typedef RayQueryStructureN<2, Real> RayQueryStructure2;

/**
 * A description of the intersection point of a ray with a structure in 3-space. Specifies the hit time, the normal at the
 * intersection point, and the index of the intersected element.
 */
typedef RayStructureIntersectionN<3, Real> RayStructureIntersection3;

/** Interface for a structure that supports ray intersection queries in 3-space. */
typedef RayQueryStructureN<3, Real> RayQueryStructure3;

} // namespace Algorithms
} // namespace Thea

#endif
