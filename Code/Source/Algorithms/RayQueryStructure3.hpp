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

#ifndef __Thea_Algorithms_RayQueryStructure3_hpp__
#define __Thea_Algorithms_RayQueryStructure3_hpp__

#include "../Common.hpp"
#include "../RayIntersectable3.hpp"

namespace Thea {
namespace Algorithms {

/**
 * A description of the intersection point of a ray with a structure. Specifies the hit time, the normal at the intersection
 * point, and the index of the intersected element.
 */
class THEA_API RayStructureIntersection3 : public RayIntersection3
{
  private:
    typedef RayIntersection3 BaseT;

    long element_index;

  public:
    /** Constructor. */
    RayStructureIntersection3(Real time_ = -1, Vector3 const * normal_ = NULL, long element_index_ = -1)
    : BaseT(time_, normal_), element_index(element_index_)
    {}

    /** Constructor. */
    RayStructureIntersection3(RayIntersection3 const & isec, long element_index_ = -1)
    : BaseT(isec), element_index(element_index_)
    {}

    /** Get the index of the intersected element. */
    long getElementIndex() const { return element_index; }

    /** Set the index of the intersected element. */
    void setElementIndex(long element_index_) { element_index = element_index_; }

}; // class RayStructureIntersection3

/**
 * Interface for a structure that supports ray intersection queries in 3-space. None of the functions are virtual, this just
 * defines a concept subclasses must implement.
 */
class THEA_API RayQueryStructure3
{
  public:
    THEA_DEF_POINTER_TYPES(RayQueryStructure3, shared_ptr, weak_ptr)

    /** Check if a ray intersects the structure in the forward direction. */
    template <typename RayIntersectionTesterT> bool rayIntersects(Ray3 const & ray, Real max_time = -1) const;

    /**
     * Get the time taken for a ray to intersect the structure, or a negative value if there was no intersection in the forward
     * direction.
     */
    template <typename RayIntersectionTesterT> Real rayIntersectionTime(Ray3 const & ray, Real max_time = -1) const;

    /**
     * Get the intersection of a ray with the structure, including the hit time, the normal at the intersection point, and the
     * index of the intersected element. A negative time is returned if there was no intersection in the forward direction. A
     * zero normal and a negative index are returned if those quantities are not known.
     */
    template <typename RayIntersectionTesterT>
    RayStructureIntersection3 rayStructureIntersection(Ray3 const & ray, Real max_time = -1) const;

}; // class RayQueryStructure3

} // namespace Algorithms
} // namespace Thea

#endif
