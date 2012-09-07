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

#ifndef __Thea_AxisAlignedBox3_hpp__
#define __Thea_AxisAlignedBox3_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"

namespace Thea {

/** A 3-dimensional axis-aligned box. */
template <typename T>
class /* THEA_API */ AxisAlignedBoxN<3, T> : public Internal::AxisAlignedBoxNBase<3, T>
{
  private:
    typedef Internal::AxisAlignedBoxNBase<3, T> BaseT;

  public:
    typedef typename BaseT::VectorT VectorT;

    /** Default constructor, creates a null box. */
    AxisAlignedBoxN() : BaseT() {}

    /** Constructor. Sets the box to be a single point. */
    AxisAlignedBoxN(VectorT const & v) : BaseT(v) {}

    /** Constructor. Sets the extents of the box. */
    AxisAlignedBoxN(VectorT const & lo_, VectorT const & hi_) : BaseT(lo_, hi_) {}

    /**
     * Transform the box and return a new axis-aligned box which tightly encloses the result.
     *
     * Algorithm taken from the Ogre source code, http://www.ogre3d.org.
     */
    template <typename TransformT> AxisAlignedBoxN transformAndBound(TransformT const & tr) const
    {
      Vector3 const & lo_ = this->getLow();
      Vector3 const & hi_ = this->getHigh();

      // We sequentially compute the corners in the following order:
      // 0, 6, 5, 1, 2, 4, 7, 3
      // This sequence allows us to only change one member at a time to get at all corners. For each one, we transform it and
      // merge the resulting point.

      // min min min
      Vector3 current_corner = lo_;
      AxisAlignedBoxN result = AxisAlignedBoxN(tr * current_corner);

      // min min max
      current_corner.z() = hi_.z();
      result.merge(tr * current_corner);

      // min max max
      current_corner.y() = hi_.y();
      result.merge(tr * current_corner);

      // min max min
      current_corner.z() = lo_.z();
      result.merge(tr * current_corner);

      // max max min
      current_corner.x() = hi_.x();
      result.merge(tr * current_corner);

      // max max max
      current_corner.z() = hi_.z();
      result.merge(tr * current_corner);

      // max min max
      current_corner.y() = lo_.y();
      result.merge(tr * current_corner);

      // max min min
      current_corner.z() = lo_.z();
      result.merge(tr * current_corner);

      return result;
    }

    /**
     * Get edge number \a i of the box, where i is between 0 and 11.
     *
     * @param i Index of edge, between 0 and 11 inclusive.
     * @param start Used to return the starting point of the edge.
     * @param end Used to return the endpoint of the edge.
     */
    void getEdge(int i, Vector3 & start, Vector3 & end) const
    {
      static int const INDICES[12][6] = {
        { 0, 0, 0, 0, 0, 1 },
        { 0, 0, 0, 0, 1, 0 },
        { 0, 0, 0, 1, 0, 0 },

        { 1, 1, 0, 1, 1, 1 },
        { 1, 1, 0, 1, 0, 0 },
        { 1, 1, 0, 0, 1, 0 },

        { 1, 0, 1, 0, 0, 1 },
        { 1, 0, 1, 1, 1, 1 },
        { 1, 0, 1, 1, 0, 0 },

        { 0, 1, 1, 1, 1, 1 },
        { 0, 1, 1, 0, 0, 1 },
        { 0, 1, 1, 0, 1, 0 }
      };

      Vector3 const * v[2] = { &this->getLow(), &this->getHigh() };
      start  =  Vector3((*v[INDICES[i][0]])[0], (*v[INDICES[i][1]])[1], (*v[INDICES[i][2]])[2]);
      end    =  Vector3((*v[INDICES[i][3]])[0], (*v[INDICES[i][4]])[1], (*v[INDICES[i][5]])[2]);
    }

}; // class AxisAlignedBoxN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API AxisAlignedBoxN<3, Real>;
#endif

/** The default axis-aligned box class in 3-dimensional real space. */
typedef AxisAlignedBoxN<3, Real> AxisAlignedBox3;

} // namespace Thea

#endif
