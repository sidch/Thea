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
//============== License for line-AAB distance calculation ===================
//
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)
//
//============================================================================

#ifndef __Thea_AxisAlignedBox3_hpp__
#define __Thea_AxisAlignedBox3_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"

namespace Thea {

namespace Internal {

// From https://www.geometrictools.com/GTEngine/Include/Mathematics/GteDistLine3AlignedBox3.h
template <typename T>
class LineBoxDistance3
{
  public:
    typedef VectorN<3, T> VectorT;

    struct Result
    {
      T sqrDistance;
      T lineParameter;
    };

    // Compute the distance and closest point between a line and an
    // axis-aligned box whose center is the origin.  On input, 'point' is the
    // line origin and 'direction' is the line direction.  On output, 'point'
    // is the point on the box closest to the line.  The 'direction' is
    // non-const to allow transforming the problem into the first octant.
    static void DoQuery(VectorT & point, VectorT & direction,
                        VectorT const & boxExtent, Result & result);

  private:
    static void Face(int i0, int i1, int i2, VectorT & pnt,
                     VectorT const & dir, VectorT const & PmE,
                     VectorT const & boxExtent, Result & result);

    static void CaseNoZeros(VectorT & pnt, VectorT const & dir,
                            VectorT const & boxExtent, Result & result);

    static void Case0(int i0, int i1, int i2, VectorT & pnt,
                      VectorT const & dir, VectorT const & boxExtent,
                      Result & result);

    static void Case00(int i0, int i1, int i2, VectorT & pnt,
                       VectorT const & dir, VectorT const & boxExtent,
                       Result & result);

    static void Case000(VectorT & pnt, VectorT const & boxExtent,
                        Result & result);

}; // class LineBoxDistance3

} // namespace Internal

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
    void getEdge(int i, VectorT & start, VectorT & end) const
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

      VectorT const * v[2] = { &this->getLow(), &this->getHigh() };
      start  =  VectorT((*v[INDICES[i][0]])[0], (*v[INDICES[i][1]])[1], (*v[INDICES[i][2]])[2]);
      end    =  VectorT((*v[INDICES[i][3]])[0], (*v[INDICES[i][4]])[1], (*v[INDICES[i][5]])[2]);
    }

    using BaseT::distance;
    using BaseT::squaredDistance;

    T distance(LineN<3, T> const & line) const { return std::sqrt(squaredDistance(line)); }

    T squaredDistance(LineN<3, T> const & line, VectorT * this_pt = NULL, VectorT * line_pt = NULL) const
    {
      // Translate the line and box so that the box has center at the origin.
      VectorT boxCenter = this->getCenter(), boxExtent = 0.5f * this->getExtent();
      VectorT point = line.getPoint() - boxCenter;
      VectorT direction = line.getDirection();

      typename Internal::LineBoxDistance3<T>::Result result;
      Internal::LineBoxDistance3<T>::DoQuery(point, direction, boxExtent, result);

      if (this_pt) *this_pt = boxCenter + point;
      if (line_pt) *line_pt = line.getPoint() + result.lineParameter * line.getDirection();

      return result.sqrDistance;
    }

    T distance(LineSegmentN<3, T> const & seg) const { return std::sqrt(squaredDistance(seg)); }

    T squaredDistance(LineSegmentN<3, T> const & seg, VectorT * this_pt = NULL, VectorT * seg_pt = NULL) const
    {
      // Translate the line and box so that the box has center at the origin.
      VectorT boxCenter = this->getCenter(), boxExtent = 0.5f * this->getExtent();
      VectorT e0 = seg.getEndpoint(0);
      VectorT point = e0 - boxCenter;

      // FIXME: Does DoQuery actually need a unit length direction vector?
      T seg_len = seg.length();
      VectorT seg_dir = (seg_len >= Math::eps<T>() ? seg.getDirection() / seg_len : VectorT::unitX());
      VectorT direction = seg_dir;  // can be overwritten by DoQuery

      typename Internal::LineBoxDistance3<T>::Result result;
      Internal::LineBoxDistance3<T>::DoQuery(point, direction, boxExtent, result);

      if (result.lineParameter >= 0 && result.lineParameter <= seg_len)
      {
        if (this_pt) *this_pt = boxCenter + point;
        if (seg_pt)  *seg_pt  = e0 + result.lineParameter * seg_dir;
      }
      else
      {
        VectorT e1 = seg.getEndpoint(1);
        VectorT c0 = this->closestPoint(e0);
        VectorT c1 = this->closestPoint(e1);
        Real sqdist0 = (c0 - e0).squaredLength();
        Real sqdist1 = (c1 - e1).squaredLength();
        if (sqdist0 < sqdist1)
        {
          if (this_pt) *this_pt = c0;
          if (seg_pt)  *seg_pt  = e0;
          result.sqrDistance = sqdist0;
        }
        else
        {
          if (this_pt) *this_pt = c1;
          if (seg_pt)  *seg_pt  = e1;
          result.sqrDistance = sqdist1;
        }
      }

      return result.sqrDistance;
    }

    T distance(RayN<3, T> const & ray) const { return std::sqrt(squaredDistance(ray)); }

    T squaredDistance(RayN<3, T> const & ray, VectorT * this_pt = NULL, VectorT * ray_pt = NULL) const
    {
      // Translate the line and box so that the box has center at the origin.
      VectorT boxCenter = this->getCenter(), boxExtent = 0.5f * this->getExtent();
      VectorT point = ray.getOrigin() - boxCenter;
      VectorT ray_dir = ray.getDirection().unit();  // FIXME: Does DoQuery actually need a unit length direction vector?
      VectorT direction = ray_dir;  // can be overwritten by DoQuery

      typename Internal::LineBoxDistance3<T>::Result result;
      Internal::LineBoxDistance3<T>::DoQuery(point, direction, boxExtent, result);

      if (result.lineParameter >= 0)
      {
        if (this_pt) *this_pt = boxCenter + point;
        if (ray_pt)  *ray_pt  = ray.getOrigin() + result.lineParameter * ray_dir;
      }
      else
      {
        VectorT c = this->closestPoint(ray.getOrigin());
        result.sqrDistance = (c - ray.getOrigin()).squaredLength();

        if (this_pt) *this_pt = c;
        if (ray_pt)  *ray_pt  = ray.getOrigin();
      }

      return result.sqrDistance;
    }

}; // class AxisAlignedBoxN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API AxisAlignedBoxN<3, Real>;
#endif

/** The default axis-aligned box class in 3-dimensional real space. */
typedef AxisAlignedBoxN<3, Real> AxisAlignedBox3;

namespace Internal {

//=============================================================================================================================
//
// Line-to-axis-aligned-box distance calculations in 3D. Note that "extent" below refers to half the extent of the box above.
//
//=============================================================================================================================

template <typename T>
void
LineBoxDistance3<T>::DoQuery(VectorT & point, VectorT & direction,
                             VectorT const & boxExtent, Result & result)
{
  result.sqrDistance = (T)0;
  result.lineParameter = (T)0;

  // Apply reflections so that direction vector has nonnegative components.
  bool reflect[3];
  for (int i = 0; i < 3; ++i)
  {
    if (direction[i] < (T)0)
    {
      point[i] = -point[i];
      direction[i] = -direction[i];
      reflect[i] = true;
    }
    else
    {
      reflect[i] = false;
    }
  }

  if (direction[0] > (T)0)
  {
    if (direction[1] > (T)0)
    {
      if (direction[2] > (T)0)  // (+,+,+)
      {
        CaseNoZeros(point, direction, boxExtent, result);
      }
      else  // (+,+,0)
      {
        Case0(0, 1, 2, point, direction, boxExtent, result);
      }
    }
    else
    {
      if (direction[2] > (T)0)  // (+,0,+)
      {
        Case0(0, 2, 1, point, direction, boxExtent, result);
      }
      else  // (+,0,0)
      {
        Case00(0, 1, 2, point, direction, boxExtent, result);
      }
    }
  }
  else
  {
    if (direction[1] > (T)0)
    {
      if (direction[2] > (T)0)  // (0,+,+)
      {
        Case0(1, 2, 0, point, direction, boxExtent, result);
      }
      else  // (0,+,0)
      {
        Case00(1, 0, 2, point, direction, boxExtent, result);
      }
    }
    else
    {
      if (direction[2] > (T)0)  // (0,0,+)
      {
        Case00(2, 0, 1, point, direction, boxExtent, result);
      }
      else  // (0,0,0)
      {
        Case000(point, boxExtent, result);
      }
    }
  }

  // Undo the reflections applied previously.
  for (int i = 0; i < 3; ++i)
  {
    if (reflect[i])
    {
      point[i] = -point[i];
    }
  }
}

template <typename T>
void
LineBoxDistance3<T>::Face(int i0, int i1,
                          int i2, VectorT & pnt, VectorT const & dir,
                          VectorT const & PmE, VectorT const & boxExtent, Result & result)
{
  VectorT PpE;
  T lenSqr, inv, tmp, param, t, delta;

  PpE[i1] = pnt[i1] + boxExtent[i1];
  PpE[i2] = pnt[i2] + boxExtent[i2];
  if (dir[i0] * PpE[i1] >= dir[i1] * PmE[i0])
  {
    if (dir[i0] * PpE[i2] >= dir[i2] * PmE[i0])
    {
      // v[i1] >= -e[i1], v[i2] >= -e[i2] (distance = 0)
      pnt[i0] = boxExtent[i0];
      inv = ((T)1) / dir[i0];
      pnt[i1] -= dir[i1] * PmE[i0] * inv;
      pnt[i2] -= dir[i2] * PmE[i0] * inv;
      result.lineParameter = -PmE[i0] * inv;
    }
    else
    {
      // v[i1] >= -e[i1], v[i2] < -e[i2]
      lenSqr = dir[i0] * dir[i0] + dir[i2] * dir[i2];
      tmp = lenSqr * PpE[i1] - dir[i1] * (dir[i0] * PmE[i0] +
                                          dir[i2] * PpE[i2]);
      if (tmp <= ((T)2)*lenSqr * boxExtent[i1])
      {
        t = tmp / lenSqr;
        lenSqr += dir[i1] * dir[i1];
        tmp = PpE[i1] - t;
        delta = dir[i0] * PmE[i0] + dir[i1] * tmp + dir[i2] * PpE[i2];
        param = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + tmp * tmp +
                              PpE[i2] * PpE[i2] + delta * param;

        result.lineParameter = param;
        pnt[i0] = boxExtent[i0];
        pnt[i1] = t - boxExtent[i1];
        pnt[i2] = -boxExtent[i2];
      }
      else
      {
        lenSqr += dir[i1] * dir[i1];
        delta = dir[i0] * PmE[i0] + dir[i1] * PmE[i1] + dir[i2] * PpE[i2];
        param = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + PmE[i1] * PmE[i1]
                              + PpE[i2] * PpE[i2] + delta * param;

        result.lineParameter = param;
        pnt[i0] = boxExtent[i0];
        pnt[i1] = boxExtent[i1];
        pnt[i2] = -boxExtent[i2];
      }
    }
  }
  else
  {
    if (dir[i0] * PpE[i2] >= dir[i2] * PmE[i0])
    {
      // v[i1] < -e[i1], v[i2] >= -e[i2]
      lenSqr = dir[i0] * dir[i0] + dir[i1] * dir[i1];
      tmp = lenSqr * PpE[i2] - dir[i2] * (dir[i0] * PmE[i0] +
                                          dir[i1] * PpE[i1]);
      if (tmp <= ((T)2)*lenSqr * boxExtent[i2])
      {
        t = tmp / lenSqr;
        lenSqr += dir[i2] * dir[i2];
        tmp = PpE[i2] - t;
        delta = dir[i0] * PmE[i0] + dir[i1] * PpE[i1] + dir[i2] * tmp;
        param = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] +
                              tmp * tmp + delta * param;

        result.lineParameter = param;
        pnt[i0] = boxExtent[i0];
        pnt[i1] = -boxExtent[i1];
        pnt[i2] = t - boxExtent[i2];
      }
      else
      {
        lenSqr += dir[i2] * dir[i2];
        delta = dir[i0] * PmE[i0] + dir[i1] * PpE[i1] + dir[i2] * PmE[i2];
        param = -delta / lenSqr;
        result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] +
                              PmE[i2] * PmE[i2] + delta * param;

        result.lineParameter = param;
        pnt[i0] = boxExtent[i0];
        pnt[i1] = -boxExtent[i1];
        pnt[i2] = boxExtent[i2];
      }
    }
    else
    {
      // v[i1] < -e[i1], v[i2] < -e[i2]
      lenSqr = dir[i0] * dir[i0] + dir[i2] * dir[i2];
      tmp = lenSqr * PpE[i1] - dir[i1] * (dir[i0] * PmE[i0] +
                                          dir[i2] * PpE[i2]);
      if (tmp >= (T)0)
      {
        // v[i1]-edge is closest
        if (tmp <= ((T)2)*lenSqr * boxExtent[i1])
        {
          t = tmp / lenSqr;
          lenSqr += dir[i1] * dir[i1];
          tmp = PpE[i1] - t;
          delta = dir[i0] * PmE[i0] + dir[i1] * tmp + dir[i2] * PpE[i2];
          param = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + tmp * tmp +
                                PpE[i2] * PpE[i2] + delta * param;

          result.lineParameter = param;
          pnt[i0] = boxExtent[i0];
          pnt[i1] = t - boxExtent[i1];
          pnt[i2] = -boxExtent[i2];
        }
        else
        {
          lenSqr += dir[i1] * dir[i1];
          delta = dir[i0] * PmE[i0] + dir[i1] * PmE[i1]
                  + dir[i2] * PpE[i2];
          param = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + PmE[i1] * PmE[i1]
                                + PpE[i2] * PpE[i2] + delta * param;

          result.lineParameter = param;
          pnt[i0] = boxExtent[i0];
          pnt[i1] = boxExtent[i1];
          pnt[i2] = -boxExtent[i2];
        }
        return;
      }

      lenSqr = dir[i0] * dir[i0] + dir[i1] * dir[i1];
      tmp = lenSqr * PpE[i2] - dir[i2] * (dir[i0] * PmE[i0] +
                                          dir[i1] * PpE[i1]);
      if (tmp >= (T)0)
      {
        // v[i2]-edge is closest
        if (tmp <= ((T)2)*lenSqr * boxExtent[i2])
        {
          t = tmp / lenSqr;
          lenSqr += dir[i2] * dir[i2];
          tmp = PpE[i2] - t;
          delta = dir[i0] * PmE[i0] + dir[i1] * PpE[i1] + dir[i2] * tmp;
          param = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1] +
                                tmp * tmp + delta * param;

          result.lineParameter = param;
          pnt[i0] = boxExtent[i0];
          pnt[i1] = -boxExtent[i1];
          pnt[i2] = t - boxExtent[i2];
        }
        else
        {
          lenSqr += dir[i2] * dir[i2];
          delta = dir[i0] * PmE[i0] + dir[i1] * PpE[i1]
                  + dir[i2] * PmE[i2];
          param = -delta / lenSqr;
          result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1]
                                + PmE[i2] * PmE[i2] + delta * param;

          result.lineParameter = param;
          pnt[i0] = boxExtent[i0];
          pnt[i1] = -boxExtent[i1];
          pnt[i2] = boxExtent[i2];
        }
        return;
      }

      // (v[i1],v[i2])-corner is closest
      lenSqr += dir[i2] * dir[i2];
      delta = dir[i0] * PmE[i0] + dir[i1] * PpE[i1] + dir[i2] * PpE[i2];
      param = -delta / lenSqr;
      result.sqrDistance += PmE[i0] * PmE[i0] + PpE[i1] * PpE[i1]
                            + PpE[i2] * PpE[i2] + delta * param;

      result.lineParameter = param;
      pnt[i0] = boxExtent[i0];
      pnt[i1] = -boxExtent[i1];
      pnt[i2] = -boxExtent[i2];
    }
  }
}

template <typename T>
void
LineBoxDistance3<T>::CaseNoZeros(VectorT & pnt, VectorT const & dir,
                                 VectorT const & boxExtent, Result & result)
{
  VectorT PmE = pnt - boxExtent;
  T prodDxPy = dir[0] * PmE[1];
  T prodDyPx = dir[1] * PmE[0];
  T prodDzPx, prodDxPz, prodDzPy, prodDyPz;

  if (prodDyPx >= prodDxPy)
  {
    prodDzPx = dir[2] * PmE[0];
    prodDxPz = dir[0] * PmE[2];
    if (prodDzPx >= prodDxPz)
    {
      // line intersects x = e0
      Face(0, 1, 2, pnt, dir, PmE, boxExtent, result);
    }
    else
    {
      // line intersects z = e2
      Face(2, 0, 1, pnt, dir, PmE, boxExtent, result);
    }
  }
  else
  {
    prodDzPy = dir[2] * PmE[1];
    prodDyPz = dir[1] * PmE[2];
    if (prodDzPy >= prodDyPz)
    {
      // line intersects y = e1
      Face(1, 2, 0, pnt, dir, PmE, boxExtent, result);
    }
    else
    {
      // line intersects z = e2
      Face(2, 0, 1, pnt, dir, PmE, boxExtent, result);
    }
  }
}

template <typename T>
void
LineBoxDistance3<T>::Case0(int i0, int i1,
                           int i2, VectorT & pnt, VectorT const & dir,
                           VectorT const & boxExtent, Result & result)
{
  T PmE0 = pnt[i0] - boxExtent[i0];
  T PmE1 = pnt[i1] - boxExtent[i1];
  T prod0 = dir[i1] * PmE0;
  T prod1 = dir[i0] * PmE1;
  T delta, invLSqr, inv;

  if (prod0 >= prod1)
  {
    // line intersects P[i0] = e[i0]
    pnt[i0] = boxExtent[i0];

    T PpE1 = pnt[i1] + boxExtent[i1];
    delta = prod0 - dir[i0] * PpE1;
    if (delta >= (T)0)
    {
      invLSqr = ((T)1) / (dir[i0] * dir[i0] + dir[i1] * dir[i1]);
      result.sqrDistance += delta * delta * invLSqr;
      pnt[i1] = -boxExtent[i1];
      result.lineParameter = -(dir[i0] * PmE0 + dir[i1] * PpE1) * invLSqr;
    }
    else
    {
      inv = ((T)1) / dir[i0];
      pnt[i1] -= prod0 * inv;
      result.lineParameter = -PmE0 * inv;
    }
  }
  else
  {
    // line intersects P[i1] = e[i1]
    pnt[i1] = boxExtent[i1];

    T PpE0 = pnt[i0] + boxExtent[i0];
    delta = prod1 - dir[i1] * PpE0;
    if (delta >= (T)0)
    {
      invLSqr = ((T)1) / (dir[i0] * dir[i0] + dir[i1] * dir[i1]);
      result.sqrDistance += delta * delta * invLSqr;
      pnt[i0] = -boxExtent[i0];
      result.lineParameter = -(dir[i0] * PpE0 + dir[i1] * PmE1) * invLSqr;
    }
    else
    {
      inv = ((T)1) / dir[i1];
      pnt[i0] -= prod1 * inv;
      result.lineParameter = -PmE1 * inv;
    }
  }

  if (pnt[i2] < -boxExtent[i2])
  {
    delta = pnt[i2] + boxExtent[i2];
    result.sqrDistance += delta * delta;
    pnt[i2] = -boxExtent[i2];
  }
  else if (pnt[i2] > boxExtent[i2])
  {
    delta = pnt[i2] - boxExtent[i2];
    result.sqrDistance += delta * delta;
    pnt[i2] = boxExtent[i2];
  }
}

template <typename T>
void
LineBoxDistance3<T>::Case00(int i0, int i1,
                            int i2, VectorT & pnt, VectorT const & dir,
                            VectorT const & boxExtent, Result & result)
{
  T delta;

  result.lineParameter = (boxExtent[i0] - pnt[i0]) / dir[i0];

  pnt[i0] = boxExtent[i0];

  if (pnt[i1] < -boxExtent[i1])
  {
    delta = pnt[i1] + boxExtent[i1];
    result.sqrDistance += delta * delta;
    pnt[i1] = -boxExtent[i1];
  }
  else if (pnt[i1] > boxExtent[i1])
  {
    delta = pnt[i1] - boxExtent[i1];
    result.sqrDistance += delta * delta;
    pnt[i1] = boxExtent[i1];
  }

  if (pnt[i2] < -boxExtent[i2])
  {
    delta = pnt[i2] + boxExtent[i2];
    result.sqrDistance += delta * delta;
    pnt[i2] = -boxExtent[i2];
  }
  else if (pnt[i2] > boxExtent[i2])
  {
    delta = pnt[i2] - boxExtent[i2];
    result.sqrDistance += delta * delta;
    pnt[i2] = boxExtent[i2];
  }
}

template <typename T>
void
LineBoxDistance3<T>::Case000(VectorT & pnt, VectorT const & boxExtent, Result & result)
{
  T delta;

  if (pnt[0] < -boxExtent[0])
  {
    delta = pnt[0] + boxExtent[0];
    result.sqrDistance += delta * delta;
    pnt[0] = -boxExtent[0];
  }
  else if (pnt[0] > boxExtent[0])
  {
    delta = pnt[0] - boxExtent[0];
    result.sqrDistance += delta * delta;
    pnt[0] = boxExtent[0];
  }

  if (pnt[1] < -boxExtent[1])
  {
    delta = pnt[1] + boxExtent[1];
    result.sqrDistance += delta * delta;
    pnt[1] = -boxExtent[1];
  }
  else if (pnt[1] > boxExtent[1])
  {
    delta = pnt[1] - boxExtent[1];
    result.sqrDistance += delta * delta;
    pnt[1] = boxExtent[1];
  }

  if (pnt[2] < -boxExtent[2])
  {
    delta = pnt[2] + boxExtent[2];
    result.sqrDistance += delta * delta;
    pnt[2] = -boxExtent[2];
  }
  else if (pnt[2] > boxExtent[2])
  {
    delta = pnt[2] - boxExtent[2];
    result.sqrDistance += delta * delta;
    pnt[2] = boxExtent[2];
  }
}

} // namespace Internal

} // namespace Thea

#endif
