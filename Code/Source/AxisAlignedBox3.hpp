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

#ifndef __Thea_AxisAlignedBox3_hpp__
#define __Thea_AxisAlignedBox3_hpp__

#include "Common.hpp"
#include "AxisAlignedBoxN.hpp"
#include <array>
#include <utility>

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
      VectorT const & lo_ = this->getLow();
      VectorT const & hi_ = this->getHigh();

      // We sequentially compute the corners in the following order:
      // 0, 6, 5, 1, 2, 4, 7, 3
      // This sequence allows us to only change one member at a time to get at all corners. For each one, we transform it and
      // merge the resulting point.

      // min min min
      VectorT current_corner = lo_;
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
     * Get the endpoint indices of edge number \a i of the box, where i is between 0 and 11.
     *
     * @param i Index of edge, between 0 and 11 inclusive.
     *
     * @return A pair containing the vertex indices of the edge. The indices can be mapped to vertex positions using
     *   getVertex().
     *
     * @see getVertex(), getFace()
     */
    std::pair<uintx, uintx> const & getEdge(intx i) const
    {
      static std::pair<uintx, uintx> const EDGES[12] = {
        { 0b000, 0b100 },
        { 0b000, 0b010 },
        { 0b000, 0b001 },

        { 0b011, 0b111 },
        { 0b011, 0b001 },
        { 0b011, 0b010 },

        { 0b101, 0b100 },
        { 0b101, 0b111 },
        { 0b101, 0b001 },

        { 0b110, 0b111 },
        { 0b110, 0b100 },
        { 0b110, 0b010 }
      };

      theaAssertM(i >= 0 && i < 12, "AxisAlignedBox3: Edge index out of bounds");

      return EDGES[i];
    }

    /**
     * Get the vertex indices of face number \a i of the box, where i is between 0 and 5.
     *
     * @param i Index of face, between 0 and 5 inclusive.
     *
     * @return A 4-tuple containing the vertex indices of the face, in counter-clockwise order when viewing the box outside-in.
     *   The indices can be mapped to vertex positions using getVertex().
     *
     * @see getVertex(), getEdge(), getFaceNormal()
     */
    std::array<uintx, 4> const & getFace(intx i) const
    {
      static std::array<uintx, 4> const FACES[6] = {
        { 0, 4, 6, 2 },  // -X
        { 1, 3, 7, 5 },  // +X

        { 0, 1, 5, 4 },  // -Y
        { 2, 6, 7, 3 },  // +Y

        { 0, 2, 3, 1 },  // -Z
        { 4, 5, 7, 6 }   // +Z
      };

      theaAssertM(i >= 0 && i < 6, "AxisAlignedBox3: Face index out of bounds");

      return FACES[i];
    }

    /**
     * Get the outward unit normal vector for face number \a i of the box, where i is between 0 and 5.
     *
     * @param i Index of face, between 0 and 5 inclusive.
     *
     * @see getFace()
     */
    VectorT const & getFaceNormal(intx i) const
    {
      // Synced with the face order in getFace()
      static VectorT const NORMALS[6] = {
        -VectorT::UnitX(),
         VectorT::UnitX(),

        -VectorT::UnitY(),
         VectorT::UnitY(),

        -VectorT::UnitZ(),
         VectorT::UnitZ()
      };

      theaAssertM(i >= 0 && i < 6, "AxisAlignedBox3: Face index out of bounds");

      return NORMALS[i];
    }

}; // class AxisAlignedBoxN<3, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API AxisAlignedBoxN<3, Real>;
#endif

/** The default axis-aligned box class in 3-dimensional real space. */
typedef AxisAlignedBoxN<3, Real> AxisAlignedBox3;

} // namespace Thea

#endif
