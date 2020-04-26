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

#ifndef __Thea_RigidTransform2_hpp__
#define __Thea_RigidTransform2_hpp__

#include "Common.hpp"
#include "AffineTransform2.hpp"
#include "RigidTransformN.hpp"

namespace Thea {

/**
 * A rigid transformation in 2-space, consisting of a rotation followed by a translation.
 *
 * @note While this is technically an affine transform, it restricts enough functionality to make a separate implementation
 * preferable. It can be trivially converted to an AffineTransform2 using the implicit conversion operator or toAffine().
 */
template <typename T>
class /* THEA_API */ RigidTransformN<2, T> : public Internal::RigidTransformNBase<2, T>
{
  private:
    typedef Internal::RigidTransformNBase<2, T> BaseT;

  public:
    typedef typename BaseT::AffineTransformT  AffineTransformT;
    typedef typename BaseT::VectorT           VectorT;
    typedef typename BaseT::MatrixT           MatrixT;

    /** Default constructor. Constructs the identity transform. */
    RigidTransformN() {}

    using BaseT::translation;

    /** Construct a translation. */
    static RigidTransformN translation(T const & tx, T const & ty)
    {
      return BaseT::translation(VectorT(tx, ty));
    }

    /** Construct a rotation specified by an angle (in radians) around the origin. */
    static RigidTransformN rotation(Real radians)
    {
      return BaseT::_fromAffine(AffineTransformT::rotation(radians));
    }

}; // class RigidTransformN<2, T>

#ifdef THEA_EXPORT_INSTANTIATION
  template class THEA_API RigidTransformN<2, Real>;
#endif

/** The default rigid transform class in real 2-space. */
typedef RigidTransformN<2, Real> RigidTransform2;

} // namespace Thea

#endif
