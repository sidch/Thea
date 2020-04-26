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

#ifndef __Thea_Algorithms_PointTraitsN_hpp__
#define __Thea_Algorithms_PointTraitsN_hpp__

#include "../Common.hpp"
#include "../MatVec.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

//=============================================================================================================================
// Flag that should be set to true for objects equivalent to single points
//=============================================================================================================================

/**
 * Has boolean member <code>value = true</code> if <code>T</code> can be identified with a single point in n-D space, else
 * <code>value = false</code>. Unless you specialize the class to set the value to true, it is false by default.
 *
 * @see TraitsN
 */
template <typename T, int N, typename Enable = void>
class /* THEA_API */ IsPointN
{
  public:
    static bool const value = false;
};

// Partial specialization for const and pointer types
template <typename T, int N> class IsPointN<T const, N>  { public: static bool const value = IsPointN<T, N>::value; };
template <typename T, int N> class IsPointN<T *, N>      { public: static bool const value = IsPointN<T, N>::value; };

// Specialization for Vector
template <typename T, int N>
class /* THEA_API */ IsPointN< T, N, typename std::enable_if< std::is_base_of< Eigen::DenseBase<T>, T >::value
                                                           && T::IsVectorAtCompileTime != 0
                                                           && T::SizeAtCompileTime == N >::type >
{
  public:
    static bool const value = true;
};

/** Same as IsPointN (no need to specialize it separately), except false for const or pointer types. */
template <typename T, int N>
class /* THEA_API */ IsRawPointN
{
  public:
    static bool const value = IsPointN<T, N>::value
                           && !std::is_const<T>::value
                           && !std::is_pointer<T>::value;
};

/** Same as IsPointN (no need to specialize it separately), except false for pointer types. */
template <typename T, int N>
class /* THEA_API */ IsNonReferencedPointN
{
  public:
    static bool const value = IsPointN<T, N>::value
                           && !std::is_pointer<T>::value;
};

//=============================================================================================================================
// Traits class to get position of a point-like object
//=============================================================================================================================

/**
 * Traits for an object which can be identified with a single point in N-space.
 *
 * @see IsPointN
 */
template <typename PointT, int N, typename ScalarT = Real>
class /* THEA_API */ PointTraitsN
{
  public:
    typedef Vector<N, ScalarT> VectorT;  ///< A vector in N-space.

    /** Get the position of the "point". The default implementation is for Vector and types implicitly convertible to it. */
    static VectorT getPosition(PointT const & p) { return p; }

}; // class PointTraitsN

// Partial specialization of PointTraitsN for const types
template <typename PointT, int N, typename ScalarT>
class /* THEA_API */ PointTraitsN<PointT const, N, ScalarT>
{
  public:
    typedef Vector<N, ScalarT> VectorT;

    static VectorT getPosition(PointT const & p) { return PointTraitsN<PointT, N, ScalarT>::getPosition(p); }

}; // class PointTraitsN<PointT const, N, ScalarT>

// Partial specialization of PointTraitsN for pointer types
template <typename PointT, int N, typename ScalarT>
class /* THEA_API */ PointTraitsN<PointT *, N, ScalarT>
{
  public:
    typedef Vector<N, ScalarT> VectorT;

    static VectorT getPosition(PointT * p) { return PointTraitsN<PointT, N, ScalarT>::getPosition(*p); }

}; // class PointTraitsN<PointT *, N, ScalarT>

} // namespace Algorithms
} // namespace Thea

#endif
