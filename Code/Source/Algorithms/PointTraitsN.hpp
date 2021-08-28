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
 * Has boolean member <code>value = true</code> if <code>T</code> can be identified with a single point in N-D space, else
 * <code>value = false</code>. Unless you specialize the class to set the value to true, it is false by default.
 *
 * @see PointTraitsN
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
// Traits class to get the position of a point-like object
//=============================================================================================================================

/**
 * Traits for an object which can be identified with a single point in N-space. If IsPointN<T>::value is false, then the
 * getPosition() function returns some arbitrary value by default. The default implementation for IsPointN<T>::value == true
 * handles Vector and types implicitly convertible to it. Specialize as needed for other classes.
 *
 * @see IsPointN
 */
template <typename T, int N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ PointTraitsN
{
  public:
    /** Get the position of the point. */
    static Vector<N, ScalarT> getPosition(T const & t) { return t; }

}; // class PointTraitsN

// Default specialization for all classes T with IsPointN<T>::value == false
template <typename T, int N, typename ScalarT>
class PointTraitsN< T, N, ScalarT, typename std::enable_if< !IsPointN<T, N>::value >::type >
{
  public:
    static Vector<N, ScalarT> getPosition(T const & t) { return Vector<N, ScalarT>::Zero(); }

}; // class PointTraitsN<T, N, ScalarT, IsPointN<T, N>::value == false>

// Partial specialization of PointTraitsN for const types
template <typename T, int N, typename ScalarT>
class /* THEA_API */ PointTraitsN< T const, N, ScalarT, typename std::enable_if< IsPointN<T, N>::value >::type >
{
  public:
    static Vector<N, ScalarT> getPosition(T const & t) { return PointTraitsN<T, N, ScalarT>::getPosition(t); }

}; // class PointTraitsN<T const, N, ScalarT>

// Partial specialization of PointTraitsN for pointer types
template <typename T, int N, typename ScalarT>
class /* THEA_API */ PointTraitsN< T *, N, ScalarT, typename std::enable_if< IsPointN<T, N>::value >::type >
{
  public:
    static Vector<N, ScalarT> getPosition(T * t) { return PointTraitsN<T, N, ScalarT>::getPosition(*t); }

}; // class PointTraitsN<T *, N, ScalarT>

} // namespace Algorithms
} // namespace Thea

#endif
