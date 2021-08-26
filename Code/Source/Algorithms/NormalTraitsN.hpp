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
// First version: 2021
//
//============================================================================

#ifndef __Thea_Algorithms_NormalTraitsN_hpp__
#define __Thea_Algorithms_NormalTraitsN_hpp__

#include "../Common.hpp"
#include "../MatVec.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

//=============================================================================================================================
// Flag that should be set to true for objects that have normal vectors
//=============================================================================================================================

/**
 * Has boolean member <code>value = true</code> if <code>T</code> has an N-D normal vector, else <code>value = false</code>.
 * Unless you specialize the class to set the value to true, it is false by default.
 *
 * @note The value is <i>false</i> for vectors -- they can <i>represent</i> normal vectors, but do not themselves <i>possess</i>
 *   normals.
 *
 * @see NormalTraitsN
 */
template <typename T, int N, typename Enable = void>
class /* THEA_API */ HasNormalN
{
  public:
    static bool const value = false;
};

// Partial specialization for const and pointer types
template <typename T, int N> class HasNormalN<T const, N>  { public: static bool const value = HasNormalN<T, N>::value; };
template <typename T, int N> class HasNormalN<T *, N>      { public: static bool const value = HasNormalN<T, N>::value; };

//=============================================================================================================================
// Traits class to get the normal of an object
//=============================================================================================================================

/**
 * Traits for an object which has an N-D normal vector. If HasNormalN<T>::value is false, then the getNormal() function returns
 * some arbitrary value by default. The default implementation for HasNormalN<T>::value == true assumes the type T has a
 * <tt>T::getNormal()</tt> member function. Specialize as needed for other classes.
 *
 * @see HasNormalN
 */
template <typename T, int N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ NormalTraitsN
{
  public:
    /** Get the normal vector of the object. */
    static Vector<N, ScalarT> getNormal(T const & t) { return t.getNormal(); }

}; // class NormalTraitsN

// Default specialization for all classes T with HasNormalN<T>::value == false
template <typename T, int N, typename ScalarT>
class NormalTraitsN< T, N, ScalarT, typename std::enable_if< !HasNormalN<T, N>::value >::type >
{
  public:
    static Vector<N, ScalarT> getNormal(T const & t) { return Vector<N, ScalarT>::Zero(); }

}; // class NormalTraitsN<T, N, ScalarT, HasNormalN<T, N>::value == false>

// Partial specialization of NormalTraitsN for const types
template <typename T, int N, typename ScalarT>
class /* THEA_API */ NormalTraitsN<T const, N, ScalarT>
{
  public:
    static Vector<N, ScalarT> getNormal(T const & t) { return NormalTraitsN<T, N, ScalarT>::getNormal(t); }

}; // class NormalTraitsN<T const, N, ScalarT>

// Partial specialization of NormalTraitsN for pointer types
template <typename T, int N, typename ScalarT>
class /* THEA_API */ NormalTraitsN<T *, N, ScalarT>
{
  public:
    static Vector<N, ScalarT> getNormal(T * t) { return NormalTraitsN<T, N, ScalarT>::getNormal(*t); }

}; // class NormalTraitsN<T *, N, ScalarT>

} // namespace Algorithms
} // namespace Thea

#endif
