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

#ifndef __Thea_VectorN_hpp__
#define __Thea_VectorN_hpp__

#include "Common.hpp"
#include "Algorithms/FastCopy.hpp"
#include "Math.hpp"
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <algorithm>
#include <iterator>
#include <limits>
#include <sstream>

namespace Thea {

// Forward declarations
template <long N, typename T> class VectorN;
template <long M, long N, typename T> class MatrixMN;

/**
 * Namespace for internal classes and functions.
 *
 * @note The members of this namespace are <b>INTERNAL</b>! Don't use them directly.
 */
namespace Internal {

// Check that a given type is a valid scalar.
template <typename S, typename T> struct ScalarCheck
{
  static bool const value = boost::is_arithmetic<S>::value || boost::is_same<S, T>::value;

}; // struct ScalarCheck

/**
 * <b>[Internal]</b> Base class for fixed-size N-dimensional vectors, where N is any <b>positive</b> (non-zero) integer and T is
 * a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <long N, typename T>
class /* THEA_DLL_LOCAL */ VectorNBase
{
  public:
    typedef VectorN<N, T>                     VectorT;               ///< N-dimensional vector type.
    typedef T                                 Value;                 ///< Type of values stored in the vector.
    typedef T                                 value_type;            ///< Type of values stored in the vector (STL convention).
    typedef T *                               Iterator;              ///< Forward iterator through elements.
    typedef T const *                         ConstIterator;         ///< Forward const iterator through elements.
    typedef std::reverse_iterator<T *>        ReverseIterator;       ///< Reverse iterator through elements.
    typedef std::reverse_iterator<T const *>  ConstReverseIterator;  ///< Reverse const iterator through elements.

    THEA_DEF_POINTER_TYPES(VectorT, shared_ptr, weak_ptr)

    /** Default constructor (does not initialize anything). */
    VectorNBase() {}

    /** Initialize all components to a single value. */
    explicit VectorNBase(T const & fill_value) { fill(fill_value); }

    /** Copy constructor. */
    template <typename U> VectorNBase(VectorNBase<N, U> const & src)
    {
      Algorithms::fastCopy(&src[0], &src[0] + N, &values[0]);
    }

    /** Initialize from a column matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorNBase(MatrixMN<N, 1, U> const & src)
    {
      Algorithms::fastCopy(&src(0, 0), &src(0, 0) + N, &values[0]);
    }

    /** Initialize from a row matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorNBase(MatrixMN<1, N, U> const & src)
    {
      Algorithms::fastCopy(&src(0, 0), &src(0, 0) + N, &values[0]);
    }

    /** Convert the vector to a column matrix (not defined unless MatrixMN.hpp is included). */
    MatrixMN<N, 1, T> toColumnMatrix() const;

    /** Convert the vector to a row matrix (not defined unless MatrixMN.hpp is included). */
    MatrixMN<1, N, T> toRowMatrix() const;

    /** Set all elements of the vector to the same value. */
    void fill(T const & fill_value)
    {
      for (long i = 0; i < N; ++i)
        values[i] = fill_value;
    }

    /** Get an iterator pointing to the beginning of the vector. */
    Iterator begin() { return &values[0]; }

    /** Get a const iterator pointing to the beginning of the vector. */
    ConstIterator begin() const { return &values[0]; }

    /** Get an iterator pointing to the end of the vector. */
    Iterator end() { return &values[0] + N; }

    /** Get a const iterator pointing to the end of the vector. */
    ConstIterator end() const { return &values[0] + N; }

    /** Get a reverse iterator pointing to the beginning of the reversed vector. */
    ReverseIterator rbegin() { return ReverseIterator(end()); }

    /** Get a const iterator pointing to the beginning of the vector. */
    ConstReverseIterator rbegin() const { return ConstReverseIterator(end()); }

    /** Get an iterator pointing to the end of the vector. */
    ReverseIterator rend() { return ReverseIterator(begin()); }

    /** Get a const iterator pointing to the end of the vector. */
    ConstReverseIterator rend() const { return ConstReverseIterator(begin()); }

    /** Access an element of the vector immutably. */
    template <typename IndexT> T const & operator[](IndexT i) const { return values[i]; }

    /** Access an element of the vector mutably. */
    template <typename IndexT> T & operator[](IndexT i) { return values[i]; }

    /** Access an element of the vector immutably. */
    template <typename IndexT> T const & operator()(IndexT i) const { return values[i]; }

    /** Access an element of the vector mutably. */
    template <typename IndexT> T & operator()(IndexT i) { return values[i]; }

    /** Equality test. */
    bool operator==(VectorT const & rhs) const
    {
      for (long i = 0; i < N; ++i)
        if (values[i] != rhs[i])
          return false;

      return true;
    }

    /** Inequality test. */
    bool operator!=(VectorT const & rhs) const
    {
      return !operator==(rhs);
    }

    /** Less-than operator (lexicographical comparison in order of increasing index). */
    bool operator<(VectorT const & rhs) const
    {
      return std::lexicographical_compare(begin(), end(), rhs.begin(), rhs.end());
    }

    /** Greater-than operator (lexicographical comparison in order of increasing index). */
    bool operator>(VectorT const & rhs) const
    {
      return std::lexicographical_compare(rhs.begin(), rhs.end(), begin(), end());
    }

    /** Less-than-or-equal-to operator (lexicographical comparison in order of increasing index). */
    bool operator<=(VectorT const & rhs) const
    {
      return !(*this > rhs);
    }

    /** Greater-than-or-equal-to operator (lexicographical comparison in order of increasing index). */
    bool operator>=(VectorT const & rhs) const
    {
      return !(rhs > *this);
    }

    /** Negation. */
    VectorT operator-() const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = -values[i];

      return result;
    }

    /** Addition (per-component). */
    VectorT operator+(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = values[i] + rhs[i];

      return result;
    }

    /** Subtraction (per-component). */
    VectorT operator-(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = values[i] - rhs[i];

      return result;
    }

    /** Multiplication (per-component). */
    VectorT operator*(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = values[i] * rhs[i];

      return result;
    }

    /** Multiply by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type operator*(S const & s) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = static_cast<T>(values[i] * s);

      return result;
    }

    /** Division (per-component). */
    VectorT operator/(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = values[i] / rhs[i];

      return result;
    }

    /** Divide by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type operator/(S const & s) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = static_cast<T>(values[i] / s);

      return result;
    }

    /** Add-and-assign. */
    VectorT & operator+=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        values[i] += rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Subtract-and-assign. */
    VectorT & operator-=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        values[i] -= rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Multiply-and-assign. */
    VectorT & operator*=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        values[i] *= rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Multiply in-place by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type & operator*=(S const & s)
    {
      for (long i = 0; i < N; ++i)
        values[i] = static_cast<T>(values[i] * s);

      return *static_cast<VectorT *>(this);
    }

    /** Divide-and-assign. */
    VectorT & operator/=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        values[i] /= rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Divide in-place by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type & operator/=(S const & s)
    {
      for (long i = 0; i < N; ++i)
        values[i] = static_cast<T>(values[i] / s);

      return *static_cast<VectorT *>(this);
    }

    /** Dot product. */
    T dot(VectorT const & rhs) const
    {
      T result = values[0] * rhs[0];  // safe if class precondition N > 0 is met
      for (long i = 1; i < N; ++i)
        result += (values[i] * rhs[i]);

      return result;
    }

    /**
     * Construct a matrix as the outer product <code>uv^T</code> of this vector u and a second vector \a v. The element (i, j)
     * of the matrix is u[i] * v[j]. This function is not defined unless MatrixMN.hpp is included.
     *
     * @note The same effect can be achieved (with two extra memory copies) by casting the vectors to a column and a row matrix
     *   -- toColumnMatrix() and toRowMatrix() -- respectively, and then multiplying them.
     *
     * @param v The second vector in the product.
     */
    MatrixMN<N, N, T> outerProduct(VectorT const & v) const;

    /**
     * Get the minimum component in the vector (according to signed comparison). To compare absolute values, use minAbs()
     * instead.
     */
    T const & min() const { return values[minAxis()]; }

    /**
     * Get the (signed) component with the minimum absolute value.
     *
     * @see min()
     */
    T const & minAbs() const { return values[minAbsAxis()]; }

    /**
     * Get the maximum component in the vector (according to signed comparison). To compare absolute values, use maxAbs()
     * instead.
     */
    T const & max() const { return values[maxAxis()]; }

    /**
     * Get the (signed) component with the maximum absolute value.
     *
     * @see max()
     */
    T const & maxAbs() const { return values[maxAbsAxis()]; }

    /**
     * Get the index of the axis of the vector with the minimum coordinate (according to signed comparison). To compare absolute
     * values, use minAbsAxis() instead.
     */
    long minAxis() const
    {
      long min_axis = 0;
      for (long i = 1; i < N; ++i)
      {
        if (values[i] < values[min_axis])
          min_axis = i;
      }

      return min_axis;
    }

    /**
     * Get the index of the axis of the vector with the numerically smallest coordinate (the one with the smallest absolute
     * value).
     *
     * @see minAxis()
     */
    long minAbsAxis() const
    {
      long min_axis = 0;
      T abs_pc = std::abs(values[0]);
      for (long i = 1; i < N; ++i)
      {
        T abs_i = std::abs(values[i]);
        if (abs_i < abs_pc)
        {
          min_axis = i;
          abs_pc = abs_i;
        }
      }

      return min_axis;
    }

    /**
     * Get the index of the axis of the vector with the maximum coordinate (according to signed comparison). To compare absolute
     * values, use maxAbsAxis() instead.
     */
    long maxAxis() const
    {
      long max_axis = 0;
      for (long i = 1; i < N; ++i)
      {
        if (values[i] > values[max_axis])
          max_axis = i;
      }

      return max_axis;
    }

    /**
     * Get the index of the axis of the vector with the numerically largest coordinate (the one with the largest absolute
     * value).
     *
     * @see maxAxis()
     */
    long maxAbsAxis() const
    {
      long max_axis = 0;
      T abs_pc = std::abs(values[0]);
      for (long i = 1; i < N; ++i)
      {
        T abs_i = std::abs(values[i]);
        if (abs_i > abs_pc)
        {
          max_axis = i;
          abs_pc = abs_i;
        }
      }

      return max_axis;
    }

    /** Return a vector containing the component-wise minima of this vector and another. */
    VectorT min(VectorT const & other) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = std::min(values[i], other[i]);

      return result;
    }

    /** Return a vector containing the component-wise maxima of this vector and another. */
    VectorT max(VectorT const & other) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = std::max(values[i], other[i]);

      return result;
    }

    /** Get the square of the L2 length of the vector. */
    T squaredLength() const { return this->dot(*static_cast<VectorT const *>(this)); }

    /** Get the length (L2-norm) of the vector. */
    T length() const { return static_cast<T>(std::sqrt(squaredLength())); }

    /** Get the length (L2-norm) of the vector, using a fast approximation to the square root function. */
    T fastLength() const { return static_cast<T>(Math::fastSqrt(squaredLength())); }

    /** Get a unit vector along the same direction. */
    VectorT unit() const
    {
      T len = length();
      if (std::abs(len) < 32 * std::numeric_limits<T>::min())
        return zero();
      else
        return *this / len;
    }

    /** Get a unit vector along the same direction, using a fast approximation to the reciprocal of the square root function. */
    VectorT fastUnit() const
    {
      return *this * Math::fastRsq(squaredLength());
    }

    /** Normalize the vector to have unit length. */
    void unitize()
    {
      *this = unit();
    }

    /** Normalize the vector to have unit length, using a fast approximation to the reciprocal of the square root function. */
    void fastUnitize()
    {
      *this *= Math::fastRsq(squaredLength());
    }

    /** Get a textual representation of the vector. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << '(' << values[0];

      for (long i = 1; i < N; ++i)
        oss << ", " << values[i];

      oss << ')';
      return oss.str();
    }

    /** Get a vector containing only zeroes. */
    static VectorT const & zero() { static VectorT const z(0); return z; }

  private:
    T values[N];  ///< Vector values.

}; // class VectorNBase

} // namespace Internal

/**
 * Fixed-size N-dimensional vectors, where N is any <b>positive</b> (non-zero) integer and T is a field. Implemented as a
 * contiguous array of N values of type T.
 */
template <long N, typename T>
class /* THEA_API */ VectorN : public Internal::VectorNBase<N, T>
{
  private:
    typedef Internal::VectorNBase<N, T> BaseT;

  public:
    /** Default constructor (does not initialize anything). */
    VectorN() {}

    /** Initialize all components to a single value. */
    explicit VectorN(T const & fill_value) : BaseT(fill_value) {}

    /** Copy constructor. */
    template <typename U> VectorN(VectorN<N, U> const & src) : BaseT(src) {}

    /** Initialize from a column matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<N, 1, U> const & src) : BaseT(src) {}

    /** Initialize from a row matrix (not defined unless MatrixMN.hpp is included). */
    template <typename U> explicit VectorN(MatrixMN<1, N, U> const & src) : BaseT(src) {}

}; // class VectorN

/** Pre-multiply an N-dimensional vector by a scalar. */
template <typename S, long N, typename T>
typename boost::enable_if< Internal::ScalarCheck<S, T>, VectorN<N, T> >::type
operator*(S const & s, VectorN<N, T> const & v)
{
  return v * s;
}

/** Pipe a textual representation of a N-dimensional vector to a <code>std::ostream</code>. */
template <long N, typename T>
std::ostream &
operator<<(std::ostream & os, VectorN<N, T> const & v)
{
  return os << v.toString();
}

} // namespace Thea

#include "Vector2.hpp"
#include "Vector3.hpp"
#include "Vector4.hpp"

#endif
