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
#include "Math.hpp"
#include <boost/array.hpp>
#include <limits>
#include <sstream>

namespace Thea {

// Forward declarations
template <long N, typename T> class VectorN;

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
 * a field. Implemented as a wrapper for <code>boost::array</code>, with the same interface plus arithmetic operators.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <long N, typename T>
class /* THEA_DLL_LOCAL */ VectorNBase : public boost::array<T, N>
{
  private:
    typedef typename boost::array<T, N> BaseT;

  public:
    typedef VectorN<N, T>                           VectorT;     ///< N-dimensional vector type.
    typedef T                                       Value;       ///< Type of values stored in the vector.
    typedef T                                       value_type;  ///< Type of values stored in the vector (STL convention).
    typedef typename BaseT::iterator                Iterator;              ///< Forward iterator through elements.
    typedef typename BaseT::const_iterator          ConstIterator;         ///< Forward const iterator through elements.
    typedef typename BaseT::reverse_iterator        ReverseIterator;       ///< Reverse iterator through elements.
    typedef typename BaseT::const_reverse_iterator  ConstReverseIterator;  ///< Reverse const iterator through elements.

    THEA_DEF_POINTER_TYPES(VectorT, shared_ptr, weak_ptr)

    /** Default constructor (does not initialize anything). */
    VectorNBase() {}

    /** Initialize all components to a single value. */
    explicit VectorNBase(T const & fill_value) { BaseT::assign(fill_value); }

    /** Copy constructor. */
    template <typename U> VectorNBase(VectorNBase<N, U> const & src) { BaseT::operator=(src); }

    /** Equality test. */
    bool operator==(VectorT const & rhs) const
    {
      for (long i = 0; i < N; ++i)
        if ((*this)[i] != rhs[i])
          return false;

      return true;
    }

    /** Inequality test. */
    bool operator!=(VectorT const & rhs) const
    {
      return !operator==(rhs);
    }

    /** Negation. */
    VectorT operator-() const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = -(*this)[i];

      return result;
    }

    /** Addition (per-component). */
    VectorT operator+(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = (*this)[i] + rhs[i];

      return result;
    }

    /** Subtraction (per-component). */
    VectorT operator-(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = (*this)[i] - rhs[i];

      return result;
    }

    /** Multiplication (per-component). */
    VectorT operator*(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = (*this)[i] * rhs[i];

      return result;
    }

    /** Multiply by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type operator*(S const & s) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = static_cast<T>((*this)[i] * s);

      return result;
    }

    /** Division (per-component). */
    VectorT operator/(VectorT const & rhs) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = (*this)[i] / rhs[i];

      return result;
    }

    /** Divide by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type operator/(S const & s) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = static_cast<T>((*this)[i] / s);

      return result;
    }

    /** Add-and-assign. */
    VectorT & operator+=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        (*this)[i] += rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Subtract-and-assign. */
    VectorT & operator-=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        (*this)[i] -= rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Multiply-and-assign. */
    VectorT & operator*=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        (*this)[i] *= rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Multiply in-place by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type & operator*=(S const & s)
    {
      for (long i = 0; i < N; ++i)
        (*this)[i] = static_cast<T>((*this)[i] * s);

      return *static_cast<VectorT *>(this);
    }

    /** Divide-and-assign. */
    VectorT & operator/=(VectorT const & rhs)
    {
      for (long i = 0; i < N; ++i)
        (*this)[i] /= rhs[i];

      return *static_cast<VectorT *>(this);
    }

    /** Divide in-place by a scalar. */
    template <typename S> typename boost::enable_if< ScalarCheck<S, T>, VectorT >::type & operator/=(S const & s)
    {
      for (long i = 0; i < N; ++i)
        (*this)[i] = static_cast<T>((*this)[i] / s);

      return *static_cast<VectorT *>(this);
    }

    /** Dot product. */
    T dot(VectorT const & rhs) const
    {
      T result = (*this)[0] * rhs[0];  // safe if class precondition N > 0 is met
      for (long i = 1; i < N; ++i)
        result += ((*this)[i] * rhs[i]);

      return result;
    }

    /**
     * Get the minimum component in the vector (according to signed comparison). To compare absolute values, use minAbs()
     * instead.
     */
    T const & min() const { return (*this)[minAxis()]; }

    /**
     * Get the (signed) component with the minimum absolute value.
     *
     * @see min()
     */
    T const & minAbs() const { return (*this)[minAbsAxis()]; }

    /**
     * Get the maximum component in the vector (according to signed comparison). To compare absolute values, use maxAbs()
     * instead.
     */
    T const & max() const { return (*this)[maxAxis()]; }

    /**
     * Get the (signed) component with the maximum absolute value.
     *
     * @see max()
     */
    T const & maxAbs() const { return (*this)[maxAbsAxis()]; }

    /**
     * Get the index of the axis of the vector with the minimum coordinate (according to signed comparison). To compare absolute
     * values, use minAbsAxis() instead.
     */
    long minAxis() const
    {
      long min_axis = 0;
      for (long i = 1; i < N; ++i)
      {
        if ((*this)[i] < (*this)[min_axis])
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
      T abs_pc = std::abs((*this)[0]);
      for (long i = 1; i < N; ++i)
      {
        T abs_i = std::abs((*this)[i]);
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
        if ((*this)[i] > (*this)[max_axis])
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
      T abs_pc = std::abs((*this)[0]);
      for (long i = 1; i < N; ++i)
      {
        T abs_i = std::abs((*this)[i]);
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
        result[i] = std::min((*this)[i], other[i]);

      return result;
    }

    /** Return a vector containing the component-wise maxima of this vector and another. */
    VectorT max(VectorT const & other) const
    {
      VectorT result;
      for (long i = 0; i < N; ++i)
        result[i] = std::max((*this)[i], other[i]);

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
      oss << '(' << (*this)[0];

      for (long i = 1; i < N; ++i)
        oss << ", " << (*this)[i];

      oss << ')';
      return oss.str();
    }

    /** Get a vector containing only zeroes. */
    static VectorT const & zero() { static VectorT const z(0); return z; }

}; // class VectorNBase

} // namespace Internal

/**
 * Fixed-size N-dimensional vectors, where N is any <b>positive</b> (non-zero) integer and T is a field. Implemented as a
 * wrapper for <code>boost::array</code>, with the same interface plus arithmetic operators.
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
