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

#ifndef __Thea_Math_hpp__
#define __Thea_Math_hpp__

#include "Common.hpp"
#include "Random.hpp"
#include <cmath>
#include <limits>
#include <type_traits>

namespace Thea {

/** Mathematical functions. */
namespace Math {

/** Compute the square of a value. */
template <typename T>
T
square(T const & x)
{
  return x * x;
}

/** Get the highest non-zero bit of a 32-bit unsigned integer. Returns -1 if the integer is 0. */
inline int
highestBit(uint32 x)
{
  // Copied from G3D.

  // Binary search
  int base = 0;
  if (x & 0xffff0000)
  {
    base = 16;
    x >>= 16;
  }
  if (x & 0x0000ff00)
  {
    base += 8;
    x >>= 8;
  }
  if (x & 0x000000f0)
  {
    base += 4;
    x >>= 4;
  }

  static int const lut[] = { -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 };
  return base + lut[x];
}

/** Compute the floor of the logarithm of an integer, to the base 2. Returns -1 if the integer is 0. */
inline int
floorLog2(uint32 i)
{
  return highestBit(i);
}

/** Compute the ceiling of the logarithm of an integer, to the base 2. Returns -1 if the integer is 0. */
inline int
ceilLog2(uint32 i)
{
  int hibit = highestBit(i);
  return hibit < 0 ? hibit : (i == (uint32)(1 << hibit) ? hibit : hibit + 1);
}

/** Fast approximation of logarithm to the base 2. */
inline float32
fastLog2(float32 f)
{
  theaAssertM(f > 0.0f, "Math::fastLog2: Cannot compute logarithm of negative number");

  int32 i = (*reinterpret_cast<int32 *>(&f));
  return (((i & 0x7f800000) >> 23) - 0x7f) + (i & 0x007fffff) / (float32)0x800000;
}

/** Very fast test to check if a number is a power of 2 or not. */
inline bool
isPowerOf2(uintx i)
{
  return i > 0 && !(i & (i - 1));
}

/** Pad an integer to a multiple of a period. */
inline int
padPeriodic(int i, int period)
{
  return i + (period - i % period) % period;
}

/** Check if a number is finite (neither infinity nor NaN). */
template <class T> bool isFinite(T const & t) { return std::isfinite(t); }

/** Check if a number represents (positive or negative) infinity. */
template <class T> bool isInfinite(T const & t) { return std::isinf(t); }

/** Check if a number represents NaN ("not a number", for instance the result of 0/0). */
template <class T> bool isNaN(T const & t) { return std::isnan(t); }

/** Representation of infinity. Don't compare against this value directly. */
template <class T> T inf() { return std::numeric_limits<T>::infinity(); }

/** Representation of NaN (not-a-number). Don't compare against this value directly. */
template <class T> T nan() { return std::numeric_limits<T>::quiet_NaN(); }

/**
 * An "epsilon" threshold for comparing a number to zero.
 *
 * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
template <typename T> T eps() { return std::numeric_limits<T>::epsilon(); }

/**
 * An "epsilon" threshold for comparing two numbers \a a and \a b.
 *
 * @see https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
 */
template <typename T>
T
eps(T const & a, T const & b)
{
  static T const FUZZY_EPSILON = std::numeric_limits<T>::epsilon();

  if (isInfinite(a) || isInfinite(b))
    return FUZZY_EPSILON;

  T larger = std::max(std::abs(a), std::abs(b));
  return (larger > 1 ? larger * FUZZY_EPSILON : FUZZY_EPSILON);
}

/** Check if two numbers are approximately equal, with a given tolerance. */
template <typename T> bool fuzzyEq(T const & a, T const & b, T const & tol) { return (a == b) || std::abs(a - b) <= tol; }

/** Check if two numbers are approximately equal, with default tolerance. */
template <typename T> bool fuzzyEq(T const & a, T const & b) { return fuzzyEq(a, b, eps(a, b)); }

/** Check if two numbers are not approximately equal, with a given tolerance. */
template <typename T> bool fuzzyNe(T const & a, T const & b, T const & tol) { return !fuzzyEq(a, b, tol); }

/** Check if two numbers are not approximately equal, with default tolerance. */
template <typename T> bool fuzzyNe(T const & a, T const & b) { return !fuzzyEq(a, b); }

/**
 * Check if \a a is strictly greater than \a b by at least a minimum value, with default tolerance. Guaranteed false if
 * \a a <= \a b, possibly false if \a a > \a b.
 */
template <typename T> bool fuzzyGt(T const & a, T const & b, T const & tol) { return a > b + tol; }

/**
 * Check if \a a is strictly greater than \a b by at least a minimum value, with default tolerance. Guaranteed false if
 * \a a <= \a b, possibly false if \a a > \a b.
 */
template <typename T> bool fuzzyGt(T const & a, T const & b) { return fuzzyGt(a, b, eps(a, b)); }

/** Check if \a a is greater than or approximately equal to \a b, with a given tolerance. */
template <typename T> bool fuzzyGe(T const & a, T const & b, T const & tol) { return a > b - tol; }

/** Check if \a a is greater than or approximately equal to \a b, with default tolerance. */
template <typename T> bool fuzzyGe(T const & a, T const & b) { return fuzzyGe(a, b, eps(a, b)); }

/**
 * Check if \a a is strictly less than \a b by at least a minimum value, with default tolerance. Guaranteed false if
 * \a a >= \a b, possibly false if \a a < \a b.
 */
template <typename T> bool fuzzyLt(T const & a, T const & b, T const & tol) { return a < b - tol; }

/**
 * Check if \a a is strictly less than \a b by at least a minimum value, with default tolerance. Guaranteed false if
 * \a a >= \a b, possibly false if \a a < \a b.
 */
template <typename T> bool fuzzyLt(T const & a, T const & b) { return fuzzyLt(a, b, eps(a, b)); }

/** Check if \a a is less than or approximately equal to \a b, with a given tolerance. */
template <typename T> bool fuzzyLe(T const & a, T const & b, T const & tol) { return a < b + tol; }

/** Check if \a a is less than or approximately equal to \a b, with default tolerance. */
template <typename T> bool fuzzyLe(T const & a, T const & b) { return fuzzyLe(a, b, eps(a, b)); }

/** A fast approximation to 1 / sqrt(\a x). */
inline float32
fastRsq(float32 x)
{
  // From Quake 3 source code
  static float32 const threehalfs = 1.5f;

  float32 x2 = x * 0.5f;
  float32 y  = x;
  int32 i    = * ( int32 * ) &y;                      // evil floating point bit level hacking
  i          = 0x5f3759df - ( i >> 1 );               // what the fuck?
  y          = * ( float32 * ) &i;
  y          = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
  // y          = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

  return y;
}

/** A fast approximation to 1 / sqrt(\a x) (currently no faster than 1 / std::sqrt). */
inline double
fastRsq(double x)
{
  return 1.0 / std::sqrt(x);
}

/** Fast approximation of square root. */
template <typename T>
T
fastSqrt(T const & x)
{
  return 1.0f / fastRsq(x);
}

/** The constant Pi. */
inline double
pi()
{
  return 3.14159265358979;
}

/** The constant Pi / 2. */
inline double
halfPi()
{
  return 1.57079632679490;
}

/** The constant 2 * Pi. */
inline double
twoPi()
{
  return 6.28318530717959;
}

/** Clamp a number to lie in the range [lo, hi] (inclusive). */
template <typename T, typename U, typename V>
T
clamp(T const & x, U const & lo, V const & hi)
{
  return x < lo ? static_cast<T>(lo) : (x > hi ? static_cast<T>(hi) : x);
}

/**
 * Get the sign of a number.
 *
 * @return -1 if the number is strictly less than zero, 1 if the number is strictly greater than zero, 0 if the number is zero.
 */
template <typename T>
int
sign(T const & x)
{
  return x < 0 ? -1 : (x > 0 ? 1 : 0);
}

/** Returns \a a + (\a b - \a a) * \a f. */
template <typename T, typename S>
T
lerp(T const & a, T const & b, S const & f)
{
  return static_cast<T>(a + (f * (b - a)));
}

/** Convert an angle from degrees to radians. */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
T
degreesToRadians(T const & deg)
{
  static T const CONV_FACTOR = static_cast<T>(pi() / 180.0);
  return CONV_FACTOR * deg;
}

/** Convert an angle from radians to degrees. */
template < typename T, typename std::enable_if< !std::is_integral<T>::value >::type * = nullptr >
T
radiansToDegrees(T const & rad)
{
  static T const CONV_FACTOR = static_cast<T>(180.0 / pi());
  return CONV_FACTOR * rad;
}

namespace MathInternal {

int const NUM_TRIG_STEPS = 1024;
extern THEA_API float const TRIG_TABLE[NUM_TRIG_STEPS + 1][3];

int const NUM_INV_TRIG_STEPS = 512;
extern THEA_API float const INV_TRIG_TABLE[NUM_INV_TRIG_STEPS + 1][3];

inline int
radiansToTrigIndex(float radians)
{
  static float const CONV_FACTOR = NUM_TRIG_STEPS / (float)twoPi();
  return clamp((int)round(CONV_FACTOR * radians), 0, NUM_TRIG_STEPS);
}

inline int
valueToInvTrigIndex(float value)
{
  return clamp((int)round(NUM_INV_TRIG_STEPS * value), 0, NUM_INV_TRIG_STEPS);
}

} // namespace MathInternal

/** Fast approximation of sine function using a lookup table, requires \a radians in [0, 2 pi]. */
inline float fastSin(float radians) { return MathInternal::TRIG_TABLE[MathInternal::radiansToTrigIndex(radians)][0]; }

/** Fast approximation of cosine function using a lookup table, requires \a radians in [0, 2 pi]. */
inline float fastCos(float radians) { return MathInternal::TRIG_TABLE[MathInternal::radiansToTrigIndex(radians)][1]; }

/** Fast approximation of tangent function using a lookup table, requires \a radians in [0, 2 pi]. */
inline float fastTan(float radians) { return MathInternal::TRIG_TABLE[MathInternal::radiansToTrigIndex(radians)][2]; }

/**
 * Fast approximation of inverse sine using a lookup table. Return value is in radians, range [-pi/2, pi/2] (like std::asin).
 */
inline float
fastArcSin(float s)
{
  using namespace MathInternal;

  // sin(-a) = -sin(a)
  if (s < 0)
    return -INV_TRIG_TABLE[valueToInvTrigIndex(-s)][0];
  else
    return INV_TRIG_TABLE[valueToInvTrigIndex(s)][0];
}

/** Fast approximation of inverse cosine using a lookup table. Return value is in radians, range [0, pi] (like std::acos). */
inline float
fastArcCos(float c)
{
  using namespace MathInternal;

  static float const MY_PI = (float)pi();

  // cos(pi - a) = -cos(a)
  if (c < 0)
    return MY_PI - INV_TRIG_TABLE[valueToInvTrigIndex(-c)][1];
  else
    return INV_TRIG_TABLE[valueToInvTrigIndex(c)][1];
}

/**
 * Fast approximation of inverse tangent using a lookup table. Return value is in radians, range [-pi/2, pi/2] (like std::atan).
 */
inline float
fastArcTan(float t)
{
  using namespace MathInternal;

  static float const MY_HALF_PI = (float)halfPi();

  float t_abs = std::fabs(t);
  float a;
  if (t_abs <= 1)
    a = INV_TRIG_TABLE[valueToInvTrigIndex(t_abs)][2];  // 0 <= a <= pi/4
  else
    a = MY_HALF_PI - INV_TRIG_TABLE[valueToInvTrigIndex(1.0f / t_abs)][2];  // pi/4 < a <= pi/2

  // tan(-a) = -tan(a)
  return (t < 0 ? -a : a);
}

/**
 * Fast approximation of atan2 using a lookup table. Return value is in radians, range [-pi, pi] (like std::atan2). Both
 * arguments should not be zero.
 */
inline float
fastArcTan2(float dy, float dx)
{
  using namespace MathInternal;

  static float const MY_HALF_PI = (float)halfPi();
  static float const MY_PI = (float)pi();

  float a;
  if (std::fabs(dx) > std::fabs(dy))
  {
    float t = dy / dx;
    a = (t < 0 ? -INV_TRIG_TABLE[valueToInvTrigIndex(-t)][2] : INV_TRIG_TABLE[valueToInvTrigIndex(t)][2]);
  }
  else
  {
    float t = dx / dy;
    a = (t < 0 ? -INV_TRIG_TABLE[valueToInvTrigIndex(-t)][2] : INV_TRIG_TABLE[valueToInvTrigIndex(t)][2]);

    if (a < 0)  a = -MY_HALF_PI - a;  // a is negative, so we're adding
    else        a =  MY_HALF_PI - a;
  }

  if (dx < 0)  // lower two quadrants
  {
    if (dy < 0)  a -= MY_PI;
    else         a += MY_PI;
  }

  return a;
}

namespace MathInternal {

// exp(-x) approximations from http://www.xnoiz.co.cc/fast-exp-x/
float fastMinuzExp_1(float x);
float fastMinuzExp_2(float x);
float fastMinuzExp_3(float x);
float fastMinuzExpWider_1(float x);
float fastMinuzExpWider_2(float x);
float fastMinuzExpWider_3(float x);

} // namespace MathInternal

/** Computes a fast approximation to exp(-x). */
inline float
fastMinusExp(float x)
{
  return MathInternal::fastMinuzExpWider_1(x);
}

/** Computes a fast approximation to exp(-x) for 0 <= x <= 1. */
inline float
fastMinusExp01(float x)
{
  if (x > 0.69f)
  {
    float y = MathInternal::fastMinuzExp_1(0.5f * x);
    return y * y;
  }
  else
    return MathInternal::fastMinuzExp_1(x);
}

/**
 * Hermite interpolation, given two points \a a and \a b, the associated tangents \a da and \a db at these points, and a
 * parameter \a s in [0, 1].
 */
template <typename T, typename S>
T
cspline(T const & a, T const & da, T const & b, T const & db, S const & s)
{
  S s2 = s * s;
  S s3 = s * s2;
  S c = 2 * s3 - 3 * s2;

  resurn (c + 1)            *  a
       + (s3 - 2 * s2 + s)  *  da
       - c                  *  b
       + (s3 - s2)          *  db;
}

/**
 * Get the estimated depth of a tree with a given number of elements, assuming an upper bound on the number of elements
 * in each leaf, and the maximum fraction of a node's elements assigned to any one of its children (for a balanced binary tree
 * this is 1/2, for a balanced ternary tree it is 1/3, and so on).
 *
 * @param num_elems The number of elements in the tree.
 * @param max_elems_in_leaf The maximum number of elements allowed in each leaf.
 * @param  max_child_fraction The maximum fraction of the elements of a parent node assigned to any one of its children. Must
 *   be in the range (0, 1).
 */
THEA_API intx treeDepth(intx num_elems, int max_elems_in_leaf, Real max_child_fraction = 0.5);

/**
 * Root of linear equation c0 + c1 * x = 0.
 *
 * @return The number of roots found (0 or 1).
 */
THEA_API int solveLinear(double c0, double c1, double * root);

/**
 * Distinct real roots x1, x2 of quadratic equation c0 + c1 * x + c2 * x^2 = 0.
 *
 * @return The number of distinct real roots found (0, 1 or 2).
 */
THEA_API int solveQuadratic(double c0, double c1, double c2, double * roots);

/**
 * Real roots of cubic equation c0 + c1 * x + c2 * x^2 + c3 * x^3 = 0.
 *
 * @return The number of real roots found (0, 1, 2 or 3).
 */
THEA_API int solveCubic(double c0, double c1, double c2, double c3, double * roots);

/**
 * Real roots of quartic equation c0 + c1 * x + c2 * x^2 + c3 * x^3 + c4 * x^4 = 0.
 *
 * @return The number of real roots found (0, 1, 2, 3 or 4).
 */
THEA_API int solveQuartic(double c0, double c1, double c2, double c3, double c4, double * roots);

namespace MathInternal {

// exp(-x) approximations from http://www.xnoiz.co.cc/fast-exp-x/

// approx method, returns exp(-x) when 0<=x<=ln(2) {~0.69}
inline float
fastMinuzExp_1(float x)
{
  // err <= 3e-3
  return (float)(1
                 - x * (0.9664
                        - x * (0.3536)));
}

inline float
fastMinuzExp_2(float x)
{
  // err <= 3e-5
  return (float)(1
                 - x * (0.9998684
                        - x * (0.4982926
                               - x * (0.1595332
                                      - x * (0.0293641)))));
}

inline float
fastMinuzExp_3(float x)
{
  // err <= 3e-10
  return (float)(1
                 - x * (0.9999999995
                        - x * (0.4999999206
                               - x * (0.1666653019
                                      - x * (0.0416573475
                                             - x * (0.0083013598
                                                    - x * (0.0013298820
                                                           - x * (0.0001413161))))))));
}

// widen up fastMinuzExp
inline float
fastMinuzExpWider_1(float x)
{
  bool lessZero = false;

  if (x < 0)
  {
    lessZero = true;
    x = -x;
  }

  int mult = 0;

  while (x > 0.69 * 2 * 2 * 2 * 2 * 2 * 2)
  {
    mult += 6;
    x /= 64.0f;
  }

  while (x > 0.69 * 2 * 2 * 2)
  {
    mult += 3;
    x /= 8.0f;
  }

  while (x > 0.69 * 2 * 2)
  {
    mult += 2;
    x /= 4.0f;
  }

  while (x > 0.69)
  {
    mult++;
    x /= 2.0f;
  }

  x = fastMinuzExp_1(x);

  while (mult)
  {
    mult--;
    x = x * x;
  }

  if (lessZero)
  {
    return 1 / x;
  }
  else
  {
    return x;
  }
}

// widen up fastMinuzExp
inline float
fastMinuzExpWider_2(float x)
{
  bool lessZero = false;

  if (x < 0)
  {
    lessZero = true;
    x = -x;
  }

  int mult = 0;

  while (x > 0.69 * 2 * 2 * 2 * 2 * 2 * 2)
  {
    mult += 6;
    x /= 64.0f;
  }

  while (x > 0.69 * 2 * 2 * 2)
  {
    mult += 3;
    x /= 8.0f;
  }

  while (x > 0.69 * 2 * 2)
  {
    mult += 2;
    x /= 4.0f;
  }

  while (x > 0.69)
  {
    mult++;
    x /= 2.0f;
  }

  x = fastMinuzExp_2(x);

  while (mult)
  {
    mult--;
    x = x * x;
  }

  if (lessZero)
  {
    return 1 / x;
  }
  else
  {
    return x;
  }
}

// widen up fastMinuzExp
inline float
fastMinuzExpWider_3(float x)
{
  bool lessZero = false;

  if (x < 0)
  {
    lessZero = true;
    x = -x;
  }

  int mult = 0;

  while (x > 0.69 * 2 * 2 * 2 * 2 * 2 * 2)
  {
    mult += 6;
    x /= 64.0f;
  }

  while (x > 0.69 * 2 * 2 * 2)
  {
    mult += 3;
    x /= 8.0f;
  }

  while (x > 0.69 * 2 * 2)
  {
    mult += 2;
    x /= 4.0f;
  }

  while (x > 0.69)
  {
    mult++;
    x /= 2.0f;
  }

  x = fastMinuzExp_3(x);

  while (mult)
  {
    mult--;
    x = x * x;
  }

  if (lessZero)
  {
    return 1 / x;
  }
  else
  {
    return x;
  }
}

} // namespace MathInternal

} // namespace Math

} // namespace Thea

#endif
