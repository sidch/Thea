//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

/*
 ORIGINAL HEADER

 @file Random.cpp

 @maintainer Morgan McGuire, http://graphics.cs.williams.edu

 @created 2009-01-02
 @edited  2009-03-29

 Copyright 2000-2009, Morgan McGuire.
 All rights reserved.
*/

#include "Random.hpp"
#include "Array.hpp"
#include "Math.hpp"
#include <ctime>
#include <limits>

#ifdef THEA_WINDOWS
#  include <windows.h>
#else
#  include <unistd.h>
#endif

namespace Thea {

Random::Random(void * x)
: state(NULL), m_threadsafe(false)
{
  (void)x;
}

Random::Random(uint32 seed, bool threadsafe)
: m_threadsafe(threadsafe)
{
  uint32 const X = 1812433253UL;
  state = new uint32[N];
  state[0] = seed;

  for (index = 1; index < (int32)N; ++index)
  {
    state[index] = X * (state[index - 1] ^ (state[index - 1] >> 30)) + index;
  }
}

Random::~Random()
{
  delete [] state;
  state = NULL;
}

uint32
Random::bits()
{
  // See http://en.wikipedia.org/wiki/Mersenne_twister

  // Make a local copy of the index variable to ensure that it is not out of bounds
  int32 localIndex = index;

  // Automatically checks for index < 0 if corrupted by unsynchronized threads.
  if ((uint32)localIndex >= (uint32)N)
  {
    generate();
    localIndex = 0;
  }

  // Increment the global index.  It may go out of bounds on multiple threads, but the above check ensures that the array index
  // actually used never goes out of bounds. It doesn't matter if we grab the same array index twice on two threads, since the
  // distribution of random numbers will still be uniform.
  ++index;

  // Return the next random in the sequence
  uint32 r = state[localIndex];

  // Temper the result
  r ^=  r >> U;
  r ^= (r << S) & B;
  r ^= (r << T) & C;
  r ^=  r >> L;
  return r;
}

void
Random::generate()
{
  // Lower R bits
  static uint32 const LOWER_MASK = (1LU << R) - 1;

  // Upper (32 - R) bits
  static uint32 const UPPER_MASK = 0xFFFFFFFF << R;
  static uint32 const mag01[2] = { 0UL, (uint32)A };

  if (m_threadsafe)
    lock.lock();

  // First N - M bits
  for (unsigned int i = 0; i < N - M; ++i)
  {
    uint32 x = (state[i] & UPPER_MASK) | (state[i + 1] & LOWER_MASK);
    state[i] = state[i + M] ^ (x >> 1) ^ mag01[x & 1];
  }

  // Remaining bits
  for (unsigned int i = N - M + 1; i < N - 1; ++i)
  {
    uint32 x = (state[i] & UPPER_MASK) | (state[i + 1] & LOWER_MASK);
    state[i] = state[i + (M - N)] ^ (x >> 1) ^ mag01[x & 1];
  }

  uint32 y = (state[N - 1] & UPPER_MASK) | (state[0] & LOWER_MASK);
  state[N - 1] = state[M - 1] ^ (y >> 1) ^ mag01[y & 1];
  index = 0;

  if (m_threadsafe)
    lock.unlock();
}

bool
Random::coinToss()
{
  // Probably not as fast as generating just the next bit but will do for now
  return (bool)(bits() & 0x01);
}

int32
Random::integer()
{
  return (int32)(bits() & MAX_INTEGER);
}

int32
Random::integer(int32 low, int32 high)
{
  int32 r = (int32)std::floor(low + (high - low + 1) * (double)bits() / 0xFFFFFFFFUL);

  // There is a *very small* chance of generating a number larger than high.
  if (r > high)
  {
    return high;
  }
  else
  {
    return r;
  }
}

Real
Random::gaussian(Real mean, Real stddev)
{
  // Using Box-Mueller method from http://www.taygeta.com/random/gaussian.html
  // Modified to specify standard deviation and mean of distribution
  Real w, x1, x2;

  // Loop until w is less than 1 so that log(w) is negative
  do
  {
    x1 = uniform(-1.0f, 1.0f);
    x2 = uniform(-1.0f, 1.0f);
    w = (Real)(Math::square(x1) + Math::square(x2));
  }
  while (w > 1.0f);

  // Transform to gassian distribution. Multiply by the variance \a sigma (stddev^2) and add the \a mean.
  return x2 * (Real)Math::square(stddev) * std::sqrt((-2.0f * std::log(w) ) / w) + mean;
}

void
Random::cosHemi(Real & x, Real & y, Real & z)
{
  Real const e1 = uniform01();
  Real const e2 = uniform01();

  // Jensen's method
  Real const sin_theta = std::sqrt(1.0f - e1);
  Real const cos_theta = std::sqrt(e1);
  Real const phi = 6.28318531f * e2;

  x = std::cos(phi) * sin_theta;
  y = std::sin(phi) * sin_theta;
  z = cos_theta;

  // We could also use Malley's method (pbrt p.657), since they are the same cost:
  //
  //  r = sqrt(e1);
  //  t = 2*pi*e2;
  //  x = cos(t)*r;
  //  y = sin(t)*r;
  //  z = sqrt(1.0 - x*x + y*y);
}

void
Random::cosPowHemi(Real const k, Real & x, Real & y, Real & z)
{
  Real const e1 = uniform01();
  Real const e2 = uniform01();

  Real const cos_theta = std::pow(e1, 1.0f / (k + 1.0f));
  Real const sin_theta = std::sqrt(1.0f - Math::square(cos_theta));
  Real const phi = 6.28318531f * e2;

  x = std::cos(phi) * sin_theta;
  y = std::sin(phi) * sin_theta;
  z = cos_theta;
}

void
Random::hemi(Real & x, Real & y, Real & z)
{
  sphere(x, y, z);
  z = std::abs(z);
}

void
Random::sphere(Real & x, Real & y, Real & z)
{
  // Adapted from: Robert E. Knop, "Algorithm 381: Random Vectors Uniform in Solid Angle", CACM, 13, p326, 1970.
  // It uses the fact that the projection of the uniform 2-sphere distribution onto any axis (in this case the Z axis) is
  // uniform.

  Real m2;  // squared magnitude
  do
  {
    x = 2 * uniform01() - 1;
    y = 2 * uniform01() - 1;
    m2 = x * x + y * y;

  } while (m2 > 1 || m2 < 1e-10);

  // The squared length happens to be uniformly distributed in [0, 1]. Since the projection of the 2-sphere distribution onto
  // the Z axis is uniform (can be derived from the observation that the area of a spherical cap is linear in its height), we
  // can use the squared length, rescaled to [-1, 1], as the Z value (latitude).
  z = 2 * m2 - 1;

  // x and y now locate the longitude of the point after scaling to lie on the sphere
  Real s = 2 * std::sqrt(1 - m2);  // this factor ensures x^2 + y^2 + z^2 = 1
  x *= s;
  y *= s;
}

void
Random::integers(int32 lo, int32 hi, int32 m, int32 * selected)
{
  // The current algorithm does no extra work to sort
  sortedIntegers(lo, hi, m, selected);
}

namespace RandomInternal {

// Get a very large random integer
int32
bigRand(Random & r)
{
  static int32 const BIG_HALF = Random::MAX_INTEGER / 2;
  return BIG_HALF + r.integer() % BIG_HALF;  // guaranteed to not exceed MAX_INTEGER
}

} // namespace RandomInternal

void
Random::sortedIntegers(int32 lo, int32 hi, int32 m, int32 * selected)
{
  int32 remaining = m;
  for (int32 i = lo; i <= hi; ++i)
  {
    // Select m of remaining n - i
    int32 r = RandomInternal::bigRand(*this);
    if ((r % (hi - i + 1)) < remaining)
    {
      selected[m - remaining] = i;
      remaining--;

      if (remaining <= 0)
        return;
    }
  }
}

namespace RandomInternal {

// http://burtleburtle.net/bob/hash/doobs.html
unsigned long
mix(unsigned long a, unsigned long b, unsigned long c)
{
  a=a-b;  a=a-c;  a=a^(c >> 13);
  b=b-c;  b=b-a;  b=b^(a << 8);
  c=c-a;  c=c-b;  c=c^(b >> 13);
  a=a-b;  a=a-c;  a=a^(c >> 12);
  b=b-c;  b=b-a;  b=b^(a << 16);
  c=c-a;  c=c-b;  c=c^(b >> 5);
  a=a-b;  a=a-c;  a=a^(c >> 3);
  b=b-c;  b=b-a;  b=b^(a << 10);
  c=c-a;  c=c-b;  c=c^(b >> 15);
  return c;
}

} // namespace RandomInternal

uint32
Random::getRandomSeed()
{
  // http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
#ifdef THEA_WINDOWS
  return (uint32)RandomInternal::mix((unsigned long)std::clock(), (unsigned long)std::time(NULL),
                                     (unsigned long)GetCurrentProcessId());
#else
  return (uint32)RandomInternal::mix((unsigned long)std::clock(), (unsigned long)std::time(NULL), (unsigned long)getpid());
#endif
}

} // namespace Thea
