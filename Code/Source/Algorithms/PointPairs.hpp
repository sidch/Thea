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

#ifndef __Thea_Algorithms_PointPairs_hpp__
#define __Thea_Algorithms_PointPairs_hpp__

#include "../Common.hpp"
#include "../VectorN.hpp"
#include <utility>

namespace Thea {
namespace Algorithms {

/** A pair of points in 1-space (single-precision) */
typedef std::pair<float, float> FloatPair;

/** A pair of points in 1-space (double-precision) */
typedef std::pair<double, double> DoublePair;

/** A pair of points in 1-space (default precision) */
typedef std::pair<Real, Real> RealPair;

/** A pair of points in 2-space. */
typedef std::pair<Vector2, Vector2> PointPair2;

/** A pair of points in 3-space. */
typedef std::pair<Vector3, Vector3> PointPair3;

/** A pair of points in 4-space. */
typedef std::pair<Vector4, Vector4> PointPair4;

/** A pair of points in n-dimensional space. */
template <size_t N, typename T>
class /* THEA_API */ PointPairN : public std::pair< VectorN<N, T>, VectorN<N, T> >
{
  private:
    typedef std::pair< VectorN<N, T>, VectorN<N, T> > BaseType;

  public:
    /** Default constructor. */
    PointPairN() : BaseType() {}

    /** Initializing constructor. */
    PointPairN(VectorN<N, T> const & a, VectorN<N, T> const & b) : BaseType(a, b) {}

}; // class PointPairN

} // namespace Algorithms
} // namespace Thea

#endif
