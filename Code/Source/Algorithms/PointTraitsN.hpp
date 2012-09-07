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

#ifndef __Thea_Algorithms_PointTraitsN_hpp__
#define __Thea_Algorithms_PointTraitsN_hpp__

#include "../Common.hpp"
#include "../VectorN.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Has boolean member <code>value = true</code> if <code>T</code> can be identified with a single point in n-D space, else
 * <code>value = false</code>. Unless you specialize the class to set the value to true, it is false by default.
 *
 * @see PointTraitsN
 */
template <typename PointT, long N>
class /* THEA_API */ IsPointN
{
  public:
    static bool const value = false;
};

// Partial specialization of PointTraitsN for const types
template <typename PointT, long N>
class /* THEA_API */ IsPointN<PointT const, N>
{
  public:
    static bool const value = IsPointN<PointT, N>::value;
};

// Specialization for VectorN
template <long N, typename ScalarT>
class /* THEA_API */ IsPointN< VectorN<N, ScalarT>, N >
{
  public:
    static bool const value = true;
};

/**
 * Traits for an object which can be identified with a single point in N-space.
 *
 * @see IsPointN
 */
template <typename PointT, long N, typename ScalarT = Real>
class /* THEA_API */ PointTraitsN
{
  public:
    typedef VectorN<N, ScalarT> VectorT;  ///< A vector in N-space.

    /** Get the position of the "point". The default implementation is for Vector and types implicitly convertible to it. */
    static VectorT getPosition(PointT const & p) { return p; }

}; // class PointTraitsN

// Partial specialization of PointTraitsN for const types
template <typename PointT, long N, typename ScalarT>
class /* THEA_API */ PointTraitsN<PointT const, N, ScalarT>
{
  public:
    typedef VectorN<N, ScalarT> VectorT;

    static VectorT getPosition(PointT const & p) { return PointTraitsN<PointT, N, ScalarT>::getPosition(p); }

}; // class PointTraitsN<PointT const, N, ScalarT>

// Partial specialization of PointTraitsN for pointer types
template <typename PointT, long N, typename ScalarT>
class /* THEA_API */ PointTraitsN<PointT *, N, ScalarT>
{
  public:
    typedef VectorN<N, ScalarT> VectorT;

    static VectorT getPosition(PointT * p) { return PointTraitsN<PointT, N, ScalarT>::getPosition(*p); }

}; // class PointTraitsN<PointT *, N, ScalarT>

} // namespace Algorithms
} // namespace Thea

#endif
