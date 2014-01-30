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

#ifndef __Thea_Algorithms_IntersectionTester_Transform_hpp__
#define __Thea_Algorithms_IntersectionTester_Transform_hpp__

#include "TransformedObject.hpp"
#include "Transformer.hpp"

namespace Thea {
namespace Algorithms {

// Default specializations
template <typename ObjectA, typename TransformA, typename ObjectB, typename TransformB, long N, typename T>
struct IntersectionTesterImpl< TransformedObject<ObjectA, TransformA>, TransformedObject<ObjectB, TransformB>, N, T >
{
  typedef TransformedObject<ObjectA, TransformA> TA;
  typedef TransformedObject<ObjectB, TransformB> TB;

  static bool intersects(TA const & a, TB const & b)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return IntersectionTester::intersects<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()),
                                                    Transformer::transform<N, T>(b.getObject(), b.getTransform()));
      else
        return IntersectionTester::intersects<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()),
                                                    b.getObject());
    }
    else
    {
      if (b.hasTransform())
        return IntersectionTester::intersects<N, T>(a.getObject(), Transformer::transform<N, T>(b.getObject(),
                                                    b.getTransform()));
      else
        return IntersectionTester::intersects<N, T>(a.getObject(), b.getObject());
    }
  }
};

template <typename ObjectA, typename TransformA, typename B, long N, typename T>
struct IntersectionTesterImpl< TransformedObject<ObjectA, TransformA>, B, N, T >
{
  typedef TransformedObject<ObjectA, TransformA> TA;

  static bool intersects(TA const & a, B const & b)
  {
    if (a.hasTransform())
      return IntersectionTester::intersects<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b);
    else
      return IntersectionTester::intersects<N, T>(a.getObject(), b);
  }
};

template <typename A, typename ObjectB, typename TransformB, long N, typename T>
struct IntersectionTesterImpl< A, TransformedObject<ObjectB, TransformB>, N, T >
{
  typedef TransformedObject<ObjectB, TransformB> TB;

  static bool intersects(A const & a, TB const & b)
  {
    if (b.hasTransform())
      return IntersectionTester::intersects<N, T>(a, Transformer::transform<N, T>(b.getObject(), b.getTransform()));
    else
      return IntersectionTester::intersects<N, T>(a, b.getObject());
  }
};

} // namespace Algorithms
} // namespace Thea

#endif
