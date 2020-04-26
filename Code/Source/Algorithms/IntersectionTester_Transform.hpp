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

#ifndef __Thea_Algorithms_IntersectionTester_Transform_hpp__
#define __Thea_Algorithms_IntersectionTester_Transform_hpp__

#include "TransformedObject.hpp"
#include "Transformer.hpp"

namespace Thea {
namespace Algorithms {

// Default specializations
template <typename ObjectA, typename TransformA, typename ObjectB, typename TransformB, int N, typename T>
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

template <typename ObjectA, typename TransformA, typename B, int N, typename T>
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

template <typename A, typename ObjectB, typename TransformB, int N, typename T>
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
