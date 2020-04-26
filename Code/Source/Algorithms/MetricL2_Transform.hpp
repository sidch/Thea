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

#ifndef __Thea_Algorithms_MetricL2_Transform_hpp__
#define __Thea_Algorithms_MetricL2_Transform_hpp__

#include "TransformedObject.hpp"
#include "Transformer.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

// Default to hard-transforming each shape, in the absence of anything smarter
template <typename ObjectA, typename TransformA, typename ObjectB, typename TransformB, int N, typename T>
struct MetricL2Impl< TransformedObject<ObjectA, TransformA>, TransformedObject<ObjectB, TransformB>, N, T >
{
  typedef TransformedObject<ObjectA, TransformA> TA;
  typedef TransformedObject<ObjectB, TransformB> TB;

  static T distance(TA const & a, TB const & b)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return MetricL2::distance<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()),
                                                 Transformer::transform<N, T>(b.getObject(), b.getTransform()));
      else
        return MetricL2::distance<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b.getObject());
    }
    else
    {
      if (b.hasTransform())
        return MetricL2::distance<N, T>(a.getObject(), Transformer::transform<N, T>(b.getObject(), b.getTransform()));
      else
        return MetricL2::distance<N, T>(a.getObject(), b.getObject());
    }
  }

  static T monotoneApproxDistance(TA const & a, TB const & b)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return MetricL2::monotoneApproxDistance<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()),
                                                      Transformer::transform<N, T>(b.getObject(), b.getTransform()));
      else
        return MetricL2::monotoneApproxDistance<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()),
                                                      b.getObject());
    }
    else
    {
      if (b.hasTransform())
        return MetricL2::monotoneApproxDistance<N, T>(a.getObject(),
                                                      Transformer::transform<N, T>(b.getObject(), b.getTransform()));
      else
        return MetricL2::monotoneApproxDistance<N, T>(a.getObject(), b.getObject());
    }
  }

  static T closestPoints(TA const & a, TB const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return MetricL2::closestPoints<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()),
                                             Transformer::transform<N, T>(b.getObject(), b.getTransform()), cpa, cpb);
      else
        return MetricL2::closestPoints<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b.getObject(),
                                             cpa, cpb);
    }
    else
    {
      if (b.hasTransform())
        return MetricL2::closestPoints<N, T>(a.getObject(), Transformer::transform<N, T>(b.getObject(), b.getTransform()),
                                             cpa, cpb);
      else
        return MetricL2::closestPoints<N, T>(a.getObject(), b.getObject(), cpa, cpb);
    }
  }
};

template <typename ObjectA, typename TransformA, typename B, int N, typename T>
struct MetricL2Impl< TransformedObject<ObjectA, TransformA>, B, N, T,
                     typename std::enable_if< !std::is_pointer<B>::value
                                           && !MetricL2Internal::TransformedObjectCheck<B>::value >::type>
{
  typedef TransformedObject<ObjectA, TransformA> TA;

  static T distance(TA const & a, B const & b)
  {
    if (a.hasTransform())
      return MetricL2::distance<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b);
    else
      return MetricL2::distance<N, T>(a.getObject(), b);
  }

  static T monotoneApproxDistance(TA const & a, B const & b)
  {
    if (a.hasTransform())
      return MetricL2::monotoneApproxDistance<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b);
    else
      return MetricL2::monotoneApproxDistance<N, T>(a.getObject(), b);
  }

  static T closestPoints(TA const & a, B const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  {
    if (a.hasTransform())
      return MetricL2::closestPoints<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b, cpa, cpb);
    else
      return MetricL2::closestPoints<N, T>(a.getObject(), b, cpa, cpb);
  }
};

template <typename A, typename ObjectB, typename TransformB, int N, typename T>
struct MetricL2Impl< A, TransformedObject<ObjectB, TransformB>, N, T,
                     typename std::enable_if< !std::is_pointer<A>::value
                                           && !MetricL2Internal::TransformedObjectCheck<A>::value >::type>
{
  typedef TransformedObject<ObjectB, TransformB> TB;

  static T distance(A const & a, TB const & b)
  {
    if (b.hasTransform())
      return MetricL2::distance<N, T>(a, Transformer::transform<N, T>(b.getObject(), b.getTransform()));
    else
      return MetricL2::distance<N, T>(a, b.getObject());
  }

  static T monotoneApproxDistance(A const & a, TB const & b)
  {
    if (b.hasTransform())
      return MetricL2::monotoneApproxDistance<N, T>(a, Transformer::transform<N, T>(b.getObject(), b.getTransform()));
    else
      return MetricL2::monotoneApproxDistance<N, T>(a, b.getObject());
  }

  static T closestPoints(A const & a, TB const & b, Vector<N, T> & cpa, Vector<N, T> & cpb)
  {
    if (b.hasTransform())
      return MetricL2::closestPoints<N, T>(a, Transformer::transform<N, T>(b.getObject(), b.getTransform()), cpa, cpb);
    else
      return MetricL2::closestPoints<N, T>(a, b.getObject(), cpa, cpb);
  }
};

} // namespace Algorithms
} // namespace Thea

#endif
