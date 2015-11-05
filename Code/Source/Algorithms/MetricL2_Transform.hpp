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

#ifndef __Thea_Algorithms_MetricL2_Transform_hpp__
#define __Thea_Algorithms_MetricL2_Transform_hpp__

#include "TransformedObject.hpp"
#include "Transformer.hpp"
#include <boost/type_traits/is_pointer.hpp>
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

// Default to hard-transforming each shape, in the absence of anything smarter
template <typename ObjectA, typename TransformA, typename ObjectB, typename TransformB, long N, typename T>
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

  static T closestPoints(TA const & a, TB const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
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

template <typename ObjectA, typename TransformA, typename B, long N, typename T>
struct MetricL2Impl< TransformedObject<ObjectA, TransformA>, B, N, T,
                     typename boost::enable_if_c< !boost::is_pointer<B>::value
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

  static T closestPoints(TA const & a, B const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
  {
    if (a.hasTransform())
      return MetricL2::closestPoints<N, T>(Transformer::transform<N, T>(a.getObject(), a.getTransform()), b, cpa, cpb);
    else
      return MetricL2::closestPoints<N, T>(a.getObject(), b, cpa, cpb);
  }
};

template <typename A, typename ObjectB, typename TransformB, long N, typename T>
struct MetricL2Impl< A, TransformedObject<ObjectB, TransformB>, N, T,
                     typename boost::enable_if_c< !boost::is_pointer<A>::value
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

  static T closestPoints(A const & a, TB const & b, VectorN<N, T> & cpa, VectorN<N, T> & cpb)
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
