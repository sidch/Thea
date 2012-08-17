//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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

namespace Thea {
namespace Algorithms {

// Default to hard-transforming each shape, assuming 3-space, in the absence of anything smarter
template <typename ObjectA, typename TransformA, typename ObjectB, typename TransformB>
struct /* THEA_API */ MetricL2Impl< TransformedObject<ObjectA, TransformA>, TransformedObject<ObjectB, TransformB> >
{
  typedef TransformedObject<ObjectA, TransformA> TA;
  typedef TransformedObject<ObjectB, TransformB> TB;

  static double distance(TA const & a, TB const & b)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return MetricL2::distance(Transformer::transform(a.getObject(), a.getTransform()),
                                  Transformer::transform(b.getObject(), b.getTransform()));
      else
        return MetricL2::distance(Transformer::transform(a.getObject(), a.getTransform()), b.getObject());
    }
    else
    {
      if (b.hasTransform())
        return MetricL2::distance(a.getObject(), Transformer::transform(b.getObject(), b.getTransform()));
      else
        return MetricL2::distance(a.getObject(), b.getObject());
    }
  }

  static double monotoneApproxDistance(TA const & a, TB const & b)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return MetricL2::monotoneApproxDistance(Transformer::transform(a.getObject(), a.getTransform()),
                                                Transformer::transform(b.getObject(), b.getTransform()));
      else
        return MetricL2::monotoneApproxDistance(Transformer::transform(a.getObject(), a.getTransform()), b.getObject());
    }
    else
    {
      if (b.hasTransform())
        return MetricL2::monotoneApproxDistance(a.getObject(), Transformer::transform(b.getObject(), b.getTransform()));
      else
        return MetricL2::monotoneApproxDistance(a.getObject(), b.getObject());
    }
  }

  static double closestPoints(TA const & a, TB const & b, Vector3 & cpa, Vector3 & cpb)
  {
    if (a.hasTransform())
    {
      if (b.hasTransform())
        return MetricL2::closestPoints(Transformer::transform(a.getObject(), a.getTransform()),
                                       Transformer::transform(b.getObject(), b.getTransform()), cpa, cpb);
      else
        return MetricL2::closestPoints(Transformer::transform(a.getObject(), a.getTransform()), b.getObject(), cpa, cpb);
    }
    else
    {
      if (b.hasTransform())
        return MetricL2::closestPoints(a.getObject(), Transformer::transform(b.getObject(), b.getTransform()), cpa, cpb);
      else
        return MetricL2::closestPoints(a.getObject(), b.getObject(), cpa, cpb);
    }
  }
};

template <typename ObjectA, typename TransformA, typename B>
struct /* THEA_API */ MetricL2Impl< TransformedObject<ObjectA, TransformA>, B >
{
  typedef TransformedObject<ObjectA, TransformA> TA;

  static double distance(TA const & a, B const & b)
  {
    if (a.hasTransform())
      return MetricL2::distance(Transformer::transform(a.getObject(), a.getTransform()), b);
    else
      return MetricL2::distance(a.getObject(), b);
  }

  static double monotoneApproxDistance(TA const & a, B const & b)
  {
    if (a.hasTransform())
      return MetricL2::monotoneApproxDistance(Transformer::transform(a.getObject(), a.getTransform()), b);
    else
      return MetricL2::monotoneApproxDistance(a.getObject(), b);
  }

  static double closestPoints(TA const & a, B const & b, Vector3 & cpa, Vector3 & cpb)
  {
    if (a.hasTransform())
      return MetricL2::closestPoints(Transformer::transform(a.getObject(), a.getTransform()), b, cpa, cpb);
    else
      return MetricL2::closestPoints(a.getObject(), b, cpa, cpb);
  }
};

template <typename ObjectA, typename TransformA, typename B>
struct /* THEA_API */ MetricL2Impl< TransformedObject<ObjectA, TransformA>, B * >
{
  typedef TransformedObject<ObjectA, TransformA> TA;

  static double distance(TA const & a, B const * b)
  {
    if (a.hasTransform())
      return MetricL2::distance(Transformer::transform(a.getObject(), a.getTransform()), *b);
    else
      return MetricL2::distance(a.getObject(), *b);
  }

  static double monotoneApproxDistance(TA const & a, B const * b)
  {
    if (a.hasTransform())
      return MetricL2::monotoneApproxDistance(Transformer::transform(a.getObject(), a.getTransform()), *b);
    else
      return MetricL2::monotoneApproxDistance(a.getObject(), *b);
  }

  static double closestPoints(TA const & a, B const * b, Vector3 & cpa, Vector3 & cpb)
  {
    if (a.hasTransform())
      return MetricL2::closestPoints(Transformer::transform(a.getObject(), a.getTransform()), *b, cpa, cpb);
    else
      return MetricL2::closestPoints(a.getObject(), *b, cpa, cpb);
  }
};

template <typename A, typename ObjectB, typename TransformB>
struct /* THEA_API */ MetricL2Impl< A, TransformedObject<ObjectB, TransformB> >
{
  typedef TransformedObject<ObjectB, TransformB> TB;

  static double distance(A const & a, TB const & b)
  {
    if (b.hasTransform())
      return MetricL2::distance(a, Transformer::transform(b.getObject(), b.getTransform()));
    else
      return MetricL2::distance(a, b.getObject());
  }

  static double monotoneApproxDistance(A const & a, TB const & b)
  {
    if (b.hasTransform())
      return MetricL2::monotoneApproxDistance(a, Transformer::transform(b.getObject(), b.getTransform()));
    else
      return MetricL2::monotoneApproxDistance(a, b.getObject());
  }

  static double closestPoints(A const & a, TB const & b, Vector3 & cpa, Vector3 & cpb)
  {
    if (b.hasTransform())
      return MetricL2::closestPoints(a, Transformer::transform(b.getObject(), b.getTransform()), cpa, cpb);
    else
      return MetricL2::closestPoints(a, b.getObject(), cpa, cpb);
  }
};

template <typename A, typename ObjectB, typename TransformB>
struct /* THEA_API */ MetricL2Impl< A *, TransformedObject<ObjectB, TransformB> >
{
  typedef TransformedObject<ObjectB, TransformB> TB;

  static double distance(A const * a, TB const & b)
  {
    if (b.hasTransform())
      return MetricL2::distance(*a, Transformer::transform(b.getObject(), b.getTransform()));
    else
      return MetricL2::distance(*a, b.getObject());
  }

  static double monotoneApproxDistance(A const * a, TB const & b)
  {
    if (b.hasTransform())
      return MetricL2::monotoneApproxDistance(*a, Transformer::transform(b.getObject(), b.getTransform()));
    else
      return MetricL2::monotoneApproxDistance(*a, b.getObject());
  }

  static double closestPoints(A const * a, TB const & b, Vector3 & cpa, Vector3 & cpb)
  {
    if (b.hasTransform())
      return MetricL2::closestPoints(*a, Transformer::transform(b.getObject(), b.getTransform()), cpa, cpb);
    else
      return MetricL2::closestPoints(*a, b.getObject(), cpa, cpb);
  }
};

} // namespace Algorithms
} // namespace Thea

#endif
