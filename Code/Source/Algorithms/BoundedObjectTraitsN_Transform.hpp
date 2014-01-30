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

#ifndef __Thea_Algorithms_BoundedObjectTraitsN_Transform_hpp__
#define __Thea_Algorithms_BoundedObjectTraitsN_Transform_hpp__

#include "TransformedObject.hpp"

namespace Thea {
namespace Algorithms {

template <typename ObjectT, typename TransformT, long N, typename ScalarT>
class /* THEA_API */ BoundedObjectTraitsN< TransformedObject<ObjectT, TransformT>, N, ScalarT >
{
  public:
    typedef TransformedObject<ObjectT, TransformT> TO;

    template <typename RangeT> static void getBounds(TO const & t, RangeT & bounds)
    {
      BoundedObjectTraitsN<ObjectT, N, ScalarT>::getBounds(t.getObject(), bounds);
      bounds = bounds.transformAndBound(t.getTransform());
    }

    static VectorN<N, ScalarT> getCenter(TO const & t)
    {
      AxisAlignedBoxN<N, ScalarT> bounds;
      getBounds(t, bounds);
      return bounds.getCenter();
    }

    static ScalarT getHigh(TO const & t, long coord)
    {
      AxisAlignedBoxN<N, ScalarT> bounds;
      getBounds(t, bounds);
      return bounds.getHigh()[coord];
    }

    static ScalarT getLow(TO const & t, long coord)
    {
      AxisAlignedBoxN<N, ScalarT> bounds;
      getBounds(t, bounds);
      return bounds.getLow()[coord];
    }

}; // class BoundedObjectTraitsN

} // namespace Algorithms
} // namespace Thea

#endif
