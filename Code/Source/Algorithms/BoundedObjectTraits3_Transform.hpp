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

#ifndef __Thea_Algorithms_BoundedObjectTraits3_Transform_hpp__
#define __Thea_Algorithms_BoundedObjectTraits3_Transform_hpp__

#include "TransformedObject.hpp"

namespace Thea {
namespace Algorithms {

template <typename ObjectT, typename TransformT>
class /* THEA_API */ BoundedObjectTraits3< TransformedObject<ObjectT, TransformT> >
{
  public:
    typedef TransformedObject<ObjectT, TransformT> TO;

    template <typename RangeT> static void getBounds(TO const & t, RangeT & bounds)
    {
      BoundedObjectTraits3<ObjectT>::getBounds(t.getObject(), bounds);
      bounds = bounds.transformAndBound(t.getTransform());
    }

    static Vector3 getCenter(TO const & t)
    {
      AxisAlignedBox3 bounds;
      getBounds(t, bounds);
      return bounds.getCenter();
    }

    static Real getHigh(TO const & t, int coord)
    {
      AxisAlignedBox3 bounds;
      getBounds(t, bounds);
      return bounds.getHigh()[coord];
    }

    static Real getLow(TO const & t, int coord)
    {
      AxisAlignedBox3 bounds;
      getBounds(t, bounds);
      return bounds.getLow()[coord];
    }

}; // class BoundedObjectTraits3

} // namespace Algorithms
} // namespace Thea

#endif
