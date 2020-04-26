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

#ifndef __Thea_Algorithms_BoundedTraitsN_Transform_hpp__
#define __Thea_Algorithms_BoundedTraitsN_Transform_hpp__

#include "TransformedObject.hpp"

namespace Thea {
namespace Algorithms {

template <typename ObjectT, typename TransformT, int N>
class IsBoundedN< TransformedObject<ObjectT, TransformT>, N >
{
  public:
    static bool const value = IsBoundedN<ObjectT, N>::value;
};

template <typename ObjectT, typename TransformT, int N, typename ScalarT>
class /* THEA_API */ BoundedTraitsN< TransformedObject<ObjectT, TransformT>, N, ScalarT >
{
  public:
    typedef TransformedObject<ObjectT, TransformT> TO;

    template <typename RangeT> static void getBounds(TO const & t, RangeT & bounds)
    {
      BoundedTraitsN<ObjectT, N, ScalarT>::getBounds(t.getObject(), bounds);
      bounds = bounds.transformAndBound(t.getTransform());
    }

    static Vector<N, ScalarT> getCenter(TO const & t)
    {
      AxisAlignedBoxN<N, ScalarT> bounds;
      getBounds(t, bounds);
      return bounds.getCenter();
    }

    static ScalarT getHigh(TO const & t, intx coord)
    {
      AxisAlignedBoxN<N, ScalarT> bounds;
      getBounds(t, bounds);
      return bounds.getHigh()[coord];
    }

    static ScalarT getLow(TO const & t, intx coord)
    {
      AxisAlignedBoxN<N, ScalarT> bounds;
      getBounds(t, bounds);
      return bounds.getLow()[coord];
    }

}; // class BoundedTraitsN

} // namespace Algorithms
} // namespace Thea

#endif
