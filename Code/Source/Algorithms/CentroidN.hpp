//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2012, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_CentroidN_hpp__
#define __Thea_Algorithms_CentroidN_hpp__

#include "../Common.hpp"
#include "../Math.hpp"
#include "../VectorN.hpp"
#include "IteratorModifiers.hpp"
#include "PointTraitsN.hpp"

namespace Thea {
namespace Algorithms {

/** Finding the centroid of N-dimensional data. */
template <typename T, long N, typename ScalarT = Real, typename Enable = void>
class /* THEA_API */ CentroidN
{
  public:
    typedef VectorN<N, ScalarT> VectorT;  ///< N-D vector used to represent points.

    /**
     * Unweighted centroid of a set of N-D objects. InputIterator must dereference to type T.
     *
     * @param begin The first object in the set.
     * @param end One position beyond the last object in the set.
     *
     * @return The centroid of the objects (or the origin if no objects were supplied).
     */
    template <typename InputIterator> static VectorT compute(InputIterator begin, InputIterator end);

    /**
     * Weighted centroid of a set of N-D objects. ObjectInputIterator must dereference to type T and WeightInputIterator to a
     * scalar type. The objects and weights sequences must correspond.
     *
     * @param objects_begin The first object in the set.
     * @param objects_end One position beyond the last object in the set.
     * @param weights_begin The first entry in the weights list.
     *
     * @return The weighted centroid of the objects (or the origin if the weights sum to zero).
     */
    template <typename ObjectInputIterator, typename WeightInputIterator>
    static VectorT compute(ObjectInputIterator objects_begin, ObjectInputIterator objects_end,
                           WeightInputIterator weights_begin);

}; // class CentroidN

// Centroid of objects passed as pointers.
template <typename T, long N, typename ScalarT>
class /* THEA_API */ CentroidN<T *, N, ScalarT>
{
  public:
    typedef VectorN<N, ScalarT> VectorT;

    template <typename InputIterator> static VectorT compute(InputIterator begin, InputIterator end)
    {
      return CentroidN<T, N, ScalarT>::compute(PtrToRefIterator<T, InputIterator>(begin),
                                               PtrToRefIterator<T, InputIterator>(end));
    }

    template <typename ObjectInputIterator, typename WeightInputIterator>
    static VectorT compute(ObjectInputIterator objects_begin, ObjectInputIterator objects_end,
                           WeightInputIterator weights_begin)
    {
      return CentroidN<T, N, ScalarT>::compute(PtrToRefIterator<T, ObjectInputIterator>(objects_begin),
                                               PtrToRefIterator<T, ObjectInputIterator>(objects_end), weights_begin);
    }

}; // class CentroidN<T *>

// Centroid of objects that map to single points in N-space.
template <typename T, long N, typename ScalarT>
class /* THEA_API */ CentroidN<T, N, ScalarT, typename boost::enable_if< IsNonReferencedPointN<T, N> >::type>
{
  public:
    typedef VectorN<N, ScalarT> VectorT;

    template <typename InputIterator> static VectorT compute(InputIterator begin, InputIterator end)
    {
      VectorT sum_points = VectorT::zero();
      long num_points = 0;
      for (InputIterator i = begin; i != end; ++i, ++num_points)
        sum_points += PointTraitsN<T, N, ScalarT>::getPosition(*i);

      return num_points > 0 ? sum_points / num_points : VectorT::zero();
    }

    template <typename ObjectInputIterator, typename WeightInputIterator>
    static VectorT compute(ObjectInputIterator objects_begin, ObjectInputIterator objects_end,
                           WeightInputIterator weights_begin)
    {
      VectorT sum_points = VectorT::zero();
      double sum_weights = 0;
      WeightInputIterator wi = weights_begin;
      for (ObjectInputIterator pi = objects_begin; pi != objects_end; ++pi, ++wi)
      {
        sum_points += PointTraitsN<T, N, ScalarT>::getPosition(*pi);
        sum_weights += static_cast<double>(*wi);
      }

      return Math::fuzzyEq(sum_weights, 0.0) ? VectorT::zero() : sum_points / sum_weights;
    }

}; // class CentroidN<PointN>

} // namespace Algorithms
} // namespace Thea

#endif
