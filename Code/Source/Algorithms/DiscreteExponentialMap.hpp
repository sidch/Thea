//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Cornell University
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

#ifndef __Thea_Algorithms_DiscreteExponentialMap_hpp__
#define __Thea_Algorithms_DiscreteExponentialMap_hpp__

#include "../Common.hpp"
#include "../CoordinateFrame3.hpp"
#include "../Noncopyable.hpp"
#include "../UnorderedMap.hpp"
#include "../Vector3.hpp"

namespace Thea {
namespace Algorithms {

namespace DiscreteExponentialMapInternal { class Impl; }

// Forward declarations
class SampleGraph;

/**
 * Compute the discrete exponential map of a surface around a point. Based (roughly) on:
 *
 * Ryan Schmidt, Cindy Grimm and Brian Wyvill, "Interactive Decal Compositing with Discrete Exponential Maps", ACM Transactions
 * on %Graphics, 25(3), July 2006, pp. 605-613.
 */
class DiscreteExponentialMap : private Noncopyable
{
  public:
    /** %Options for computing the discrete exponential map of a surface around a point. */
    class Options
    {
      public:
        /**
         * Select if the predecessor's parameters should be computed as a weighted average of nearby points or not. If selected,
         * this can slightly slow down the parametrization.
         */
        Options & setBlendUpwind(bool blend_upwind_) { blend_upwind = blend_upwind_; return *this; }

        /** Check if the predecessor's parameters should be computed as a weighted average of nearby points or not. */
        bool getBlendUpwind() const { return blend_upwind; }

        /** Select if parameters should be normalized to [-1, 1] or not. */
        Options & setNormalize(bool normalize_) { normalize = normalize_; return *this; }

        /** Check if parameters should be normalized to [-1, 1] or not. */
        bool getNormalize() const { return normalize; }

        /** Construct with default values. */
        Options() : blend_upwind(true), normalize(true) {}

        /** Get a set of options with default values. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        bool blend_upwind;  ///< Compute predecessor's parameters as a weighted average of nearby points.
        bool normalize;  ///< Normalize parameters to [-1, 1].

    }; // class Options

    typedef TheaUnorderedMap<long, Vector2> ParameterMap;  ///< Map from sample indices to parameters.

    /** Constructor. */
    DiscreteExponentialMap(Options const & options = Options::defaults());

    /** Destructor. */
    ~DiscreteExponentialMap();

    /**
     * Compute the exponential map parametrization of a surface, centered at the sample with index \a origin_index, with a basis
     * (\a u_axis, \a v_axis) in the tangent plane, and up to a geodesic radius of \a radius. The surface is represented as an
     * adjacency graph of surface samples. The sample coordinates can be recovered using getParameters().
     *
     * @see getParameters(), getParameterMap()
     */
    void parametrize(SampleGraph const & sample_graph, long origin_index, Vector3 const & u_axis, Vector3 const & v_axis,
                     Real radius);

    /**
     * Get the exponential map parameters of a particular sample.
     *
     * @param sample_index Index of the query sample.
     * @param has_parameters Set to false if the sample was not parametrized (outside radius), else true.
     *
     * @see parametrize()
     */
    Vector2 getParameters(long sample_index, bool & has_parameters) const;

    /**
     * Get the full map from sample indices to parameters. Samples that were not parametrized (outside radius) are not included.
     *
     * @see parametrize()
     */
    ParameterMap const & getParameterMap() const;

    /** Get the coordinate frame centered at the parametrization origin and having X and Y axes in the tangent plane. */
    CoordinateFrame3 getTangentFrame() const;

    /** Get the bounding radius of the region which is parametrized. */
    Real getRadius() const;

    /** Clear any existing parametrization. */
    void clear();

  private:
    DiscreteExponentialMapInternal::Impl * impl;

}; // class DiscreteExponentialMap

} // namespace Algorithms
} // namespace Thea

#endif
