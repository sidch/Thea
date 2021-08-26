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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_DiscreteExponentialMap_hpp__
#define __Thea_Algorithms_DiscreteExponentialMap_hpp__

#include "../Common.hpp"
#include "../CoordinateFrame3.hpp"
#include "../MatVec.hpp"
#include "../Noncopyable.hpp"
#include "../UnorderedMap.hpp"

namespace Thea {
namespace Algorithms {

namespace DiscreteExponentialMapInternal { class Impl; }

// Forward declarations
class PointCloud3;

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
        bool blendUpwind() const { return blend_upwind; }

        /** Select if parameters should be normalized to [-1, 1] or not. */
        Options & setNormalize(bool normalize_) { do_normalize = normalize_; return *this; }

        /** Check if parameters should be normalized to [-1, 1] or not. */
        bool normalize() const { return do_normalize; }

        /** Construct with default values. */
        Options() : blend_upwind(true), do_normalize(true) {}

        /** Get a set of options with default values. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        bool blend_upwind;  ///< Compute predecessor's parameters as a weighted average of nearby points.
        bool do_normalize;  ///< Normalize parameters to [-1, 1].

    }; // class Options

    typedef UnorderedMap<intx, Vector2> ParameterMap;  ///< Map from sample indices to parameters.

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
    void parametrize(PointCloud3 const & surf, intx origin_index, Vector3 const & u_axis, Vector3 const & v_axis, Real radius);

    /**
     * Get the exponential map parameters of a particular sample.
     *
     * @param sample_index Index of the query sample.
     * @param has_parameters Set to false if the sample was not parametrized (outside radius), else true.
     *
     * @see param
     etrize()
     */
    Vector2 getParameters(intx sample_index, bool & has_parameters) const;

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
