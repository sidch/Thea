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
// First version: 2013
//
//============================================================================

#ifndef __Thea_Algorithms_ImageFeatures_ShapeContext_hpp__
#define __Thea_Algorithms_ImageFeatures_ShapeContext_hpp__

#include "../../Common.hpp"
#include "../../Array.hpp"

namespace Thea {

// Forward declarations
class Image;

namespace Algorithms {

/** Namespace for classes that compute features on an image. */
namespace ImageFeatures {

namespace ShapeContextInternal {

struct QuadTree;

} // namespace ShapeContextInternal

/**
 * Computes the shape context at pixels of a contour image. If the image has more than one channel, the luminance channel is
 * used. Follows the algorithm of:
 *
 * Mori, Belongie and Malik, "Shape contexts enable efficient retrieval of similar shapes", CVPR 2001.
 */
class THEA_API ShapeContext
{
  public:
    /** Constructor. */
    ShapeContext(Image const & image);

    /** Destructor. */
    ~ShapeContext();

    /**
     * Computes the shape context at each pixel of a contour image.
     *
     * @param num_radial_bins Number of divisions in the radial direction.
     * @param num_polar_bins Number of divisions in the angular direction.
     * @param values Computed shape contexts, in a num_pixels * \a num_radial_bins * \a num_polar_bins array (num_pixels is
     *   major dimension for packing).
     * @param ignore_empty_pixels Does not compute shape contexts for pixels that are empty (zero luminance, i.e. black). The
     *   corresponding entries in \a values are set to zero.
     * @param max_radius Limits the area of the context for a pixel to this radius (negative for default).
     */
    void compute(intx num_radial_bins, intx num_polar_bins, Array<Real> & values, bool ignore_empty_pixels = true,
                 Real max_radius = -1) const;

    /**
     * Computes the shape context at a single pixel of a contour image.
     *
     * @param row Row of pixel.
     * @param col Column of pixel.
     * @param num_radial_bins Number of divisions in the radial direction.
     * @param num_polar_bins Number of divisions in the angular direction.
     * @param values Computed shape context, in a \a num_radial_bins * \a num_polar_bins array (\a num_radial_bins is major
     *   dimension for packing).
     * @param max_radius Limits the area of the context to this radius (negative for default).
     */
    void compute(int row, int col, intx num_radial_bins, intx num_polar_bins, Array<Real> & values, Real max_radius = -1)
         const;

  private:
    ShapeContextInternal::QuadTree * qtree;

}; // class ShapeContext

} // namespace ImageFeatures
} // namespace Algorithms

} // namespace Thea

#endif
