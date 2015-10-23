//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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
    void compute(long num_radial_bins, long num_polar_bins, TheaArray<Real> & values, bool ignore_empty_pixels = true,
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
    void compute(int row, int col, long num_radial_bins, long num_polar_bins, TheaArray<Real> & values, Real max_radius = -1)
         const;

  private:
    ShapeContextInternal::QuadTree * qtree;

}; // class ShapeContext

} // namespace ImageFeatures
} // namespace Algorithms

} // namespace Thea

#endif
