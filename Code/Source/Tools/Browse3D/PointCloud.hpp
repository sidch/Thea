//============================================================================
//
// This file is part of the Browse3D project.
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

#ifndef __Browse3D_PointCloud_hpp__
#define __Browse3D_PointCloud_hpp__

#include "Common.hpp"
#include "GraphicsWidget.hpp"
#include "../../AxisAlignedBox3.hpp"

namespace Thea {
namespace Graphics {

class VARArea;
class VAR;

} // namespace Graphics
} // namespace Thea

namespace Browse3D {

/** The model manipulated by the user. */
class PointCloud : public virtual NamedObject, public GraphicsWidget
{
  private:
    /** A point on the surface. */
    struct Point
    {
      Point() {}
      Point(Vector3 const & p_, Vector3 const & n_ = Vector3::zero()) : p(p_), n(n_) {}

      Vector3 p;
      Vector3 n;
    };

  public:
    THEA_DEF_POINTER_TYPES(PointCloud, shared_ptr, weak_ptr)

    /** Constructor. */
    PointCloud(std::string const & path = "", std::string const & features_path = "");

    /** Destructor. */
    ~PointCloud();

    /** Is the model empty? */
    bool isEmpty() const { return points.empty(); }

    /** Clear the point cloud. */
    void clear();

    /**
     * Load the point cloud from a disk file.
     *
     * @return True if the model was successfully loaded, else false.
     */
    bool load(std::string const & path, std::string const & features_path = "");

    /** Load features from a file. */
    bool loadFeatures(std::string const & filename_);

    /** Explicitly set the color of each point. Overrides any colors derived from features. */
    bool setPointColors(TheaArray<ColorRGBA> const & colors_);

    AxisAlignedBox3 const & getBounds() const;

    void updateBounds();

    void uploadToGraphicsSystem(Graphics::RenderSystem & render_system);

    void draw(Graphics::RenderSystem & render_system,
              Graphics::RenderOptions const & options = Graphics::RenderOptions::defaults()) const;

  private:
    /** Invalidate the bounding box of the point cloud. */
    void invalidateBounds();

    /**
     * Get the path to the file in which features are stored, given the path to the point cloud and optionally a file/directory
     * containing the features.
     */
    std::string getDefaultFeaturesFilename(std::string const & filename, std::string const & features_path) const;

    /** Reconstruct an approximate surface from the point cloud. */
    void reconstructSurface(Graphics::RenderSystem & render_system,
                            Graphics::RenderOptions const & options = Graphics::RenderOptions::defaults()) const;

    /** Get the color of a point. */
    ColorRGBA getColor(size_t point_index) const;

    TheaArray<Point> points;

    bool has_normals;
    bool normals_are_normalized;
    AxisAlignedBox3 bounds;

    TheaArray<ColorRGBA> colors;
    TheaArray< TheaArray<Real> > features;

    bool has_graph;
    TheaArray< TheaArray<long> > graph;

    int changed_buffers;           ///< A bitwise OR of the flags of the buffers that have changed.
    Graphics::VARArea * var_area;  ///< GPU buffer area.
    Graphics::VAR * vertices_var;  ///< GPU buffer for vertex positions.
    Graphics::VAR * colors_var;    ///< GPU buffer for vertex colors.

}; // class PointCloud

} // namespace Browse3D

#endif
