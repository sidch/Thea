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
// First version: 2012
//
//============================================================================

#ifndef __Browse3D_PointCloud_hpp__
#define __Browse3D_PointCloud_hpp__

#include "Common.hpp"
#include "GraphicsWidget.hpp"
#include "../../AxisAlignedBox3.hpp"

namespace Thea {
namespace Graphics {

class IBufferPool;
class IBuffer;

} // namespace Graphics
} // namespace Thea

namespace Browse3D {

/** The model manipulated by the user. */
class PointCloud : public NamedObject, public GraphicsWidget
{
  private:
    /** A point on the surface. */
    struct Point
    {
      Point() {}
      Point(Vector3 const & p_, Vector3 const & n_ = Vector3::Zero()) : p(p_), n(n_) {}

      Vector3 p;
      Vector3 n;
    };

  public:
    THEA_DECL_SMART_POINTERS(PointCloud)

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
    bool setPointColors(Array<ColorRgba> const & colors_);

    AxisAlignedBox3 const & getBounds() const;

    void updateBounds();

    int8 draw(Graphics::IRenderSystem * render_system, Graphics::IRenderOptions const * options = nullptr) const;

  private:
    /** Invalidate the bounding box of the point cloud. */
    void invalidateBounds();

    /**
     * Get the path to the file in which features are stored, given the path to the point cloud and optionally a file/directory
     * containing the features.
     */
    std::string getDefaultFeaturesFilename(std::string const & filename, std::string const & features_path) const;

    /** Reconstruct an approximate surface from the point cloud. */
    void reconstructSurface(Graphics::IRenderSystem * render_system, Graphics::IRenderOptions const * options = nullptr) const;

    /** Get the color of a point. */
    ColorRgba getColor(size_t point_index) const;

    /** Upload graphics buffers etc to GPU. */
    bool uploadToGraphicsSystem(Graphics::IRenderSystem & render_system);

    Array<Point> points;

    bool has_normals;
    bool normals_are_normalized;
    AxisAlignedBox3 bounds;

    Array<ColorRgba> colors;
    Array< Array<Real> > features;

    bool has_graph;
    Array< Array<intx> > graph;

    int changed_buffers;                   ///< A bitwise OR of the flags of the buffers that have changed.
    Graphics::IBufferPool * buf_pool;      ///< GPU buffer area.
    Graphics::IBuffer     * vertices_buf;  ///< GPU buffer for vertex positions.
    Graphics::IBuffer     * colors_buf;    ///< GPU buffer for vertex colors.

}; // class PointCloud

} // namespace Browse3D

#endif
