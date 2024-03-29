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

#ifndef __Thea_Graphics_Camera_hpp__
#define __Thea_Graphics_Camera_hpp__

#include "../Common.hpp"
#include "../NamedObject.hpp"
#include "../RigidTransform3.hpp"
#include "../Ray3.hpp"
#include "../Serializable.hpp"
#include "IRenderSystem.hpp"

namespace Thea {
namespace Graphics {

/**
 * A camera for viewing a scene.
 *
 * The transformation M from a point v in world space to its image M.v in the camera's projection space is:
 * \code
 * getProjectionTransform() * getWorldToCameraTransform()
 * \endcode
 */
class THEA_API Camera : public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(Camera)

    /** Type of projection used by the camera (enum class). */
    struct THEA_API ProjectionType
    {
      /** Supported values. */
      enum Value
      {
        ORTHOGRAPHIC,  ///< Orthographic projection.
        PERSPECTIVE    ///< Perspective projection.
      };

      THEA_ENUM_CLASS_BODY(ProjectionType)

      THEA_ENUM_CLASS_STRINGS_BEGIN(ProjectionType)
        THEA_ENUM_CLASS_STRING(ORTHOGRAPHIC,  "orthographic")
        THEA_ENUM_CLASS_STRING(PERSPECTIVE,   "perspective")
      THEA_ENUM_CLASS_STRINGS_END(ProjectionType)
    };

    /** The direction in which projected Y coordinates increase, relative to the up vector of the camera frame. */
    struct THEA_API ProjectedYDirection
    {
      /** Supported values. */
      enum Value
      {
        UP,   ///< Projected Y coordinates increase bottom to top.
        DOWN  ///< Projected Y coordinates increase top to bottom.
      };

      THEA_ENUM_CLASS_BODY(ProjectedYDirection)

      THEA_ENUM_CLASS_STRINGS_BEGIN(ProjectedYDirection)
        THEA_ENUM_CLASS_STRING(UP,    "up")
        THEA_ENUM_CLASS_STRING(DOWN,  "down")
      THEA_ENUM_CLASS_STRINGS_END(ProjectedYDirection)
    };

    /**
     * Default constructor. No fields are initialized, this constructor is provided only for compatibility with standard
     * containers.
     */
    Camera() {}

    /**
     * Initializing constructor. \a near_dist_ and \a far_dist_ are <b>positive</b> distances to the near and far clipping
     * planes.
     */
    Camera(CoordinateFrame3 const & frame_, ProjectionType projection_type_, Real left_, Real right_, Real bottom_, Real top_,
           Real near_dist_, Real far_dist_, ProjectedYDirection proj_y_dir_);

    /** Set all the parameters of the camera in one go. */
    void set(CoordinateFrame3 const & frame_, ProjectionType projection_type_, Real left_, Real right_, Real bottom_, Real top_,
             Real near_dist_, Real far_dist_, ProjectedYDirection proj_y_dir_);

    /** Set the coordinate frame. The camera looks along the negative local Z-axis of the frame. */
    void setFrame(CoordinateFrame3 const & frame_);

    /**
     * Set the projection parameters. \a near_dist_ and \a far_dist_ are <b>positive</b> distances to the near and far clipping
     * planes.
     */
    void setProjection(ProjectionType projection_type_, Real left_, Real right_, Real bottom_, Real top_,
                       Real near_dist_, Real far_dist_, ProjectedYDirection proj_y_dir_);

    /**
     * Infer the projection parameters from a projection matrix. The matrix should have been generated exactly as
     * getProjectionTransform() would do it, i.e. using Math::orthogonalProjection() or Math::perspectiveProjection().
     */
    bool setProjection(Matrix4 const & m);

    /** Get the coordinate frame of the camera. The camera looks along the negative local Z-axis of the frame. */
    CoordinateFrame3 const & getFrame() const { return frame; }

    /** Get the world-space position of the camera (the eye position). */
    Vector3 getPosition() const { return frame.getTranslation(); }

    /** Get a unit vector (in world space) in the direction the camera is facing. This is the camera's negative Z axis. */
    Vector3 getLookDirection() const { return frame.lookVector(); }

    /** Get a unit vector (in world space) in the up direction of the camera. This is the camera's Y axis. */
    Vector3 getUpDirection() const { return frame.upVector(); }

    /** Get a unit vector (in world space) in the direction to the right of the camera. This is the camera's X axis. */
    Vector3 getRightDirection() const { return frame.rightVector(); }

    /** Get the projection type of the camera. */
    ProjectionType getProjectionType() const { return projection_type; }

    /** Get the coordinate of the left margin of the viewport, in camera space. */
    Real getLeftMargin() const { return left; }

    /** Get the coordinate of the right margin of the viewport, in camera space. */
    Real getRightMargin() const { return right; }

    /** Get the coordinate of the bottom margin of the viewport, in camera space. */
    Real getBottomMargin() const { return bottom; }

    /** Get the coordinate of the bottom margin of the viewport, in camera space. */
    Real getTopMargin() const { return top; }

    /** Get the distance (always <b>positive</b>) to the near clipping plane, in camera space. */
    Real getNearDistance() const { return near_dist; }

    /** Get the distance (always <b>positive</b>) to the far clipping plane, in camera space. */
    Real getFarDistance() const { return far_dist; }

    /** Get the distance between the left and right margins of the viewport, in camera space. */
    Real getWidth() const { return std::abs(right - left); }

    /** Get the distance between the bottom and top margins of the viewport, in camera space. */
    Real getHeight() const { return std::abs(top - bottom); }

    /** Get the distance between the near and far clipping planes, in camera space. */
    Real getDepthRange() const { return std::abs(far_dist - near_dist); }

    /** Get the direction in which projected Y coordinates increase. */
    ProjectedYDirection getProjectedYDirection() const { return proj_y_dir; }

    /** Set the projection type of the camera. */
    void setProjectionType(ProjectionType projection_type_)
    {
      if (projection_type_ != projection_type)
      {
        projection_type = projection_type_;
        proj_changed = true;
      }
    }

    /** Set the coordinate of the left margin of the viewport, in camera space. */
    void setLeftMargin(Real left_) { left = left_; proj_changed = true; }

    /** Set the coordinate of the right margin of the viewport, in camera space. */
    void setRightMargin(Real right_) { right = right_; proj_changed = true; }

    /** Set the coordinate of the bottom margin of the viewport, in camera space. */
    void setBottomMargin(Real bottom_) { bottom = bottom_; proj_changed = true; }

    /** Set the coordinate of the bottom margin of the viewport, in camera space. */
    void setTopMargin(Real top_) { top = top_; proj_changed = true; }

    /** Get the distance (always <b>positive</b>) to the near clipping plane, in camera space. */
    void setNearDistance(Real near_dist_) { near_dist = near_dist_; proj_changed = true; }

    /** Set the distance (always <b>positive</b>) to the far clipping plane, in camera space. */
    void setFarDistance(Real far_dist_) { far_dist = far_dist_; proj_changed = true; }

    /** Set the direction in which projected Y coordinates increase. */
    void setProjectedYDirection(ProjectedYDirection proj_y_dir_)
    {
      if (proj_y_dir_ != proj_y_dir)
      {
        proj_y_dir = proj_y_dir_;
        proj_changed = true;
      }
    }

    /** Invert the direction in which projected Y coordinates increase. */
    void flipProjectedYDirection()
    {
      if (proj_y_dir == ProjectedYDirection::UP)
        proj_y_dir = ProjectedYDirection::DOWN;
      else
        proj_y_dir = ProjectedYDirection::UP;

      proj_changed = true;
    }

    /** Get the transformation from camera coordinates to world coordinates. */
    RigidTransform3 getCameraToWorldTransform() const { return RigidTransform3(frame); }

    /** Get the transformation from world coordinates to camera coordinates. */
    RigidTransform3 getWorldToCameraTransform() const { return RigidTransform3(frame).inverse(); }

    /** Get the transformation from camera coordinates to projection space coordinates. */
    Matrix4 getProjectionTransform() const;

    /** Get the transformation from projection space coordinates to camera coordinates. */
    Matrix4 getInverseProjectionTransform() const;

    /**
     * Project a point from world space to projection space.
     *
     * @see unproject(), computePickRay()
     */
    Vector3 project(Vector3 const & world_space_point) const
    {
      return Math::hmul(getProjectionTransform(), frame.pointToObjectSpace(world_space_point));
    }

    /**
     * Project a point in homogeneous coordinates from world space to projection space.
     *
     * @see unproject(), computePickRay()
     */
    Vector4 project(Vector4 const & world_space_point) const
    {
      return getProjectionTransform() * (frame.inverse().homogeneous() * world_space_point);
    }

    /**
     * Unproject a point from projection space to world space. If you only know the 2D screen coordinates of the point and not
     * its projected depth, use computePickRay().
     *
     * @see project(), computePickRay()
     */
    Vector3 unproject(Vector3 const & projection_space_point) const
    {
      return frame.pointToWorldSpace(Math::hmul(getInverseProjectionTransform(), projection_space_point));
    }

    /**
     * Unproject a point in homogeneous coordinates from projection space to world space. If you only know the 2D screen
     * coordinates of the point and not its projected depth, use computePickRay().
     *
     * @see project(), computePickRay()
     */
    Vector4 unproject(Vector4 const & projection_space_point) const
    {
      return frame.homogeneous() * (getInverseProjectionTransform() * projection_space_point);
    }

    /** Check if the view is orthographic or not. */
    bool isOrthographic() const { return projection_type == ProjectionType::ORTHOGRAPHIC; }

    /** Check if the view is perspective or not. */
    bool isPerspective() const { return projection_type == ProjectionType::PERSPECTIVE; }

    /**
     * Compute the world-space ray entering the screen at a given position, expressed in normalized coordinates ([-1, 1]^2).
     *
     * The normalized coordinates correspond to the X and Y coordinates of the camera's projection space. In particular, the
     * screen Y coordinate is assumed to increase in the same direction as the projected Y coordinate of the camera.
     */
    Ray3 computePickRay(Vector2 const & screen_pos) const;

    /** Set the camera as the current viewing camera on a rendersystem. */
    void makeCurrent(IRenderSystem * render_system) const;

    /** Get the name of the camera. */
    std::string getName() const
    {
      return (projection_type == ProjectionType::ORTHOGRAPHIC) ? "Orthographic Camera" : "Perspective Camera";
    }

    /** Get a string describing the camera. */
    std::string toString() const;

    void read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false);
    void write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false) const;

  private:
    /** Update the cached projection transform and its inverse. */
    void updateCachedProjectionTransform() const;

    CoordinateFrame3 frame;
    ProjectionType projection_type;
    Real left, right, bottom, top, near_dist, far_dist;
    ProjectedYDirection proj_y_dir;

    mutable bool proj_changed;
    mutable Matrix4 cached_proj_transform;
    mutable Matrix4 cached_inv_proj_transform;

}; // class Camera

} // namespace Graphics
} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Graphics::Camera)

#endif
