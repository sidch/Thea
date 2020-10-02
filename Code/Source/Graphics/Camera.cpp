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

#include "Camera.hpp"
#include "../MatrixWrapper.hpp"

THEA_INSTANTIATE_SMART_POINTERS(Thea::Graphics::Camera)

namespace Thea {
namespace Graphics {

Camera::Camera(CoordinateFrame3 const & frame_, ProjectionType projection_type_, Real left_, Real right_, Real bottom_,
               Real top_, Real near_dist_, Real far_dist_, ProjectedYDirection proj_y_dir_)
: frame(frame_), projection_type(projection_type_), left(left_), right(right_), bottom(bottom_), top(top_),
  near_dist(near_dist_), far_dist(far_dist_), proj_y_dir(proj_y_dir_), proj_changed(true)
{}

void
Camera::set(CoordinateFrame3 const & frame_, ProjectionType projection_type_, Real left_, Real right_, Real bottom_,
            Real top_, Real near_dist_, Real far_dist_, ProjectedYDirection proj_y_dir_)
{
  setFrame(frame_);
  setProjection(projection_type_, left_, right_, bottom_, top_, near_dist_, far_dist_, proj_y_dir_);
}

void
Camera::setFrame(CoordinateFrame3 const & frame_)
{
  frame = frame_;
}

void
Camera::setProjection(ProjectionType projection_type_, Real left_, Real right_, Real bottom_, Real top_, Real near_dist_,
                      Real far_dist_, ProjectedYDirection proj_y_dir_)
{
  projection_type = projection_type_;
  left = left_;
  right = right_;
  bottom = bottom_;
  top = top_;
  near_dist = near_dist_;
  far_dist = far_dist_;
  proj_y_dir = proj_y_dir_;

  proj_changed = true;
}

Matrix4
Camera::getProjectionTransform() const
{
  updateCachedProjectionTransform();
  return cached_proj_transform;
}

Matrix4
Camera::getInverseProjectionTransform() const
{
  updateCachedProjectionTransform();
  return cached_inv_proj_transform;
}

void
Camera::updateCachedProjectionTransform() const
{
  if (proj_changed)
  {
    cached_proj_transform = (projection_type == ProjectionType::ORTHOGRAPHIC)
        ? Math::orthogonalProjection(left, right, bottom, top, near_dist, far_dist,
                                     proj_y_dir == ProjectedYDirection::UP ? true : false)
        : Math::perspectiveProjection(left, right, bottom, top, near_dist, far_dist,
                                      proj_y_dir == ProjectedYDirection::UP ? true : false);

    cached_inv_proj_transform = cached_proj_transform.inverse();

    proj_changed = false;
  }
}

Ray3
Camera::computePickRay(Vector2 const & screen_pos) const
{
  Vector2 p = 0.5f * (screen_pos + Vector2(1, 1));
  Ray3 view_ray(Vector3::Zero(),
                Vector3(left + p.x() * (right - left), bottom + p.y() * (top - bottom), -near_dist).normalized());
  return view_ray.toWorldSpace(frame);
}

void
Camera::makeCurrent(IRenderSystem * render_system) const
{
  alwaysAssertM(render_system, "Camera: Cannot set current camera on null rendersystem");

  Matrix4 proj = getProjectionTransform();
  Matrix4 modelview = getWorldToCameraTransform().homogeneous();
  MatrixWrapper wp(&proj), wmv(&modelview);

  render_system->setMatrixMode(IRenderSystem::MatrixMode::PROJECTION);
  render_system->setMatrix(&wp);

  render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW);
  render_system->setMatrix(&wmv);
}

std::string
Camera::toString() const
{
  std::ostringstream oss;
  oss << "Frame = " << frame.toString()
      << ", ProjectionType = " << ((projection_type == ProjectionType::ORTHOGRAPHIC) ? "Orthographic" : "Perspective")
      << ", Left = "   << left   << ", Right = " << right
      << ", Bottom = " << bottom << ", Top = "   << top
      << ", NearDist = " << near_dist << ", FarDist = " << far_dist
      << ", Projected Y increases " << (proj_y_dir == ProjectedYDirection::UP ? "upwards" : "downwards");
  return oss.str();
}

void
Camera::read(BinaryInputStream & input, Codec const & codec, bool read_block_header)
{
  if (read_block_header)
    input.skip(Codec::BlockHeader::SERIALIZED_LENGTH);  // header is not needed

  { BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);

    frame = input.readCoordinateFrame3();

    uint8 pt = input.readUInt8();
    projection_type = (pt == 0 ? ProjectionType::ORTHOGRAPHIC : ProjectionType::PERSPECTIVE);

    left       =  static_cast<Real>(input.readFloat32());
    right      =  static_cast<Real>(input.readFloat32());
    bottom     =  static_cast<Real>(input.readFloat32());
    top        =  static_cast<Real>(input.readFloat32());
    near_dist  =  static_cast<Real>(input.readFloat32());
    far_dist   =  static_cast<Real>(input.readFloat32());

    uint8 pyd = input.readUInt8();
    proj_y_dir = (pyd == 0 ? ProjectedYDirection::UP : ProjectedYDirection::DOWN);
  }
}

void
Camera::write(BinaryOutputStream & output, Codec const & codec, bool write_block_header) const
{
  Codec::BlockHeader header("Camera");
  if (write_block_header)
    header.markAndSkip(output);

  { BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);

    output.writeCoordinateFrame3(frame);
    output.writeUInt8(projection_type == ProjectionType::ORTHOGRAPHIC ? 0 : 1);
    output.writeFloat32(static_cast<float32>(left));
    output.writeFloat32(static_cast<float32>(right));
    output.writeFloat32(static_cast<float32>(bottom));
    output.writeFloat32(static_cast<float32>(top));
    output.writeFloat32(static_cast<float32>(near_dist));
    output.writeFloat32(static_cast<float32>(far_dist));
    output.writeUInt8(proj_y_dir == ProjectedYDirection::UP ? 0 : 1);
  }

  if (write_block_header)
    header.calcAndWrite(output);
}

} // namespace Graphics
} // namespace Thea
