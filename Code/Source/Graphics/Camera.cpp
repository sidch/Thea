//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#include "Camera.hpp"

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
        ? Matrix4::orthogonalProjection(left, right, bottom, top, near_dist, far_dist,
                                        proj_y_dir == ProjectedYDirection::UP ? true : false)
        : Matrix4::perspectiveProjection(left, right, bottom, top, near_dist, far_dist,
                                         proj_y_dir == ProjectedYDirection::UP ? true : false);

    cached_inv_proj_transform = cached_proj_transform.inverse();

    proj_changed = false;
  }
}

Ray3
Camera::computePickRay(Vector2 const & screen_pos) const
{
  Vector2 p = 0.5f * (screen_pos + Vector2(1, 1));
  Ray3 view_ray(Vector3::zero(), Vector3(left + p.x() * (right - left), bottom + p.y() * (top - bottom), -near_dist).unit());
  return view_ray.toWorldSpace(frame);
}

std::string
Camera::toString() const
{
  std::ostringstream oss;
  oss << "Frame = [R: "  << frame.getRotation().toString() << ", T: " << frame.getTranslation().toString() << ']'
      << ", ProjectionType = " << ((projection_type == ProjectionType::ORTHOGRAPHIC) ? "Orthographic" : "Perspective")
      << ", Left = "   << left   << ", Right = " << right
      << ", Bottom = " << bottom << ", Top = "   << top
      << ", NearDist = " << near_dist << ", FarDist = " << far_dist
      << ", Projected Y increases " << (proj_y_dir == ProjectedYDirection::UP ? "upwards" : "downwards");
  return oss.str();
}

void
Camera::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  output.setEndianness(getEndianness());

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

void
Camera::deserialize(BinaryInputStream & input, Codec const & codec)
{
  input.setEndianness(getEndianness());

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

} // namespace Graphics
} // namespace Thea
