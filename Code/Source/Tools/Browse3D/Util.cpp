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
// First version: 2011
//
//============================================================================

#include "Util.hpp"
#include "Math.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/IRenderSystem.hpp"
#include "../../Ball3.hpp"
#include "../../FileSystem.hpp"
#include "../../Image.hpp"
#include <wx/gdicmn.h>
#include <cmath>
#include <functional>
#include <limits>

namespace Browse3D {

// Polar axis (center to pole) is w
void
drawHemisphere(Graphics::IRenderSystem & render_system, Vector3 const & center, Vector3 const & u, Vector3 const & v,
               Vector3 const & w, Real radius, int num_longitude_steps, int num_latitude_steps)
{
  using namespace Graphics;

  alwaysAssertM(num_longitude_steps >= 3, "A hemisphere must be drawn with at least 3 longitudinal steps");
  alwaysAssertM(num_latitude_steps >= 1, "A hemisphere must be drawn with at least 3 latitudinal steps");

  float longitude_conv_factor = (float)Math::twoPi() / num_longitude_steps;
  float latitude_conv_factor = (float)Math::halfPi() / num_latitude_steps;

  float prev_s_lat = 0;
  float prev_c_lat = 1;

  Vector3 offset_dir, offset;

  for (int i = 1; i <= num_latitude_steps; ++i)
  {
    float latitude = latitude_conv_factor * i;
    float s_lat = Math::fastSin(latitude);
    float c_lat = Math::fastCos(latitude);

    if (i > 1)
    {
      render_system.beginPrimitive(IRenderSystem::Primitive::QUAD_STRIP);

        for (int j = 0; j <= num_longitude_steps; ++j)
        {
          float longitude = longitude_conv_factor * j;
          float s_lng = Math::fastSin(longitude);
          float c_lng = Math::fastCos(longitude);

          offset_dir = (c_lng * prev_s_lat) * u + (s_lng * prev_s_lat) * v + prev_c_lat * w;
          offset = radius * offset_dir;
          render_system.sendNormal(offset_dir.data());
          render_system.sendVertex(3, Vector3(center + offset).data());

          offset_dir = (c_lng * s_lat) * u + (s_lng * s_lat) * v + c_lat * w;
          offset = radius * offset_dir;
          render_system.sendNormal(offset_dir.data());
          render_system.sendVertex(3, Vector3(center + offset).data());
        }

      render_system.endPrimitive();
    }
    else
    {
      render_system.beginPrimitive(IRenderSystem::Primitive::TRIANGLE_FAN);

        render_system.sendNormal(w.data());
        render_system.sendVertex(3, Vector3(center + radius * w).data());

        for (int j = 0; j <= num_longitude_steps; ++j)
        {
          float longitude = longitude_conv_factor * j;
          float s_lng = Math::fastSin(longitude);
          float c_lng = Math::fastCos(longitude);

          offset_dir = (c_lng * s_lat) * u + (s_lng * s_lat) * v + c_lat * w;
          offset = radius * offset_dir;
          render_system.sendNormal(offset_dir.data());
          render_system.sendVertex(3, Vector3(center + offset).data());
        }

      render_system.endPrimitive();
    }

    prev_s_lat = s_lat;
    prev_c_lat = c_lat;
  }
}

void
drawSphere(Graphics::IRenderSystem & render_system, Vector3 const & center, Real radius, int num_steps)
{
  int num_latitude_steps = (int)std::ceil(num_steps / 4.0);

  drawHemisphere(render_system, center, Vector3::UnitX(), Vector3::UnitY(),  Vector3::UnitZ(), radius, num_steps,
                 num_latitude_steps);
  drawHemisphere(render_system, center,  Vector3::UnitX(), Vector3::UnitY(), -Vector3::UnitZ(), radius, num_steps,
                 num_latitude_steps);
}

void
drawCylinder(Graphics::IRenderSystem & render_system, Vector3 const & base_center, Vector3 const & top_center, Vector3 const & u,
             Vector3 const & v, Real radius, int num_steps)
{
  using namespace Graphics;

  alwaysAssertM(num_steps >= 3, "A cylinder must be drawn with at least 3 steps");

  float conv_factor = (float)Math::twoPi() / num_steps;

  render_system.beginPrimitive(IRenderSystem::Primitive::QUAD_STRIP);

    for (int i = 0; i <= num_steps; ++i)
    {
      float angle = conv_factor * i;
      float s = Math::fastSin(angle);
      float c = Math::fastCos(angle);
      Vector3 offset_dir = s * u + c * v;
      Vector3 offset = radius * offset_dir;
      render_system.sendNormal(offset_dir.data());
      render_system.sendVertex(3, Vector3(base_center + offset).data());
      render_system.sendVertex(3, Vector3(top_center + offset).data());
    }

  render_system.endPrimitive();
}

void
drawCapsule(Graphics::IRenderSystem & render_system, Vector3 const & base_center, Vector3 const & top_center, Real radius,
            int num_steps)
{
  Vector3 w = (top_center - base_center).normalized();
  Matrix3 basis = Math::orthonormalBasis(w);
  Vector3 u = basis.col(0);
  Vector3 v = basis.col(1);

  // Sides
  drawCylinder(render_system, base_center, top_center, u, v, radius, num_steps);

  // Base and top
  int num_latitude_steps = std::ceil(num_steps / 4.0);
  drawHemisphere(render_system, base_center, u, v, -w, radius, num_steps, num_latitude_steps);
  drawHemisphere(render_system, top_center,  u, v,  w, radius, num_steps, num_latitude_steps);
}

void
drawTorus(Graphics::IRenderSystem & render_system, Vector3 const & center, Vector3 const & u, Vector3 const & v, Real radius,
          Real width, int num_major_steps, int num_minor_steps, bool alternate_colors, ColorRgba const & color1,
          ColorRgba const & color2)
{
  using namespace Graphics;

  float major_conv_factor = (float)Math::twoPi() / num_major_steps;
  float minor_conv_factor = (float)Math::twoPi() / num_minor_steps;

  Vector3 w = u.cross(v);

  if (alternate_colors)
    render_system.pushColorFlags();

  render_system.beginPrimitive(IRenderSystem::Primitive::QUAD_STRIP);

    for (int i = 0; i < num_minor_steps; ++i)
    {
      float minor_angle = minor_conv_factor * i;
      float s_min = Math::fastSin(minor_angle);
      float c_min = Math::fastCos(minor_angle);

      float s_min_offset = Math::fastSin(minor_conv_factor * (i + 1)) - s_min;
      float c_min_offset = Math::fastCos(minor_conv_factor * (i + 1)) - c_min;

      for (int j = 0; j <= num_major_steps; ++j)
      {
        float major_angle = major_conv_factor * j;
        float s_maj = Math::fastSin(major_angle);
        float c_maj = Math::fastCos(major_angle);

        Vector3 p(c_maj * (radius + width * c_min),
                  s_maj * (radius + width * c_min),
                  width * s_min);
        Vector3 n(c_maj * c_min, s_maj * c_min, s_min);

        if (alternate_colors)
          render_system.setColor(j & 0x01 ? color2.data() : color1.data());

        render_system.sendNormal(Vector3(n.x() * u + n.y() * v + n.z() * w).data());
        render_system.sendVertex(3, Vector3(center + p.x() * u + p.y() * v + p.z() * w).data());

        p += Vector3(c_maj * width * c_min_offset, s_maj * width * c_min_offset, width * s_min_offset);
        n += Vector3(c_maj * c_min_offset, s_maj * c_min_offset, s_min_offset);
        render_system.sendNormal(Vector3(n.x() * u + n.y() * v + n.z() * w).data());
        render_system.sendVertex(3, Vector3(center + p.x() * u + p.y() * v + p.z() * w).data());
      }
    }

  render_system.endPrimitive();

  if (alternate_colors)
    render_system.popColorFlags();
}

Ray3
computePickRay(wxRealPoint const & p, Graphics::Camera const & camera, int width, int height)
{
  Vector2 screen_pos = Vector2(2 * p.x / (Real)width - 1, 1 - 2 * p.y / (Real)height);
  return camera.computePickRay(screen_pos);
}

Vector3
dragToTranslation(wxPoint const & start, wxPoint const & end, int width, int height, Graphics::Camera const & camera,
                  Real object_distance)
{
  // Remember pixel coordinates increase top to bottom, so diff.y() is downwards in camera space

  wxPoint diff = end - start;

  Real d_scale = object_distance / camera.getNearDistance();
  Real horz_scale = d_scale * (camera.getRightMargin() - camera.getLeftMargin()) / (Real)width;
  Real vert_scale = d_scale * (camera.getTopMargin() - camera.getBottomMargin()) / (Real)height;

  return diff.x * horz_scale * camera.getRightDirection()
       - diff.y * vert_scale * camera.getUpDirection();
}

Matrix3
dragToRotation(wxPoint const & start, wxPoint const & end, int width, int height, Graphics::Camera const & camera)
{
  // Remember pixel coordinates increase top to bottom, so diff.y() is downwards in camera space

  wxPoint diff = end - start;
  if (diff.x == 0 && diff.y == 0) return Matrix3::Identity();

  Vector3 axis = (Real)diff.y * camera.getRightDirection() + (Real)diff.x * camera.getUpDirection();

  static Real const ROT_SPEED = 5;
  int size = width < height ? width : height;
  Real angle = ROT_SPEED * Vector2(diff.x, -diff.y).norm() / size;

  return Math::rotationAxisAngle(axis, angle);
}

Matrix3
dragToRotationAroundAxis(wxPoint const & start, Vector3 const & start_pick, wxPoint const & end, Vector3 const & axis,
                         Vector3 const & center, int width, int height, Graphics::Camera const & camera)
{
  // Remember pixel coordinates increase top to bottom, so diff.y() is downwards in camera space

  Vector3 axis_dir = axis.normalized();
  Vector3 offset = start_pick - center;
  Vector3 perp_offset = offset - offset.dot(axis_dir) * axis_dir;
  Vector3 tangent = axis.cross(perp_offset).normalized();
  Vector2 proj_tangent = (camera.getWorldToCameraTransform().getRotation() * tangent).head<2>().normalized();

  wxPoint diff = end - start;
  if (diff.x == 0 && diff.y == 0) return Matrix3::Identity();

  static Real const ROT_SPEED = 0.02f;
  Real angle = ROT_SPEED * (diff.x * proj_tangent.x() - diff.y * proj_tangent.y());

  return Math::rotationAxisAngle(axis, angle);
}

Matrix3
dragToJoystickRotation(wxPoint const & start, wxPoint const & end, Vector3 const & center, Real offset, int width, int height,
                       Graphics::Camera const & camera)
{
  // Remember pixel coordinates increase top to bottom, so diff.y() is downwards in camera space

  Ray3 rays[2] = { computePickRay(start, camera, width, height), computePickRay(end, camera, width, height) };

  Vector3 dirs[2];
  for (int i = 0; i < 2; ++i)
  {
    Ball3 ball(center, offset);
    Real t = ball.rayIntersectionTime(rays[i]);
    Vector3 p = (t >= 0 ? rays[i].getPoint(t) : rays[i].closestPoint(center));

    dirs[i] = (p - center).normalized();
  }

  return Math::rotationArc(dirs[0], dirs[1], false);
}

Real
dragToScale(wxPoint const & start, wxPoint const & end, int width, int height, Graphics::Camera const & camera)
{
  // Remember pixel coordinates increase top to bottom, so diff.y() is downwards in camera space

  wxPoint diff = end - start;

  static Real const SCALE_INC_SPEED = 4;
  static Real const SCALE_DEC_SPEED = 2;
  static Real const MIN_SCALE = 0.25;

  if (diff.y < 0)  // drag up => increase
    return std::max(1 - SCALE_INC_SPEED * diff.y / (Real)height, MIN_SCALE);
  else  // drag down => decrease
    return std::max(1 - SCALE_DEC_SPEED * diff.y / (Real)height, MIN_SCALE);
}

int const NUM_PALETTE_COLORS = 24;
ColorRgba COLOR_PALETTE[NUM_PALETTE_COLORS] = {
  ColorRgb::fromARGB(0x298edb),
  ColorRgb::fromARGB(0x982411),
  ColorRgb::fromARGB(0x6d4e25),
  ColorRgb::fromARGB(0x1b5043),
  ColorRgb::fromARGB(0x6e7662),
  ColorRgb::fromARGB(0xa08b00),
  ColorRgb::fromARGB(0x58427b),
  ColorRgb::fromARGB(0x1d2f5b),
  ColorRgb::fromARGB(0xac5e34),
  ColorRgb::fromARGB(0x804055),
  ColorRgb::fromARGB(0x6d7a00),
  ColorRgb::fromARGB(0x572e2c),

  // Invert each color above
  ColorRgb::fromARGB(~0x298edb),
  ColorRgb::fromARGB(~0x982411),
  ColorRgb::fromARGB(~0x6d4e25),
  ColorRgb::fromARGB(~0x1b5043),
  ColorRgb::fromARGB(~0x6e7662),
  ColorRgb::fromARGB(~0xa08b00),
  ColorRgb::fromARGB(~0x58427b),
  ColorRgb::fromARGB(~0x1d2f5b),
  ColorRgb::fromARGB(~0xac5e34),
  ColorRgb::fromARGB(~0x804055),
  ColorRgb::fromARGB(~0x6d7a00),
  ColorRgb::fromARGB(~0x572e2c)
};

int
numPaletteColors()
{
  return NUM_PALETTE_COLORS;
}

ColorRgba const &
getPaletteColor(intx i)
{
  intx n = numPaletteColors();
  intx index;
  if (i < 0)
    index = (n - ((-i) % n)) % n;
  else
    index = i % n;

  return COLOR_PALETTE[index];
}

ColorRgba
getLabelColor(std::string const & label)
{
  std::hash<std::string> hasher;
  Random rnd((uint32)hasher(label));

  return ColorRgba(rnd.uniform01(), rnd.uniform01(), rnd.uniform01(), 1);
}

bool
loadImage(Image & image, std::string const & path)
{
  if (!FileSystem::fileExists(path)) return false;

  image.load(path);
  if (image.isValid())
  {
    // This used to be necessary on Linux. No longer, since we've fixed TextureFormat::fromImageType() to detect and handle BGR
    // images.
    // fixChannelOrdering(image);

    return true;
  }

  return false;
}

void
fixChannelOrdering(Image & image)
{
  // Make sure the red channel is at pixel[0] and the blue channel at pixel[2]
  int bytes_pp = 0;
  switch (image.getType())
  {
    case Image::Type::RGB_8U: bytes_pp = 3; break;
    case Image::Type::RGBA_8U: bytes_pp = 4; break;
    default: break;
  }

  if (bytes_pp > 0)
  {
    for (int r = 0; r < image.getHeight(); ++r)
    {
      uint8 * pixels = (uint8 *)image.getScanLine(r);
      for (int c = 0; c < image.getWidth(); ++c)
      {
        uint8 * pixel = pixels + bytes_pp * c;

        uint8 red    =  pixel[Image::Channel::RED];
        uint8 green  =  pixel[Image::Channel::GREEN];
        uint8 blue   =  pixel[Image::Channel::BLUE];

        pixel[0] = red;
        pixel[1] = green;
        pixel[2] = blue;
      }
    }
  }
}

} // namespace Browse3D
