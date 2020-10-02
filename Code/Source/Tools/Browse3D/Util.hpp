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

#ifndef __Browse3D_Util_hpp__
#define __Browse3D_Util_hpp__

#include "Common.hpp"
#include "../../Colors.hpp"
#include "../../Ray3.hpp"

namespace Thea {

class Image;

namespace Graphics {

class IRenderSystem;
class Camera;

} // namespace Graphics

} // namespace Thea

namespace Browse3D {

// Draw a sphere.
void drawSphere(Graphics::IRenderSystem & render_system, Vector3 const & center, Real radius, int num_steps = 16);

// Draw a capsule (cylinder with edges capped by hemispheres).
void drawCapsule(Graphics::IRenderSystem & render_system, Vector3 const & base_center, Vector3 const & top_center, Real radius,
                 int num_steps = 16);

// Draw a torus with the given center and primary axes given by the unit vectors u, v.
void drawTorus(Graphics::IRenderSystem & render_system, Vector3 const & center, Vector3 const & u, Vector3 const & v,
               Real radius, Real width, int num_major_steps = 16, int num_minor_steps = 8, bool alternate_dark_light = false,
               ColorRgba const & color1 = ColorRgba(1, 0, 0, 1), ColorRgba const & color2 = ColorRgba(0, 1, 0, 1));

// Get the number of colors in the standard palette.
int numPaletteColors();

// Get the i-th color in the standard palette.
ColorRgba const & getPaletteColor(intx i);

// Map a label to a color.
ColorRgba getLabelColor(std::string const & label);

// Compute a picking ray, given a screen point and a camera.
Ray3 computePickRay(wxRealPoint const & p, Graphics::Camera const & camera, int width, int height);

// Map mouse drags to transforms.
Vector3 dragToTranslation(wxPoint const & start, wxPoint const & end, int width, int height, Graphics::Camera const & camera,
                          Real object_distance);
Matrix3 dragToRotation(wxPoint const & start, wxPoint const & end, int width, int height, Graphics::Camera const & camera);
Matrix3 dragToRotationAroundAxis(wxPoint const & start, Vector3 const & start_pick, wxPoint const & end, Vector3 const & axis,
                                 Vector3 const & center, int width, int height, Graphics::Camera const & camera);
Matrix3 dragToJoystickRotation(wxPoint const & start, wxPoint const & end, Vector3 const & center, Real offset, int width,
                               int height, Graphics::Camera const & camera);
Real dragToScale(wxPoint const & start, wxPoint const & end, int width, int height, Graphics::Camera const & camera);

// Load an image from a file. Fixes channel ordering if necessary.
bool loadImage(Image & image, std::string const & path);

// Make sure image channels are in order RGB[A]. Currently works only for RGB_8U and RGBA_8U.
void fixChannelOrdering(Image & image);

} // namespace Browse3D

#endif
