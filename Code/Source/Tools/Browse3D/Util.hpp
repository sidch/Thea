//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Browse3D_Util_hpp__
#define __Browse3D_Util_hpp__

#include "Common.hpp"
#include "../../Colors.hpp"
#include "../../Ray3.hpp"

namespace Thea {

class Image;

namespace Graphics {

class RenderSystem;
class Camera;

} // namespace Graphics

} // namespace Thea

namespace Browse3D {

// Draw a sphere.
void drawSphere(Graphics::RenderSystem & render_system, Vector3 const & center, Real radius, int num_steps = 16);

// Draw a capsule (cylinder with edges capped by hemispheres).
void drawCapsule(Graphics::RenderSystem & render_system, Vector3 const & base_center, Vector3 const & top_center, Real radius,
                 int num_steps = 16);

// Draw a torus with the given center and primary axes given by the unit vectors u, v.
void drawTorus(Graphics::RenderSystem & render_system, Vector3 const & center, Vector3 const & u, Vector3 const & v,
               Real radius, Real width, int num_major_steps = 16, int num_minor_steps = 8, bool alternate_dark_light = false,
               ColorRGBA const & color1 = ColorRGBA(1, 0, 0, 1), ColorRGBA const & color2 = ColorRGBA(0, 1, 0, 1));

// Get the number of colors in the standard palette.
int numPaletteColors();

// Get the i-th color in the standard palette.
ColorRGB const & getPaletteColor(long i);

// Map a label to a color.
ColorRGB getLabelColor(std::string const & label);

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
