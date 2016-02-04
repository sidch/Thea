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

#ifndef __Thea_Graphics_RenderOptions_hpp__
#define __Thea_Graphics_RenderOptions_hpp__

#include "../Common.hpp"
#include "../Colors.hpp"
#include "Camera.hpp"

namespace Thea {
namespace Graphics {

/** %Options controlling the display of a DrawableObject. */
class THEA_API RenderOptions
{
  private:
    bool       send_normals;
    bool       send_colors;
    bool       send_texcoords;
    bool       use_vertex_normals;
    bool       use_vertex_data;
    bool       draw_faces;
    bool       draw_edges;
    bool       override_edge_color;
    ColorRGBA  edge_color;
    Camera     viewing_camera;

  public:
    /** Default constructor. */
    RenderOptions()
    : send_normals(true),
      send_colors(true),
      send_texcoords(false),
      use_vertex_normals(true),
      use_vertex_data(true),
      draw_faces(true),
      draw_edges(false),
      override_edge_color(false),
      edge_color(ColorRGB::white()),
      viewing_camera()
    {}

    /** Get the default set of options. */
    static RenderOptions const & defaults() { static RenderOptions def; return def; }

    /** Send colors to the rendersystem? */
    bool sendColors() const { return send_colors; }

    /** Send colors to the rendersystem? */
    bool & sendColors() { return send_colors; }

    /** Send normals to the rendersystem? */
    bool sendNormals() const { return send_normals; }

    /** Send normals to the rendersystem? */
    bool & sendNormals() { return send_normals; }

    /** Send texture coordinates to the rendersystem? */
    bool sendTexCoords() const { return send_texcoords; }

    /** Send texture coordinates to the rendersystem? */
    bool & sendTexCoords() { return send_texcoords; }

    /** Use vertex normals instead of face normals (for smooth shading)? */
    bool useVertexNormals() const { return use_vertex_normals; }

    /** Use vertex normals instead of face normals (for smooth shading)? */
    bool & useVertexNormals() { return use_vertex_normals; }

    /** Use data at vertices instead of faces? Does <b>not</b> apply to normals, see useVertexNormals(). */
    bool useVertexData() const { return use_vertex_data; }

    /** Use data at vertices instead of faces? Does <b>not</b> apply to normals, see useVertexNormals(). */
    bool & useVertexData() { return use_vertex_data; }

    /** Draw polygon faces? */
    bool drawFaces() const { return draw_faces; }

    /** Draw polygon faces? */
    bool & drawFaces() { return draw_faces; }

    /** Draw polygon edges? */
    bool drawEdges() const { return draw_edges; }

    /** Draw polygon edges? */
    bool & drawEdges() { return draw_edges; }

    /** Override edge-specific colors with the value of edgeColor() when drawing edges? */
    bool overrideEdgeColor() const { return override_edge_color; }

    /** Override edge-specific colors with the value of edgeColor() when drawing edges? */
    bool & overrideEdgeColor() { return override_edge_color; }

    /** Color for drawing edges when overrideEdgeColor() is true. */
    ColorRGBA const & edgeColor() const { return edge_color; }

    /** Color for drawing edges when overrideEdgeColor() is true. */
    ColorRGBA & edgeColor() { return edge_color; }

    /** The camera being used for rendering. */
    Camera const & viewingCamera() const { return viewing_camera; }

    /** The camera being used for rendering. */
    Camera & viewingCamera() { return viewing_camera; }

}; // class RenderOptions

} // namespace Graphics
} // namespace Thea

#endif
