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
#include "Camera.hpp"

namespace Thea {
namespace Graphics {

/** %Options controlling the display of a DrawableObject. */
class THEA_API RenderOptions
{
  private:
    bool     _send_normals;
    bool     _send_colors;
    bool     _send_texcoords;
    bool     _use_vertex_data;
    bool     _draw_faces;
    bool     _draw_edges;
    bool     _override_edge_color;
    Color4   _edge_color;
    Camera   _viewing_camera;

  public:
    /** Default constructor. */
    RenderOptions()
    : _send_normals(true),
      _send_colors(true),
      _send_texcoords(false),
      _use_vertex_data(true),
      _draw_faces(true),
      _draw_edges(false),
      _override_edge_color(false),
      _edge_color(Color3::white()),
      _viewing_camera()
    {}

    /** Get the default set of options. */
    static RenderOptions const & defaults() { static RenderOptions def; return def; }

    /** Send colors to the rendersystem? */
    bool sendColors() const { return _send_colors; }

    /** Send colors to the rendersystem? */
    bool & sendColors() { return _send_colors; }

    /** Send normals to the rendersystem? */
    bool sendNormals() const { return _send_normals; }

    /** Send normals to the rendersystem? */
    bool & sendNormals() { return _send_normals; }

    /** Send texture coordinates to the rendersystem? */
    bool sendTexCoords() const { return _send_texcoords; }

    /** Send texture coordinates to the rendersystem? */
    bool & sendTexCoords() { return _send_texcoords; }

    /** Use data at vertices instead of faces? */
    bool useVertexData() const { return _use_vertex_data; }

    /** Use data at vertices instead of faces? */
    bool & useVertexData() { return _use_vertex_data; }

    /** Draw polygon faces? */
    bool drawFaces() const { return _draw_faces; }

    /** Draw polygon faces? */
    bool & drawFaces() { return _draw_faces; }

    /** Draw polygon edges? */
    bool drawEdges() const { return _draw_edges; }

    /** Draw polygon edges? */
    bool & drawEdges() { return _draw_edges; }

    /** Override edge-specific colors with the value of edgeColor() when drawing edges? */
    bool overrideEdgeColor() const { return _override_edge_color; }

    /** Override edge-specific colors with the value of edgeColor() when drawing edges? */
    bool & overrideEdgeColor() { return _override_edge_color; }

    /** Color for drawing edges when overrideEdgeColor() is true. */
    Color4 const & edgeColor() const { return _edge_color; }

    /** Color for drawing edges when overrideEdgeColor() is true. */
    Color4 & edgeColor() { return _edge_color; }

    /** The camera being used for rendering. */
    Camera const & viewingCamera() const { return _viewing_camera; }

    /** The camera being used for rendering. */
    Camera & viewingCamera() { return _viewing_camera; }

}; // class RenderOptions

} // namespace Thea
} // namespace Graphics

#endif
