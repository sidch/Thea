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

namespace Thea {
namespace Graphics {

/** Abstract base class for options controlling the display of a Drawable, safe to pass across DLL boundaries. */
class THEA_API AbstractRenderOptions
{
  public:
    /** Destructor. */
    virtual ~AbstractRenderOptions() = 0;

    /** Send colors to the rendersystem? */
    virtual int8 sendColors() const = 0;

    /** Set whether colors will be sent to the rendersystem. */
    virtual AbstractRenderOptions & setSendColors(int8 value) = 0;

    /** Send normals to the rendersystem? */
    virtual int8 sendNormals() const = 0;

    /** Set whether normals will be sent to the rendersystem. */
    virtual AbstractRenderOptions & setSendNormals(int8 value) = 0;

    /** Send texture coordinates to the rendersystem? */
    virtual int8 sendTexCoords() const = 0;

    /** Set whether texture coordinates will be sent to the rendersystem. */
    virtual AbstractRenderOptions & setSendTexCoords(int8 value) = 0;

    /** Use vertex normals instead of face normals (for smooth shading)? */
    virtual int8 useVertexNormals() const = 0;

    /** Set whether vertex normals be used instead of face normals (for smooth shading). */
    virtual AbstractRenderOptions & setUseVertexNormals(int8 value) = 0;

    /** Use data at vertices instead of faces? Does <b>not</b> apply to normals, see useVertexNormals(). */
    virtual int8 useVertexData() const = 0;

    /**
     * Set whether data at vertices will be used instead of faces. Does <b>not</b> apply to normals, see setUseVertexNormals().
     */
    virtual AbstractRenderOptions & setUseVertexData(int8 value) = 0;

    /** Draw mesh faces? */
    virtual int8 drawFaces() const = 0;

    /** Set whether mesh faces will be drawn. */
    virtual AbstractRenderOptions & setDrawFaces(int8 value) = 0;

    /** Draw mesh edges? */
    virtual int8 drawEdges() const = 0;

    /** Set whether mesh edges will be drawn. */
    virtual AbstractRenderOptions & setDrawEdges(int8 value) = 0;

    /** Override edge-specific colors with the value of edgeColor() when drawing edges? */
    virtual int8 overrideEdgeColor() const = 0;

    /** Set whether edge-specific colors will be overridden with the value of edgeColor() when drawing edges. */
    virtual AbstractRenderOptions & setOverrideEdgeColor(int8 value) = 0;

    /** Color for drawing edges when overrideEdgeColor() is true, stored as an RGBA quadruplet. */
    virtual Real const * edgeColor() const = 0;

    /** Set the color for drawing edges when overrideEdgeColor() is true, passed as an RGBA quadruplet. */
    virtual AbstractRenderOptions & setEdgeColor(Real const * rgba) = 0;

}; // class AbstractRenderOptions

// Pure virtual destructor should have a body
// http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
inline AbstractRenderOptions::~AbstractRenderOptions() {}

/** %Options controlling the display of a Drawable. */
class THEA_API RenderOptions : public AbstractRenderOptions
{
  private:
    int8       send_normals;
    int8       send_colors;
    int8       send_texcoords;
    int8       use_vertex_normals;
    int8       use_vertex_data;
    int8       draw_faces;
    int8       draw_edges;
    int8       override_edge_color;
    ColorRGBA  edge_color;

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
      edge_color(ColorRGB::white())
    {}

    /** Copy from base class object. */
    RenderOptions(AbstractRenderOptions const & rhs) { *this = rhs; }

    /** Assign from base class object. */
    RenderOptions & operator=(AbstractRenderOptions const & rhs)
    {
      send_normals         =  rhs.sendNormals();
      send_colors          =  rhs.sendColors();
      send_texcoords       =  rhs.sendTexCoords();
      use_vertex_normals   =  rhs.useVertexNormals();
      use_vertex_data      =  rhs.useVertexData();
      draw_faces           =  rhs.drawFaces();
      draw_edges           =  rhs.drawEdges();
      override_edge_color  =  rhs.overrideEdgeColor();
      edge_color           =  ColorRGBA(rhs.edgeColor());

      return *this;
    }

    static RenderOptions const & defaults() { static RenderOptions def; return def; }

    int8 sendColors() const { return send_colors; }
    AbstractRenderOptions & setSendColors(int8 value) { send_colors = value; return *this; }
    int8 sendNormals() const { return send_normals; }
    AbstractRenderOptions & setSendNormals(int8 value) { send_normals = value; return *this; }
    int8 sendTexCoords() const { return send_texcoords; }
    AbstractRenderOptions & setSendTexCoords(int8 value) { send_texcoords = value; return *this; }
    int8 useVertexNormals() const { return use_vertex_normals; }
    AbstractRenderOptions & setUseVertexNormals(int8 value) { use_vertex_normals = value; return *this; }
    int8 useVertexData() const { return use_vertex_data; }
    AbstractRenderOptions & setUseVertexData(int8 value) { use_vertex_data = value; return *this; }
    int8 drawFaces() const { return draw_faces; }
    AbstractRenderOptions & setDrawFaces(int8 value) { draw_faces = value; return *this; }
    int8 drawEdges() const { return draw_edges; }
    AbstractRenderOptions & setDrawEdges(int8 value) { draw_edges = value; return *this; }
    int8 overrideEdgeColor() const { return override_edge_color; }
    AbstractRenderOptions & setOverrideEdgeColor(int8 value) { override_edge_color = value; return *this; }
    Real const * edgeColor() const { return edge_color.data(); }
    AbstractRenderOptions & setEdgeColor(Real const * rgba) { edge_color = ColorRGBA(rgba); return *this; }

}; // class RenderOptions

} // namespace Graphics
} // namespace Thea

#endif
