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

#ifndef __Thea_Graphics_RenderOptions_hpp__
#define __Thea_Graphics_RenderOptions_hpp__

#include "../Common.hpp"
#include "../Colors.hpp"

namespace Thea {
namespace Graphics {

/** Interface for options controlling the display of an IDrawable, safe to pass across DLL boundaries. */
class THEA_API IRenderOptions
{
  public:
    /** Destructor. */
    virtual ~IRenderOptions() = 0;

    /** Send colors to the rendersystem? */
    virtual int8 THEA_ICALL sendColors() const = 0;

    /** Set whether colors will be sent to the rendersystem. */
    virtual IRenderOptions & THEA_ICALL setSendColors(int8 value) = 0;

    /** Send normals to the rendersystem? */
    virtual int8 THEA_ICALL sendNormals() const = 0;

    /** Set whether normals will be sent to the rendersystem. */
    virtual IRenderOptions & THEA_ICALL setSendNormals(int8 value) = 0;

    /** Send texture coordinates to the rendersystem? */
    virtual int8 THEA_ICALL sendTexCoords() const = 0;

    /** Set whether texture coordinates will be sent to the rendersystem. */
    virtual IRenderOptions & THEA_ICALL setSendTexCoords(int8 value) = 0;

    /** Use vertex normals instead of face normals (for smooth shading)? */
    virtual int8 THEA_ICALL useVertexNormals() const = 0;

    /** Set whether vertex normals be used instead of face normals (for smooth shading). */
    virtual IRenderOptions & THEA_ICALL setUseVertexNormals(int8 value) = 0;

    /** Use data at vertices instead of faces? Does <b>not</b> apply to normals, see useVertexNormals(). */
    virtual int8 THEA_ICALL useVertexData() const = 0;

    /**
     * Set whether data at vertices will be used instead of faces. Does <b>not</b> apply to normals, see setUseVertexNormals().
     */
    virtual IRenderOptions & THEA_ICALL setUseVertexData(int8 value) = 0;

    /** Draw mesh faces? */
    virtual int8 THEA_ICALL drawFaces() const = 0;

    /** Set whether mesh faces will be drawn. */
    virtual IRenderOptions & THEA_ICALL setDrawFaces(int8 value) = 0;

    /** Draw mesh edges? */
    virtual int8 THEA_ICALL drawEdges() const = 0;

    /** Set whether mesh edges will be drawn. */
    virtual IRenderOptions & THEA_ICALL setDrawEdges(int8 value) = 0;

    /** Override edge-specific colors with the value of edgeColor() when drawing edges? */
    virtual int8 THEA_ICALL overrideEdgeColor() const = 0;

    /** Set whether edge-specific colors will be overridden with the value of edgeColor() when drawing edges. */
    virtual IRenderOptions & THEA_ICALL setOverrideEdgeColor(int8 value) = 0;

    /** Color for drawing edges when overrideEdgeColor() is true, stored as an RGBA quadruplet. */
    virtual Real const * THEA_ICALL edgeColor() const = 0;

    /** Set the color for drawing edges when overrideEdgeColor() is true, passed as an RGBA quadruplet. */
    virtual IRenderOptions & THEA_ICALL setEdgeColor(Real const * rgba) = 0;

}; // class IRenderOptions

// Pure virtual destructor should have a body
// http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
inline IRenderOptions::~IRenderOptions() {}

/** %Options controlling the display of a IDrawable. */
class THEA_API RenderOptions : public virtual IRenderOptions
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
    ColorRgba  edge_color;

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
      edge_color(ColorRgb::white())
    {}

    /** Copy from base class object. */
    RenderOptions(IRenderOptions const & rhs) { *this = rhs; }

    /** Assign from base class object. */
    RenderOptions & operator=(IRenderOptions const & rhs)
    {
      send_normals         =  rhs.sendNormals();
      send_colors          =  rhs.sendColors();
      send_texcoords       =  rhs.sendTexCoords();
      use_vertex_normals   =  rhs.useVertexNormals();
      use_vertex_data      =  rhs.useVertexData();
      draw_faces           =  rhs.drawFaces();
      draw_edges           =  rhs.drawEdges();
      override_edge_color  =  rhs.overrideEdgeColor();
      edge_color           =  ColorRgba(rhs.edgeColor());

      return *this;
    }

    /** Get the default set of options. */
    static RenderOptions const * defaults() { static RenderOptions def; return &def; }

    int8 THEA_ICALL sendColors() const { return send_colors; }
    IRenderOptions & THEA_ICALL setSendColors(int8 value) { send_colors = value; return *this; }
    int8 THEA_ICALL sendNormals() const { return send_normals; }
    IRenderOptions & THEA_ICALL setSendNormals(int8 value) { send_normals = value; return *this; }
    int8 THEA_ICALL sendTexCoords() const { return send_texcoords; }
    IRenderOptions & THEA_ICALL setSendTexCoords(int8 value) { send_texcoords = value; return *this; }
    int8 THEA_ICALL useVertexNormals() const { return use_vertex_normals; }
    IRenderOptions & THEA_ICALL setUseVertexNormals(int8 value) { use_vertex_normals = value; return *this; }
    int8 THEA_ICALL useVertexData() const { return use_vertex_data; }
    IRenderOptions & THEA_ICALL setUseVertexData(int8 value) { use_vertex_data = value; return *this; }
    int8 THEA_ICALL drawFaces() const { return draw_faces; }
    IRenderOptions & THEA_ICALL setDrawFaces(int8 value) { draw_faces = value; return *this; }
    int8 THEA_ICALL drawEdges() const { return draw_edges; }
    IRenderOptions & THEA_ICALL setDrawEdges(int8 value) { draw_edges = value; return *this; }
    int8 THEA_ICALL overrideEdgeColor() const { return override_edge_color; }
    IRenderOptions & THEA_ICALL setOverrideEdgeColor(int8 value) { override_edge_color = value; return *this; }
    Real const * THEA_ICALL edgeColor() const { return edge_color.data(); }
    IRenderOptions & THEA_ICALL setEdgeColor(Real const * rgba) { edge_color = ColorRgba(rgba); return *this; }

}; // class RenderOptions

} // namespace Graphics
} // namespace Thea

#endif
