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

#ifndef __Browse3D_ModelDisplay_hpp__
#define __Browse3D_ModelDisplay_hpp__

#include "Common.hpp"
#include "../../AffineTransform3.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/RenderOptions.hpp"
#include "../../Plugins/GL/GLHeaders.hpp"
#include <wx/glcanvas.h>

namespace Thea {
namespace Graphics {

class RenderSystem;
class Texture;
class Shader;

} // namespace Graphics
} // namespace Thea

namespace Browse3D {

class Model;

/** An OpenGL widget to display and interact with a model. */
class ModelDisplay : public wxGLCanvas
{
  private:
    typedef Graphics::Camera Camera;
    typedef Graphics::RenderOptions RenderOptions;

    /** The interaction mode the widget is in (enum class). */
    struct Mode
    {
      enum Value
      {
        DEFAULT,    ///< Default state, no interaction.
        EDIT_VIEW,  ///< The view is being changed.
        EDIT_MODEL  ///< The model is being changed.
      };

      THEA_ENUM_CLASS_BODY(Mode)
    };

    /** How the view is being modified (enum class). */
    struct ViewEditMode
    {
      enum Value
      {
        DEFAULT,      ///< Default state, no modification.
        PAN,          ///< The view is being panned.
        ROTATE,       ///< The view is being freely rotated (inverse of mouse-look).
        ROTATE_ROLL,  ///< The view is being rolled around the viewing direction.
        ZOOM          ///< The view is being zoomed in or out.
      };

      THEA_ENUM_CLASS_BODY(ViewEditMode)
    };

  public:
    /** Constructs the widget as a viewer for a given model. */
    explicit ModelDisplay(wxWindow * parent, Model * model_);

    /** Destructor. */
    ~ModelDisplay();

    /** Set the reference to the model, which may be null. */
    void setModel(Model * model_);

    /** Get the viewing camera. */
    Graphics::Camera const & getCamera() const { return camera; }

    /** Get the point the camera is looking at. */
    Vector3 const & getCameraLookAt() const { return camera_look_at; }

    /**
     * Compute a ray that emanates from the viewer's eye and passes through a given pixel on the widget. The (roughly) inverse
     * operation is project().
     */
    Ray3 computePickRay(wxRealPoint const & p) const;

    /** Project a 3D point to the viewing plane. The (roughly) inverse operation is computePickRay(). */
    wxRealPoint project(Vector3 const & p) const;

    /** Check if flat shading is on/off. */
    bool flatShading() const { return !render_opts.useVertexNormals(); }

    /** Check if two-sided lighting is on/off. */
    bool twoSided() const;

    /** Set two-sided lighting on/off. */
    void setTwoSided(bool value);

    /** Set flat shading on/off. */
    void setFlatShading(bool value);

    /** Save a screenshot to a file. If the path is empty, a default path is used. */
    void saveScreenshot(std::string path = "") const;

    //=========================================================================================================================
    // GUI callbacks
    //=========================================================================================================================

    /** Adjust the view to fit the current model. */
    void fitViewToModel(wxEvent & event = DUMMY_EVENT);

    /** Called when the model geometry changes. */
    void modelGeometryChanged(wxEvent & event = DUMMY_EVENT);

    /** Called when the model needs to be redrawn. */
    void modelNeedsRedraw(wxEvent & event = DUMMY_EVENT);

    /** Render shaded polygons, without edges. */
    void renderShaded(wxEvent & event = DUMMY_EVENT);

    /** Render only polygon edges. */
    void renderWireframe(wxEvent & event = DUMMY_EVENT);

    /** Render shaded polygons, with edges in a different color. */
    void renderShadedWireframe(wxEvent & event = DUMMY_EVENT);

    /** Set two-sided lighting on/off. */
    void setTwoSided(wxCommandEvent & event);

    /** Set flat shading on/off. */
    void setFlatShading(wxCommandEvent & event);

    /** Save a screenshot to a file.. */
    void saveScreenshot(wxEvent & event = DUMMY_EVENT);

    /** Called to render display. */
    void paintGL(wxPaintEvent & event);

    /** Called when widget is resized. */
    void resize(wxSizeEvent & event);

    /** Called when a key is pressed. */
    void keyPressEvent(wxKeyEvent & event);

    /** Called when a mouse button is pressed. */
    void mousePressEvent(wxMouseEvent & event);

    /** Called when the mouse is moved. */
    void mouseMoveEvent(wxMouseEvent & event);

    /** Called when a mouse button is released. */
    void mouseReleaseEvent(wxMouseEvent & event);

    /** Called when the mouse wheel is turned. */
    void wheelEvent(wxMouseEvent & event);

  private:
    /** Update the viewing camera to fit the current model. */
    void updateCameraFromModel();

    /** Update the camera view frame to fit the current model. */
    void updateCameraFrameFromModel();

    /** Update the camera projection to fit the current model. */
    void updateCameraProjection();

    /** Get the distance to the center of the model. */
    Real getModelDistance() const;

    /** Drag the view window with the mouse. */
    void panView(wxMouseEvent & event);

    /** Rotate the view with the mouse. */
    void rotateView(wxMouseEvent & event);

    /** Roll the view with the mouse around the viewing direction. */
    void rollView(wxMouseEvent & event);

    /** Zoom the view in or out by dragging with the mouse. */
    void zoomView(wxMouseEvent & event);

    /** Zoom the view in or out by rotating the mouse wheel. */
    void zoomViewWheel(wxMouseEvent & event);

    /** Multiply the current view transform by an increment. */
    void incrementViewTransform(AffineTransform3 const & tr);

    /** Draw a corner icon showing the coordinate axes. */
    void drawAxes(Graphics::RenderSystem & rs);

    /** Draw the background image. */
    void drawBackground(Graphics::RenderSystem & rs);

    /** Get the width of the widget. */
    int width() const { return GetSize().GetWidth(); }

    /** Get the height of the widget. */
    int height() const { return GetSize().GetHeight(); }

    Model * model;

    Camera camera;
    Vector3 camera_look_at;
    RenderOptions render_opts;

    Mode mode;
    ViewEditMode view_edit_mode;
    wxPoint view_drag_start;
    wxPoint last_cursor;

    Graphics::Texture * background_texture;
    Graphics::Shader * background_shader;

    wxGLContext * context;

}; // class ModelDisplay

} // namespace Browse3D

#endif
