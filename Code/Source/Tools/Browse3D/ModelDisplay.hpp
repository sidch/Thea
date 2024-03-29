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

#ifndef __Browse3D_ModelDisplay_hpp__
#define __Browse3D_ModelDisplay_hpp__

#include "Common.hpp"
#include "../../AffineTransform3.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/RenderOptions.hpp"
#include "../../Plugins/GL/GlHeaders.hpp"
#include <wx/glcanvas.h>

namespace Thea {
namespace Graphics {

class IRenderSystem;
class ITexture;
class IShader;

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
    bool flatShaded() const;

    /** Check if two-sided lighting is on/off. */
    bool twoSided() const;

    /** Set two-sided lighting on/off. */
    void setTwoSided(bool value);

    /** Set flat shading on/off. */
    void setFlatShaded(bool value);

    /** Save a screenshot to a file. If the path is empty, a default path is used. */
    void saveScreenshot(std::string path = "") const;

    //=========================================================================================================================
    // GUI callbacks
    //=========================================================================================================================

    /** Adjust the view to fit the current model. */
    void fitViewToModel(wxEvent & event = DUMMY_EVENT);

    /** Print the camera parameters. */
    void printCamera(wxEvent & event = DUMMY_EVENT);

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
    void setFlatShaded(wxCommandEvent & event);

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
    void drawAxes(Graphics::IRenderSystem & rs);

    /** Draw the background image. */
    void drawBackground(Graphics::IRenderSystem & rs);

    /** Get the width of the drawable area of the widget. */
    int width() const { return GetClientSize().GetWidth(); }

    /** Get the height of the drawable area of the widget. */
    int height() const { return GetClientSize().GetHeight(); }

    Model * model;

    Camera camera;
    Vector3 camera_look_at;
    RenderOptions render_opts;

    Mode mode;
    ViewEditMode view_edit_mode;
    wxPoint view_drag_start;
    wxPoint last_cursor;

    Graphics::ITexture * background_texture;
    Graphics::IShader * background_shader;

    wxGLContext * context;

}; // class ModelDisplay

} // namespace Browse3D

#endif
