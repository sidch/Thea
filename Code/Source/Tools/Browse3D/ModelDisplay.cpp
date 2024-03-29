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

#include "ModelDisplay.hpp"
#include "App.hpp"
#include "GraphicsWidget.hpp"
#include "MainWindow.hpp"
#include "Math.hpp"
#include "Model.hpp"
#include "Util.hpp"
#include "../../Application.hpp"
#include "../../Colors.hpp"
#include "../../Image.hpp"
#include "../../Ray3.hpp"
#include "../../Graphics/IRenderSystem.hpp"
#include "../../Graphics/IShader.hpp"
#include "../../Graphics/ITexture.hpp"
#include <wx/datetime.h>
#include <wx/dcclient.h>
#include <wx/image.h>
#include <wx/stdpaths.h>

namespace Browse3D {

#define NEW_WXGLCANVAS (wxMAJOR_VERSION > 3 || (wxMAJOR_VERSION == 3 && wxMINOR_VERSION >= 1))

#if NEW_WXGLCANVAS
wxGLAttributes
getGlAttributes()
{
  wxGLAttributes attribs;
  attribs.PlatformDefaults().RGBA().Depth(16).SampleBuffers(1).Samplers(4).DoubleBuffer().EndList();
  return attribs;
}
#else
int const *
getGlAttributes()
{
  static int attribs[32];
  int n = 0;
  attribs[n++] = WX_GL_RGBA;
  attribs[n++] = WX_GL_DEPTH_SIZE;
  attribs[n++] = 16;
  attribs[n++] = WX_GL_SAMPLE_BUFFERS;
  attribs[n++] = 1;
  attribs[n++] = WX_GL_SAMPLES;
  attribs[n++] = 4;
  attribs[n++] = WX_GL_DOUBLEBUFFER;
  attribs[n] = 0; // terminate the list

  return (int const *)attribs;
}
#endif

ModelDisplay::ModelDisplay(wxWindow * parent, Model * model_)
#if NEW_WXGLCANVAS
: wxGLCanvas(parent, getGlAttributes()),
#else
: wxGLCanvas(parent, wxID_ANY, getGlAttributes()),
#endif
  model(model_),
  camera(CoordinateFrame3::identity(),
         Camera::ProjectionType::PERSPECTIVE, -1, 1, -1, 1, 1, 3, Camera::ProjectedYDirection::UP),
  camera_look_at(0, 0, -1),
  render_opts(*RenderOptions::defaults()),
  mode(Mode::DEFAULT),
  view_edit_mode(ViewEditMode::DEFAULT),
  background_texture(nullptr),
  background_shader(nullptr),
  context(nullptr)
{
  alwaysAssertM(model, "ModelDisplay: Can't create a display without a valid model");

//   setFocusPolicy(Qt::StrongFocus);
//   setMouseTracking(true);

  context = new wxGLContext(this);

  static ColorRgba const DEFAULT_EDGE_COLOR(0.15f, 0.25f, 0.5f, 1.0f);
  render_opts.setSendNormals(true)
             .setSendColors(false)
             .setSendTexCoords(false)
             .setDrawFaces(true)
             .setDrawEdges(true)
             .setOverrideEdgeColor(true)
             .setEdgeColor(DEFAULT_EDGE_COLOR.data());

  Bind(wxEVT_PAINT, &ModelDisplay::paintGL, this);
  Bind(wxEVT_SIZE, &ModelDisplay::resize, this);

  Bind(wxEVT_LEFT_DOWN, &ModelDisplay::mousePressEvent, this);
  Bind(wxEVT_RIGHT_DOWN, &ModelDisplay::mousePressEvent, this);
  Bind(wxEVT_MIDDLE_DOWN, &ModelDisplay::mousePressEvent, this);

  Bind(wxEVT_LEFT_UP, &ModelDisplay::mouseReleaseEvent, this);
  Bind(wxEVT_RIGHT_UP, &ModelDisplay::mouseReleaseEvent, this);
  Bind(wxEVT_MIDDLE_UP, &ModelDisplay::mouseReleaseEvent, this);

  Bind(wxEVT_MOTION, &ModelDisplay::mouseMoveEvent, this);
  Bind(wxEVT_MOUSEWHEEL, &ModelDisplay::wheelEvent, this);

  model->registerDisplay(this);
}

ModelDisplay::~ModelDisplay()
{
  if (model)
    model->deregisterDisplay(this);
}

void
ModelDisplay::setModel(Model * model_)
{
  if (model)
    model->deregisterDisplay(this);

  model = model_;

  if (model)
    model->registerDisplay(this);
}

Ray3
ModelDisplay::computePickRay(wxRealPoint const & p) const
{
  return Browse3D::computePickRay(p, camera, width(), height());
}

wxRealPoint
ModelDisplay::project(Vector3 const & p) const
{
  Vector3 proj_p = camera.project(p);
  return wxRealPoint(0.5 * (proj_p.x() + 1) * width(), 0.5 * (proj_p.y() + 1) * height());
}

void
ModelDisplay::fitViewToModel(wxEvent & event)
{
  updateCameraFromModel();
  Refresh();
}

void
ModelDisplay::updateCameraFromModel()
{
  if (!model || model->empty()) return;

  model->updateBounds();
  updateCameraFrameFromModel();
  updateCameraProjection();
}

void
ModelDisplay::updateCameraFrameFromModel()
{
  if (!model || model->empty())
  {
    camera.setFrame(CoordinateFrame3::identity());
    return;
  }

  camera_look_at = model->getTransformedBounds().getCenter();
  Real model_scale = model->getTransformedBounds().getExtent().norm();
  Real camera_separation = model_scale > 1.0e-3f ? 2 * model_scale : 1.0e-3f;

  // Maintain current orientation
  Vector3 dir = camera.getLookDirection();
  Vector3 up  = camera.getUpDirection();
  camera.setFrame(CoordinateFrame3::fromViewFrame(camera_look_at - camera_separation * dir,  // eye
                                                  camera_look_at,                            // look-at
                                                  up));                                      // up
}

void
ModelDisplay::updateCameraProjection()
{
  if (!model || model->empty())
  {
    camera.setProjection(Camera::ProjectionType::PERSPECTIVE, -1, 1, -1, 1, 1, 3, Camera::ProjectedYDirection::UP);
    return;
  }

  static Real const HALF_WIDTH = 0.5;
  int w = width(), h = height();
  Real hw = 0, hh = 0;
  if (h > w)
  {
    Real aspect_ratio = h / (Real)w;
    hw = HALF_WIDTH;
    hh = aspect_ratio * HALF_WIDTH;
  }
  else
  {
    Real aspect_ratio = w / (Real)h;
    hw = aspect_ratio * HALF_WIDTH;
    hh = HALF_WIDTH;
  }

  Vector3 center = model->getTransformedBounds().getCenter();
  Vector3 dir = camera.getLookDirection();
  Real center_dist = std::max((center - camera.getPosition()).dot(dir), (Real)0);
  Real scale = model->getTransformedBounds().getExtent().norm();

  Real near_dist = std::max(center_dist - 1.1f * scale, 0.01f * scale);
  Real far_dist  = center_dist + 2 * scale;

  hw *= (0.5f * near_dist);
  hh *= (0.5f * near_dist);

  camera.setProjection(Camera::ProjectionType::PERSPECTIVE, -hw, hw, -hh, hh, near_dist, far_dist,
                       Camera::ProjectedYDirection::UP);
}

void
ModelDisplay::printCamera(wxEvent & event)
{
  THEA_CONSOLE << "Camera is: " << camera.toString();
  THEA_CONSOLE << "Viewing matrix (world to camera) is: " << toString(camera.getWorldToCameraTransform().homogeneous());
  THEA_CONSOLE << "Projection matrix (camera to projection) is: " << toString(camera.getProjectionTransform());
}

void
ModelDisplay::modelGeometryChanged(wxEvent & event)
{
  fitViewToModel(event);
}

void
ModelDisplay::modelNeedsRedraw(wxEvent & event)
{
  Refresh();
}

void
ModelDisplay::renderShaded(wxEvent & event)
{
  render_opts.setDrawFaces(true).setDrawEdges(false);
  Refresh();

  THEA_CONSOLE << "Rendering shaded faces";
}

void
ModelDisplay::renderWireframe(wxEvent & event)
{
  render_opts.setDrawFaces(false).setDrawEdges(true);
  Refresh();

  THEA_CONSOLE << "Rendering wireframe";
}

void
ModelDisplay::renderShadedWireframe(wxEvent & event)
{
  render_opts.setDrawFaces(true).setDrawEdges(true);
  Refresh();

  THEA_CONSOLE << "Rendering shaded faces with wireframe edges";
}

bool
ModelDisplay::twoSided() const
{
  return GraphicsWidget::isTwoSided();
}

void
ModelDisplay::setTwoSided(bool value)
{
  if (GraphicsWidget::isTwoSided() != value)
  {
    GraphicsWidget::setTwoSided(value);
    Refresh();

    THEA_CONSOLE << "Two-sided lighting = " << value;
  }
}

void
ModelDisplay::setTwoSided(wxCommandEvent & event)
{
  setTwoSided(event.IsChecked());
}

bool
ModelDisplay::flatShaded() const
{
  return GraphicsWidget::isFlatShaded();
}

void
ModelDisplay::setFlatShaded(bool value)
{
  if (GraphicsWidget::isFlatShaded() != value)
  {
    GraphicsWidget::setFlatShaded(value);
    Refresh();

    THEA_CONSOLE << "Flat shading = " << value;
  }
}

void
ModelDisplay::setFlatShaded(wxCommandEvent & event)
{
  setFlatShaded(event.IsChecked());
}

void
ModelDisplay::paintGL(wxPaintEvent & event)
{
  using namespace Graphics;

  if (!model) return;
  if (!app().hasRenderSystem()) return;

  SetCurrent(*context);
  wxPaintDC(this);

#ifndef THEA_WINDOWS
  glEnable(GL_MULTISAMPLE);
#endif

  IRenderSystem & rs = *app().getRenderSystem();

  // wxGLCanvas uses physical pixels, unlike the rest of the framework which mostly uses logical pixels
  wxSize physical_size = GetClientSize() * GetContentScaleFactor();
  rs.setViewport(0, 0, physical_size.GetWidth(), physical_size.GetHeight());

  rs.setClearColor(app().options().bg_color.data());
  rs.clear();

  rs.setColorWrite(true, true, true, true);
  rs.setDepthWrite(true);
  rs.setDepthTest(IRenderSystem::DepthTest::LESS);
  camera.makeCurrent(&rs);

  if (!app().options().bg_plain)
    drawBackground(rs);

  model->draw(&rs, &render_opts);

  intx num_overlays = app().getMainWindow()->numOverlays();
  Model const * const * overlays = app().getMainWindow()->getOverlays();
  for (intx i = 0; i < num_overlays; ++i)
    overlays[i]->draw(&rs, &render_opts);

  if (!app().options().no_axes)
    drawAxes(rs);

  glFlush();
  SwapBuffers();
}

void
ModelDisplay::resize(wxSizeEvent & event)
{
  fitViewToModel();
}

void
ModelDisplay::drawBackground(Graphics::IRenderSystem & rs)
{
  // THEA_CONSOLE << "drawBackground";

  using namespace Graphics;

  if (!background_texture)
  {
    static std::string const BG_IMAGE = "background.jpg";
    Image bg_img;
    if (loadImage(bg_img, Application::getResourcePath("Images/" + BG_IMAGE)))
    {
      TextureOptions opts;
      opts.setInterpolateMode(ITextureOptions::InterpolateMode::BILINEAR_NO_MIPMAP);

      background_texture = rs.createTexture("Background texture", &bg_img, TextureFormat::RGBA8(), ITexture::Dimension::DIM_2D,
                                            &opts);
    }

    if (background_texture)
     THEA_CONSOLE << "Loaded " << background_texture->getWidth() << "x" << background_texture->getHeight()
                  << " background image " << BG_IMAGE;
    else
    {
      THEA_CONSOLE << "Couldn't load background texture";
      return;  // just use whatever we had earlier as a fallback
    }
  }

  if (!background_shader)
  {
    background_shader = rs.createShader("Background shader");

    if (!background_shader->attachModuleFromFile(IShader::ModuleType::VERTEX,
                                                 Application::getResourcePath("Materials/FlatTextureVert.glsl").c_str())
     || !background_shader->attachModuleFromFile(IShader::ModuleType::FRAGMENT,
                                                 Application::getResourcePath("Materials/FlatTextureFrag.glsl").c_str()))
      return;

    background_shader->setUniform("texture", background_texture);
  }

  rs.pushShader();
    rs.setShader(background_shader);

    // FIXME: For some strange reason calling glActiveTextureARB (as done by IRenderSystem::sendTexCoord) within the
    // glBegin/glEnd block below triggers a GL invalid operation error. Hence we'll just use the plain vanilla GL calls instead
    // of the IRenderSystem wrappers. (Update: this might have been fixed by the recent fix to GlTexture.)
    glPushAttrib(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
      glDepthFunc(GL_ALWAYS);
      glDepthMask(GL_FALSE);
      glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

      glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();
      glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();

        glBegin(GL_QUADS);
           glTexCoord2f(0, 0); glVertex3f(-1, -1, 0.5f);
           glTexCoord2f(1, 0); glVertex3f( 1, -1, 0.5f);
           glTexCoord2f(1, 1); glVertex3f( 1,  1, 0.5f);
           glTexCoord2f(0, 1); glVertex3f(-1,  1, 0.5f);
        glEnd();

      glMatrixMode(GL_PROJECTION); glPopMatrix();
      glMatrixMode(GL_MODELVIEW); glPopMatrix();
    glPopAttrib();

  rs.popShader();
}

void
ModelDisplay::drawAxes(Graphics::IRenderSystem & rs)
{
  // THEA_CONSOLE << "drawAxes";

  // TODO: Make the axes and labels occlude each other properly depending on orientation.

  using namespace Graphics;

  Matrix3 rot = camera.getFrame().getRotation();

  Vector2 x_axis = rot.row(0).head<2>();
  Vector2 y_axis = rot.row(1).head<2>();
  Vector2 z_axis = rot.row(2).head<2>();

  Vector2 aspect_ratio_compensation = width() < height() ? Vector2(1.0, width() / (Real)height())
                                                         : Vector2(height() / (Real)width(), 1.0);

  static Real const ARROW_LENGTH = 0.1f;
  Vector2 arrow_origin = Vector2(-1, -1) + (2.0f * ARROW_LENGTH) * aspect_ratio_compensation;
  Vector2 x_arrow_tip = arrow_origin + ARROW_LENGTH * aspect_ratio_compensation.cwiseProduct(x_axis);
  Vector2 y_arrow_tip = arrow_origin + ARROW_LENGTH * aspect_ratio_compensation.cwiseProduct(y_axis);
  Vector2 z_arrow_tip = arrow_origin + ARROW_LENGTH * aspect_ratio_compensation.cwiseProduct(z_axis);

  rs.pushColorFlags();
  rs.pushDepthFlags();
  rs.setDepthTest(IRenderSystem::DepthTest::ALWAYS_PASS);

  rs.pushShader();
  rs.setShader(nullptr);

    rs.setMatrixMode(IRenderSystem::MatrixMode::PROJECTION); rs.pushMatrix(); rs.setIdentityMatrix();
    rs.setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); rs.pushMatrix(); rs.setIdentityMatrix();

      rs.beginPrimitive(IRenderSystem::Primitive::LINES);

        rs.setColor(ColorRgba::red().data());
        rs.sendVertex(arrow_origin[0], arrow_origin[1]);
        rs.sendVertex(x_arrow_tip[0],  x_arrow_tip[1]);

        rs.setColor(ColorRgba::green().data());
        rs.sendVertex(arrow_origin[0], arrow_origin[1]);
        rs.sendVertex(y_arrow_tip[0],  y_arrow_tip[1]);

        rs.setColor(ColorRgba::blue().data());
        rs.sendVertex(arrow_origin[0], arrow_origin[1]);
        rs.sendVertex(z_arrow_tip[0],  z_arrow_tip[1]);

      rs.endPrimitive();

    rs.setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); rs.popMatrix();
    rs.setMatrixMode(IRenderSystem::MatrixMode::PROJECTION); rs.popMatrix();

// #define DRAW_AXIS_LABELS
#ifdef DRAW_AXIS_LABELS

  // For some reason drawing text disables this and doesn't restore it properly
#ifndef THEA_WINDOWS
  bool multisample = (glIsEnabled(GL_MULTISAMPLE) == GL_TRUE);
#endif

//   FIXME
//   static QFont label_font("Helvetica", 12);
//   static QFontMetrics fm(label_font);
//   static QSize half_x_label_size = 0.5 * fm.tightBoundingRect("X").size();
//   static QSize half_y_label_size = 0.5 * fm.tightBoundingRect("Y").size();
//   static QSize half_z_label_size = 0.5 * fm.tightBoundingRect("Z").size();
//
//   static Real const TEXT_OFFSET = 0.5;
//   Vector2 x_text_center = x_arrow_tip + TEXT_OFFSET * (x_arrow_tip - arrow_origin);
//   Vector2 y_text_center = y_arrow_tip + TEXT_OFFSET * (y_arrow_tip - arrow_origin);
//   Vector2 z_text_center = z_arrow_tip + TEXT_OFFSET * (z_arrow_tip - arrow_origin);
//
//   wxPoint x_text_bottom_left((int)std::round((0.5 * (1 + x_text_center.x())) * width()) - half_x_label_size.width(),
//                              (int)std::round((0.5 * (1 - x_text_center.y())) * height()) + half_x_label_size.height());
//
//   wxPoint y_text_bottom_left((int)std::round((0.5 * (1 + y_text_center.x())) * width()) - half_y_label_size.width(),
//                              (int)std::round((0.5 * (1 - y_text_center.y())) * height()) + half_y_label_size.height());
//
//   wxPoint z_text_bottom_left((int)std::round((0.5 * (1 + z_text_center.x())) * width()) - half_z_label_size.width(),
//                              (int)std::round((0.5 * (1 - z_text_center.y())) * height()) + half_z_label_size.height());
//
//   rs.setColor(ColorRgb::red());
//   renderText(x_text_bottom_left.x(), x_text_bottom_left.y(), "X", label_font);
//
//   rs.setColor(ColorRgb::green());
//   renderText(y_text_bottom_left.x(), y_text_bottom_left.y(), "Y", label_font);
//
//   rs.setColor(ColorRgb::blue());
//   renderText(z_text_bottom_left.x(), z_text_bottom_left.y(), "Z", label_font);

#ifndef THEA_WINDOWS
  if (multisample) glEnable(GL_MULTISAMPLE);  // for some reason drawing text disables this and doesn't restore it properly
#endif

#endif

  rs.popShader();
  rs.popDepthFlags();
  rs.popColorFlags();
}

void
ModelDisplay::keyPressEvent(wxKeyEvent & event)
{
}

void
ModelDisplay::mousePressEvent(wxMouseEvent & event)
{
  if (!model) return;

  last_cursor = event.GetPosition();

  bool left    =  event.LeftDown();
  bool right   =  event.RightDown();
  bool middle  =  event.MiddleDown();

  bool shift  =  event.ShiftDown();
  bool ctrl   =  event.ControlDown();
  bool alt    =  event.AltDown();
  bool no_mods     =  !(shift || ctrl || alt);
  bool shift_only  =   shift && !ctrl && !alt;
  bool ctrl_only   =  !shift &&  ctrl && !alt;
  bool alt_only    =  !shift && !ctrl &&  alt;

  bool pick_points = app().getMainWindow()->pickPoints();
  bool pick_segments = app().getMainWindow()->pickSegments();
  bool edit_model = (pick_points || pick_segments);

  if (edit_model && (right || (ctrl && left)))
  {
    Ray3 ray = computePickRay(event.GetPosition());

    bool success = false;
    if (pick_points)
      success = (model->pick(ray) >= 0);
    else if (pick_segments)
      success = (model->togglePickMesh(ray, shift) >= 0);

    if (success)
    {
      model->mousePressEvent(event);
      if (!event.ShouldPropagate())
      {
        mode = Mode::EDIT_MODEL;
        return;
      }
    }
  }

  bool pan_view     =  (no_mods && right) || middle;
  bool rotate_view  =  (no_mods && left) || (alt_only && left);
  bool roll_view    =  (ctrl_only && left);
  bool zoom_view    =  (shift_only && left);

  if (pan_view || rotate_view || roll_view || zoom_view)
  {
    mode = Mode::EDIT_VIEW;

    if (pan_view)
      view_edit_mode = ViewEditMode::PAN;
    else if (rotate_view)
      view_edit_mode = ViewEditMode::ROTATE;
    else if (roll_view)
      view_edit_mode = ViewEditMode::ROTATE_ROLL;
    else  // zoom_view
      view_edit_mode = ViewEditMode::ZOOM;

    view_drag_start = event.GetPosition();

    event.StopPropagation();
    Refresh();
    return;
  }
}

void
ModelDisplay::mouseMoveEvent(wxMouseEvent & event)
{
  // THEA_CONSOLE << "mouseMoveEvent";

  if (!model) return;

  if (mode == Mode::EDIT_VIEW)
  {
    switch (view_edit_mode)
    {
      case ViewEditMode::PAN:
        panView(event);
        event.StopPropagation();
        break;

      case ViewEditMode::ROTATE:
        rotateView(event);
        event.StopPropagation();
        break;

      case ViewEditMode::ROTATE_ROLL:
        rollView(event);
        event.StopPropagation();
        break;

      case ViewEditMode::ZOOM:
        zoomView(event);
        event.StopPropagation();
        break;

      default: break;
    }
  }
  else  // register all mouse movement over the model, not just in EDIT_MODEL mode
  {
    model->mouseMoveEvent(event);
  }

  last_cursor = event.GetPosition();
}

void
ModelDisplay::incrementViewTransform(AffineTransform3 const & tr)
{
  // THEA_CONSOLE << "incrementViewTransform";

  AffineTransform3 inv_vt = tr.inverse();
  CoordinateFrame3 const & old_cframe = camera.getFrame();
  camera.setFrame(CoordinateFrame3(RigidTransform3::_fromAffine(AffineTransform3(inv_vt.getLinear() * old_cframe.getRotation(),
                                                                                 inv_vt * old_cframe.getTranslation()))));

  if (view_edit_mode == ViewEditMode::PAN)
    camera_look_at = inv_vt * camera_look_at;

  updateCameraProjection();  // needed because it adjusts to the new separation from the model

  Refresh();
}

void
ModelDisplay::mouseReleaseEvent(wxMouseEvent & event)
{
  // THEA_CONSOLE << "mouseReleaseEvent";

  if (!model) return;

  if (mode == Mode::EDIT_VIEW)
  {
    if (view_edit_mode != ViewEditMode::DEFAULT)
    {
      view_edit_mode = ViewEditMode::DEFAULT;
      mode = Mode::DEFAULT;

      event.StopPropagation();
      Refresh();
    }
  }
  else if (mode == Mode::EDIT_MODEL)
  {
    model->mouseReleaseEvent(event);
  }

  last_cursor = event.GetPosition();
}

void
ModelDisplay::wheelEvent(wxMouseEvent & event)
{
  // THEA_CONSOLE << "wheelEvent";

  mode = Mode::EDIT_VIEW;
  view_edit_mode = ViewEditMode::ZOOM;

    zoomViewWheel(event);

  mode = Mode::DEFAULT;
  view_edit_mode = ViewEditMode::DEFAULT;
}

Real
ModelDisplay::getModelDistance() const
{
  // THEA_CONSOLE << "getModelDistance";

  if (!model) return 0.0;

  return model ? (model->getTransformedBounds().getCenter() - camera.getPosition()).norm()
               : (camera_look_at - camera.getPosition()).norm();
}

void
ModelDisplay::panView(wxMouseEvent & event)
{
  // THEA_CONSOLE << "panView";

  if (!model) return;

  Vector3 trn = dragToTranslation(last_cursor, event.GetPosition(), width(), height(), camera, getModelDistance());
  incrementViewTransform(AffineTransform3(Matrix3::Identity(), trn));
}

void
ModelDisplay::rotateView(wxMouseEvent & event)
{
  // THEA_CONSOLE << "rotateView";

  if (!model) return;

  Matrix3 rot = dragToRotation(last_cursor, event.GetPosition(), width(), height(), camera);
  Vector3 trn = camera_look_at - rot * camera_look_at;
  incrementViewTransform(AffineTransform3(rot, trn));
}

void
ModelDisplay::rollView(wxMouseEvent & event)
{
  // THEA_CONSOLE << "rollView";

  if (!model) return;

  Vector2 center(0.5f * width(), 0.5f * height());
  Vector2 u = Vector2(last_cursor.x, last_cursor.y) - center;
  Vector2 v = Vector2(event.GetPosition().x, event.GetPosition().y) - center;
  Real s = u.x() * v.y() - u.y() * v.x();
  Real c = u.dot(v);
  Real angle = Math::fastArcTan2(s, c);

  Matrix3 rot = Math::rotationAxisAngle(camera.getLookDirection(), angle);
  Vector3 trn = camera_look_at - rot * camera_look_at;
  incrementViewTransform(AffineTransform3(rot, trn));
}

namespace ModelDisplayInternal {

AffineTransform3
zoomTransform(Real zoom, Real camera_separation, Vector3 const & look_dir)
{
  static Real const MIN_ZOOM = 0.25;
  Vector3 trn = (1.0f / std::max(zoom, MIN_ZOOM) - 1) * camera_separation * look_dir;
  return AffineTransform3(Matrix3::Identity(), trn);
}

} // namespace ModelDisplayInternal

void
ModelDisplay::zoomView(wxMouseEvent & event)
{
  // THEA_CONSOLE << "zoomView";

  using namespace ModelDisplayInternal;

  if (!model) return;

  Real zoom = dragToScale(last_cursor, event.GetPosition(), width(), height(), camera);

  // Zoom in at mouse position, zoom out at screen center
  Vector3 dir;
  if (zoom > 1)
    dir = computePickRay(view_drag_start).getDirection().normalized();
  else
    dir = camera.getLookDirection();

  Real camera_separation = (camera_look_at - camera.getPosition()).norm();
  incrementViewTransform(zoomTransform(zoom, camera_separation, dir));
}

void
ModelDisplay::zoomViewWheel(wxMouseEvent & event)
{
  // THEA_CONSOLE << "zoomView";

  using namespace ModelDisplayInternal;

  if (!model) return;

  static Real const ZOOM_PER_DEGREE = 2.0f / 180.0f;
  Real num_degrees = event.GetWheelRotation() / 8;  // delta is in eighths of a degree
  Real zoom = 1 + ZOOM_PER_DEGREE * num_degrees;

  // Zoom in at mouse position, zoom out at screen center
  Vector3 dir;
  if (zoom > 1)
    dir = computePickRay(event.GetPosition()).getDirection().normalized();
  else
    dir = camera.getLookDirection();

  Real camera_separation = (camera_look_at - camera.getPosition()).norm();
  incrementViewTransform(zoomTransform(zoom, camera_separation, dir));
}

void
ModelDisplay::saveScreenshot(std::string path) const
{
  if (!model) return;

  if (path.empty())
  {
    Model const * model = app().getMainWindow()->getModel();
    std::string prefix = model ? FilePath::baseName(model->getPath()) : "Browse3D-Screenshot";
    path = FilePath::concat(wxStandardPaths::Get().GetDocumentsDir().ToStdString(),
                            prefix + wxDateTime::Now().Format("-%Y-%m-%d-%H-%M-%S").ToStdString() + ".png");
  }

  // Following method adapted from
  // http://www.organicvectory.com/index.php?option=com_content&view=article&id=73:save-a-screenshot

  SetCurrent(*context);

  // Retrieve size of frame buffer
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  // Allocate memory for the pixel data
  int w = viewport[2];
  int h = viewport[3];
  wxImage img(w, h, false);

  // Are we in single or double-buffered mode?
  GLboolean double_buffer;
  glGetBooleanv(GL_DOUBLEBUFFER, &double_buffer);

  // Retrieve framebuffer pixels
  glPushAttrib(GL_PIXEL_MODE_BIT);
  glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(double_buffer == GL_TRUE ? GL_BACK : GL_FRONT);  // this is the default, but be explicit in case it was
                                                                  // modified elsewhere
    glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3], GL_RGB, GL_UNSIGNED_BYTE, img.GetData());
  glPopClientAttrib();
  glPopAttrib();

  // Image is flipped vertically
  img = img.Mirror(false);

  // Save image
  if (img.SaveFile(path, wxBITMAP_TYPE_PNG))
    THEA_CONSOLE << "Saved screenshot to " << path;
  else
    THEA_ERROR << "Could not save screenshot to " << path;
}

void
ModelDisplay::saveScreenshot(wxEvent & event)
{
  // This function is a callback, so it can't be const else it doesn't match the Bind signature
  const_cast<ModelDisplay *>(this)->saveScreenshot("");
}

} // namespace Browse3D
