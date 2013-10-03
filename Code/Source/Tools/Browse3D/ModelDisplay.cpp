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

#include "ModelDisplay.hpp"
#include "App.hpp"
#include "MainWindow.hpp"
#include "Math.hpp"
#include "Model.hpp"
#include "Util.hpp"
#include "../../Application.hpp"
#include "../../Colors.hpp"
#include "../../Image.hpp"
#include "../../Ray3.hpp"
#include "../../Graphics/RenderSystem.hpp"
#include "../../Graphics/Shader.hpp"
#include "../../Graphics/Texture.hpp"
#include <QDateTime>
#include <QDir>
#include <QFont>
#include <QFontMetrics>
#include <QImage>
#include <QMouseEvent>

namespace Browse3D {

ModelDisplay::ModelDisplay(QWidget * parent, Model * model_)
: QGLWidget(parent),
  model(model_),
  camera(CoordinateFrame3(), Camera::ProjectionType::PERSPECTIVE, -1, 1, -1, 1, 0, 1, Camera::ProjectedYDirection::UP),
  camera_look_at(0, 0, -1),
  render_opts(RenderOptions::defaults()),
  mode(Mode::DEFAULT),
  view_edit_mode(ViewEditMode::DEFAULT),
  background_texture(NULL),
  background_shader(NULL)
{
  alwaysAssertM(model, "ModelDisplay: Can't create a display without a valid model");

  setFocusPolicy(Qt::StrongFocus);
  setMouseTracking(true);

  render_opts.useVertexData() = true;
  render_opts.sendNormals() = true;
  render_opts.sendColors() = true;
  render_opts.sendTexCoords() = false;
  render_opts.drawEdges() = true;
  render_opts.overrideEdgeColor() = true;
  render_opts.edgeColor() = ColorRGB(0.15f, 0.25f, 0.5f);

  connect(model, SIGNAL(geometryChanged(Model const *)), this, SLOT(modelGeometryChanged()));
  connect(model, SIGNAL(needsRedraw(Model const *)), this, SLOT(update()));
}

void
ModelDisplay::initializeGL()
{
#ifndef THEA_WINDOWS
  glEnable(GL_MULTISAMPLE);
#endif
}

Ray3
ModelDisplay::computePickRay(QPointF const & p) const
{
  return Browse3D::computePickRay(p, camera, width(), height());
}

QPointF
ModelDisplay::project(Vector3 const & p) const
{
  Vector3 proj_p = camera.project(p);
  return QPointF(0.5 * (p.x() + 1) * width(), 0.5 * (p.y() + 1) * height());
}

void
ModelDisplay::fitViewToModel()
{
  updateCameraFromModel();
  update();
}

void
ModelDisplay::updateCameraFromModel()
{
  model->updateBounds();
  updateCameraFrameFromModel();
  updateCameraProjection();
}

void
ModelDisplay::updateCameraFrameFromModel()
{
  camera_look_at = model->getBounds().getCenter();
  Real model_scale = model->getBounds().getExtent().fastLength();
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

  Real camera_separation = (camera_look_at - camera.getPosition()).length();
  Real near_dist = 0.01 * camera_separation;
  Real far_dist  = 1000 * camera_separation;

  hw *= (0.5f * near_dist);
  hh *= (0.5f * near_dist);

  camera.setProjection(Camera::ProjectionType::PERSPECTIVE, -hw, hw, -hh, hh, near_dist, far_dist,
                       Camera::ProjectedYDirection::UP);
}

void
ModelDisplay::modelGeometryChanged()
{
  fitViewToModel();
}

void
ModelDisplay::renderShaded()
{
  render_opts.drawFaces() = true;
  render_opts.drawEdges() = false;
  update();

  qDebug() << "Rendering shaded faces";
}

void
ModelDisplay::renderWireframe()
{
  render_opts.drawFaces() = false;
  render_opts.drawEdges() = true;
  update();

  qDebug() << "Rendering wireframe";
}

void
ModelDisplay::renderShadedWireframe()
{
  render_opts.drawFaces() = true;
  render_opts.drawEdges() = true;
  update();

  qDebug() << "Rendering shaded faces with wireframe edges";
}

void
ModelDisplay::resizeGL(int w, int h)
{
  glViewport(0, 0, w, h);
  fitViewToModel();
  update();
}

void
ModelDisplay::paintGL()
{
  using namespace Graphics;

  if (!app().hasRenderSystem()) return;

  RenderSystem & rs = *app().getRenderSystem();

  rs.setColorClearValue(app().options().bg_color);
  rs.clear();

  rs.setColorWrite(true, true, true, true);
  rs.setCamera(camera);

  if (!app().options().bg_plain)
    drawBackground(rs);

  drawAxes(rs);

  model->draw(rs, render_opts);

#ifdef THEA_OSX
  swapBuffers();  // for some reason this is necessary even if auto buffer swap is on
#endif
}

void
ModelDisplay::drawBackground(Graphics::RenderSystem & rs)
{
  // qDebug() << "drawBackground";

  using namespace Graphics;

  if (!background_texture)
  {
    static std::string const BG_IMAGE = "wood.jpg";
    Image bg_img;
    if (loadImage(bg_img, toQString(Application::getFullResourcePath("Images/" + BG_IMAGE))))
    {
      Texture::Options opts = Texture::Options::defaults();
      opts.interpolateMode = Texture::InterpolateMode::BILINEAR_NO_MIPMAP;

      background_texture = rs.createTexture("Background texture", bg_img, Texture::Format::RGBA8(), Texture::Dimension::DIM_2D,
                                            opts);
    }

    if (background_texture)
     qDebug().nospace() << "Loaded " << background_texture->getWidth() << "x" << background_texture->getHeight()
                        << " background image " << BG_IMAGE;
    else
    {
      qDebug() << "Couldn't load background texture";
      return;  // just use whatever we had earlier as a fallback
    }
  }

  if (!background_shader)
  {
    background_shader = rs.createShader("Background shader");

    background_shader->attachModuleFromFile(Shader::ModuleType::VERTEX,
                                            Application::getFullResourcePath("Materials/FlatTextureVert.glsl"));
    background_shader->attachModuleFromFile(Shader::ModuleType::FRAGMENT,
                                            Application::getFullResourcePath("Materials/FlatTextureFrag.glsl"));

    background_shader->setUniform("texture", background_texture);
  }

  rs.pushShader();
    rs.setShader(background_shader);

    // FIXME: For some strange reason calling glActiveTextureARB (as done by RenderSystem::sendTexCoord) within the
    // glBegin/glEnd block below triggers a GL invalid operation error. Hence we'll just use the plain vanilla GL calls instead
    // of the RenderSystem wrappers. (Update: this might have been fixed by the recent fix to GLTexture.)
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
ModelDisplay::drawAxes(Graphics::RenderSystem & rs)
{
  // qDebug() << "drawAxes";

  // TODO: Make the axes and labels occlude each other properly depending on orientation.

  using namespace Graphics;

  Matrix3 rot = camera.getFrame().getRotation();

  Vector2 x_axis = rot.getRow(0).xy();
  Vector2 y_axis = rot.getRow(1).xy();
  Vector2 z_axis = rot.getRow(2).xy();

  Vector2 aspect_ratio_compensation = width() < height() ? Vector2(1.0, width() / (Real)height())
                                                         : Vector2(height() / (Real)width(), 1.0);

  static Real const ARROW_LENGTH = 0.1f;
  Vector2 arrow_origin = Vector2(-1, -1) + (2.0f * ARROW_LENGTH) * aspect_ratio_compensation;
  Vector2 x_arrow_tip = arrow_origin + ARROW_LENGTH * aspect_ratio_compensation * x_axis;
  Vector2 y_arrow_tip = arrow_origin + ARROW_LENGTH * aspect_ratio_compensation * y_axis;
  Vector2 z_arrow_tip = arrow_origin + ARROW_LENGTH * aspect_ratio_compensation * z_axis;

  rs.pushColorFlags();
  rs.pushDepthFlags();
  rs.setDepthTest(RenderSystem::DepthTest::ALWAYS_PASS);

  rs.pushShader();
  rs.setShader(NULL);

    rs.setMatrixMode(RenderSystem::MatrixMode::PROJECTION); rs.pushMatrix(); rs.setIdentityMatrix();
    rs.setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); rs.pushMatrix(); rs.setIdentityMatrix();

      rs.beginPrimitive(RenderSystem::Primitive::LINES);

        rs.setColor(ColorRGB::red());
        rs.sendVertex(arrow_origin);
        rs.sendVertex(x_arrow_tip);

        rs.setColor(ColorRGB::green());
        rs.sendVertex(arrow_origin);
        rs.sendVertex(y_arrow_tip);

        rs.setColor(ColorRGB::blue());
        rs.sendVertex(arrow_origin);
        rs.sendVertex(z_arrow_tip);

      rs.endPrimitive();

    rs.setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); rs.popMatrix();
    rs.setMatrixMode(RenderSystem::MatrixMode::PROJECTION); rs.popMatrix();

#define DRAW_AXIS_LABELS
#ifdef DRAW_AXIS_LABELS

  // For some reason drawing text disables this and doesn't restore it properly
#ifndef THEA_WINDOWS
  bool multisample = (glIsEnabled(GL_MULTISAMPLE) == GL_TRUE);
#endif

  static QFont label_font("Helvetica", 12);
  static QFontMetrics fm(label_font);
  static QSize half_x_label_size = 0.5 * fm.tightBoundingRect("X").size();
  static QSize half_y_label_size = 0.5 * fm.tightBoundingRect("Y").size();
  static QSize half_z_label_size = 0.5 * fm.tightBoundingRect("Z").size();

  static Real const TEXT_OFFSET = 0.5;
  Vector2 x_text_center = x_arrow_tip + TEXT_OFFSET * (x_arrow_tip - arrow_origin);
  Vector2 y_text_center = y_arrow_tip + TEXT_OFFSET * (y_arrow_tip - arrow_origin);
  Vector2 z_text_center = z_arrow_tip + TEXT_OFFSET * (z_arrow_tip - arrow_origin);

  QPoint x_text_bottom_left((int)Math::round((0.5 * (1 + x_text_center.x())) * width()) - half_x_label_size.width(),
                            (int)Math::round((0.5 * (1 - x_text_center.y())) * height()) + half_x_label_size.height());

  QPoint y_text_bottom_left((int)Math::round((0.5 * (1 + y_text_center.x())) * width()) - half_y_label_size.width(),
                            (int)Math::round((0.5 * (1 - y_text_center.y())) * height()) + half_y_label_size.height());

  QPoint z_text_bottom_left((int)Math::round((0.5 * (1 + z_text_center.x())) * width()) - half_z_label_size.width(),
                            (int)Math::round((0.5 * (1 - z_text_center.y())) * height()) + half_z_label_size.height());

  rs.setColor(ColorRGB::red());
  renderText(x_text_bottom_left.x(), x_text_bottom_left.y(), "X", label_font);

  rs.setColor(ColorRGB::green());
  renderText(y_text_bottom_left.x(), y_text_bottom_left.y(), "Y", label_font);

  rs.setColor(ColorRGB::blue());
  renderText(z_text_bottom_left.x(), z_text_bottom_left.y(), "Z", label_font);

#ifndef THEA_WINDOWS
  if (multisample) glEnable(GL_MULTISAMPLE);  // for some reason drawing text disables this and doesn't restore it properly
#endif

#endif

  rs.popShader();
  rs.popDepthFlags();
  rs.popColorFlags();
}

void
ModelDisplay::keyPressEvent(QKeyEvent * event)
{
}

void
ModelDisplay::mousePressEvent(QMouseEvent * event)
{
  // qDebug() << "mousePressEvent";

  last_cursor = event->pos();

  Qt::MouseButton button = event->button();
  bool left    =  (button == Qt::LeftButton);
  bool right   =  (button == Qt::RightButton);
  bool middle  =  (button == Qt::MidButton);

  Qt::KeyboardModifiers kbmods = event->modifiers();
  bool shift  =  kbmods.testFlag(Qt::ShiftModifier);
  bool ctrl   =  kbmods.testFlag(Qt::ControlModifier);
  bool alt    =  kbmods.testFlag(Qt::AltModifier);
  bool no_mods     =  !(shift || ctrl || alt);
  bool shift_only  =   shift && !ctrl && !alt;
  bool ctrl_only   =  !shift &&  ctrl && !alt;
  bool alt_only    =  !shift && !ctrl &&  alt;

  bool edit_model = app().getMainWindow()->pickPoints();

  if (edit_model && (right || (ctrl && left)))
  {
    Ray3 ray = computePickRay(event->pos());
    if (model->pick(ray) >= 0)
    {
      model->mousePressEvent(event);
      if (event->isAccepted())
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

    view_drag_start = event->pos();

    event->accept();
    update();
    return;
  }

  event->ignore();
}

void
ModelDisplay::mouseMoveEvent(QMouseEvent * event)
{
  // qDebug() << "mouseMoveEvent";

  if (mode == Mode::EDIT_VIEW)
  {
    switch (view_edit_mode)
    {
      case ViewEditMode::PAN:
        panView(event);
        event->accept();
        break;

      case ViewEditMode::ROTATE:
        rotateView(event);
        event->accept();
        break;

      case ViewEditMode::ROTATE_ROLL:
        rollView(event);
        event->accept();
        break;

      case ViewEditMode::ZOOM:
        zoomView(event);
        event->accept();
        break;

      default: break;
    }
  }
  else  // register all mouse movement over the model, not just in EDIT_MODEL mode
  {
    model->mouseMoveEvent(event);
  }

  last_cursor = event->pos();
}

void
ModelDisplay::incrementViewTransform(AffineTransform3 const & tr)
{
  // qDebug() << "incrementViewTransform";

  AffineTransform3 inv_vt = tr.inverse();
  CoordinateFrame3 const & old_cframe = camera.getFrame();
  camera.setFrame(CoordinateFrame3(RigidTransform3::_fromAffine(AffineTransform3(inv_vt.getLinear() * old_cframe.getRotation(),
                                                                                 inv_vt * old_cframe.getTranslation()))));

  if (view_edit_mode == ViewEditMode::PAN)
    camera_look_at = inv_vt * camera_look_at;

  update();
}

void
ModelDisplay::mouseReleaseEvent(QMouseEvent * event)
{
  // qDebug() << "mouseReleaseEvent";

  if (mode == Mode::EDIT_VIEW)
  {
    if (view_edit_mode != ViewEditMode::DEFAULT)
    {
      view_edit_mode = ViewEditMode::DEFAULT;
      mode = Mode::DEFAULT;

      event->accept();
      update();
    }
  }
  else if (mode == Mode::EDIT_MODEL)
  {
    model->mouseReleaseEvent(event);
  }

  last_cursor = event->pos();
}

void
ModelDisplay::wheelEvent(QWheelEvent * event)
{
  // qDebug() << "wheelEvent";

  mode = Mode::EDIT_VIEW;
  view_edit_mode = ViewEditMode::ZOOM;

    zoomView(event);

  mode = Mode::DEFAULT;
  view_edit_mode = ViewEditMode::DEFAULT;
}

Real
ModelDisplay::getModelDistance() const
{
  // qDebug() << "getModelDistance";

  return model ? (model->getBounds().getCenter() - camera.getPosition()).length()
               : (camera_look_at - camera.getPosition()).length();
}

void
ModelDisplay::panView(QMouseEvent * event)
{
  // qDebug() << "panView";

  Vector3 trn = dragToTranslation(last_cursor, event->pos(), width(), height(), camera, getModelDistance());
  incrementViewTransform(AffineTransform3(Matrix3::identity(), trn));
}

void
ModelDisplay::rotateView(QMouseEvent * event)
{
  // qDebug() << "rotateView";

  Matrix3 rot = dragToRotation(last_cursor, event->pos(), width(), height(), camera);
  Vector3 trn = camera_look_at - rot * camera_look_at;
  incrementViewTransform(AffineTransform3(rot, trn));
}

void
ModelDisplay::rollView(QMouseEvent * event)
{
  // qDebug() << "rollView";

  Vector2 center(0.5f * width(), 0.5f * height());
  Vector2 u = Vector2(last_cursor.x(), last_cursor.y()) - center;
  Vector2 v = Vector2(event->pos().x(), event->pos().y()) - center;
  Real s = u.x() * v.y() - u.y() * v.x();
  Real c = u.dot(v);
  Real angle = Math::fastArcTan2(s, c);

  Matrix3 rot = Matrix3::rotationAxisAngle(camera.getLookDirection(), angle);
  Vector3 trn = camera_look_at - rot * camera_look_at;
  incrementViewTransform(AffineTransform3(rot, trn));
}

namespace ModelDisplayInternal {

AffineTransform3
zoomTransform(Real zoom, Real camera_separation, Vector3 const & look_dir)
{
  static Real const MIN_ZOOM = 0.25;
  Vector3 trn = (1.0f / std::max(zoom, MIN_ZOOM) - 1) * camera_separation * look_dir;
  return AffineTransform3(Matrix3::identity(), trn);
}

} // namespace ModelDisplayInternal

void
ModelDisplay::zoomView(QMouseEvent * event)
{
  // qDebug() << "zoomView";

  using namespace ModelDisplayInternal;

  Real zoom = dragToScale(last_cursor, event->pos(), width(), height(), camera);

  // Zoom in at mouse position, zoom out at screen center
  Vector3 dir;
  if (zoom > 1)
    dir = computePickRay(view_drag_start).getDirection().unit();
  else
    dir = camera.getLookDirection();

  Real camera_separation = (camera_look_at - camera.getPosition()).length();
  incrementViewTransform(zoomTransform(zoom, camera_separation, dir));
}

void
ModelDisplay::zoomView(QWheelEvent * event)
{
  // qDebug() << "zoomView";

  using namespace ModelDisplayInternal;

  static Real const ZOOM_PER_DEGREE = 2.0f / 180.0f;
  Real num_degrees = event->delta() / 8;  // delta is in eighths of a degree
  Real zoom = 1 + ZOOM_PER_DEGREE * num_degrees;

  // Zoom in at mouse position, zoom out at screen center
  Vector3 dir;
  if (zoom > 1)
    dir = computePickRay(event->pos()).getDirection().unit();
  else
    dir = camera.getLookDirection();

  Real camera_separation = (camera_look_at - camera.getPosition()).length();
  incrementViewTransform(zoomTransform(zoom, camera_separation, dir));
}

void
ModelDisplay::saveScreenshot(QString path)
{
  if (path.isEmpty())
  {
    Model const * model = app().getMainWindow()->getModel();
    QString prefix = model ? QFileInfo(model->getFilename()).baseName() : "Browse3D-Screenshot";
    path = getFullPath(QDir::homePath(), prefix + QDateTime::currentDateTime().toString("-yyyy-MM-dd-hh-mm-ss") + ".png");
  }

  QImage img = grabFrameBuffer();
  if (img.save(path))
    THEA_CONSOLE << "Saved screenshot to " << path;
  else
    THEA_ERROR << "Could not save screenshot to " << path;
}

} // namespace Browse3D
