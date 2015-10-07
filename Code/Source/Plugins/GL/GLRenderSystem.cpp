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

#include "GLRenderSystem.hpp"
#include "GLCaps.hpp"
#include "GLHeaders.hpp"
#include "GLTexture.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

// Convert a RenderSystem primitive to the corresponding OpenGL enum.
static GLenum
RenderSystem__primitiveToGLenum(RenderSystem::Primitive primitive)
{
  switch (primitive)
  {
    case RenderSystem::Primitive::POINTS:         return GL_POINTS;
    case RenderSystem::Primitive::LINES:          return GL_LINES;
    case RenderSystem::Primitive::LINE_STRIP:     return GL_LINE_STRIP;
    case RenderSystem::Primitive::LINE_LOOP:      return GL_LINE_LOOP;
    case RenderSystem::Primitive::TRIANGLES:      return GL_TRIANGLES;
    case RenderSystem::Primitive::TRIANGLE_STRIP: return GL_TRIANGLE_STRIP;
    case RenderSystem::Primitive::TRIANGLE_FAN:   return GL_TRIANGLE_FAN;
    case RenderSystem::Primitive::QUADS:          return GL_QUADS;
    case RenderSystem::Primitive::QUAD_STRIP:     return GL_QUAD_STRIP;
    case RenderSystem::Primitive::POLYGON:        return GL_POLYGON;
  }

  throw Error("GLRenderSystem: Unknown primitive");
}

// Convert a RenderSystem depth test to the corresponding OpenGL enum.
static GLenum
RenderSystem__depthTestToGLenum(RenderSystem::DepthTest depth_test)
{
  switch (depth_test)
  {
    case RenderSystem::DepthTest::GREATER     : return GL_GREATER;
    case RenderSystem::DepthTest::LESS        : return GL_LESS;
    case RenderSystem::DepthTest::GEQUAL      : return GL_GEQUAL;
    case RenderSystem::DepthTest::LEQUAL      : return GL_LEQUAL;
    case RenderSystem::DepthTest::NOTEQUAL    : return GL_NOTEQUAL;
    case RenderSystem::DepthTest::EQUAL       : return GL_EQUAL;
    case RenderSystem::DepthTest::ALWAYS_PASS : return GL_ALWAYS;
    case RenderSystem::DepthTest::NEVER_PASS  : return GL_NEVER;
  }

  throw Error("GLRenderSystem: Unknown depth test");
}

GLRenderSystem::GLRenderSystem(char const * name_)
: name(name_), current_framebuffer(NULL), current_shader(NULL)
{
  GLCaps::init();

  if (!THEA_GL_SUPPORTS(ARB_texture_non_power_of_two))
    THEA_WARNING << getName() << ": Non-power-of-two textures are not supported";

  if (!THEA_GL_SUPPORTS(ARB_vertex_buffer_object))
    THEA_WARNING << getName() << ": Vertex/index buffers in GPU memory are not supported";

  if (!THEA_GL_SUPPORTS(ARB_shader_objects))
    THEA_WARNING << getName() << ": Programmable shaders are not supported";

  if (!THEA_GL_SUPPORTS(EXT_framebuffer_object))
    THEA_WARNING << getName() << ": Framebuffer objects are not supported";

  glDisable(GL_LIGHTING);
  setCullFace(CullFace::NONE);
  setDepthTest(DepthTest::LESS);
  THEA_CHECK_GL_OK
}

GLRenderSystem::~GLRenderSystem()
{
  for (FramebufferSet::const_iterator fi = created_framebuffers.begin(); fi != created_framebuffers.end(); ++fi)
    delete *fi;

  created_framebuffers.clear();

  for (TextureSet::const_iterator ti = created_textures.begin(); ti != created_textures.end(); ++ti)
    delete *ti;

  created_textures.clear();

  for (ShaderSet::const_iterator si = created_shaders.begin(); si != created_shaders.end(); ++si)
    delete *si;

  created_shaders.clear();

  for (VARAreaSet::const_iterator vi = created_varareas.begin(); vi != created_varareas.end(); ++vi)
    delete *vi;

  created_varareas.clear();
}

char const *
GLRenderSystem::describeSystem() const
{
  desc = GLCaps::describeSystem();
  return desc.c_str();
}

Framebuffer *
GLRenderSystem::createFramebuffer(char const * name_)
{
  GLFramebuffer * fb = new GLFramebuffer(this, name_);
  if (fb)
    created_framebuffers.insert(fb);

  return fb;
}

void
GLRenderSystem::destroyFramebuffer(Framebuffer * framebuffer)
{
  if (!framebuffer)
    return;

  if (created_framebuffers.erase(dynamic_cast<GLFramebuffer *>(framebuffer)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy framebuffer '" << framebuffer->getName()
               << "' which was not created using this rendersystem";
    return;
  }

  delete framebuffer;
}

Shader *
GLRenderSystem::createShader(char const * name_)
{
  GLShader * shader = new GLShader(this, name_);
  if (shader)
    created_shaders.insert(shader);

  return shader;
}

void
GLRenderSystem::destroyShader(Shader * shader)
{
  if (!shader)
    return;

  if (created_shaders.erase(dynamic_cast<GLShader *>(shader)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy shader '" << shader->getName()
               << "' which was not created using this rendersystem";
    return;
  }

  delete shader;
}

Texture *
GLRenderSystem::createTexture(char const * name_, int width, int height, int depth, Texture::Format const * desired_format,
                              Texture::Dimension dimension, Texture::Options const & options)
{
  GLTexture * tex = new GLTexture(this, name_, width, height, depth, desired_format, dimension, options);
  if (tex)
    created_textures.insert(tex);

  return tex;
}

Texture *
GLRenderSystem::createTexture(char const * name_, AbstractImage const & image, Texture::Format const * desired_format,
                              Texture::Dimension dimension, Texture::Options const & options)
{
  GLTexture * tex = new GLTexture(this, name_, image, desired_format, dimension, options);
  if (tex)
    created_textures.insert(tex);

  return tex;
}

Texture *
GLRenderSystem::createTexture(char const * name_, AbstractImage const * images[6], Texture::Format const * desired_format,
                              Texture::Options const & options)
{
  GLTexture * tex = new GLTexture(this, name_, images, desired_format, options);
  if (tex)
    created_textures.insert(tex);

  return tex;
}

void
GLRenderSystem::destroyTexture(Texture * texture)
{
  if (!texture)
    return;

  if (created_textures.erase(dynamic_cast<GLTexture *>(texture)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy texture '" << texture->getName()
               << "' which was not created using this rendersystem";
    return;
  }

  delete texture;
}

VARArea *
GLRenderSystem::createVARArea(char const * name_, long num_bytes, VARArea::Usage usage, bool gpu_memory)
{
  GLVARArea * vararea = new GLVARArea(this, name_, num_bytes, usage, gpu_memory);
  if (vararea)
    created_varareas.insert(vararea);

  return vararea;
}

void
GLRenderSystem::destroyVARArea(VARArea * area)
{
  if (!area)
    return;

  if (created_varareas.erase(dynamic_cast<GLVARArea *>(area)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy VAR area '" << area->getName()
               << "' which was not created using this rendersystem";
    return;
  }

  delete area;
}

void
GLRenderSystem::pushFramebuffer()
{
  glPushAttrib(GL_VIEWPORT_BIT | GL_COLOR_BUFFER_BIT);
  THEA_CHECK_GL_OK

  framebuffer_stack.push(current_framebuffer);
}

void
GLRenderSystem::setFramebuffer(Framebuffer * framebuffer)
{
  if (framebuffer)
  {
    GLFramebuffer * glfb = dynamic_cast<GLFramebuffer *>(framebuffer);
    debugAssertM(glfb, std::string(getName()) + ": Attempt to use a non-OpenGL framebuffer with a GL rendersystem");

    if (glfb != current_framebuffer)
    {
      glfb->use();
      current_framebuffer = glfb;
    }
  }
  else
  {
    if (current_framebuffer)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      current_framebuffer = NULL;
      THEA_CHECK_GL_OK
    }
  }
}

Framebuffer const *
GLRenderSystem::getFramebuffer() const
{
  return current_framebuffer;
}

Framebuffer *
GLRenderSystem::getFramebuffer()
{
  return current_framebuffer;
}

void
GLRenderSystem::popFramebuffer()
{
  alwaysAssertM(!framebuffer_stack.empty(), std::string(getName()) + ": No framebuffer to pop");

  Framebuffer * fb = framebuffer_stack.top();
  framebuffer_stack.pop();
  setFramebuffer(fb);

  glPopAttrib();
  THEA_CHECK_GL_OK
}

void
GLRenderSystem::pushShader()
{
  shader_stack.push(current_shader);
  pushTextures();  // since binding a shader can overwrite current texture bindings
}

void
GLRenderSystem::setShader(Shader * shader)
{
  alwaysAssertM(THEA_GL_SUPPORTS(ARB_shader_objects),
                std::string(getName()) + ": This OpenGL installation does not support shader objects");

  if (shader)
  {
    GLShader * glshader = dynamic_cast<GLShader *>(shader);
    debugAssertM(glshader, std::string(getName()) + ": Attempt to use a non-OpenGL shader with a GL rendersystem");

    if (glshader != current_shader)
    {
      glshader->use();
      current_shader = glshader;
    }
  }
  else
  {
    if (current_shader)
    {
      glUseProgramObjectARB(0);
      current_shader = NULL;
      THEA_CHECK_GL_OK
    }
  }
}

Shader const *
GLRenderSystem::getShader() const
{
  return current_shader;
}

Shader *
GLRenderSystem::getShader()
{
  return current_shader;
}

void
GLRenderSystem::popShader()
{
  debugAssertM(!shader_stack.empty(), std::string(getName()) + ": push/popShader calls not matched");

  popTextures();  // must be called before binding shader below

  setShader(shader_stack.top());
  shader_stack.pop();
}

void
GLRenderSystem::pushTextures()
{
  glPushAttrib(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
}

void
GLRenderSystem::setTexture(int texunit, Texture * texture)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
  else if (texunit != 0)
    throw Error(std::string(getName()) + ": Non-zero texture unit specified but multitexturing is not supported by OpenGL");

  if (texture)
  {
    GLTexture * gltex = dynamic_cast<GLTexture *>(texture);
    debugAssertM(gltex, std::string(getName()) + ": Attempt to use a non-OpenGL texture with a GL rendersystem");

    GLenum target = gltex->getGLTarget();
    GLuint id     = gltex->getGLID();

    glEnable(target);
    glBindTexture(target, id);
  }
  else
  {
    // Disable the texture unit
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_3D);
    glDisable(GL_TEXTURE_CUBE_MAP_ARB);
    glDisable(GL_TEXTURE_RECTANGLE_ARB);
  }

  THEA_CHECK_GL_OK
}

void
GLRenderSystem::popTextures()
{
  glPopAttrib();
}

void
GLRenderSystem::beginIndexedPrimitives()
{
  buffer_stack.push(current_buffer_state);
  current_buffer_state = BufferState();

  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
}

void
GLRenderSystem::setVertexAreaFromVAR(GLVAR const & v)
{
  alwaysAssertM(!current_buffer_state.vertex_area || (v.getArea() == current_buffer_state.vertex_area),
      std::string(getName()) + ": All vertex arrays used within a single begin/endIndexedPrimitives block must share the same VARArea");

  if (v.getArea() != current_buffer_state.vertex_area)
  {
    current_buffer_state.vertex_area = const_cast<GLVAR &>(v).getArea();
    if (current_buffer_state.vertex_area->inGPUMemory())
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, current_buffer_state.vertex_area->getGLBuffer());
  }
}

void
GLRenderSystem::setIndexAreaFromVAR(GLVAR const & v)
{
  alwaysAssertM(!current_buffer_state.index_area || (v.getArea() == current_buffer_state.index_area),
                std::string(getName())
              + ": All index arrays used within a single begin/endIndexedPrimitives block must share the same VARArea");

  if (v.getArea() != current_buffer_state.index_area)
  {
    current_buffer_state.index_area = const_cast<GLVAR &>(v).getArea();
    if (current_buffer_state.index_area->inGPUMemory())
      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, current_buffer_state.index_area->getGLBuffer());
  }
}

void
GLRenderSystem::setVertexArray(VAR const * vertices)
{
  if (vertices)
  {
    assert(vertices->isValid());

    GLVAR const & g = dynamic_cast<GLVAR const &>(*vertices);
    assert(g.getGLType() != GL_UNSIGNED_BYTE && g.getGLType() != GL_UNSIGNED_SHORT && g.getGLType() != GL_UNSIGNED_INT);

    setVertexAreaFromVAR(g);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(g.getNumComponents(), g.getGLType(), g.getElementSize(), g.getBasePointer());
  }
  else
    glDisableClientState(GL_VERTEX_ARRAY);
}

void
GLRenderSystem::setColorArray(VAR const * colors)
{
  if (colors)
  {
    assert(colors->isValid());
    GLVAR const & g = dynamic_cast<GLVAR const &>(*colors);

    setVertexAreaFromVAR(g);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(g.getNumComponents(), g.getGLType(), g.getElementSize(), g.getBasePointer());
  }
  else
    glDisableClientState(GL_COLOR_ARRAY);
}

void
GLRenderSystem::setTexCoordArray(int texunit, VAR const * texcoords)
{
  if (texcoords)
  {
    assert(texcoords->isValid());
    debugAssertM(THEA_GL_SUPPORTS(ARB_multitexture) || (texunit == 0),
                 std::string(getName()) + ": Graphics card does not support multitexture");

    GLVAR const & g = dynamic_cast<GLVAR const &>(*texcoords);

    if (THEA_GL_SUPPORTS(ARB_multitexture))
      glClientActiveTextureARB(GL_TEXTURE0_ARB + texunit);

    setVertexAreaFromVAR(g);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer(g.getNumComponents(), g.getGLType(), g.getElementSize(), g.getBasePointer());

    if (THEA_GL_SUPPORTS(ARB_multitexture))
      glClientActiveTextureARB(GL_TEXTURE0_ARB);
  }
  else
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

void
GLRenderSystem::setNormalArray(VAR const * normals)
{
  if (normals)
  {
    assert(normals->isValid());

    GLVAR const & g = dynamic_cast<GLVAR const &>(*normals);
    assert(g.getNumComponents() == 3);
    assert(g.getGLType() != GL_UNSIGNED_BYTE && g.getGLType() != GL_UNSIGNED_SHORT && g.getGLType() != GL_UNSIGNED_INT);

    setVertexAreaFromVAR(g);
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(g.getGLType(), g.getElementSize(), g.getBasePointer());
  }
  else
    glDisableClientState(GL_NORMAL_ARRAY);
}

void
GLRenderSystem::setIndexArray(VAR const * indices)
{
  if (indices)
  {
    assert(indices->isValid());

    GLVAR const & g = dynamic_cast<GLVAR const &>(*indices);
    assert(g.getNumComponents() == 1);
    assert(g.getGLTarget() == GL_ELEMENT_ARRAY_BUFFER_ARB);

    setIndexAreaFromVAR(g);
    glEnableClientState(GL_VERTEX_ARRAY);
    current_buffer_state.index_var = g;
  }
  else
  {
    current_buffer_state.index_var = GLVAR();  // an invalid buffer
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
  }
}

void
GLRenderSystem::endIndexedPrimitives()
{
  debugAssertM(!buffer_stack.empty(), std::string(getName()) + ": begin/endIndexedPrimitives calls not matched");

  glPopClientAttrib();

  current_buffer_state = buffer_stack.top();
  buffer_stack.pop();
}

RenderSystem::MatrixMode
GLRenderSystem::getMatrixMode() const
{
  GLint gl_mode;
  glGetIntegerv(GL_MATRIX_MODE, &gl_mode);
  switch (gl_mode)
  {
    case GL_MODELVIEW:   return MatrixMode::MODELVIEW;
    case GL_PROJECTION:  return MatrixMode::PROJECTION;
    case GL_TEXTURE:     return MatrixMode::TEXTURE;
    case GL_COLOR:       return MatrixMode::COLOR;
    default:             throw Error(std::string(getName()) + "Unknown matrix mode");
  }
}

void
GLRenderSystem::setMatrixMode(MatrixMode mode)
{
  switch (mode)
  {
    case MatrixMode::MODELVIEW:   glMatrixMode(GL_MODELVIEW); break;
    case MatrixMode::PROJECTION:  glMatrixMode(GL_PROJECTION); break;
    case MatrixMode::TEXTURE:     glMatrixMode(GL_TEXTURE); break;
    case MatrixMode::COLOR:       glMatrixMode(GL_COLOR); break;
    default:                      throw Error(std::string(getName()) + "Unsupported matrix mode");
  }
}

void
GLRenderSystem::pushMatrix()
{
  glPushMatrix();
}

Matrix4
GLRenderSystem::getMatrix(MatrixMode mode) const
{
  float f[16];
  glGetFloatv((mode == MatrixMode::MODELVIEW ? GL_MODELVIEW_MATRIX : GL_PROJECTION_MATRIX), f);

  THEA_CHECK_GL_OK

  return Matrix4(f[0], f[4], f[ 8], f[12],
                 f[1], f[5], f[ 9], f[13],
                 f[2], f[6], f[10], f[14],
                 f[3], f[7], f[11], f[15]);
}

void
GLRenderSystem::setMatrix(Matrix4 const & m)
{
  GLfloat f[16];
  m.getElementsColumnMajor(f);
  glLoadMatrixf(f);
}

void
GLRenderSystem::setIdentityMatrix()
{
  glLoadIdentity();
}

void
GLRenderSystem::multMatrix(Matrix4 const & m)
{
  GLfloat f[16];
  m.getElementsColumnMajor(f);
  glMultMatrixf(f);
}

void
GLRenderSystem::popMatrix()
{
  glPopMatrix();
}

void
GLRenderSystem::sendIndices(Primitive primitive, long num_indices, uint8 const * indices)
{
  glDrawElements(RenderSystem__primitiveToGLenum(primitive), (GLsizei)num_indices, GL_UNSIGNED_BYTE, indices);
}

void
GLRenderSystem::sendIndices(Primitive primitive, long num_indices, uint16 const * indices)
{
  glDrawElements(RenderSystem__primitiveToGLenum(primitive), (GLsizei)num_indices, GL_UNSIGNED_SHORT, indices);
}

void
GLRenderSystem::sendIndices(Primitive primitive, long num_indices, uint32 const * indices)
{
  glDrawElements(RenderSystem__primitiveToGLenum(primitive), (GLsizei)num_indices, GL_UNSIGNED_INT, indices);
}

void
GLRenderSystem::sendSequentialIndices(Primitive primitive, int first_index, int num_indices)
{
  glDrawArrays(RenderSystem__primitiveToGLenum(primitive), first_index, (GLsizei)num_indices);
}

void
GLRenderSystem::sendIndicesFromArray(Primitive primitive, long offset, long num_indices)
{
  alwaysAssertM(current_buffer_state.index_var.isValid(), std::string(getName()) + ": No valid index array set");

  uint8 * ptr = static_cast<uint8 *>(current_buffer_state.index_var.getBasePointer());
  int elem_size = current_buffer_state.index_var.getElementSize();

  glDrawElements(RenderSystem__primitiveToGLenum(primitive),
                 (GLsizei)num_indices,
                 current_buffer_state.index_var.getGLType(),
                 ptr + offset * elem_size);
}

void
GLRenderSystem::beginPrimitive(Primitive primitive)
{
  glBegin(RenderSystem__primitiveToGLenum(primitive));
}

void
GLRenderSystem::sendVertex(Vector2 const & vertex)
{
  glVertex2f(vertex.x(), vertex.y());
}

void
GLRenderSystem::sendVertex(float x, float y)
{
  glVertex2f(x, y);
}

void
GLRenderSystem::sendVertex(double x, double y)
{
  glVertex2d(x, y);
}

void
GLRenderSystem::sendVertex(Vector3 const & vertex)
{
  glVertex3f(vertex.x(), vertex.y(), vertex.z());
}

void
GLRenderSystem::sendVertex(float x, float y, float z)
{
  glVertex3f(x, y, z);
}

void
GLRenderSystem::sendVertex(double x, double y, double z)
{
  glVertex3d(x, y, z);
}

void
GLRenderSystem::sendVertex(Vector4 const & vertex)
{
  glVertex4f(vertex.x(), vertex.y(), vertex.z(), vertex.w());
}

void
GLRenderSystem::sendVertex(float x, float y, float z, float w)
{
  glVertex4f(x, y, z, w);
}

void
GLRenderSystem::sendVertex(double x, double y, double z, double w)
{
  glVertex4d(x, y, z, w);
}

void
GLRenderSystem::sendNormal(Vector3 const & normal)
{
  glNormal3f(normal.x(), normal.y(), normal.z());
}

void
GLRenderSystem::sendNormal(float x, float y, float z)
{
  glNormal3f(x, y, z);
}

void
GLRenderSystem::sendNormal(double x, double y, double z)
{
  glNormal3d(x, y, z);
}

void
GLRenderSystem::sendTexCoord(int texunit, float texcoord)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord1fARB(GL_TEXTURE0_ARB + texunit, texcoord);
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord1f(texcoord);
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, double texcoord)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord1dARB(GL_TEXTURE0_ARB + texunit, texcoord);
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord1d(texcoord);
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, Vector2 const & texcoord)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord2fARB(GL_TEXTURE0_ARB + texunit, (float)texcoord.x(), (float)texcoord.y());
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord2f((float)texcoord.x(), (float)texcoord.y());
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, float x, float y)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord2fARB(GL_TEXTURE0_ARB + texunit, x, y);
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord2f(x, y);
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, double x, double y)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord2dARB(GL_TEXTURE0_ARB + texunit, x, y);
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord2d(x, y);
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, Vector3 const & texcoord)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord3fARB(GL_TEXTURE0_ARB + texunit, (float)texcoord.x(), (float)texcoord.y(), (float)texcoord.z());
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord3f((float)texcoord.x(), (float)texcoord.y(), (float)texcoord.z());
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, float x, float y, float z)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord3fARB(GL_TEXTURE0_ARB + texunit, x, y, z);
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord3f(x, y, z);
  }
}

void
GLRenderSystem::sendTexCoord(int texunit, double x, double y, double z)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glMultiTexCoord3dARB(GL_TEXTURE0_ARB + texunit, x, y, z);
  else
  {
    debugAssertM(texunit == 0, std::string(getName()) + ": Multitexturing not supported, texture unit must be zero");
    glTexCoord3d(x, y, z);
  }
}

void
GLRenderSystem::endPrimitive()
{
  glEnd();
}

void
GLRenderSystem::pushState()
{
  pushFramebuffer();
  pushShader();
  pushTextures();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);
}

void
GLRenderSystem::pushColorFlags()
{
  glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT);
}

void
GLRenderSystem::pushDepthFlags()
{
  glPushAttrib(GL_DEPTH_BUFFER_BIT);
}

void
GLRenderSystem::pushStencilFlags()
{
  glPushAttrib(GL_STENCIL_BUFFER_BIT);
}

void
GLRenderSystem::pushShapeFlags()
{
  glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
}

void
GLRenderSystem::setColorWrite(bool red, bool green, bool blue, bool alpha)
{
  glColorMask(red, green, blue, alpha);
}

void
GLRenderSystem::setDepthWrite(bool value)
{
  glDepthMask(value);
}

void
GLRenderSystem::setStencilWrite(uint32 mask)
{
  glStencilMask((GLuint)mask);
}

void
GLRenderSystem::setColor(ColorRGB const & value)
{
  glColor3f((GLfloat)value.r(), (GLfloat)value.g(), (GLfloat)value.b());
}

void
GLRenderSystem::setColor(ColorRGBA const & value)
{
  glColor4f((GLfloat)value.r(), (GLfloat)value.g(), (GLfloat)value.b(), (GLfloat)value.a());
}

void
GLRenderSystem::setColorClearValue(ColorRGB const & value)
{
  glClearColor((GLclampf)value.r(), (GLclampf)value.g(), (GLclampf)value.b(), 1.0f);
}

void
GLRenderSystem::setColorClearValue(ColorRGBA const & value)
{
  glClearColor((GLclampf)value.r(), (GLclampf)value.g(), (GLclampf)value.b(), (GLclampf)value.a());
}

void
GLRenderSystem::setDepthClearValue(Real value)
{
  glClearDepth((GLclampd)value);
}

void
GLRenderSystem::setStencilClearValue(int value)
{
  glClearStencil(value);
}

void
GLRenderSystem::clear()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
}

void
GLRenderSystem::clear(bool color, bool depth, bool stencil)
{
  glClear((color ? GL_COLOR_BUFFER_BIT : 0) | (depth ? GL_DEPTH_BUFFER_BIT : 0) | (stencil ? GL_STENCIL_BUFFER_BIT : 0));
}

void
GLRenderSystem::setDepthTest(DepthTest test)
{
  if (test == DepthTest::ALWAYS_PASS)
  {
    glDisable(GL_DEPTH_TEST);
    glDepthFunc(GL_ALWAYS);  // just to be safe
    return;
  }
  else
    glEnable(GL_DEPTH_TEST);

  glDepthFunc(RenderSystem__depthTestToGLenum(test));
}

void
GLRenderSystem::setCullFace(CullFace cull)
{
  if (cull == CullFace::NONE)
  {
    glDisable(GL_CULL_FACE);
    return;
  }
  else
  {
    glEnable(GL_CULL_FACE);
    glCullFace((cull == CullFace::FRONT) ? GL_FRONT : GL_BACK);
  }
}

void
GLRenderSystem::setPolygonOffset(bool enable, Real offset)
{
  if (enable)
  {
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(offset, offset);
  }
  else
    glDisable(GL_POLYGON_OFFSET_FILL);
}

void
GLRenderSystem::setPointSize(Real size)
{
  glEnable(GL_POINT_SMOOTH);
  glPointSize((GLfloat)size);
}

void
GLRenderSystem::popColorFlags()
{
  glPopAttrib();
}

void
GLRenderSystem::popDepthFlags()
{
  glPopAttrib();
}

void
GLRenderSystem::popStencilFlags()
{
  glPopAttrib();
}

void
GLRenderSystem::popShapeFlags()
{
  glPopAttrib();
}

void
GLRenderSystem::popState()
{
  glPopClientAttrib();
  glPopAttrib();
  popTextures();
  popShader();
  popFramebuffer();
}

void
GLRenderSystem::finishAllOperations()
{
  glFlush();
}

RenderSystem * GLRenderSystemFactory::singleton          =  NULL;
bool           GLRenderSystemFactory::singleton_created  =  false;

GLRenderSystemFactory::~GLRenderSystemFactory()
{
  delete singleton;
}

RenderSystem *
GLRenderSystemFactory::createRenderSystem(char const * name)
{
  if (singleton)
  {
    THEA_ERROR << "GLRenderSystemFactory: Only one OpenGL rendersystem can be created per process";
    return NULL;
  }

  try
  {
    singleton = new GLRenderSystem(name);
  }
  THEA_STANDARD_CATCH_BLOCKS(return NULL;, ERROR, "%s", "Could not create new OpenGL rendersystem")

  singleton_created = true;
  return singleton;
}

void
GLRenderSystemFactory::destroyRenderSystem(RenderSystem * render_system)
{
  alwaysAssertM(singleton == render_system,
                "GLRenderSystemFactory: Trying to destroy rendersystem not created with this factory");

  destroyAllRenderSystems();  // there's only one rendersystem
}

void
GLRenderSystemFactory::destroyAllRenderSystems()
{
  delete singleton;
  singleton = NULL;
  /* singleton_created = false; */  // Disabled until we figure out how to completely clean up the previous rendersystem.
}

} // namespace GL
} // namespace Graphics
} // namespace Thea
