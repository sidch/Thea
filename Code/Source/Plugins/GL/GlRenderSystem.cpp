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

#include "GlRenderSystem.hpp"
#include "GlCaps.hpp"
#include "GlHeaders.hpp"
#include "GlTexture.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

namespace GlRenderSystemInternal {

static_assert(sizeof(Real) == sizeof(float32) || sizeof(Real) == sizeof(float64),
              "GlRenderSystem: Real number type must be either float32 or float64");

// Convert an IRenderSystem primitive to the corresponding OpenGL enum.
static GLenum
primitiveToGLenum(int32 primitive)
{
  switch (primitive)
  {
    case IRenderSystem::Primitive::POINTS:         return GL_POINTS;
    case IRenderSystem::Primitive::LINES:          return GL_LINES;
    case IRenderSystem::Primitive::LINE_STRIP:     return GL_LINE_STRIP;
    case IRenderSystem::Primitive::LINE_LOOP:      return GL_LINE_LOOP;
    case IRenderSystem::Primitive::TRIANGLES:      return GL_TRIANGLES;
    case IRenderSystem::Primitive::TRIANGLE_STRIP: return GL_TRIANGLE_STRIP;
    case IRenderSystem::Primitive::TRIANGLE_FAN:   return GL_TRIANGLE_FAN;
    case IRenderSystem::Primitive::QUADS:          return GL_QUADS;
    case IRenderSystem::Primitive::QUAD_STRIP:     return GL_QUAD_STRIP;
    case IRenderSystem::Primitive::POLYGON:        return GL_POLYGON;
    default: THEA_ERROR << "GlRenderSystem: Unknown primitive"; return GL_INVALID_ENUM;
  }
}

// Convert a IRenderSystem depth test to the corresponding OpenGL enum.
static GLenum
depthTestToGLenum(int32 depth_test)
{
  switch (depth_test)
  {
    case IRenderSystem::DepthTest::GREATER     : return GL_GREATER;
    case IRenderSystem::DepthTest::LESS        : return GL_LESS;
    case IRenderSystem::DepthTest::GEQUAL      : return GL_GEQUAL;
    case IRenderSystem::DepthTest::LEQUAL      : return GL_LEQUAL;
    case IRenderSystem::DepthTest::NOTEQUAL    : return GL_NOTEQUAL;
    case IRenderSystem::DepthTest::EQUAL       : return GL_EQUAL;
    case IRenderSystem::DepthTest::ALWAYS_PASS : return GL_ALWAYS;
    case IRenderSystem::DepthTest::NEVER_PASS  : return GL_NEVER;
    default: THEA_ERROR << "GlRenderSystem: Unknown depth test"; return GL_INVALID_ENUM;
  }
}

} // namespace GlRenderSystemInternal

GlRenderSystem::GlRenderSystem(char const * name_)
: name(name_), current_framebuffer(nullptr), current_shader(nullptr)
{
  GlCaps::init();

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

GlRenderSystem::~GlRenderSystem()
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

  for (BufferPoolSet::const_iterator vi = created_bufpools.begin(); vi != created_bufpools.end(); ++vi)
    delete *vi;

  created_bufpools.clear();
}

char const *
GlRenderSystem::describeSystem() const
{
  desc = GlCaps::describeSystem();
  return desc.c_str();
}

char const *
GlRenderSystem::getErrorString() const
{
  GLenum err_code;
  if ((err_code = glGetError()) != GL_NO_ERROR)
    return theaGlErrorString(err_code);
  else
    return nullptr;
}

IFramebuffer *
GlRenderSystem::createFramebuffer(char const * name_)
{
  GlFramebuffer * fb = nullptr;
  try
  {
    fb = new GlFramebuffer(this, name_);
    if (fb)
      created_framebuffers.insert(fb);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s: Could not create new OpenGL framebuffer", getName())

  return fb;
}

int8
GlRenderSystem::destroyFramebuffer(IFramebuffer * framebuffer)
{
  if (!framebuffer)
    return true;

  if (created_framebuffers.erase(dynamic_cast<GlFramebuffer *>(framebuffer)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy framebuffer '" << framebuffer->getName()
               << "' which was not created using this rendersystem";
    return false;
  }

  delete framebuffer;

  return true;
}

IShader *
GlRenderSystem::createShader(char const * name_)
{
  GlShader * shader = nullptr;
  try
  {
    shader = new GlShader(this, name_);
    if (shader)
      created_shaders.insert(shader);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s: Could not create new OpenGL shader", getName())

  return shader;
}

int8
GlRenderSystem::destroyShader(IShader * shader)
{
  if (!shader)
    return true;

  if (created_shaders.erase(dynamic_cast<GlShader *>(shader)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy shader '" << shader->getName()
               << "' which was not created using this rendersystem";
    return false;
  }

  delete shader;

  return true;
}

ITexture *
GlRenderSystem::createTexture(char const * name_, int64 width, int64 height, int64 depth,
                              ITexture::Format const * desired_format, int32 dimension, ITexture::Options const * options)
{
  GlTexture * tex = nullptr;
  try
  {
    tex = new GlTexture(this, name_, width, height, depth, desired_format, dimension, options);
    if (tex)
      created_textures.insert(tex);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s: Could not create new OpenGL texture", getName())

  return tex;
}

ITexture *
GlRenderSystem::createTexture(char const * name_, IImage const * image, ITexture::Format const * desired_format,
                              int32 dimension, ITexture::Options const * options)
{
  GlTexture * tex = nullptr;
  try
  {
    tex = new GlTexture(this, name_, image, desired_format, dimension, options);
    if (tex)
      created_textures.insert(tex);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s: Could not create new OpenGL texture", getName())

  return tex;
}

ITexture *
GlRenderSystem::createTexture(char const * name_, IImage const * images[6], ITexture::Format const * desired_format,
                              ITexture::Options const * options)
{
  GlTexture * tex = nullptr;
  try
  {
    GlTexture * tex = new GlTexture(this, name_, images, desired_format, options);
    if (tex)
      created_textures.insert(tex);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s: Could not create new OpenGL cubemap texture", getName())

  return tex;
}

int8
GlRenderSystem::destroyTexture(ITexture * texture)
{
  if (!texture)
    return true;

  if (created_textures.erase(dynamic_cast<GlTexture *>(texture)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy texture '" << texture->getName()
               << "' which was not created using this rendersystem";
    return false;
  }

  delete texture;

  return true;
}

IBufferPool *
GlRenderSystem::createBufferPool(char const * name_, int64 num_bytes, int32 usage, int8 gpu_memory)
{
  GlBufferPool * bufpool = nullptr;
  try
  {
    bufpool = new GlBufferPool(this, name_, num_bytes, usage, gpu_memory);
    if (bufpool)
      created_bufpools.insert(bufpool);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s: Could not create new OpenGL buffer pool", getName())

  return bufpool;
}

int8
GlRenderSystem::destroyBufferPool(IBufferPool * pool)
{
  if (!pool)
    return true;

  if (created_bufpools.erase(dynamic_cast<GlBufferPool *>(pool)) < 1)
  {
    THEA_ERROR << getName() << ": Attempting to destroy buffer pool '" << pool->getName()
               << "' which was not created using this rendersystem";
    return false;
  }

  delete pool;

  return true;
}

int8
GlRenderSystem::pushFramebuffer()
{
  glPushAttrib(GL_VIEWPORT_BIT | GL_COLOR_BUFFER_BIT);
  THEA_GL_OK_OR_RETURN(0)

  framebuffer_stack.push(current_framebuffer);

  return true;
}

int8
GlRenderSystem::setFramebuffer(IFramebuffer * framebuffer)
{
  if (framebuffer)
  {
    GlFramebuffer * glfb = dynamic_cast<GlFramebuffer *>(framebuffer);
    if (!glfb) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL framebuffer with a GL rendersystem"; return false; }

    if (glfb != current_framebuffer)
    {
      if (!glfb->use()) return false;
      current_framebuffer = glfb;
    }
  }
  else
  {
    if (current_framebuffer)
    {
      glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
      current_framebuffer = nullptr;
      THEA_GL_OK_OR_RETURN(0)
    }
  }

  return true;
}

IFramebuffer const *
GlRenderSystem::getFramebuffer() const
{
  return current_framebuffer;
}

IFramebuffer *
GlRenderSystem::getFramebuffer()
{
  return current_framebuffer;
}

int8
GlRenderSystem::popFramebuffer()
{
  if (framebuffer_stack.empty()) { THEA_ERROR << getName() << ": No framebuffer to pop"; return false; }

  if (!setFramebuffer(framebuffer_stack.top()))
    return false;

  framebuffer_stack.pop();
  glPopAttrib();
  THEA_GL_OK_OR_RETURN(0)

  return true;
}

int8
GlRenderSystem::pushShader()
{
  shader_stack.push(current_shader);
  return pushTextures();  // since binding a shader can overwrite current texture bindings
}

int8
GlRenderSystem::setShader(IShader * shader)
{
  if (!THEA_GL_SUPPORTS(ARB_shader_objects))
  { THEA_ERROR << getName() << ": This OpenGL installation does not support shader objects"; return false; }

  if (shader)
  {
    GlShader * glshader = dynamic_cast<GlShader *>(shader);
    if (!glshader)
    { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL shader with a GL rendersystem"; return false; }

    if (glshader != current_shader)
    {
      if (!glshader->use()) return false;
      current_shader = glshader;
    }
  }
  else
  {
    if (current_shader)
    {
      glUseProgramObjectARB(0);
      current_shader = nullptr;
      THEA_GL_OK_OR_RETURN(0)
    }
  }

  return true;
}

IShader const *
GlRenderSystem::getShader() const
{
  return current_shader;
}

IShader *
GlRenderSystem::getShader()
{
  return current_shader;
}

int8
GlRenderSystem::popShader()
{
  if (shader_stack.empty()) { THEA_ERROR << getName() << ": push/popShader calls not matched"; return false; }

  if (!popTextures()) return false;  // must be called before binding shader below

  if (!setShader(shader_stack.top())) return false;
  shader_stack.pop();

  return true;
}

int8
GlRenderSystem::pushTextures()
{
  glPushAttrib(GL_TEXTURE_BIT | GL_ENABLE_BIT);  // Can we do without ENABLE_BIT? The doc is unclear.
  return true;
}

int8
GlRenderSystem::setTexture(int32 texunit, ITexture * texture)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
    glActiveTextureARB(GL_TEXTURE0_ARB + texunit);
  else if (texunit != 0)
  { THEA_ERROR << getName() << ": Non-zero texture unit specified but multitexturing isn't supported by OpenGL"; return false; }

  if (texture)
  {
    GlTexture * gltex = dynamic_cast<GlTexture *>(texture);
    if (!gltex) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL texture with a GL rendersystem"; return false; }

    GLenum target = gltex->getGlTarget();
    GLuint id     = gltex->getGlId();

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

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlRenderSystem::popTextures()
{
  glPopAttrib();
  return true;
}

int32
GlRenderSystem::getMatrixMode() const
{
  GLint gl_mode;
  glGetIntegerv(GL_MATRIX_MODE, &gl_mode);
  switch (gl_mode)
  {
    case GL_MODELVIEW:   return MatrixMode::MODELVIEW;
    case GL_PROJECTION:  return MatrixMode::PROJECTION;
    case GL_TEXTURE:     return MatrixMode::TEXTURE;
    case GL_COLOR:       return MatrixMode::COLOR;
    default:             THEA_ERROR << getName() << ": Unknown matrix mode"; return -1;
  }
}

int8
GlRenderSystem::setMatrixMode(int32 mode)
{
  switch (mode)
  {
    case MatrixMode::MODELVIEW:   glMatrixMode(GL_MODELVIEW); return true;
    case MatrixMode::PROJECTION:  glMatrixMode(GL_PROJECTION); return true;
    case MatrixMode::TEXTURE:     glMatrixMode(GL_TEXTURE); return true;
    case MatrixMode::COLOR:       glMatrixMode(GL_COLOR); return true;
  }

  THEA_ERROR << getName() << ": Unknown matrix mode";
  return false;
}

int8
GlRenderSystem::pushMatrix()
{
  glPushMatrix();
  return true;
}

int8
GlRenderSystem::pushViewMatrices()
{
  GLint gl_mode;
  glGetIntegerv(GL_MATRIX_MODE, &gl_mode);

  glMatrixMode(GL_MODELVIEW); glPushMatrix();
  glMatrixMode(GL_PROJECTION); glPushMatrix();

  glMatrixMode(gl_mode);

  return true;
}

int8
GlRenderSystem::getMatrix(int32 mode, IDenseMatrix<Real> * m) const
{
  if (!m) { THEA_ERROR << getName() << ": No output matrix supplied for getMatrix()"; return false; }
  if (m->rows() != 4 || m->cols() != 4) { THEA_ERROR << getName() << ": getMatrix() output matrix isn't 4x4"; return false; }

  float32 f[16];
  glGetFloatv((mode == MatrixMode::MODELVIEW ? GL_MODELVIEW_MATRIX : GL_PROJECTION_MATRIX), f);

  THEA_GL_OK_OR_RETURN(0)

  if (m->isRowMajor())
  {
    auto wm = Math::mapTo< Matrix<4, 4, Real, MatrixLayout::ROW_MAJOR> >(*m);
    wm << f[0], f[4], f[ 8], f[12],
          f[1], f[5], f[ 9], f[13],
          f[2], f[6], f[10], f[14],
          f[3], f[7], f[11], f[15];
  }
  else
  {
    auto wm = Math::mapTo< Matrix<4, 4, Real, MatrixLayout::COLUMN_MAJOR> >(*m);
    wm << f[0], f[4], f[ 8], f[12],
          f[1], f[5], f[ 9], f[13],
          f[2], f[6], f[10], f[14],
          f[3], f[7], f[11], f[15];
  }

  return true;
}

int8
GlRenderSystem::setMatrix(IDenseMatrix<Real> const * m)
{
  if (!m) { THEA_ERROR << getName() << ": Can't set null matrix"; return false; }
  if (m->rows() != 4 || m->cols() != 4) { THEA_ERROR << getName() << ": Can't set non-4x4 matrix"; return false; }

  GLfloat f[16];
  GlInternal::toArray(*m, f);
  glLoadMatrixf(f);

  return true;
}

int8
GlRenderSystem::setIdentityMatrix()
{
  glLoadIdentity();
  return true;
}

int8
GlRenderSystem::multMatrix(IDenseMatrix<Real> const * m)
{
  if (!m) { THEA_ERROR << getName() << ": Can't multiply by null matrix"; return false; }
  if (m->rows() != 4 || m->cols() != 4) { THEA_ERROR << getName() << ": Can't multiply by non-4x4 matrix"; return false; }

  GLfloat f[16];
  GlInternal::toArray(*m, f);
  glMultMatrixf(f);

  return true;
}

int8
GlRenderSystem::popMatrix()
{
  glPopMatrix();
  return true;
}

int8
GlRenderSystem::popViewMatrices()
{
  GLint gl_mode;
  glGetIntegerv(GL_MATRIX_MODE, &gl_mode);

  glMatrixMode(GL_PROJECTION); glPopMatrix();
  glMatrixMode(GL_MODELVIEW); glPopMatrix();

  glMatrixMode(gl_mode);

  return true;
}

int8
GlRenderSystem::beginIndexedPrimitives()
{
  buffer_stack.push(current_buffer_state);
  current_buffer_state = BufferState();

  glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);

  return true;
}

int8
GlRenderSystem::setAttributePoolFromBuffer(GlBuffer const & buf)
{
  if (current_buffer_state.attrib_pool && buf.getPool() != current_buffer_state.attrib_pool)
  {
    THEA_ERROR << getName() << ": Attribute arrays used in a begin/endIndexedPrimitives block must share the same BufferPool";
    return false;
  }

  if (buf.getPool() != current_buffer_state.attrib_pool)
  {
    current_buffer_state.attrib_pool = const_cast<GlBuffer &>(buf).getPool();
    if (current_buffer_state.attrib_pool->inGpuMemory())
    {
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, current_buffer_state.attrib_pool->getGlBuffer());
      THEA_GL_OK_OR_RETURN(0)
    }
  }

  return true;
}

int8
GlRenderSystem::setIndexPoolFromBuffer(GlBuffer const & buf)
{
  if (current_buffer_state.index_pool && buf.getPool() != current_buffer_state.index_pool)
  {
    THEA_ERROR << getName() << ": Index arrays used in a begin/endIndexedPrimitives block must share the same BufferPool";
    return false;
  }

  if (buf.getPool() != current_buffer_state.index_pool)
  {
    current_buffer_state.index_pool = const_cast<GlBuffer &>(buf).getPool();
    if (current_buffer_state.index_pool->inGpuMemory())
    {
      glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, current_buffer_state.index_pool->getGlBuffer());
      THEA_GL_OK_OR_RETURN(0)
    }
  }

  return true;
}

int8
GlRenderSystem::setVertexBuffer(IBuffer const * vertices)
{
  if (vertices)
  {
    if (!vertices->isValid()) { THEA_ERROR << getName() << ": Invalid vertex buffer"; return false; }

    GlBuffer const * g = dynamic_cast<GlBuffer const *>(vertices);
    if (!g) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL buffer with a GL rendersystem"; return false; }
    if (g->getGlType() == GL_UNSIGNED_BYTE || g->getGlType() == GL_UNSIGNED_SHORT || g->getGlType() == GL_UNSIGNED_INT)
    { THEA_ERROR << getName() << ": Buffer has unsigned integers -- did you mean to use this an index buffer?"; return false; }

    if (!setAttributePoolFromBuffer(*g)) return false;
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(g->numComponents(), g->getGlType(), g->getValueSize(), g->getBasePointer());
  }
  else
    glDisableClientState(GL_VERTEX_ARRAY);

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlRenderSystem::setColorBuffer(IBuffer const * colors)
{
  if (colors)
  {
    if (!colors->isValid()) { THEA_ERROR << getName() << ": Invalid color buffer"; return false; }

    GlBuffer const * g = dynamic_cast<GlBuffer const *>(colors);
    if (!g) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL buffer with a GL rendersystem"; return false; }

    if (!setAttributePoolFromBuffer(*g)) return false;
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(g->numComponents(), g->getGlType(), g->getValueSize(), g->getBasePointer());
  }
  else
    glDisableClientState(GL_COLOR_ARRAY);

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlRenderSystem::setTexCoordBuffer(int32 texunit, IBuffer const * texcoords)
{
  if (texcoords)
  {
    if (!THEA_GL_SUPPORTS(ARB_multitexture) && texunit != 0)
    { THEA_ERROR << getName() << ": OpenGL system does not support multitexture"; return false; }

    GlBuffer const * g = dynamic_cast<GlBuffer const *>(texcoords);
    if (!g) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL buffer with a GL rendersystem"; return false; }

    if (THEA_GL_SUPPORTS(ARB_multitexture))
      glClientActiveTextureARB(GL_TEXTURE0_ARB + texunit);

    if (!setAttributePoolFromBuffer(*g)) return false;
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    glTexCoordPointer(g->numComponents(), g->getGlType(), g->getValueSize(), g->getBasePointer());

    if (THEA_GL_SUPPORTS(ARB_multitexture))
      glClientActiveTextureARB(GL_TEXTURE0_ARB);
  }
  else
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlRenderSystem::setNormalBuffer(IBuffer const * normals)
{
  if (normals)
  {
    if (!normals->isValid()) { THEA_ERROR << getName() << ": Invalid normal buffer"; return false; }

    GlBuffer const * g = dynamic_cast<GlBuffer const *>(normals);
    if (!g) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL buffer with a GL rendersystem"; return false; }
    if (g->numComponents() != 3) { THEA_ERROR << getName() << ": Normals should be 3D"; return false; }
    if (g->getGlType() == GL_UNSIGNED_BYTE || g->getGlType() == GL_UNSIGNED_SHORT || g->getGlType() == GL_UNSIGNED_INT)
    { THEA_ERROR << getName() << ": Buffer has unsigned integers -- did you mean to use this an index buffer?"; return false; }

    if (!setAttributePoolFromBuffer(*g)) return false;
    glEnableClientState(GL_NORMAL_ARRAY);
    glNormalPointer(g->getGlType(), g->getValueSize(), g->getBasePointer());
  }
  else
    glDisableClientState(GL_NORMAL_ARRAY);

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlRenderSystem::setIndexBuffer(IBuffer const * indices)
{
  if (indices)
  {
    if (!indices->isValid()) { THEA_ERROR << getName() << ": Invalid index buffer"; return false; }

    GlBuffer const * g = dynamic_cast<GlBuffer const *>(indices);
    if (!g) { THEA_ERROR << getName() << ": Attempt to use a non-OpenGL buffer with a GL rendersystem"; return false; }
    if (g->numComponents() != 1) { THEA_ERROR << getName() << ": Indices should be 1D"; return false; }
    if (g->getGlTarget() != GL_ELEMENT_ARRAY_BUFFER_ARB)
    { THEA_ERROR << getName() << ": Index buffer target should be GL_ELEMENT_ARRAY_BUFFER_ARB"; return false; }

    if (!setIndexPoolFromBuffer(*g)) return false;
    glEnableClientState(GL_VERTEX_ARRAY);
    current_buffer_state.index_buf = *g;
  }
  else
  {
    current_buffer_state.index_buf = GlBuffer();  // an invalid buffer
    glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
  }

  THEA_GL_OK_OR_RETURN(0)
  return true;
}

int8
GlRenderSystem::sendIndices(int32 primitive, int64 num_indices, uint8 const * indices)
{
  glDrawElements(GlRenderSystemInternal::primitiveToGLenum(primitive), (GLsizei)num_indices, GL_UNSIGNED_BYTE, indices);
  return true;
}

int8
GlRenderSystem::sendIndices(int32 primitive, int64 num_indices, uint16 const * indices)
{
  glDrawElements(GlRenderSystemInternal::primitiveToGLenum(primitive), (GLsizei)num_indices, GL_UNSIGNED_SHORT, indices);
  return true;
}

int8
GlRenderSystem::sendIndices(int32 primitive, int64 num_indices, uint32 const * indices)
{
  glDrawElements(GlRenderSystemInternal::primitiveToGLenum(primitive), (GLsizei)num_indices, GL_UNSIGNED_INT, indices);
  return true;
}

int8
GlRenderSystem::sendSequentialIndices(int32 primitive, uint64 first_index, int64 num_indices)
{
  glDrawArrays(GlRenderSystemInternal::primitiveToGLenum(primitive), (GLint)first_index, (GLsizei)num_indices);
  return true;
}

int8
GlRenderSystem::sendIndicesFromBuffer(int32 primitive, uint64 offset, int64 num_indices)
{
  if (!current_buffer_state.index_buf.isValid()) { THEA_ERROR << getName() << ": No valid index array set"; return false; }

  uint8 * ptr = static_cast<uint8 *>(current_buffer_state.index_buf.getBasePointer());
  int32 elem_size = current_buffer_state.index_buf.getValueSize();

  glDrawElements(GlRenderSystemInternal::primitiveToGLenum(primitive), (GLsizei)num_indices,
                 current_buffer_state.index_buf.getGlType(), ptr + offset * elem_size);

  return true;
}

int8
GlRenderSystem::endIndexedPrimitives()
{
  if (buffer_stack.empty()) { THEA_ERROR << getName() << ": begin/endIndexedPrimitives calls not matched"; return false; }

  glPopClientAttrib();

  current_buffer_state = buffer_stack.top();
  buffer_stack.pop();

  return true;
}

int8
GlRenderSystem::beginPrimitive(int32 primitive)
{
  glBegin(GlRenderSystemInternal::primitiveToGLenum(primitive));
  return true;
}

int8
GlRenderSystem::sendVertex(int32 dims, Real const * coords)
{
  // The conditional should be totally optimized out at compile-time
  if (sizeof(Real) == 4)
  {
    GLfloat const * gl_coords = static_cast<GLfloat const *>((void const *)coords);
    switch (dims)
    {
      case 1: glVertex2f(gl_coords[0], 0); return true;
      case 2: glVertex2fv(gl_coords); return true;
      case 3: glVertex3fv(gl_coords); return true;
      case 4: glVertex4fv(gl_coords); return true;
    }
  }
  else
  {
    GLdouble const * gl_coords = static_cast<GLdouble const *>((void const *)coords);
    switch (dims)
    {
      case 1: glVertex2d(gl_coords[0], 0); return true;
      case 2: glVertex2dv(gl_coords); return true;
      case 3: glVertex3dv(gl_coords); return true;
      case 4: glVertex4dv(gl_coords); return true;
    }
  }

  THEA_ERROR << getName() << ": Unsupported vertex dimension " << dims;
  return false;
}

int8
GlRenderSystem::sendVertex(Real x, Real y)
{
  if (sizeof(Real) == 4)
    glVertex2f((GLfloat)x, (GLfloat)y);  // cast helps avoid loss of precision warning when Real == float64
  else
    glVertex2d(x, y);

  return true;
}

int8
GlRenderSystem::sendVertex(Real x, Real y, Real z)
{
  if (sizeof(Real) == 4)
    glVertex3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
  else
    glVertex3d(x, y, z);

  return true;
}

int8
GlRenderSystem::sendVertex(Real x, Real y, Real z, Real w)
{
  if (sizeof(Real) == 4)
    glVertex4f((GLfloat)x, (GLfloat)y, (GLfloat)z, (GLfloat)w);
  else
    glVertex4d(x, y, z, w);

  return true;
}

int8
GlRenderSystem::sendNormal(Real const * coords)
{
  if (sizeof(Real) == 4)
    glNormal3fv(static_cast<GLfloat const *>((void const *)coords));
  else
    glNormal3dv(static_cast<GLdouble const *>((void const *)coords));

  return true;
}

int8
GlRenderSystem::sendNormal(Real x, Real y, Real z)
{
  if (sizeof(Real) == 4)
    glNormal3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
  else
    glNormal3d(x, y, z);

  return true;
}

int8
GlRenderSystem::sendTexCoord(int32 texunit, int32 dims, Real const * coords)
{
  bool multi = THEA_GL_SUPPORTS(ARB_multitexture);
  if (!multi && texunit != 0)
  { THEA_ERROR << getName() << ": Multitexturing not supported, texture unit must be zero"; return false; }

  // The conditional should be totally optimized out at compile-time
  if (sizeof(Real) == 4)
  {
    GLfloat const * gc = static_cast<GLfloat const *>((void const *)coords);
    switch (dims)
    {
      case 1: if (multi) glMultiTexCoord1fvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord1fv(gc); return true;
      case 2: if (multi) glMultiTexCoord2fvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord2fv(gc); return true;
      case 3: if (multi) glMultiTexCoord3fvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord3fv(gc); return true;
      case 4: if (multi) glMultiTexCoord4fvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord4fv(gc); return true;
    }
  }
  else
  {
    GLdouble const * gc = static_cast<GLdouble const *>((void const *)coords);
    switch (dims)
    {
      case 1: if (multi) glMultiTexCoord1dvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord1dv(gc); return true;
      case 2: if (multi) glMultiTexCoord2dvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord2dv(gc); return true;
      case 3: if (multi) glMultiTexCoord3dvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord3dv(gc); return true;
      case 4: if (multi) glMultiTexCoord4dvARB(GL_TEXTURE0_ARB + texunit, gc); else glTexCoord4dv(gc); return true;
    }
  }

  THEA_ERROR << getName() << ": Unsupported texture coordinates dimension " << dims;
  return false;
}

int8
GlRenderSystem::sendTexCoord(int32 texunit, Real x)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
  {
    if (sizeof(Real) == 4) glMultiTexCoord1fARB(GL_TEXTURE0_ARB + texunit, (GLfloat)x);
    else                   glMultiTexCoord1dARB(GL_TEXTURE0_ARB + texunit, x);
  }
  else
  {
    if (texunit != 0) { THEA_ERROR << getName() << ": Multitexturing not supported, texture unit must be zero"; return false; }

    if (sizeof(Real) == 4) glTexCoord1f((GLfloat)x);
    else                   glTexCoord1d(x);
  }

  return true;
}

int8
GlRenderSystem::sendTexCoord(int32 texunit, Real x, Real y)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
  {
    if (sizeof(Real) == 4) glMultiTexCoord2fARB(GL_TEXTURE0_ARB + texunit, (GLfloat)x, (GLfloat)y);
    else                   glMultiTexCoord2dARB(GL_TEXTURE0_ARB + texunit, x, y);
  }
  else
  {
    if (texunit != 0) { THEA_ERROR << getName() << ": Multitexturing not supported, texture unit must be zero"; return false; }

    if (sizeof(Real) == 4) glTexCoord2f((GLfloat)x, (GLfloat)y);
    else                   glTexCoord2d(x, y);
  }

  return true;
}

int8
GlRenderSystem::sendTexCoord(int32 texunit, Real x, Real y, Real z)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
  {
    if (sizeof(Real) == 4) glMultiTexCoord3fARB(GL_TEXTURE0_ARB + texunit, (GLfloat)x, (GLfloat)y, (GLfloat)z);
    else                   glMultiTexCoord3dARB(GL_TEXTURE0_ARB + texunit, x, y, z);
  }
  else
  {
    if (texunit != 0) { THEA_ERROR << getName() << ": Multitexturing not supported, texture unit must be zero"; return false; }

    if (sizeof(Real) == 4) glTexCoord3f((GLfloat)x, (GLfloat)y, (GLfloat)z);
    else                   glTexCoord3d(x, y, z);
  }

  return true;
}

int8
GlRenderSystem::sendTexCoord(int32 texunit, Real x, Real y, Real z, Real w)
{
  if (THEA_GL_SUPPORTS(ARB_multitexture))
  {
    if (sizeof(Real) == 4) glMultiTexCoord4fARB(GL_TEXTURE0_ARB + texunit, (GLfloat)x, (GLfloat)y, (GLfloat)z, (GLfloat)w);
    else                   glMultiTexCoord4dARB(GL_TEXTURE0_ARB + texunit, x, y, z, w);
  }
  else
  {
    if (texunit != 0) { THEA_ERROR << getName() << ": Multitexturing not supported, texture unit must be zero"; return false; }

    if (sizeof(Real) == 4) glTexCoord4f((GLfloat)x, (GLfloat)y, (GLfloat)z, (GLfloat)w);
    else                   glTexCoord4d(x, y, z, w);
  }

  return true;
}

int8
GlRenderSystem::endPrimitive()
{
  glEnd();
  return true;
}

int8
GlRenderSystem::clear()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  return true;
}

int8
GlRenderSystem::clear(int8 color, int8 depth, int8 stencil)
{
  glClear((color ? GL_COLOR_BUFFER_BIT : 0) | (depth ? GL_DEPTH_BUFFER_BIT : 0) | (stencil ? GL_STENCIL_BUFFER_BIT : 0));
  return true;
}

int8
GlRenderSystem::pushState()
{
  if (!pushFramebuffer()) return false;
  if (!pushShader()) return false;
  if (!pushTextures()) return false;
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);

  return true;
}

int8
GlRenderSystem::pushColorFlags()
{
  glPushAttrib(GL_COLOR_BUFFER_BIT | GL_CURRENT_BIT);
  return true;
}

int8
GlRenderSystem::pushDepthFlags()
{
  glPushAttrib(GL_DEPTH_BUFFER_BIT);
  return true;
}

int8
GlRenderSystem::pushStencilFlags()
{
  glPushAttrib(GL_STENCIL_BUFFER_BIT);
  return true;
}

int8
GlRenderSystem::pushShapeFlags()
{
  glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_POLYGON_BIT);
  return true;
}

int8
GlRenderSystem::setColorWrite(int8 red, int8 green, int8 blue, int8 alpha)
{
  glColorMask(red, green, blue, alpha);
  return true;
}

int8
GlRenderSystem::setDepthWrite(int8 value)
{
  glDepthMask(value);
  return true;
}

int8
GlRenderSystem::setStencilWrite(uint32 mask)
{
  glStencilMask((GLuint)mask);
  return true;
}

int8
GlRenderSystem::setColor(Real const * rgba)
{
  if (sizeof(Real) == 4)
    glColor4fv(static_cast<GLfloat const *>((void const *)rgba));
  else
    glColor4dv(static_cast<GLdouble const *>((void const *)rgba));

  return true;
}

int8
GlRenderSystem::setColor(Real r, Real g, Real b, Real a)
{
  if (sizeof(Real) == 4)
    glColor4f((GLfloat)r, (GLfloat)g, (GLfloat)b, (GLfloat)a);
  else
    glColor4f(r, g, b, a);

  return true;
}

int8
GlRenderSystem::setClearColor(Real const * rgba)
{
  glClearColor((GLclampf)rgba[0], (GLclampf)rgba[1], (GLclampf)rgba[2], (GLclampf)rgba[3]);
  return true;
}

int8
GlRenderSystem::setClearColor(Real r, Real g, Real b, Real a)
{
  glClearColor((GLclampf)r, (GLclampf)g, (GLclampf)b, (GLclampf)a);
  return true;
}

int8
GlRenderSystem::setClearDepth(float64 value)
{
  glClearDepth((GLclampd)value);
  return true;
}

int8
GlRenderSystem::setClearStencil(int64 value)
{
  glClearStencil((GLint)value);
  return true;
}

int8
GlRenderSystem::setDepthTest(int32 test)
{
  if (test == DepthTest::ALWAYS_PASS)
  {
    glDisable(GL_DEPTH_TEST);
    glDepthFunc(GL_ALWAYS);  // just to be safe
  }
  else
  {
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GlRenderSystemInternal::depthTestToGLenum(test));
  }

  return true;
}

int8
GlRenderSystem::setCullFace(int32 cull)
{
  if (cull == CullFace::NONE)
    glDisable(GL_CULL_FACE);
  else
  {
    glEnable(GL_CULL_FACE);
    glCullFace((cull == CullFace::FRONT) ? GL_FRONT : GL_BACK);
  }

  return true;
}

int8
GlRenderSystem::setPolygonOffset(int8 enable, float64 offset)
{
  if (enable)
  {
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset((GLfloat)offset, (GLfloat)offset);
  }
  else
    glDisable(GL_POLYGON_OFFSET_FILL);

  return true;
}

int8
GlRenderSystem::setPolygonSmooth(int8 enable)
{
  if (enable)
    glEnable(GL_POLYGON_SMOOTH);
  else
  {
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_MULTISAMPLE);
  }

  return true;
}

int8
GlRenderSystem::setLineSmooth(int8 enable)
{
  if (enable)
    glEnable(GL_LINE_SMOOTH);
  else
  {
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_MULTISAMPLE);
  }

  return true;
}

int8
GlRenderSystem::setPointSmooth(int8 enable)
{
  if (enable)
    glEnable(GL_POINT_SMOOTH);
  else
  {
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_MULTISAMPLE);
  }

  return true;
}

int8
GlRenderSystem::setPointSize(Real size)
{
  glPointSize((GLfloat)size);
  return true;
}

int8
GlRenderSystem::popColorFlags()
{
  glPopAttrib();
  return true;
}

int8
GlRenderSystem::popDepthFlags()
{
  glPopAttrib();
  return true;
}

int8
GlRenderSystem::popStencilFlags()
{
  glPopAttrib();
  return true;
}

int8
GlRenderSystem::popShapeFlags()
{
  glPopAttrib();
  return true;
}

int8
GlRenderSystem::popState()
{
  glPopClientAttrib();
  glPopAttrib();
  if (!popTextures()) return false;
  if (!popShader()) return false;
  if (!popFramebuffer()) return false;

  return true;
}

int8
GlRenderSystem::finishAllOperations()
{
  glFlush();
  return true;
}

IRenderSystem * GlRenderSystemFactory::singleton          =  nullptr;
bool            GlRenderSystemFactory::singleton_created  =  false;

GlRenderSystemFactory::~GlRenderSystemFactory()
{
  delete singleton;
}

IRenderSystem *
GlRenderSystemFactory::createRenderSystem(char const * name)
{
  if (singleton)
  {
    THEA_ERROR << "GlRenderSystemFactory: Only one OpenGL rendersystem can be created per process";
    return nullptr;
  }

  try
  {
    singleton = new GlRenderSystem(name);
  }
  THEA_STANDARD_CATCH_BLOCKS(return nullptr;, ERROR, "%s", "Could not create new OpenGL rendersystem")

  singleton_created = true;
  return singleton;
}

int8
GlRenderSystemFactory::destroyRenderSystem(IRenderSystem * render_system)
{
  if (render_system != singleton)
  { THEA_ERROR << "GlRenderSystemFactory: Trying to destroy rendersystem not created with this factory"; return false; }

  return destroyAllRenderSystems();  // there's only one rendersystem
}

int8
GlRenderSystemFactory::destroyAllRenderSystems()
{
  delete singleton;
  singleton = nullptr;
  /* singleton_created = false; */  // Disabled until we figure out how to completely clean up the previous rendersystem.

  return true;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea
