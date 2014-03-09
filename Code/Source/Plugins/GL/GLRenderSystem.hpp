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

#ifndef __Thea_Graphics_GL_GLRenderSystem_hpp__
#define __Thea_Graphics_GL_GLRenderSystem_hpp__

#include "../../Graphics/RenderSystem.hpp"
#include "../../Stack.hpp"
#include "../../UnorderedSet.hpp"
#include "GLCommon.hpp"
#include "GLFramebuffer.hpp"
#include "GLShader.hpp"
#include "GLVAR.hpp"

namespace Thea {
namespace Graphics {
namespace GL {

/** An OpenGL rendersystem. <b>Not threadsafe</b>. */
class THEA_GL_DLL_LOCAL GLRenderSystem : public RenderSystem
{
  public:
    /** Constructor. Creates a headless context for rendering if no rendering context currently exists. */
    GLRenderSystem(char const * name_);

    /** Destructor. */
    ~GLRenderSystem();

    char const * getName() const { return name.c_str(); }

    char const * describeSystem() const;

    Framebuffer * createFramebuffer(char const * name_);
    void destroyFramebuffer(Framebuffer * framebuffer);

    Shader * createShader(char const * name_);
    void destroyShader(Shader * shader);

    Texture * createTexture(char const * name_, int width, int height, int depth,
                            Texture::Format const * desired_format, Texture::Dimension dimension = Texture::Dimension::DIM_2D,
                            Texture::Options const & options = Texture::Options::defaults());

    Texture * createTexture(char const * name_, AbstractImage const & image,
                            Texture::Format const * desired_format = Texture::Format::AUTO(),
                            Texture::Dimension dimension = Texture::Dimension::DIM_2D,
                            Texture::Options const & options = Texture::Options::defaults());

    Texture * createTexture(char const * name_, AbstractImage const * images[6],
                            Texture::Format const * desired_format = Texture::Format::AUTO(),
                            Texture::Options const & options = Texture::Options::defaults());

    void destroyTexture(Texture * texture);

    VARArea * createVARArea(char const * name_, long num_bytes, VARArea::Usage usage, bool gpu_memory = true);
    void destroyVARArea(VARArea * area);

    void pushFramebuffer();
    void setFramebuffer(Framebuffer * framebuffer);
    Framebuffer const * getFramebuffer() const;
    Framebuffer * getFramebuffer();
    void popFramebuffer();

    void pushShader();
    void setShader(Shader * shader);
    Shader const * getShader() const;
    Shader * getShader();
    void popShader();

    void pushTextures();
    void setTexture(int texunit, Texture * texture);
    void popTextures();

    MatrixMode getMatrixMode() const;
    void setMatrixMode(MatrixMode mode);
    void pushMatrix();
    Matrix4 getMatrix(MatrixMode mode) const;
    void setMatrix(Matrix4 const & m);
    void setIdentityMatrix();
    void multMatrix(Matrix4 const & m);
    void popMatrix();

    void beginIndexedPrimitives();
    void setVertexArray(VAR const * vertices);
    void setColorArray(VAR const * colors);
    void setTexCoordArray(int texunit, VAR const * texcoords);
    void setNormalArray(VAR const * normals);
    void setIndexArray(VAR const * indices);
    void sendIndices(Primitive primitive, long num_indices, uint8 const * indices);
    void sendIndices(Primitive primitive, long num_indices, uint16 const * indices);
    void sendIndices(Primitive primitive, long num_indices, uint32 const * indices);
    void sendSequentialIndices(Primitive primitive, int first_index, int num_indices);
    void sendIndicesFromArray(Primitive primitive, long offset, long num_indices);
    void endIndexedPrimitives();

    void beginPrimitive(Primitive primitive);
    void sendVertex(Vector2 const & vertex);
    void sendVertex(float x, float y);
    void sendVertex(double x, double y);
    void sendVertex(Vector3 const & vertex);
    void sendVertex(float x, float y, float z);
    void sendVertex(double x, double y, double z);
    void sendVertex(Vector4 const & vertex);
    void sendVertex(float x, float y, float z, float w);
    void sendVertex(double x, double y, double z, double w);
    void sendNormal(Vector3 const & normal);
    void sendNormal(float x, float y, float z);
    void sendNormal(double x, double y, double z);
    void sendTexCoord(int texunit, float texcoord);
    void sendTexCoord(int texunit, double texcoord);
    void sendTexCoord(int texunit, Vector2 const & texcoord);
    void sendTexCoord(int texunit, float x, float y);
    void sendTexCoord(int texunit, double x, double y);
    void sendTexCoord(int texunit, Vector3 const & texcoord);
    void sendTexCoord(int texunit, float x, float y, float z);
    void sendTexCoord(int texunit, double x, double y, double z);
    void endPrimitive();

    void pushState();
    void pushColorFlags();
    void pushDepthFlags();
    void pushStencilFlags();
    void pushShapeFlags();
    void setColorWrite(bool red, bool green, bool blue, bool alpha);
    void setDepthWrite(bool value);
    void setStencilWrite(uint32 mask);
    void setColor(ColorRGB const & value);
    void setColor(ColorRGBA const & value);
    void setColorClearValue(ColorRGB const & value);
    void setColorClearValue(ColorRGBA const & value);
    void setDepthClearValue(Real value);
    void setStencilClearValue(int value);
    void clear();
    void clear(bool color, bool depth, bool stencil);
    void setDepthTest(DepthTest test);
    void setCullFace(CullFace cull);
    void setPolygonOffset(bool enable, Real offset = 1);
    void setPointSize(Real size = 1);
    void popColorFlags();
    void popDepthFlags();
    void popStencilFlags();
    void popShapeFlags();
    void popState();

    void finishAllOperations();

  private:
    /** Stores the state of the current set of GPU buffers. */
    struct BufferState
    {
      GLVARArea * vertex_area;  ///< Vertex buffer.
      GLVARArea * index_area;  ///< Index buffer.
      GLVAR index_var;  ///< Index array.

      /** Constructor. */
      BufferState(GLVARArea * vertex_area_ = NULL, GLVARArea * index_area_ = NULL, GLVAR index_var_ = GLVAR())
      : vertex_area(vertex_area_), index_area(index_area_), index_var(index_var_)
      {}

    }; // struct BufferState

    typedef TheaStack<GLFramebuffer *>  FramebufferStack;
    typedef TheaStack<BufferState>      BufferStack;
    typedef TheaStack<GLShader *>       ShaderStack;

    typedef TheaUnorderedSet<GLFramebuffer *>  FramebufferSet;
    typedef TheaUnorderedSet<GLTexture *>      TextureSet;
    typedef TheaUnorderedSet<GLShader *>       ShaderSet;
    typedef TheaUnorderedSet<GLVARArea *>      VARAreaSet;

    /**
     * Set the current vertex area to match the specified VAR. The vertex area used within a single
     * beginIndexedPrimitives()/endIndexedPrimitives() block must remain constant.
     */
    void setVertexAreaFromVAR(GLVAR const & v);

    /**
     * Set the current index area to match the specified VAR. The vertex area used within a single
     * beginIndexedPrimitives()/endIndexedPrimitives() block must remain constant.
     */
    void setIndexAreaFromVAR(GLVAR const & v);

    std::string name;
    mutable std::string desc;

    GLFramebuffer * current_framebuffer;
    FramebufferStack framebuffer_stack;

    BufferState current_buffer_state;
    BufferStack buffer_stack;

    GLShader * current_shader;
    ShaderStack shader_stack;

    FramebufferSet  created_framebuffers;
    TextureSet      created_textures;
    ShaderSet       created_shaders;
    VARAreaSet      created_varareas;

}; // class GLRenderSystem

/** Factory for creating OpenGL rendersystems. */
class THEA_GL_DLL_LOCAL GLRenderSystemFactory : public RenderSystemFactory
{
  public:
    /** Destructor. */
    ~GLRenderSystemFactory();

    RenderSystem * createRenderSystem(char const * name);
    void destroyRenderSystem(RenderSystem * render_system);

    /** Destroy all rendersystems created with this factory. */
    void destroyAllRenderSystems();

  private:
    static RenderSystem * singleton;
    static bool singleton_created;
};

} // namespace GL
} // namespace Graphics
} // namespace Thea

#endif
