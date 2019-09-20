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

    Texture * createTexture(char const * name_, int64 width, int64 height, int64 depth,
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

    VARArea * createVARArea(char const * name_, int64 num_bytes, VARArea::Usage usage, int8 gpu_memory = true);
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
    void setTexture(int32 texunit, Texture * texture);
    void popTextures();

    MatrixMode getMatrixMode() const;
    void setMatrixMode(MatrixMode mode);
    void pushMatrix();
    void pushViewMatrices();
    Matrix4 getMatrix(MatrixMode mode) const;
    void setMatrix(Matrix4 const & m);
    void setIdentityMatrix();
    void multMatrix(Matrix4 const & m);
    void popMatrix();
    void popViewMatrices();

    void beginIndexedPrimitives();
    void setVertexArray(VAR const * vertices);
    void setColorArray(VAR const * colors);
    void setTexCoordArray(int32 texunit, VAR const * texcoords);
    void setNormalArray(VAR const * normals);
    void setIndexArray(VAR const * indices);
    void sendIndices(Primitive primitive, int64 num_indices, uint8 const * indices);
    void sendIndices(Primitive primitive, int64 num_indices, uint16 const * indices);
    void sendIndices(Primitive primitive, int64 num_indices, uint32 const * indices);
    void sendSequentialIndices(Primitive primitive, int32 first_index, int32 num_indices);
    void sendIndicesFromArray(Primitive primitive, int64 offset, int64 num_indices);
    void endIndexedPrimitives();

    void beginPrimitive(Primitive primitive);
    void sendVertex(Vector2 const & vertex);
    void sendVertex(float32 x, float32 y);
    void sendVertex(float64 x, float64 y);
    void sendVertex(Vector3 const & vertex);
    void sendVertex(float32 x, float32 y, float32 z);
    void sendVertex(float64 x, float64 y, float64 z);
    void sendVertex(Vector4 const & vertex);
    void sendVertex(float32 x, float32 y, float32 z, float32 w);
    void sendVertex(float64 x, float64 y, float64 z, float64 w);
    void sendNormal(Vector3 const & normal);
    void sendNormal(float32 x, float32 y, float32 z);
    void sendNormal(float64 x, float64 y, float64 z);
    void sendTexCoord(int32 texunit, float32 texcoord);
    void sendTexCoord(int32 texunit, float64 texcoord);
    void sendTexCoord(int32 texunit, Vector2 const & texcoord);
    void sendTexCoord(int32 texunit, float32 x, float32 y);
    void sendTexCoord(int32 texunit, float64 x, float64 y);
    void sendTexCoord(int32 texunit, Vector3 const & texcoord);
    void sendTexCoord(int32 texunit, float32 x, float32 y, float32 z);
    void sendTexCoord(int32 texunit, float64 x, float64 y, float64 z);
    void endPrimitive();

    void pushState();
    void pushColorFlags();
    void pushDepthFlags();
    void pushStencilFlags();
    void pushShapeFlags();
    void setColorWrite(int8 red, int8 green, int8 blue, int8 alpha);
    void setDepthWrite(int8 value);
    void setStencilWrite(uint32 mask);
    void setColor(ColorRGB const & value);
    void setColor(ColorRGBA const & value);
    void setColorClearValue(ColorRGB const & value);
    void setColorClearValue(ColorRGBA const & value);
    void setDepthClearValue(Real value);
    void setStencilClearValue(int32 value);
    void clear();
    void clear(int8 color, int8 depth, int8 stencil);
    void setDepthTest(DepthTest test);
    void setCullFace(CullFace cull);
    void setPolygonOffset(int8 enable, float64 offset = 1);
    void setPolygonSmooth(int8 enable);
    void setLineSmooth(int8 enable);
    void setPointSmooth(int8 enable);
    void setPointSize(float64 size = 1);
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

    typedef Stack<GLFramebuffer *>  FramebufferStack;
    typedef Stack<BufferState>      BufferStack;
    typedef Stack<GLShader *>       ShaderStack;

    typedef UnorderedSet<GLFramebuffer *>  FramebufferSet;
    typedef UnorderedSet<GLTexture *>      TextureSet;
    typedef UnorderedSet<GLShader *>       ShaderSet;
    typedef UnorderedSet<GLVARArea *>      VARAreaSet;

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
