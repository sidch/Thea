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

#ifndef __Thea_Graphics_RenderSystem_hpp__
#define __Thea_Graphics_RenderSystem_hpp__

#include "../Common.hpp"
#include "../Image.hpp"
#include "../Colors.hpp"
#include "../Map.hpp"
#include "../NamedObject.hpp"
#include "Camera.hpp"
#include "Framebuffer.hpp"
#include "Texture.hpp"
#include "Shader.hpp"
#include "VAR.hpp"
#include "VARArea.hpp"

namespace Thea {
namespace Graphics {

/**
 * An interface for a rendersystem. Should be easily implementable in both OpenGL and Direct3D.
 *
 * To create an instance of a RenderSystem, one typically loads the plugin for the relevant implementation and calls
 * RenderSystemFactory::createRenderSystem().
 *
 * If no rendering context is available when a rendersystem is constructed, the new rendersystem will be set up for offscreen
 * rendering if possible. In this case, you must explicitly create and attach a framebuffer to the rendersystem before you can
 * call any drawing functions.
 */
class THEA_API RenderSystem : public AbstractNamedObject
{
  public:
    /** Basic drawing primitives (enum class). */
    struct THEA_API Primitive
    {
      /** Supported values. */
      enum Value
      {
        POINTS,          ///< Set of points.
        LINES,           ///< Set of line segments.
        LINE_STRIP,      ///< Sequence of connected line segments.
        LINE_LOOP,       ///< Loop of line segments.
        TRIANGLES,       ///< Set of triangles.
        TRIANGLE_STRIP,  ///< Triangle strip.
        TRIANGLE_FAN,    ///< Triangle fan.
        QUADS,           ///< Set of quads.
        QUAD_STRIP,      ///< Quad strip.
        POLYGON          ///< A single polygon with an arbitrary number of edges.
      };

      THEA_ENUM_CLASS_BODY(Primitive)
    };

    /** Matrix-based transformation modes (enum class). */
    struct THEA_API MatrixMode
    {
      /** Supported values. */
      enum Value
      {
        MODELVIEW,   ///< Model-view matrix (cf. GL_MODELVIEW).
        PROJECTION,  ///< Projection matrix (cf. GL_PROJECTION).
        TEXTURE,     ///< Matrix to transform texture coordinates.
        COLOR        ///< Matrix to transform colors.
      };

      THEA_ENUM_CLASS_BODY(MatrixMode)
    };

    /** Depth tests (enum class). */
    struct THEA_API DepthTest
    {
      /** Supported values. */
      enum Value
      {
        GREATER,      ///< Accept if value is strictly greater than the threshold.
        LESS,         ///< Accept if value is strictly less than the threshold.
        GEQUAL,       ///< Accept if value is greater than or equal to the threshold.
        LEQUAL,       ///< Accept if value is less than or equal to the threshold.
        NOTEQUAL,     ///< Accept if value is not equal to the threshold.
        EQUAL,        ///< Accept if value is equal to the threshold.
        ALWAYS_PASS,  ///< Always accept.
        NEVER_PASS    ///< Never accept.
      };

      THEA_ENUM_CLASS_BODY(DepthTest)
    };

    /** Faces to be culled (enum class). */
    struct THEA_API CullFace
    {
      /** Supported values. */
      enum Value
      {
        NONE,   ///< No front/back culling.
        FRONT,  ///< Cull front faces.
        BACK    ///< Cull back faces.
      };

      THEA_ENUM_CLASS_BODY(CullFace)
    };

    /** Destructor. Frees all resources (textures, shaders, framebuffers, VAR areas etc) created using this rendersystem. */
    virtual ~RenderSystem() {}

    /**
     * Get a string describing the render system. The string is guaranteed to be valid only till the next operation on the
     * object.
     */
    virtual char const * describeSystem() const = 0;

    /** Create a new, blank framebuffer with nothing attached. The framebuffer must be destroyed using destroyFramebuffer(). */
    virtual Framebuffer * createFramebuffer(char const * name) = 0;

    /** Destroy a framebuffer created with createFramebuffer(). */
    virtual void destroyFramebuffer(Framebuffer * framebuffer) = 0;

    /** Create a new, uninitialized shader. The shader must be destroyed using destroyShader(). */
    virtual Shader * createShader(char const * name) = 0;

    /** Destroy a shader created with createShader(). */
    virtual void destroyShader(Shader * shader) = 0;

    /** Create an empty texture of the specified format and size. The texture must be destroyed using destroyTexture(). */
    virtual Texture * createTexture(char const * name, int width, int height, int depth,
                                    Texture::Format const * desired_format, Texture::Dimension dimension,
                                    Texture::Options const & options = Texture::Options::defaults()) = 0;

    /**
     * Create a texture from a pixel buffer. The dimension argument <em>cannot</em> be DIM_CUBE_MAP. The texture must be
     * destroyed using destroyTexture().
     */
    virtual Texture * createTexture(char const * name, AbstractImage const & image,
                                    Texture::Format const * desired_format = Texture::Format::AUTO(),
                                    Texture::Dimension dimension = Texture::Dimension::DIM_2D,
                                    Texture::Options const & options = Texture::Options::defaults()) = 0;

    /**
     * Create a cube-map from six pixel buffers, representing 2D images of identical format and size. The texture must be
     * destroyed using destroyTexture().
     */
    virtual Texture * createTexture(char const * name, AbstractImage const * images[6],
                                    Texture::Format const * desired_format = Texture::Format::AUTO(),
                                    Texture::Options const & options = Texture::Options::defaults()) = 0;

    /** Destroy a texture created with createTexture(). */
    virtual void destroyTexture(Texture * texture) = 0;

    /**
     * Create a new, uninitialized area for storing vertex/normal/texture-coordinate/index arrays. The area must be destroyed
     * using destroyVARArea().
     */
    virtual VARArea * createVARArea(char const * name, long num_bytes, VARArea::Usage usage, bool gpu_memory = true) = 0;

    /** Destroy a memory area created with createVARArea(). */
    virtual void destroyVARArea(VARArea * area) = 0;

    /** Save the current framebuffer by pushing it onto the stack. */
    virtual void pushFramebuffer() = 0;

    /** Set the current framebuffer. The drawing viewport is set to the framebuffer's entire area. */
    virtual void setFramebuffer(Framebuffer * framebuffer) = 0;

    /**
     * Get the current framebuffer (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual framebuffer currently in use!
     */
    virtual Framebuffer const * getFramebuffer() const = 0;

    /**
     * Get the current framebuffer (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual framebuffer currently in use!
     */
    virtual Framebuffer * getFramebuffer() = 0;

    /** Restore the last saved framebuffer from the stack. */
    virtual void popFramebuffer() = 0;

    /** Save the current shader by pushing it onto the stack. */
    virtual void pushShader() = 0;

    /** Set the current shader. */
    virtual void setShader(Shader * shader) = 0;

    /**
     * Get the current shader (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual shader currently in use!
     */
    virtual Shader const * getShader() const = 0;

    /**
     * Get the current shader (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual shader currently in use!
     */
    virtual Shader * getShader() = 0;

    /** Restore the last saved shader from the stack. */
    virtual void popShader() = 0;

    /** Save all texture bindings and related state on the stack. */
    virtual void pushTextures() = 0;

    /**
     * Bind a texture to a texture unit. Passing a null pointer disables the unit. The function signals an error if the render
     * system does not support multitexturing and the specified texture unit is non-zero.
     */
    virtual void setTexture(int texunit, Texture * texture) = 0;

    /** Restore the last saved set of texture bindings and related state from the stack. */
    virtual void popTextures() = 0;

    /** Get the current matrix mode. */
    virtual MatrixMode getMatrixMode() const = 0;

    /** Set the current matrix mode. */
    virtual void setMatrixMode(MatrixMode mode) = 0;

    /** Save the matrix of the current matrix mode by pushing it onto the stack. */
    virtual void pushMatrix() = 0;

    /** Get the current matrix of the specified matrix mode. */
    virtual Matrix4 getMatrix(MatrixMode mode) const = 0;

    /** Set the matrix of the current matrix mode. */
    virtual void setMatrix(Matrix4 const & m) = 0;

    /**
     * Set the current viewing matrices from a camera. Subsequent changes to the original camera will <b>NOT</b> affect the
     * rendersystem.
     */
    virtual void setCamera(Camera const & camera)
    {
      setMatrixMode(MatrixMode::PROJECTION);
      setMatrix(camera.getProjectionTransform());

      setMatrixMode(MatrixMode::MODELVIEW);
      setMatrix(camera.getWorldToCameraTransform().toHomMatrix());
    }

    /** Set the matrix of the current matrix mode to the identity transformation. */
    virtual void setIdentityMatrix() = 0;

    /**
     * Post-multiply the matrix of the current matrix mode by the given matrix. The given transformation will be applied
     * <b>before</b> any previous transformations.
     */
    virtual void multMatrix(Matrix4 const & m) = 0;

    /** Restore the last saved matrix of the current matrix mode from the stack. */
    virtual void popMatrix() = 0;

    /** Prepare to draw a set of indexed primitives, after saving the current array states on the stack. */
    virtual void beginIndexedPrimitives() = 0;

    /** Set the current vertex array. Passing a null VAR unbinds all vertex data. */
    virtual void setVertexArray(VAR const * vertices) = 0;

    /** Set the current color array. Passing a null VAR unbinds all color data. */
    virtual void setColorArray(VAR const * colors) = 0;

    /** Set the current texture coordinate array. Passing a null VAR unbinds all texture coordinate data. */
    virtual void setTexCoordArray(int texunit, VAR const * texcoords) = 0;

    /** Set the current normal array. Passing a null VAR unbinds all normal data. */
    virtual void setNormalArray(VAR const * normals) = 0;

    /** Set the current index array. Passing a null VAR unbinds all index data. */
    virtual void setIndexArray(VAR const * indices) = 0;

    /** Draw a set of primitives by sending a set of 8-bit indices to the rendersystem. */
    virtual void sendIndices(Primitive primitive, long num_indices, uint8 const * indices) = 0;

    /** Draw a set of primitives by sending a set of 16-bit indices to the rendersystem. */
    virtual void sendIndices(Primitive primitive, long num_indices, uint16 const * indices) = 0;

    /** Draw a set of primitives by sending a set of 32-bit indices to the rendersystem. */
    virtual void sendIndices(Primitive primitive, long num_indices, uint32 const * indices) = 0;

    /**
     * Draw a set of primitives by sending \a num_indices consecutive indices, starting from \a first_index, to the
     * rendersystem.
     */
    virtual void sendSequentialIndices(Primitive primitive, int first_index, int num_indices) = 0;

    /**
     * Draw a set of primitives by sending indices from the current index array to the rendersystem. You must call
     * setIndexArray() to set a valid index array before you call this function.
     *
     * @param primitive The type of primitive to draw.
     * @param offset The number of initial indices to skip.
     * @param num_indices The number of indices to send to the rendersystem.
     *
     * @see setIndexArray()
     */
    virtual void sendIndicesFromArray(Primitive primitive, long offset, long num_indices) = 0;

    /** Finish drawing the current set of indexed primitives and restore the last saved array states from the stack. */
    virtual void endIndexedPrimitives() = 0;

    /** Start drawing a primitive of the given type. Must be matched with endPrimitive(). */
    virtual void beginPrimitive(Primitive primitive) = 0;

    /** Send a 2-vertex to the rendersystem. */
    virtual void sendVertex(Vector2 const & vertex) = 0;

    /** Send a 2-vertex to the rendersystem. */
    virtual void sendVertex(float x, float y) { sendVertex(Vector2(x, y)); }

    /** Send a 2-vertex to the rendersystem. */
    virtual void sendVertex(double x, double y) { sendVertex(Vector2((float)x, (float)y)); }

    /** Send a 3-vertex to the rendersystem. */
    virtual void sendVertex(Vector3 const & vertex) = 0;

    /** Send a 3-vertex to the rendersystem. */
    virtual void sendVertex(float x, float y, float z) { sendVertex(Vector3(x, y, z)); }

    /** Send a 3-vertex to the rendersystem. */
    virtual void sendVertex(double x, double y, double z) { sendVertex(Vector3((float)x, (float)y, (float)z)); };

    /** Send a 4-vertex to the rendersystem. */
    virtual void sendVertex(Vector4 const & vertex) = 0;

    /** Send a 4-vertex to the rendersystem. */
    virtual void sendVertex(float x, float y, float z, float w) { sendVertex(Vector4(x, y, z, w)); };

    /** Send a 4-vertex to the rendersystem. */
    virtual void sendVertex(double x, double y, double z, double w)
    { sendVertex(Vector4((float)x, (float)y, (float)z, (float)w)); };

    /** Send a normal to the rendersystem. */
    virtual void sendNormal(Vector3 const & normal) = 0;

    /** Send a normal to the rendersystem. */
    virtual void sendNormal(float x, float y, float z) { sendNormal(Vector3(x, y, z)); }

    /** Send a normal to the rendersystem. */
    virtual void sendNormal(double x, double y, double z) { sendNormal(Vector3((float)x, (float)y, (float)z)); };

    /** Send a floating-point texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, float texcoord) = 0;

    /** Send a floating-point texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, double texcoord) { sendTexCoord(texunit, (float)texcoord); }

    /** Send a 2-component texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, Vector2 const & texcoord) = 0;

    /** Send a 2-component texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, float x, float y) { sendTexCoord(texunit, Vector2(x, y)); }

    /** Send a 2-component texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, double x, double y) { sendTexCoord(texunit, Vector2((float)x, (float)y)); }

    /** Send a 3-component texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, Vector3 const & texcoord) = 0;

    /** Send a 3-component texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, float x, float y, float z) { sendTexCoord(texunit, Vector3(x, y, z)); }

    /** Send a 3-component texture coordinate to the rendersystem. */
    virtual void sendTexCoord(int texunit, double x, double y, double z)
    { sendTexCoord(texunit, Vector3((float)x, (float)y, (float)z)); }

    /** Finish drawing the primitive started by beginPrimitive(). */
    virtual void endPrimitive() = 0;

    /** Save the current state of the rendersystem. This includes all attributes saved by other <code>push...</code> calls. */
    virtual void pushState() = 0;

    /** Save the current set of color flags (color write, clear value etc) by pushing it onto the stack. */
    virtual void pushColorFlags() = 0;

    /** Save the current set of depth flags (depth write, depth test, clear value etc) by pushing it onto the stack. */
    virtual void pushDepthFlags() = 0;

    /** Save the current set of stencil flags (stencil write, clear value etc) by pushing it onto the stack. */
    virtual void pushStencilFlags() = 0;

    /**
     * Save current set of parameters for rasterizing shapes (line width, smoothing flag for lines/points, polygon offsets
     * etc).
     */
    virtual void pushShapeFlags() = 0;

    /** Set the color write state. */
    virtual void setColorWrite(bool red, bool green, bool blue, bool alpha) = 0;

    /** Set the depth write state. */
    virtual void setDepthWrite(bool value) = 0;

    /** Set the stencil write mask. */
    virtual void setStencilWrite(uint32 mask) = 0;

    /** Set the current drawing color. */
    virtual void setColor(ColorRGB const & value) = 0;

    /** Set the current drawing color. */
    virtual void setColor(ColorRGBA const & value) = 0;

    /** Set the value to clear the color buffer with. */
    virtual void setColorClearValue(ColorRGB const & value) = 0;

    /** Set the value to clear the color buffer with. */
    virtual void setColorClearValue(ColorRGBA const & value) = 0;

    /** Set the value to clear the depth buffer with. */
    virtual void setDepthClearValue(Real value) = 0;

    /** Set the value to clear the stencil buffer with. */
    virtual void setStencilClearValue(int value) = 0;

    /** Clear color, depth and stencil buffers. */
    virtual void clear() = 0;

    /** Clear the selected buffers. */
    virtual void clear(bool color, bool depth, bool stencil) = 0;

    /** Set the depth test. */
    virtual void setDepthTest(DepthTest test) = 0;

    /** Set the faces to be culled. */
    virtual void setCullFace(CullFace cull) = 0;

    /**
     * Set the depth offset, if any, to be applied to polygon faces. The supplied value is scaled by an implementation-specific
     * offset.
     */
    virtual void setPolygonOffset(bool enable, Real offset = 1) = 0;

    /** Set the size (diameter) of rasterized points. This may be ignored if shaders are being used. */
    virtual void setPointSize(Real size = 1) = 0;

    /** Restore the last saved set of color flags from the stack. */
    virtual void popColorFlags() = 0;

    /** Restore the last saved set of depth flags from the stack. */
    virtual void popDepthFlags() = 0;

    /** Restore the last saved set of stencil flags from the stack. */
    virtual void popStencilFlags() = 0;

    /** Restore the last saved set of shape flags from the stack. */
    virtual void popShapeFlags() = 0;

    /** Restore the last saved rendersystem state from the stack. */
    virtual void popState() = 0;

    /** The function returns only when the rendersystem has completely executed all previously-issued drawing calls. */
    virtual void finishAllOperations() = 0;

}; // class RenderSystem

/** An interface for a rendersystem factory. Should be implemented and registered by each actual rendersystem. */
class THEA_API RenderSystemFactory
{
  public:
    /** Destructor. */
    virtual ~RenderSystemFactory() {}

    /**
     * Create a rendersystem with the given name. The rendersystem must be destroyed using destroyRenderSystem(). If no
     * rendering context is available, the new rendersystem will be set up for offscreen rendering if possible. In this case,
     * you must explictly create and attach a framebuffer to the rendersystem before you can call any drawing functions.
     */
    virtual RenderSystem * createRenderSystem(char const * name) = 0;

    /** Destroy a rendersystem created with createRenderSystem(). */
    virtual void destroyRenderSystem(RenderSystem * render_system) = 0;

}; // class RenderSystemFactory

/** Manages available rendersystem factories. */
class THEA_API RenderSystemManager
{
  public:
    /**
     * Install a factory for a particular rendersystem type. The factory pointer should not be null.
     *
     * @return True if the factory was successfully installed, false if a factory of the specified type (with case-insensitive
     *   matching) is already installed.
     */
    bool installFactory(std::string const & type, RenderSystemFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for rendersystem of a given type. An error is thrown if no such factory has been installed. */
    RenderSystemFactory * getFactory(std::string const & type);

  private:
    typedef TheaMap<std::string, RenderSystemFactory *> FactoryMap;  ///< Maps rendersystem types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each rendersystem type.

}; // class RenderSystemManager

} // namespace Graphics
} // namespace Thea

#endif
