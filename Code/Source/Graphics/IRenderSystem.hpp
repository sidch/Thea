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

#ifndef __Thea_Graphics_IRenderSystem_hpp__
#define __Thea_Graphics_IRenderSystem_hpp__

#include "../Common.hpp"
#include "../IDenseMatrix.hpp"
#include "../IImage.hpp"
#include "../Map.hpp"
#include "../NamedObject.hpp"
#include "IBuffer.hpp"
#include "IBufferPool.hpp"
#include "IFramebuffer.hpp"
#include "IShader.hpp"
#include "ITexture.hpp"

namespace Thea {
namespace Graphics {

/**
 * Interface for a rendersystem. Should be easily implementable in both OpenGL and Direct3D.
 *
 * To create an instance of an IRenderSystem, one typically loads the plugin for the relevant implementation and calls
 * IRenderSystemFactory::createRenderSystem().
 *
 * If no rendering context is available when a rendersystem is constructed, the new rendersystem will be set up for offscreen
 * rendering if possible. In this case, you must explicitly create and attach a framebuffer to the rendersystem before you can
 * call any drawing functions.
 */
class THEA_API IRenderSystem : public INamedObject
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
        MODELVIEW,   ///< Model-view matrix (see GL_MODELVIEW).
        PROJECTION,  ///< Projection matrix (see GL_PROJECTION).
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

    /** Destructor. Frees all resources (textures, shaders, framebuffers, buffer areas etc) created using this rendersystem. */
    virtual ~IRenderSystem() = 0;

    /**
     * Get a string describing the render system. The string is guaranteed to be valid only till the next operation on the
     * object.
     */
    virtual char const * THEA_ICALL describeSystem() const = 0;

    /**
     * Get a string describing the message corresponding to the last error. Rendersystems may choose not to implement this
     * functionality and make this method a no-op, e.g. in a multithreaded rendersystem where it would be difficult to
     * synchronize.
     */
    virtual char const * THEA_ICALL getErrorString() const = 0;

    /**
     * Create a new, blank framebuffer with nothing attached. The framebuffer must be destroyed using destroyFramebuffer().
     * Returns a null pointer on error.
     */
    virtual IFramebuffer * THEA_ICALL createFramebuffer(char const * name) = 0;

    /**
     * Destroy a framebuffer created with createFramebuffer().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL destroyFramebuffer(IFramebuffer * framebuffer) = 0;

    /**
     * Create a new, uninitialized shader. The shader must be destroyed using destroyShader(). Returns a null pointer on error.
     */
    virtual IShader * THEA_ICALL createShader(char const * name) = 0;

    /**
     * Destroy a shader created with createShader().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL destroyShader(IShader * shader) = 0;

    /**
     * Create an empty texture of the specified format and size. The texture must be destroyed using destroyTexture(). Returns
     * a null pointer on error.
     */
    virtual ITexture * THEA_ICALL createTexture(char const * name, int64 width, int64 height, int64 depth,
                                                ITexture::Format const * desired_format,
                                                int32 dimension = ITexture::Dimension::DIM_2D,
                                                ITexture::Options const * options = nullptr) = 0;

    /**
     * Create a texture from a pixel buffer. The dimension argument <em>cannot</em> be Dimension::DIM_CUBE_MAP. The texture must
     * be destroyed using destroyTexture(). Returns a null pointer on error.
     */
    virtual ITexture * THEA_ICALL createTexture(char const * name, IImage const * image,
                                                ITexture::Format const * desired_format = nullptr,
                                                int32 dimension = ITexture::Dimension::DIM_2D,
                                                ITexture::Options const * options = nullptr) = 0;

    /**
     * Create a cube-map from six pixel buffers, representing 2D images of identical format and size. The texture must be
     * destroyed using destroyTexture(). Returns a null pointer on error.
     */
    virtual ITexture * THEA_ICALL createTexture(char const * name, IImage const * images[6],
                                                ITexture::Format const * desired_format = nullptr,
                                                ITexture::Options const * options = nullptr) = 0;

    /**
     * Destroy a texture created with createTexture().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL destroyTexture(ITexture * texture) = 0;

    /**
     * Create a new, uninitialized memory pool for storing vertex/normal/texture-coordinate/index buffers. The pool must be
     * destroyed using destroyBufferPool(). \a usage should be a value from IBufferPool::Usage. Returns a null pointer on error.
     */
    virtual IBufferPool * THEA_ICALL createBufferPool(char const * name, int64 num_bytes, int32 usage,
                                                      int8 gpu_memory = true) = 0;

    /**
     * Destroy a memory area created with createBufferPool().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL destroyBufferPool(IBufferPool * pool) = 0;

    /** Save the current framebuffer by pushing it onto the stack. */
    virtual int8 THEA_ICALL pushFramebuffer() = 0;

    /**
     * Set the current framebuffer. The drawing viewport is set to the framebuffer's entire area.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setFramebuffer(IFramebuffer * framebuffer) = 0;

    /**
     * Get the current framebuffer (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual framebuffer currently in use!
     */
    virtual IFramebuffer const * THEA_ICALL getFramebuffer() const = 0;

    /**
     * Get the current framebuffer (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual framebuffer currently in use!
     */
    virtual IFramebuffer * THEA_ICALL getFramebuffer() = 0;

    /**
     * Restore the last saved framebuffer from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popFramebuffer() = 0;

    /**
     * Save the current shader by pushing it onto the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushShader() = 0;

    /**
     * Set the current shader.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setShader(IShader * shader) = 0;

    /**
     * Get the current shader (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual shader currently in use!
     */
    virtual IShader const * THEA_ICALL getShader() const = 0;

    /**
     * Get the current shader (may be null). If the rendersystem has been externally modified via direct API (e.g.
     * OpenGL/Direct3D) calls, this may not be the actual shader currently in use!
     */
    virtual IShader * THEA_ICALL getShader() = 0;

    /**
     * Restore the last saved shader from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popShader() = 0;

    /**
     * Save all texture bindings and related state on the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushTextures() = 0;

    /**
     * Bind a texture to a texture unit. Passing a null pointer disables the unit. The function signals an error if the render
     * system does not support multitexturing and the specified texture unit is non-zero.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setTexture(int32 texunit, ITexture * texture) = 0;

    /**
     * Restore the last saved set of texture bindings and related state from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popTextures() = 0;

    /** Get the current matrix mode, as a value from the MatrixMode enum. */
    virtual int32 THEA_ICALL getMatrixMode() const = 0;

    /**
     * Set the current matrix mode, as a value from the MatrixMode enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setMatrixMode(int32 mode) = 0;

    /**
     * Save the matrix of the current matrix mode by pushing it onto the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushMatrix() = 0;

    /**
     * Save modelview and projection matrices, ensuring the current matrix mode is restored at the end.
     *
     * @return True on success, false on error.
     *
     * @see popViewMatrices()
     */
    virtual int8 THEA_ICALL pushViewMatrices() = 0;

    /**
     * Get the current matrix of a matrix mode, specified as a MatrixMode enum value. The result matrix must be a 4x4 matrix of
     * Real values.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL getMatrix(int32 mode, IDenseMatrix<Real> * m) const = 0;

    /**
     * Set the matrix of the current matrix mode. \a m must be a 4x4 dense matrix.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setMatrix(IDenseMatrix<Real> const * m) = 0;

    /**
     * Set the matrix of the current matrix mode to the identity transformation.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setIdentityMatrix() = 0;

    /**
     * Post-multiply the matrix of the current matrix mode by the given matrix. The given transformation will be applied
     * <b>before</b> any previous transformations. \a m must be a 4x4 dense matrix.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL multMatrix(IDenseMatrix<Real> const * m) = 0;

    /**
     * Restore the last saved matrix of the current matrix mode from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popMatrix() = 0;

    /**
     * Restore modelview and projection matrices from the respective stacks, ensuring the current matrix mode is restored at the
     * end.
     *
     * @return True on success, false on error.
     *
     * @see pushViewMatrices().
     */
    virtual int8 THEA_ICALL popViewMatrices() = 0;

    /**
     * Prepare to draw a set of indexed primitives, after saving the current buffer states on the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL beginIndexedPrimitives() = 0;

    /**
     * Set the current vertex buffer. Passing a null buffer unbinds all vertex data.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setVertexBuffer(IBuffer const * vertices) = 0;

    /**
     * Set the current color buffer. Passing a null buffer unbinds all color data.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setColorBuffer(IBuffer const * colors) = 0;

    /**
     * Set the current texture coordinate buffer. Passing a null buffer unbinds all texture coordinate data.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setTexCoordBuffer(int32 texunit, IBuffer const * texcoords) = 0;

    /**
     * Set the current normal buffer. Passing a null buffer unbinds all normal data.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setNormalBuffer(IBuffer const * normals) = 0;

    /**
     * Set the current index buffer. Passing a null buffer unbinds all index data.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setIndexBuffer(IBuffer const * indices) = 0;

    /**
     * Draw a set of primitives by sending a set of 8-bit indices to the rendersystem. \a primitive should be a value from the
     * Primitive enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendIndices(int32 primitive, int64 num_indices, uint8 const * indices) = 0;

    /**
     * Draw a set of primitives by sending a set of 16-bit indices to the rendersystem. \a primitive should be a value from the
     * Primitive enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendIndices(int32 primitive, int64 num_indices, uint16 const * indices) = 0;

    /**
     * Draw a set of primitives by sending a set of 32-bit indices to the rendersystem. \a primitive should be a value from the
     * Primitive enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendIndices(int32 primitive, int64 num_indices, uint32 const * indices) = 0;

    /**
     * Draw a set of primitives by sending \a num_indices consecutive indices, starting from \a first_index, to the
     * rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendSequentialIndices(int32 primitive, uint64 first_index, int64 num_indices) = 0;

    /**
     * Draw a set of primitives by sending indices from the current index buffer to the rendersystem. You must call
     * setIndexBuffer() to set a valid index buffer before you call this function.
     *
     * @param primitive The type of primitive to draw, as a value from the Primitive enum.
     * @param offset The number of initial indices to skip.
     * @param num_indices The number of indices to send to the rendersystem.
     *
     * @return True on success, false on error.
     *
     * @see setIndexBuffer()
     */
    virtual int8 THEA_ICALL sendIndicesFromBuffer(int32 primitive, uint64 offset, int64 num_indices) = 0;

    /**
     * Finish drawing the current set of indexed primitives and restore the last saved array states from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL endIndexedPrimitives() = 0;

    /**
     * Start drawing a primitive of the given type, specified as a value from the Primitive enum. Must be matched with
     * endPrimitive().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL beginPrimitive(int32 primitive) = 0;

    /**
     * Send a vertex to the rendersystem. \a dims must be 1, 2, 3 or 4.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendVertex(int32 dims, Real const * coords) = 0;

    /**
     * Send a 2-vertex to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendVertex(Real x, Real y) = 0;

    /**
     * Send a 3-vertex to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendVertex(Real x, Real y, Real z) = 0;

    /**
     * Send a 4-vertex to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendVertex(Real x, Real y, Real z, Real w) = 0;

    /**
     * Send a 3D normal to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendNormal(Real const * coords) = 0;

    /**
     * Send a 3D normal to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendNormal(Real x, Real y, Real z) = 0;

    /**
     * Send texture coordinates to the rendersystem. \a dims must be 1, 2 or 3.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendTexCoord(int32 texunit, int32 dims, Real const * coords) = 0;

    /**
     * Send a 1D texture coordinate to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendTexCoord(int32 texunit, Real x) = 0;

    /**
     * Send 2D texture coordinates to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendTexCoord(int32 texunit, Real x, Real y) = 0;

    /**
     * Send 3D texture coordinates to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendTexCoord(int32 texunit, Real x, Real y, Real z) = 0;

    /**
     * Send 4D texture coordinates to the rendersystem.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL sendTexCoord(int32 texunit, Real x, Real y, Real z, Real w) = 0;

    /**
     * Finish drawing the primitive started by beginPrimitive().
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL endPrimitive() = 0;

    /**
     * Save the current state of the rendersystem. This includes all attributes saved by other <code>push...</code> calls.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushState() = 0;

    /**
     * Save the current set of color flags (color write, clear value etc) by pushing it onto the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushColorFlags() = 0;

    /**
     * Save the current set of depth flags (depth write, depth test, clear value etc) by pushing it onto the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushDepthFlags() = 0;

    /**
     * Save the current set of stencil flags (stencil write, clear value etc) by pushing it onto the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushStencilFlags() = 0;

    /**
     * Save current set of parameters for rasterizing shapes (line width, smoothing flag for lines/points, polygon offsets
     * etc).
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL pushShapeFlags() = 0;

    /**
     * Set the color write state, specifying which color channels will be written and which will be suppressed.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setColorWrite(int8 red, int8 green, int8 blue, int8 alpha) = 0;

    /**
     * Set the depth write state (off or on).
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setDepthWrite(int8 value) = 0;

    /**
     * Set the stencil write mask.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setStencilWrite(uint32 mask) = 0;

    /**
     * Set the current drawing color, as a 4-element array in RGBA order.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setColor(Real const * rgba) = 0;

    /**
     * Set the current drawing color.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setColor(Real r, Real g, Real b, Real a = 1) = 0;

    /**
     * Set the color to clear the drawing buffer with, as a 4-element array in RGBA order.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setClearColor(Real const * rgba) = 0;

    /**
     * Set the color to clear the drawing buffer with.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setClearColor(Real r, Real g, Real b, Real a = 1) = 0;

    /**
     * Set the value to clear the depth buffer with.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setClearDepth(float64 value) = 0;

    /**
     * Set the value to clear the stencil buffer with.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setClearStencil(int64 value) = 0;

    /**
     * Clear color, depth and stencil buffers.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL clear() = 0;

    /**
     * Clear the selected buffers.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL clear(int8 color, int8 depth, int8 stencil) = 0;

    /**
     * Set the depth test, as a value from the DepthTest enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setDepthTest(int32 test) = 0;

    /**
     * Set the faces to be culled, as a value from the CullFace enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setCullFace(int32 cull) = 0;

    /**
     * Set the depth offset, if any, to be applied to polygon faces. The supplied value is scaled by an implementation-specific
     * offset.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setPolygonOffset(int8 enable, float64 offset = 1) = 0;

    /**
     * Set smoothing of rasterized polygons on/off.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setPolygonSmooth(int8 enable) = 0;

    /**
     * Set smoothing of rasterized lines on/off.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setLineSmooth(int8 enable) = 0;

    /**
     * Set smoothing of rasterized points on/off.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setPointSmooth(int8 enable) = 0;

    /**
     * Set the size (diameter) of rasterized points. This may be ignored if shaders are being used.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setPointSize(Real size = 1) = 0;

    /**
     * Restore the last saved set of color flags from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popColorFlags() = 0;

    /**
     * Restore the last saved set of depth flags from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popDepthFlags() = 0;

    /**
     * Restore the last saved set of stencil flags from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popStencilFlags() = 0;

    /**
     * Restore the last saved set of shape flags from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popShapeFlags() = 0;

    /**
     * Restore the last saved rendersystem state from the stack.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL popState() = 0;

    /**
     * The function returns only when the rendersystem has completely executed all previously-issued drawing calls.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL finishAllOperations() = 0;

}; // class IRenderSystem

// Pure virtual destructor should have implementation
inline IRenderSystem::~IRenderSystem() {}

/** Interface for a rendersystem factory. Should be implemented and registered by each actual rendersystem. */
class THEA_API IRenderSystemFactory
{
  public:
    /** Destructor. */
    virtual ~IRenderSystemFactory() = 0;

    /**
     * Create a rendersystem with the given name. The rendersystem must be destroyed using destroyRenderSystem(). If no
     * rendering context is available, the new rendersystem will be set up for offscreen rendering if possible. In this case,
     * you must explictly create and attach a framebuffer to the rendersystem before you can call any drawing functions.
     */
    virtual IRenderSystem * THEA_ICALL createRenderSystem(char const * name) = 0;

    /** Destroy a rendersystem created with createRenderSystem(). */
    virtual int8 THEA_ICALL destroyRenderSystem(IRenderSystem * render_system) = 0;

}; // class IRenderSystemFactory

inline IRenderSystemFactory::~IRenderSystemFactory() {}

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
    bool installFactory(std::string const & type, IRenderSystemFactory * factory);

    /** Uninstall a factory for a particular renderystem type. The match is case-insensitive. */
    void uninstallFactory(std::string const & type);

    /** Get a factory for rendersystem of a given type. An error is thrown if no such factory has been installed. */
    IRenderSystemFactory * getFactory(std::string const & type);

  private:
    typedef Map<std::string, IRenderSystemFactory *> FactoryMap;  ///< Maps rendersystem types to factory instances.

    FactoryMap installed_factories;  ///< Set of installed factories, one for each rendersystem type.

}; // class RenderSystemManager

} // namespace Graphics
} // namespace Thea

#endif
