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

#ifndef __Thea_Graphics_Gl_GlRenderSystem_hpp__
#define __Thea_Graphics_Gl_GlRenderSystem_hpp__

#include "../../Graphics/IRenderSystem.hpp"
#include "../../Stack.hpp"
#include "../../UnorderedSet.hpp"
#include "GlCommon.hpp"
#include "GlFramebuffer.hpp"
#include "GlShader.hpp"
#include "GlBuffer.hpp"

namespace Thea {
namespace Graphics {
namespace Gl {

/** An OpenGL rendersystem. <b>Not threadsafe</b>. */
class THEA_GL_DLL_LOCAL GlRenderSystem : public virtual IRenderSystem
{
  public:
    /** Constructor. Creates a headless context for rendering if no rendering context currently exists. */
    GlRenderSystem(char const * name_);

    /** Destructor. */
    ~GlRenderSystem();

    char const * THEA_ICALL getName() const { return name.c_str(); }
    int8 THEA_ICALL setName(char const * s) { return false;  /* name is read-only */ }

    char const * THEA_ICALL describeSystem() const;
    char const * THEA_ICALL getLastError() const;
    void THEA_ICALL clearErrors();

    IFramebuffer * THEA_ICALL createFramebuffer(char const * name);
    int8 THEA_ICALL destroyFramebuffer(IFramebuffer * framebuffer);

    IShader * THEA_ICALL createShader(char const * name);
    int8 THEA_ICALL destroyShader(IShader * shader);

    ITexture * THEA_ICALL createTexture(char const * name, int64 width, int64 height, int64 depth,
                             ITexture::Format const * desired_format,
                             int32 dimension = ITexture::Dimension::DIM_2D,
                             ITexture::Options const * options = nullptr);
    ITexture * THEA_ICALL createTexture(char const * name, IImage const * image,
                             ITexture::Format const * desired_format = nullptr,
                             int32 dimension = ITexture::Dimension::DIM_2D,
                             ITexture::Options const * options = nullptr);
    ITexture * THEA_ICALL createTexture(char const * name, IImage const * images[6],
                             ITexture::Format const * desired_format = nullptr,
                             ITexture::Options const * options = nullptr);
    int8 THEA_ICALL destroyTexture(ITexture * texture);

    IBufferPool * THEA_ICALL createBufferPool(char const * name, int64 num_bytes, int32 usage, int8 gpu_memory = true);
    int8 THEA_ICALL destroyBufferPool(IBufferPool * pool);

    int8 THEA_ICALL pushFramebuffer();
    int8 THEA_ICALL setFramebuffer(IFramebuffer * framebuffer);
    IFramebuffer const * THEA_ICALL getFramebuffer() const;
    IFramebuffer * THEA_ICALL getFramebuffer();
    int8 THEA_ICALL popFramebuffer();

    int8 THEA_ICALL pushShader();
    int8 THEA_ICALL setShader(IShader * shader);
    IShader const * THEA_ICALL getShader() const;
    IShader * THEA_ICALL getShader();
    int8 THEA_ICALL popShader();

    int8 THEA_ICALL pushTextures();
    int8 THEA_ICALL setTexture(int32 texunit, ITexture * texture);
    int8 THEA_ICALL popTextures();

    int32 THEA_ICALL getMatrixMode() const;
    int8 THEA_ICALL setMatrixMode(int32 mode);
    int8 THEA_ICALL pushMatrix();
    int8 THEA_ICALL pushViewMatrices();
    int8 THEA_ICALL getMatrix(int32 mode, IDenseMatrix<Real> * m) const;
    int8 THEA_ICALL setMatrix(IDenseMatrix<Real> const * m);
    int8 THEA_ICALL setIdentityMatrix();
    int8 THEA_ICALL multMatrix(IDenseMatrix<Real> const * m);
    int8 THEA_ICALL popMatrix();
    int8 THEA_ICALL popViewMatrices();

    int8 THEA_ICALL beginIndexedPrimitives();
    int8 THEA_ICALL setVertexBuffer(IBuffer const * vertices);
    int8 THEA_ICALL setColorBuffer(IBuffer const * colors);
    int8 THEA_ICALL setTexCoordBuffer(int32 texunit, IBuffer const * texcoords);
    int8 THEA_ICALL setNormalBuffer(IBuffer const * normals);
    int8 THEA_ICALL setIndexBuffer(IBuffer const * indices);
    int8 THEA_ICALL sendIndices(int32 primitive, int64 num_indices, uint8 const * indices);
    int8 THEA_ICALL sendIndices(int32 primitive, int64 num_indices, uint16 const * indices);
    int8 THEA_ICALL sendIndices(int32 primitive, int64 num_indices, uint32 const * indices);
    int8 THEA_ICALL sendSequentialIndices(int32 primitive, uint64 first_index, int64 num_indices);
    int8 THEA_ICALL sendIndicesFromBuffer(int32 primitive, uint64 offset, int64 num_indices);
    int8 THEA_ICALL endIndexedPrimitives();

    int8 THEA_ICALL beginPrimitive(int32 primitive);
    int8 THEA_ICALL sendVertex(int32 dims, Real const * coords);
    int8 THEA_ICALL sendVertex(Real x, Real y);
    int8 THEA_ICALL sendVertex(Real x, Real y, Real z);
    int8 THEA_ICALL sendVertex(Real x, Real y, Real z, Real w);
    int8 THEA_ICALL sendNormal(Real const * coords);
    int8 THEA_ICALL sendNormal(Real x, Real y, Real z);
    int8 THEA_ICALL sendTexCoord(int32 texunit, int32 dims, Real const * coords);
    int8 THEA_ICALL sendTexCoord(int32 texunit, Real x);
    int8 THEA_ICALL sendTexCoord(int32 texunit, Real x, Real y);
    int8 THEA_ICALL sendTexCoord(int32 texunit, Real x, Real y, Real z);
    int8 THEA_ICALL sendTexCoord(int32 texunit, Real x, Real y, Real z, Real w);
    int8 THEA_ICALL endPrimitive();

    int8 THEA_ICALL clear();
    int8 THEA_ICALL clear(int8 color, int8 depth, int8 stencil);

    int8 THEA_ICALL pushState();
    int8 THEA_ICALL pushColorFlags();
    int8 THEA_ICALL pushDepthFlags();
    int8 THEA_ICALL pushStencilFlags();
    int8 THEA_ICALL pushShapeFlags();
    int8 THEA_ICALL setColorWrite(int8 red, int8 green, int8 blue, int8 alpha);
    int8 THEA_ICALL setDepthWrite(int8 value);
    int8 THEA_ICALL setStencilWrite(uint32 mask);
    int8 THEA_ICALL setColor(Real const * rgba);
    int8 THEA_ICALL setColor(Real r, Real g, Real b, Real a = 1);
    int8 THEA_ICALL setClearColor(Real const * rgba);
    int8 THEA_ICALL setClearColor(Real r, Real g, Real b, Real a = 1);
    int8 THEA_ICALL setClearDepth(float64 value);
    int8 THEA_ICALL setClearStencil(int64 value);
    int8 THEA_ICALL setDepthTest(int32 test);
    int8 THEA_ICALL setCullFace(int32 cull);
    int8 THEA_ICALL setPolygonOffset(int8 enable, float64 offset = 1);
    int8 THEA_ICALL setPolygonSmooth(int8 enable);
    int8 THEA_ICALL setLineSmooth(int8 enable);
    int8 THEA_ICALL setPointSmooth(int8 enable);
    int8 THEA_ICALL setPointSize(Real size = 1);
    int8 THEA_ICALL popColorFlags();
    int8 THEA_ICALL popDepthFlags();
    int8 THEA_ICALL popStencilFlags();
    int8 THEA_ICALL popShapeFlags();
    int8 THEA_ICALL popState();
    int8 THEA_ICALL finishAllOperations();

  private:
    /** Stores the state of the current set of GPU buffers. */
    struct BufferState
    {
      GlBufferPool * attrib_pool;  ///< Vertex buffer.
      GlBufferPool * index_pool;   ///< Index buffer.
      GlBuffer       index_buf;    ///< Index array.

      /** Constructor. */
      BufferState(GlBufferPool * attrib_pool_ = nullptr, GlBufferPool * index_pool_ = nullptr, GlBuffer index_buf_ = GlBuffer())
      : attrib_pool(attrib_pool_), index_pool(index_pool_), index_buf(index_buf_)
      {}

    }; // struct BufferState

    typedef Stack<GlFramebuffer *>  FramebufferStack;
    typedef Stack<BufferState>      BufferStack;
    typedef Stack<GlShader *>       ShaderStack;

    typedef UnorderedSet<GlFramebuffer *>  FramebufferSet;
    typedef UnorderedSet<GlTexture *>      TextureSet;
    typedef UnorderedSet<GlShader *>       ShaderSet;
    typedef UnorderedSet<GlBufferPool *>   BufferPoolSet;

    /**
     * Set the current vertex pool to match the specified buffer. The vertex pool used within a single
     * beginIndexedPrimitives()/endIndexedPrimitives() block must remain constant.
     */
    int8 setAttributePoolFromBuffer(GlBuffer const & buf);

    /**
     * Set the current index pool to match the specified buffer. The vertex pool used within a single
     * beginIndexedPrimitives()/endIndexedPrimitives() block must remain constant.
     */
    int8 setIndexPoolFromBuffer(GlBuffer const & buf);

    std::string name;
    mutable std::string desc;

    GlFramebuffer * current_framebuffer;
    FramebufferStack framebuffer_stack;

    BufferState current_buffer_state;
    BufferStack buffer_stack;

    GlShader * current_shader;
    ShaderStack shader_stack;

    FramebufferSet  created_framebuffers;
    TextureSet      created_textures;
    ShaderSet       created_shaders;
    BufferPoolSet   created_bufpools;

}; // class GlRenderSystem

/** Factory for creating OpenGL rendersystems. */
class THEA_GL_DLL_LOCAL GlRenderSystemFactory : public virtual IRenderSystemFactory
{
  public:
    /** Destructor. */
    ~GlRenderSystemFactory();

    IRenderSystem * THEA_ICALL createRenderSystem(char const * name);
    int8 THEA_ICALL destroyRenderSystem(IRenderSystem * render_system);

    /** Destroy all rendersystems created with this factory. */
    int8 destroyAllRenderSystems();

  private:
    static IRenderSystem * singleton;
    static bool singleton_created;
};

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
