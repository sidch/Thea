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

#ifndef __Thea_Graphics_GlShader_hpp__
#define __Thea_Graphics_GlShader_hpp__

#include "../../Graphics/IShader.hpp"
#include "../../Algorithms/FastCopy.hpp"
#include "../../Array.hpp"
#include "../../Map.hpp"
#include "GlCommon.hpp"
#include "GlTexture.hpp"
#include "GlHeaders.hpp"
#include <type_traits>

namespace Thea {
namespace Graphics {
namespace Gl {

// Forward declarations
class GlRenderSystem;

/** An OpenGL shader. */
class THEA_GL_DLL_LOCAL GlShader : public IShader
{
  public:
    /** Constructor. */
    GlShader(GlRenderSystem * render_system_, char const * name_);

    /** Destructor. */
    ~GlShader();

    /** Get the parent rendersystem. */
    GlRenderSystem * getRenderSystem() const { return render_system; }

    char const * THEA_ICALL getName() const { return name.c_str(); }
    int8 THEA_ICALL isComplete() const { return complete; }
    int8 THEA_ICALL attachModuleFromFile(int32 type, char const * path);
    int8 THEA_ICALL attachModuleFromString(int32 type, char const * source);
    int8 THEA_ICALL hasUniform(char const * uniform_name) const { return uniforms.find(uniform_name) != uniforms.end(); }

    int8 THEA_ICALL setUniform(char const * uniform_name, int32 value)
    { return setUniformHelper(uniform_name, value); }

    int8 THEA_ICALL setUniform(char const * uniform_name, Real value)
    { return setUniformHelper(uniform_name, value); }

    int8 THEA_ICALL setUniform(char const * uniform_name, IDenseMatrix<int32> const * value)
    { return setUniformHelper(uniform_name, value); }

    int8 THEA_ICALL setUniform(char const * uniform_name, IDenseMatrix<Real> const * value)
    { return setUniformHelper(uniform_name, value); }

    int8 THEA_ICALL setUniform(char const * uniform_name, ITexture * value);

    int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, int32 const * values)
    { return setUniformHelper(uniform_name, num_values, values); }

    int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, Real const * values)
    { return setUniformHelper(uniform_name, num_values, values); }

    int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, IDenseMatrix<int32> const * const * values)
    { return setUniformHelper(uniform_name, num_values, values); }

    int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, IDenseMatrix<Real> const * const * values)
    { return setUniformHelper(uniform_name, num_values, values); }

    int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, ITexture * const * values);

    /** Link the various modules of the shader into a single program. */
    int8 link();

    /** Use the shader for rendering. */
    int8 use();

    /** Get the OpenGL ID of the shader. */
    GLhandleARB getGlId() const { return program_id; }

  private:
    /** A value for a uniform variable. */
    struct UniformValue
    {
      GLint i_value;
      GLfloat f_value;
      Array<GLint> i_array;    // matrices are unrolled in column-major order
      Array<GLfloat> f_array;  // matrices are unrolled in column-major order
      GlTexture * texture;

      /** Get the cached value of the specified type. */
      template <typename T> T const & getValue() const;

      /** Get the cached value of the specified type. */
      template <typename T> T & getValue();

      /** Get the cached array of the specified type. */
      template <typename T> Array<T> const & getArray() const;

      /** Get the cached array of the specified type. */
      template <typename T> Array<T> & getArray();
    };

    /** Data related to an uniform variable. */
    struct UniformData
    {
      GLenum type;
      GLint size;
      GLint location;
      int32 texunit;
      bool has_value;
      UniformValue value;
      bool requires_rebind;

      /** Constructor. */
      UniformData() : has_value(false), requires_rebind(false) {}

      /** Note that the value has been changed. */
      void valueChanged() { has_value = true; requires_rebind = true; }
    };

    /** A set of uniforms read from source code. */
    typedef Map<std::string, UniformData> Uniforms;

    /** Read the list of active uniforms in the shader object. */
    int8 readActiveUniforms();

    /** Bind the user-provided uniforms to the shader object. */
    int8 bindUniforms();

    /**
     * Check if a build step (compile or link) succeeded, and signal an error by printing a message and returning false if it
     * did not.
     */
    int8 checkBuildStatus(GLhandleARB obj_id, GLenum status_field, std::string const & error_msg);

    /** Helper function for setting a single scalar uniform. */
    template <typename T> int8 setUniformHelper(char const * uniform_name, T value);

    /** Helper function for setting a single matrix or vector uniform. */
    template <typename T> int8 setUniformHelper(char const * uniform_name, IDenseMatrix<T> const * value);

    /** Helper function for setting a uniform which is an array of scalars. */
    template <typename T> int8 setUniformHelper(char const * uniform_name, int64 num_values, T const * values);

    /** Helper function for setting a uniform which is an array of matrices/vectors. */
    template <typename T>
    int8 THEA_ICALL setUniformHelper(char const * uniform_name, int64 num_values, IDenseMatrix<T> const * const * values);

    GlRenderSystem * render_system;
    std::string name;
    bool complete;
    bool linked;
    bool has_vertex_module;
    bool has_fragment_module;
    GLhandleARB program_id;
    Uniforms uniforms;

}; // class GlShader

template <> inline GLint const & GlShader::UniformValue::getValue<GLint>() const { return i_value; }
template <> inline GLfloat  const & GlShader::UniformValue::getValue<GLfloat >() const { return f_value; }

template <> inline GLint & GlShader::UniformValue::getValue<GLint>() { return i_value; }
template <> inline GLfloat  & GlShader::UniformValue::getValue<GLfloat >() { return f_value; }

template <> inline Array<GLint> const & GlShader::UniformValue::getArray<GLint>() const { return i_array; }
template <> inline Array<GLfloat>  const & GlShader::UniformValue::getArray<GLfloat >() const { return f_array; }

template <> inline Array<GLint> & GlShader::UniformValue::getArray<GLint>() { return i_array; }
template <> inline Array<GLfloat>  & GlShader::UniformValue::getArray<GLfloat >() { return f_array; }

namespace GlInternal {

template <typename U, typename V>
void
toArray(IDenseMatrix<U> const & m, V * arr)
{
  if (m.isRowMajor())
    Math::getElementsColumnMajor(Math::mapTo< MatrixX<U, MatrixLayout::ROW_MAJOR> const >(m), arr);
  else
    Math::getElementsColumnMajor(Math::mapTo< MatrixX<U, MatrixLayout::COLUMN_MAJOR> const >(m), arr);
}

} // namespace GlInternal

namespace GlShaderInternal {

inline GLenum
toCanonicalType(GLenum type)
{
  switch (type)
  {
    case GL_BOOL : return GL_INT;
    case GL_BOOL_VEC2_ARB : return GL_INT_VEC2_ARB;
    case GL_BOOL_VEC3_ARB : return GL_INT_VEC3_ARB;
    case GL_BOOL_VEC4_ARB : return GL_INT_VEC4_ARB;
    case GL_SAMPLER_1D_ARB: return GL_TEXTURE_1D;
    case GL_SAMPLER_2D_ARB: return GL_TEXTURE_2D;
    case GL_SAMPLER_3D_ARB: return GL_TEXTURE_3D;
    case GL_SAMPLER_CUBE_ARB: return GL_TEXTURE_CUBE_MAP_ARB;
    case GL_SAMPLER_2D_RECT_ARB: return GL_TEXTURE_RECTANGLE_ARB;
    case GL_SAMPLER_2D_SHADOW_ARB: return GL_TEXTURE_2D;
    case GL_SAMPLER_2D_RECT_SHADOW_ARB: return GL_TEXTURE_RECTANGLE_ARB;
    default: return type;
  }
}

#define THEA_GLSHADER_SET_UNIFORM_CHECKS(uniform_gl_type, uniform_size)                                                       \
  if (entry == uniforms.end())                                                                                                \
  {                                                                                                                           \
    THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' not found";                                                \
    return false;                                                                                                             \
  }                                                                                                                           \
                                                                                                                              \
  if (GlShaderInternal::toCanonicalType(entry->second.type) != uniform_gl_type)                                               \
  {                                                                                                                           \
    THEA_ERROR << getName() << ": Argument does not match the declared type of uniform '" << uniform_name << '\'';            \
    return false;                                                                                                             \
  }                                                                                                                           \
                                                                                                                              \
  if (entry->second.size != (GLint)uniform_size)                                                                              \
  { THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' expects size " << uniform_size; return false; }

template <typename T, typename Enable = void> struct GlType {};
template <typename T> struct GlType<T, typename std::enable_if<  std::is_integral<T>::value >::type> { typedef GLint   type; };
template <typename T> struct GlType<T, typename std::enable_if< !std::is_integral<T>::value >::type> { typedef GLfloat type; };

template <typename T>
GLenum
uniformTypeFromScalar()
{
  return std::is_integral<T>::value ? GL_INT : GL_FLOAT;
}

template <typename T>
GLenum
uniformTypeFromMatrix(IDenseMatrix<T> const & m)
{
  switch (m.rows())
  {
    case 1:
      switch (m.cols())
      {
        case 2: return std::is_integral<T>::value ? GL_INT_VEC2_ARB : GL_FLOAT_VEC2_ARB;
        case 3: return std::is_integral<T>::value ? GL_INT_VEC3_ARB : GL_FLOAT_VEC3_ARB;
        case 4: return std::is_integral<T>::value ? GL_INT_VEC4_ARB : GL_FLOAT_VEC4_ARB;
      }
      break;

    case 2:
      switch (m.cols())
      {
        case 1: return std::is_integral<T>::value ? GL_INT_VEC2_ARB : GL_FLOAT_VEC2_ARB;
        case 2: return std::is_integral<T>::value ? GL_INVALID_ENUM : GL_FLOAT_MAT2_ARB;
      }
      break;

    case 3:
      switch (m.cols())
      {
        case 1: return std::is_integral<T>::value ? GL_INT_VEC3_ARB : GL_FLOAT_VEC3_ARB;
        case 3: return std::is_integral<T>::value ? GL_INVALID_ENUM : GL_FLOAT_MAT3_ARB;
      }
      break;

    case 4:
      switch (m.cols())
      {
        case 1: return std::is_integral<T>::value ? GL_INT_VEC4_ARB : GL_FLOAT_VEC4_ARB;
        case 4: return std::is_integral<T>::value ? GL_INVALID_ENUM : GL_FLOAT_MAT4_ARB;
      }
      break;
  }

  return GL_INVALID_ENUM;
}

} // namespace GlShaderInternal

template <typename T>
int8
GlShader::setUniformHelper(char const * uniform_name, T value)
{
  typedef typename GlShaderInternal::GlType<T>::type GT;

  Uniforms::iterator entry = uniforms.find(uniform_name);
  THEA_GLSHADER_SET_UNIFORM_CHECKS(GlShaderInternal::uniformTypeFromScalar<T>(), 1)
  entry->second.value.getValue<GT>() = static_cast<GT>(value);
  entry->second.valueChanged();
  return true;
}

template <typename T>
int8
GlShader::setUniformHelper(char const * uniform_name, IDenseMatrix<T> const * value)
{
  typedef typename GlShaderInternal::GlType<T>::type GT;

  if (!value) { THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' can't be assigned a null value"; return false; }

  GLenum gl_type = GlShaderInternal::uniformTypeFromMatrix(*value);
  if (gl_type == GL_INVALID_ENUM)
  {
    THEA_ERROR << getName() << ": Matrix provided for uniform '" << uniform_name << "' doesn't correspond to a valid GLSL type";
    return false;
  }

  Uniforms::iterator entry = uniforms.find(uniform_name);
  THEA_GLSHADER_SET_UNIFORM_CHECKS(gl_type, 1)

  size_t num_elems = (size_t)(value->rows() * value->cols());
  entry->second.value.getArray<GT>().resize(num_elems);
  GlInternal::toArray(*value, entry->second.value.getArray<GT>().data());

  entry->second.valueChanged();
  return true;
}

template <typename T>
int8
GlShader::setUniformHelper(char const * uniform_name, int64 num_values, T const * values)
{
  typedef typename GlShaderInternal::GlType<T>::type GT;

  if (!values) { THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' can't be assigned a null array"; return false; }

  Uniforms::iterator entry = uniforms.find(uniform_name);
  THEA_GLSHADER_SET_UNIFORM_CHECKS(GlShaderInternal::uniformTypeFromScalar<T>(), num_values)

  entry->second.value.getArray<GT>().resize((size_t)num_values);
  Algorithms::fastCopy(values, values + num_values, entry->second.value.getArray<GT>().data());

  return true;
}

template <typename T>
int8
GlShader::setUniformHelper(char const * uniform_name, int64 num_values, IDenseMatrix<T> const * const * values)
{
  typedef typename GlShaderInternal::GlType<T>::type GT;

  if (!values) { THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' can't be assigned a null array"; return false; }
  if (num_values < 1) { THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' assigned empty array"; return false; }

  for (intx i = 0; i < num_values; ++i)
    if (!values[0])
    { THEA_ERROR << getName() << ": Uniform '" << uniform_name << "' assigned null array entry"; return false; }

  GLenum gl_type = GlShaderInternal::uniformTypeFromMatrix(*values[0]);
  if (gl_type == GL_INVALID_ENUM)
  {
    THEA_ERROR << getName() << ": Matrix provided for uniform '" << uniform_name << "' doesn't correspond to a valid GLSL type";
    return false;
  }

  for (intx i = 1; i < num_values; ++i)
    if (GlShaderInternal::uniformTypeFromMatrix(*values[1]) != gl_type)
    {
      THEA_ERROR << getName() << ": Matrices provided for uniform '" << uniform_name << "' must have identical dimensions";
      return false;
    }

  Uniforms::iterator entry = uniforms.find(uniform_name);
  THEA_GLSHADER_SET_UNIFORM_CHECKS(gl_type, num_values)

  // Pack the unrolled matrices densely in the cached array
  intx num_elems = values[0]->rows() * values[0]->cols();
  entry->second.value.getArray<GT>().resize((size_t)(num_values * num_elems));
  for (intx i = 0; i < num_values; ++i)
    GlInternal::toArray(*values[i], entry->second.value.getArray<GT>().data() + i * num_elems);

  entry->second.valueChanged();
  return true;
}

} // namespace Gl
} // namespace Graphics
} // namespace Thea

#endif
