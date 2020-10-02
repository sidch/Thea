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

#ifndef __Thea_Graphics_IShader_hpp__
#define __Thea_Graphics_IShader_hpp__

#include "../Common.hpp"
#include "../IDenseMatrix.hpp"
#include "../NamedObject.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
class ITexture;

/**
 * Abstract base class for a shader.
 *
 * @todo Make this an interface, safe for passing across shared library boundaries.
 */
class THEA_API IShader : public INamedObject
{
  public:
    /** Shader module types (enum class). */
    struct THEA_API ModuleType
    {
      /** Supported values. */
      enum Value
      {
        VERTEX = 0,  ///< Vertex shader.
        FRAGMENT,    ///< Fragment (pixel) shader.
        GEOMETRY,    ///< Geometry shader.
        NUM          ///< [Internal] Number of allowed module types.
      };

      THEA_ENUM_CLASS_BODY(ModuleType)
    };

    /** Destructor. */
    virtual ~IShader() = 0;

    /**
     * Check if the shader is ready to be used for rendering or not. Typically this requires both a vertex and a fragment
     * program to be attached.
     */
    virtual int8 THEA_ICALL isComplete() const = 0;

    /**
     * Attach a program module to the shader from a file containing its source code. \a type should be a value from the
     * ModuleType enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL attachModuleFromFile(int32 type, char const * path) = 0;

    /**
     * Attach a program module to the shader from a string containing its source code. \a type should be a value from the
     * ModuleType enum.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL attachModuleFromString(int32 type, char const * source) = 0;

    /** Check if the shader has an active uniform of the given name. */
    virtual int8 THEA_ICALL hasUniform(char const * uniform_name) const = 0;

    /**
     * Set a uniform that is a single integer.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, int32 value) = 0;

    /**
     * Set a uniform that is a single real number.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, Real value) = 0;

    /**
     * Set a uniform that is a matrix of real numbers (handles vectors, matrices and colors).
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, IDenseMatrix<Real> const * value) = 0;

    /**
     * Set a uniform that is a matrix of integers (handles vectors, matrices and colors).
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, IDenseMatrix<int32> const * value) = 0;

    /**
     * Set a uniform that is a single texture.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, ITexture * value) = 0;

    /**
     * Set a uniform that is an array of integers.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, int32 const * values) = 0;

    /**
     * Set a uniform that is an array of real numbers.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, Real const * values) = 0;

    /**
     * Set a uniform that is an array of integer matrices (handles arrays of vectors, matrices and colors).
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, IDenseMatrix<int32> const * const * values) = 0;

    /**
     * Set a uniform that is an array of real matrices (handles arrays of vectors, matrices and colors).
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, IDenseMatrix<Real> const * const * values) = 0;

    /**
     * Set a uniform that is an array of textures.
     *
     * @return True on success, false on error.
     */
    virtual int8 THEA_ICALL setUniform(char const * uniform_name, int64 num_values, ITexture * const * values) = 0;

}; // class IShader

inline IShader::~IShader() {}

} // namespace Graphics
} // namespace Thea

#endif
