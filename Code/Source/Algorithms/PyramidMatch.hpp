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

#ifndef __Thea_Algorithms_PyramidMatch_hpp__
#define __Thea_Algorithms_PyramidMatch_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Math.hpp"
#include "../Matrix.hpp"
#include "../Serializable.hpp"

namespace Thea {
namespace Algorithms {

namespace PyramidInternal {

void smoothIncrement1D(Real * dst, int dst_size, int index, Real mean, Real value);

void smoothIncrement2D(Real * dst, int dst_sx, int dst_sy, int x, int y, Real mean_x, Real mean_y, Real value);

void smoothIncrement3D(Real * dst, int dst_sx, int dst_sy, int dst_sz, int x, int y, int z, Vector3 const & mean,
                       Real value);

} // namespace PyramidInternal

/** A 1D pyramid. */
class THEA_API Pyramid1D : public Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(Pyramid1D, shared_ptr, weak_ptr)

    /** Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling. */
    template <typename T> Pyramid1D(T const * base_data, int n) { construct(base_data, n, NULL, false); }

    /** Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling. */
    template <typename T> Pyramid1D(TheaArray<T> const & base_data)
    { construct(&base_data[0], (int)base_data.size(), NULL, false); }

    /**
     * Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling with gaussian smoothing. The
     * mean of each bin should be specified relative to the limits of the bin, in the range [0, 1].
     */
    template <typename T>
    Pyramid1D(T const * base_data, int n, Real const * means, bool smooth_base_level = false)
    {
      TheaArray<Real> means_copy(means, means + n);
      construct(base_data, n, &means_copy[0], smooth_base_level);
    }

    /**
     * Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling with gaussian smoothing. The
     * mean of each bin should be specified relative to the limits of the bin, in the range [0, 1].
     */
    template <typename T>
    Pyramid1D(TheaArray<T> const & base_data, TheaArray<Real> const & means, bool smooth_base_level = false)
    {
      TheaArray<Real> means_copy = means;
      construct(&base_data[0], (int)base_data.size(), &means_copy[0], smooth_base_level);
    }

    /** Load the pyramid from a binary input stream. */
    Pyramid1D(BinaryInputStream & input) { deserialize(input); }

    /** Load the pyramid from a text input stream. */
    Pyramid1D(TextInputStream & input) { deserialize(input); }

    /** Get the number of levels in the pyramid. */
    int numLevels() const { return num_levels; }

    /** Get the number of elements in the base (input) level. */
    int baseSize() const { return levels.empty() ? 0 : (int)levels[0].size(); }

    void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const;
    void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO());
    void serialize(TextOutputStream & output, Codec const & codec = Codec_AUTO()) const;
    void deserialize(TextInputStream & input, Codec const & codec = Codec_AUTO());

  protected:
    /** Default constructor. */
    Pyramid1D() {}

    /** Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling. */
    template <typename T> void construct(T const * base_data, int n, Real * means, bool smooth_base_level)
    {
      alwaysAssertM(Math::isPowerOf2((unsigned)n), "Pyramid1D: Base size must be a power of 2");

      num_levels = Math::floorLog2((uint32)n) + 1;
      levels.resize(num_levels);

      if (smooth_base_level)
      {
        levels[0].resize((array_size_t)n, 0);
        Real * dst = &levels[0][0];

        for (int i = 0; i < n; ++i)
          PyramidInternal::smoothIncrement1D(dst, n, i, means[i], base_data[i]);
      }
      else
      {
        levels[0].reserve((array_size_t)n);
        levels[0].insert(levels[0].end(), base_data, base_data + n);
      }

      createPyramid(means);
    }

    /** Create the pyramid by recursively downsampling the base level. */
    void createPyramid(Real * means);

    /** Resize a 1D array to half its size, combining each successive pair of elements. */
    void downsample(TheaArray<Real> const & src, TheaArray<Real> & dst, Real * means) const;

    int num_levels;
    TheaArray< TheaArray<Real> > levels;

    friend class PyramidMatch;

}; // class Pyramid1D

/** A 2D pyramid. */
class THEA_API Pyramid2D : public Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(Pyramid2D, shared_ptr, weak_ptr)

    /**
     * Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. The input array is assumed to
     * be in row-major order, i.e. y is the major index and x is the minor index.
     */
    template <typename T> Pyramid2D(T const * base_data, int nx_, int ny_)
    { construct(base_data, nx_, ny_, NULL, false); }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. */
    template <typename T> Pyramid2D(Matrix<T, MatrixLayout::ROW_MAJOR> const & base_data)
    { construct(base_data.data(), base_data.numColumns(), base_data.numRows(), NULL, false); }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. */
    template <typename T, MatrixLayout::Value Layout> Pyramid2D(Matrix<T, Layout> const & base_data)
    {
      Matrix<T, MatrixLayout::ROW_MAJOR> base_data_copy(base_data);
      construct(base_data_copy.data(), base_data.numColumns(), base_data.numRows(), NULL, false);
    }

    /**
     * Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling with gaussian smoothing. The
     * input array is assumed to be in row-major order, i.e. y is the major index and x is the minor index. The mean of each bin
     * should be specified relative to the limits of the bin, in the range [0, 1] x [0, 1].
     */
    template <typename T>
    Pyramid2D(T const * base_data, int nx_, int ny_, Vector2 const * means, bool smooth_base_level = false)
    {
      TheaArray<Vector2> means_copy(means, means + nx_ * ny_);
      construct(base_data, nx_, ny_, &means_copy[0], smooth_base_level);
    }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling with gaussian smoothing. */
    template <typename T, MatrixLayout::Value MeansLayout>
    Pyramid2D(Matrix<T, MatrixLayout::ROW_MAJOR> const & base_data, Matrix<Vector2, MeansLayout> const & means,
              bool smooth_base_level = false)
    {
      Matrix<Vector2, MatrixLayout::ROW_MAJOR> means_copy(means);
      construct(base_data.data(), base_data.numColumns(), base_data.numRows(), means_copy.data(), smooth_base_level);
    }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling with gaussian smoothing. */
    template <typename T, MatrixLayout::Value BaseLayout, MatrixLayout::Value MeansLayout>
    Pyramid2D(Matrix<T, BaseLayout> const & base_data, Matrix<Vector2, MeansLayout> const & means,
              bool smooth_base_level = false)
    {
      Matrix<T, MatrixLayout::ROW_MAJOR> base_data_copy(base_data);
      Matrix<Vector2, MatrixLayout::ROW_MAJOR> means_copy(means);
      construct(base_data_copy.data(), base_data.numColumns(), base_data.numRows(), means_copy.data(), smooth_base_level);
    }

    /** Load the pyramid from a binary input stream. */
    Pyramid2D(BinaryInputStream & input) { deserialize(input); }

    /** Load the pyramid from a text input stream. */
    Pyramid2D(TextInputStream & input) { deserialize(input); }

    /** Get the number of levels in the pyramid. */
    int numLevels() const { return num_levels; }

    /** Get the number of columns in the base (input) level. */
    int baseSizeX() const { return nx; }

    /** Get the number of rows in the base (input) level. */
    int baseSizeY() const { return ny; }

    void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const;
    void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO());
    void serialize(TextOutputStream & output, Codec const & codec = Codec_AUTO()) const;
    void deserialize(TextInputStream & input, Codec const & codec = Codec_AUTO());

  protected:
    /** 2D array of scalars. */
    struct Array2D
    {
      TheaArray<Real> data;
      int sx, sy;

      Array2D() : sx(0), sy(0) {}
      Array2D(int sx_, int sy_) : data(sx_ * sy_), sx(sx_), sy(sy_) {}
      Array2D(int sx_, int sy_, Real value) : data(sx_ * sy_, value), sx(sx_), sy(sy_) {}

      template <typename T> Array2D(T const * data_, int sx_, int sy_) : data(data_, data_ + sx_ * sy_), sx(sx_), sy(sy_) {}

      Real operator()(int x, int y) const { return data[y * sx + x]; }
      Real & operator()(int x, int y) { return data[y * sx + x]; }
    };

    /** Default constructor. */
    Pyramid2D() {}

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. */
    template <typename T>
    void construct(T const * base_data, int nx_, int ny_, Vector2 * means, bool smooth_base_level)
    {
      alwaysAssertM(Math::isPowerOf2((unsigned)nx_) && Math::isPowerOf2((unsigned)ny_),
                    "Pyramid2D: Base dimensions must be powers of 2");

      nx = nx_;
      ny = ny_;
      num_levels = Math::floorLog2((uint32)std::min(nx, ny)) + 1;

      levels.resize(num_levels);

      if (smooth_base_level)
      {
        levels[0] = Array2D(nx, ny, 0);

        Real * dst = &levels[0].data[0];

        long i = 0;
        for (int y = 0; y < ny; ++y)
          for (int x = 0; x < nx; ++x, ++i)
            PyramidInternal::smoothIncrement2D(dst, nx, ny, x, y, means[i].x(), means[i].y(), base_data[i]);
      }
      else
        levels[0] = Array2D(base_data, nx, ny);

      createPyramid(means);
    }

    /** Create the pyramid by recursively downsampling the base level. */
    void createPyramid(Vector2 * means);

    /** Resize a 2D array to half its size in each specified dimension. */
    void downsample(Array2D const & src, Array2D & dst, int num_dims_to_downsample, Vector2 * means) const;

    int num_levels;
    int nx, ny;
    TheaArray<Array2D> levels;

    friend class PyramidMatch;

}; // class Pyramid2D

/** A 3D pyramid. */
class THEA_API Pyramid3D : public Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(Pyramid3D, shared_ptr, weak_ptr)

    /**
     * Construct a pyramid from the base (highest-resolution) 3D array by recursive downsampling. The input array is assumed to
     * be in z-major order, i.e. z is the major index, y is the submajor index and x is the minor index.
     */
    template <typename T> Pyramid3D(T const * base_data, int nx_, int ny_, int nz_)
    { construct(base_data, nx_, ny_, nz_, NULL, false); }

    /**
     * Construct a pyramid from the base (highest-resolution) 3D array by recursive downsampling with gaussian smoothing. The
     * input array is assumed to be in z-major order, i.e. z is the major index, y is the submajor index and x is the minor
     * index. The mean of each bin should be specified relative to the limits of the bin, in the range [0, 1] x [0, 1] x [0, 1].
     */
    template <typename T> Pyramid3D(T const * base_data, int nx_, int ny_, int nz_, Vector3 const * means,
                                    bool smooth_base_level = false)
    {
      TheaArray<Vector3> means_copy(means, means + nx_ * ny_ * nz_);
      construct(base_data, nx_, ny_, nz_, &means_copy[0], smooth_base_level);
    }

    /** Load the pyramid from a binary input stream. */
    Pyramid3D(BinaryInputStream & input) { deserialize(input); }

    /** Load the pyramid from a text input stream. */
    Pyramid3D(TextInputStream & input) { deserialize(input); }

    /** Get the number of levels in the pyramid. */
    int numLevels() const { return num_levels; }

    /** Get the number of bins in the X direction in the base (input) level. */
    int baseSizeX() const { return nx; }

    /** Get the number of bins in the Y direction in the base (input) level. */
    int baseSizeY() const { return ny; }

    /** Get the number of bins in the Z direction in the base (input) level. */
    int baseSizeZ() const { return nz; }

    void serialize(BinaryOutputStream & output, Codec const & codec = Codec_AUTO()) const;
    void deserialize(BinaryInputStream & input, Codec const & codec = Codec_AUTO());
    void serialize(TextOutputStream & output, Codec const & codec = Codec_AUTO()) const;
    void deserialize(TextInputStream & input, Codec const & codec = Codec_AUTO());

  protected:
    /** 3D array of scalars. */
    struct Array3D
    {
      TheaArray<Real> data;
      int sx, sy, sz;

      Array3D() : sx(0), sy(0), sz(0) {}
      Array3D(int sx_, int sy_, int sz_) : data(sx_ * sy_ * sz_), sx(sx_), sy(sy_), sz(sz_) {}
      Array3D(int sx_, int sy_, int sz_, Real value)
      : data(sx_ * sy_ * sz_, value), sx(sx_), sy(sy_), sz(sz_) {}

      template <typename T>
      Array3D(T const * data_, int sx_, int sy_, int sz_) : data(data_, data_ + sx_ * sy_ * sz_), sx(sx_), sy(sy_), sz(sz_) {}

      Real operator()(int x, int y, int z) const { return data[(z * sy + y) * sx + x]; }
      Real & operator()(int x, int y, int z) { return data[(z * sy + y) * sx + x]; }
    };

    /** Default constructor. */
    Pyramid3D() {}

    /** Construct a pyramid from the base (highest-resolution) 3D array by recursive downsampling. */
    template <typename T>
    void construct(T const * base_data, int nx_, int ny_, int nz_, Vector3 * means, bool smooth_base_level)
    {
      alwaysAssertM(Math::isPowerOf2((unsigned)nx_) && Math::isPowerOf2((unsigned)ny_) && Math::isPowerOf2((unsigned)nz_),
                    "Pyramid3D: Base dimensions must be powers of 2");

      nx = nx_;
      ny = ny_;
      nz = nz_;
      num_levels = Math::floorLog2((uint32)std::min(nz, std::min(nx, ny))) + 1;

      levels.resize(num_levels);

      if (smooth_base_level)
      {
        levels[0] = Array3D(nx, ny, nz, 0);

        Real * dst = &levels[0].data[0];

        int i = 0;
        for (int z = 0; z < nz; ++z)
          for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x, ++i)
              PyramidInternal::smoothIncrement3D(dst, nx, ny, nz, x, y, z, means[i], base_data[i]);
      }
      else
        levels[0] = Array3D(base_data, nx, ny, nz);

      createPyramid(means);
    }

    /** Create the pyramid by recursively downsampling the base level. */
    void createPyramid(Vector3 * means);

    /** Resize a 3D array to half its size in each specified dimension. */
    void downsample(Array3D const & src, Array3D & dst, int num_dims_to_downsample, Vector3 * means) const;

    int num_levels;
    int nx, ny, nz;
    TheaArray<Array3D> levels;

    friend class PyramidMatch;

}; // class Pyramid3D

/**
 * Compare two 1, 2 or 3-dimensional scalar arrays of the same size using a pyramid-matching technique. See Kristen Grauman and
 * Trevor Darrell, "The Pyramid Match Kernel: Discriminative Classification with Sets of Image Features," in Proc. IEEE Int'l
 * Conference on Computer Vision (ICCV), Beijing, China, Oct 2005.
 */
class THEA_API PyramidMatch
{
  public:
    /** The kernel to use for comparing pyramid levels (enum class). */
    struct Kernel
    {
      /** Supported values. */
      enum Value
      {
        HIK,          ///< Histogram intersection kernel.
        L1,           ///< Euclidean L1 norm.
        L2,           ///< Euclidean L2 norm.
        CHI_SQUARED   ///< Chi-squared norm.
      };

      THEA_ENUM_CLASS_BODY(Kernel)
    };

    /**
     * Compute the similarity of two pyramids with 1D supports. More similar pyramids will return a greater value (need not be
     * positive).
     */
    static Real similarity(Pyramid1D const & pyramid1, Pyramid1D const & pyramid2, Kernel kernel = Kernel::HIK, int levels = -1,
                           Real attenuation_factor = -1);

    /**
     * Compute the similarity of two pyramids with 2D supports. More similar pyramids will return a greater value (need not be
     * positive).
     */
    static Real similarity(Pyramid2D const & pyramid1, Pyramid2D const & pyramid2, Kernel kernel = Kernel::HIK, int levels = -1,
                           Real attenuation_factor = -1);

    /**
     * Compute the similarity of two pyramids with 3D supports. More similar pyramids will return a greater value (need not be
     * positive).
     */
    static Real similarity(Pyramid3D const & pyramid1, Pyramid3D const & pyramid2, Kernel kernel = Kernel::HIK, int levels = -1,
                           Real attenuation_factor = -1);

    /** Run some unit tests. Throws an error if a test fails. */
    static void test();

}; // class PyramidMatch

} // namespace Algorithms
} // namespace Thea

#endif
