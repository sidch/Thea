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

#ifndef __Thea_Algorithms_PyramidMatch_hpp__
#define __Thea_Algorithms_PyramidMatch_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "../Serializable.hpp"

namespace Thea {
namespace Algorithms {

namespace PyramidInternal {

void smoothIncrement1d(Real * dst, int dst_size, int index, Real mean, Real value);

void smoothIncrement2d(Real * dst, int dst_sx, int dst_sy, int x, int y, Real mean_x, Real mean_y, Real value);

void smoothIncrement3d(Real * dst, int dst_sx, int dst_sy, int dst_sz, int x, int y, int z, Vector3 const & mean,
                       Real value);

} // namespace PyramidInternal

/** A 1D pyramid. */
class THEA_API Pyramid1d : public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(Pyramid1d)

    /** Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling. */
    template <typename T> Pyramid1d(T const * base_data, int n) { construct(base_data, n, nullptr, false); }

    /** Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling. */
    template <typename T> Pyramid1d(Array<T> const & base_data)
    { construct(&base_data[0], (int)base_data.size(), nullptr, false); }

    /**
     * Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling with gaussian smoothing. The
     * mean of each bin should be specified relative to the limits of the bin, in the range [0, 1].
     */
    template <typename T>
    Pyramid1d(T const * base_data, int n, Real const * means, bool smooth_base_level = false)
    {
      Array<Real> means_copy(means, means + n);
      construct(base_data, n, &means_copy[0], smooth_base_level);
    }

    /**
     * Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling with gaussian smoothing. The
     * mean of each bin should be specified relative to the limits of the bin, in the range [0, 1].
     */
    template <typename T>
    Pyramid1d(Array<T> const & base_data, Array<Real> const & means, bool smooth_base_level = false)
    {
      Array<Real> means_copy = means;
      construct(&base_data[0], (int)base_data.size(), &means_copy[0], smooth_base_level);
    }

    /** Load the pyramid from a binary input stream. */
    Pyramid1d(BinaryInputStream & input) { read(input); }

    /** Load the pyramid from a text input stream. */
    Pyramid1d(TextInputStream & input) { read(input); }

    /** Get the number of levels in the pyramid. */
    int numLevels() const { return num_levels; }

    /** Get the number of elements in the base (input) level. */
    int baseSize() const { return levels.empty() ? 0 : (int)levels[0].size(); }

    void read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false);
    void write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false) const;
    void read(TextInputStream & input, Codec const & codec = CodecAuto());
    void write(TextOutputStream & output, Codec const & codec = CodecAuto()) const;

  protected:
    /** Default constructor. */
    Pyramid1d() {}

    /** Construct a pyramid from the base (highest-resolution) 1D array by recursive downsampling. */
    template <typename T> void construct(T const * base_data, int n, Real * means, bool smooth_base_level)
    {
      alwaysAssertM(Math::isPowerOf2((unsigned)n), "Pyramid1d: Base size must be a power of 2");

      num_levels = Math::floorLog2((uint32)n) + 1;
      levels.resize(num_levels);

      if (smooth_base_level)
      {
        levels[0].resize((size_t)n, 0);
        Real * dst = &levels[0][0];

        for (int i = 0; i < n; ++i)
          PyramidInternal::smoothIncrement1d(dst, n, i, means[i], base_data[i]);
      }
      else
      {
        levels[0].reserve((size_t)n);
        levels[0].insert(levels[0].end(), base_data, base_data + n);
      }

      createPyramid(means);
    }

    /** Create the pyramid by recursively downsampling the base level. */
    void createPyramid(Real * means);

    /** Resize a 1D array to half its size, combining each successive pair of elements. */
    void downsample(Array<Real> const & src, Array<Real> & dst, Real * means) const;

    int num_levels;
    Array< Array<Real> > levels;

    friend class PyramidMatch;

}; // class Pyramid1d

/** A 2D pyramid. */
class THEA_API Pyramid2d : public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(Pyramid2d)

    /**
     * Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. The input array is assumed to
     * be in row-major order, i.e. y is the major index and x is the minor index.
     */
    template <typename T> Pyramid2d(T const * base_data, int nx_, int ny_)
    { construct(base_data, nx_, ny_, nullptr, false); }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. */
    template <typename T> Pyramid2d(MatrixX<T, MatrixLayout::ROW_MAJOR> const & base_data)
    { construct(base_data.data(), base_data.cols(), base_data.rows(), nullptr, false); }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. */
    template <typename T, MatrixLayout::Value Layout> Pyramid2d(MatrixX<T, Layout> const & base_data)
    {
      MatrixX<T, MatrixLayout::ROW_MAJOR> base_data_copy(base_data);
      construct(base_data_copy.data(), base_data.cols(), base_data.rows(), nullptr, false);
    }

    /**
     * Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling with gaussian smoothing. The
     * input array is assumed to be in row-major order, i.e. y is the major index and x is the minor index. The mean of each bin
     * should be specified relative to the limits of the bin, in the range [0, 1] x [0, 1].
     */
    template <typename T>
    Pyramid2d(T const * base_data, int nx_, int ny_, Vector2 const * means, bool smooth_base_level = false)
    {
      Array<Vector2> means_copy(means, means + nx_ * ny_);
      construct(base_data, nx_, ny_, &means_copy[0], smooth_base_level);
    }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling with gaussian smoothing. */
    template <typename T, MatrixLayout::Value MeansLayout>
    Pyramid2d(MatrixX<T, MatrixLayout::ROW_MAJOR> const & base_data, MatrixX<Vector2, MeansLayout> const & means,
              bool smooth_base_level = false)
    {
      MatrixX<Vector2, MatrixLayout::ROW_MAJOR> means_copy(means);
      construct(base_data.data(), base_data.cols(), base_data.rows(), means_copy.data(), smooth_base_level);
    }

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling with gaussian smoothing. */
    template <typename T, MatrixLayout::Value BaseLayout, MatrixLayout::Value MeansLayout>
    Pyramid2d(MatrixX<T, BaseLayout> const & base_data, MatrixX<Vector2, MeansLayout> const & means,
              bool smooth_base_level = false)
    {
      MatrixX<T, MatrixLayout::ROW_MAJOR> base_data_copy(base_data);
      MatrixX<Vector2, MatrixLayout::ROW_MAJOR> means_copy(means);
      construct(base_data_copy.data(), base_data.cols(), base_data.rows(), means_copy.data(), smooth_base_level);
    }

    /** Load the pyramid from a binary input stream. */
    Pyramid2d(BinaryInputStream & input) { read(input); }

    /** Load the pyramid from a text input stream. */
    Pyramid2d(TextInputStream & input) { read(input); }

    /** Get the number of levels in the pyramid. */
    int numLevels() const { return num_levels; }

    /** Get the number of columns in the base (input) level. */
    int baseSizeX() const { return nx; }

    /** Get the number of rows in the base (input) level. */
    int baseSizeY() const { return ny; }

    void read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false);
    void write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false) const;
    void read(TextInputStream & input, Codec const & codec = CodecAuto());
    void write(TextOutputStream & output, Codec const & codec = CodecAuto()) const;

  protected:
    /** 2D array of scalars. */
    struct Array2d
    {
      Array<Real> data;
      int sx, sy;

      Array2d() : sx(0), sy(0) {}
      Array2d(int sx_, int sy_) : data(sx_ * sy_), sx(sx_), sy(sy_) {}
      Array2d(int sx_, int sy_, Real value) : data(sx_ * sy_, value), sx(sx_), sy(sy_) {}

      template <typename T> Array2d(T const * data_, int sx_, int sy_) : data(data_, data_ + sx_ * sy_), sx(sx_), sy(sy_) {}

      Real operator()(int x, int y) const { return data[y * sx + x]; }
      Real & operator()(int x, int y) { return data[y * sx + x]; }
    };

    /** Default constructor. */
    Pyramid2d() {}

    /** Construct a pyramid from the base (highest-resolution) 2D array by recursive downsampling. */
    template <typename T>
    void construct(T const * base_data, int nx_, int ny_, Vector2 * means, bool smooth_base_level)
    {
      alwaysAssertM(Math::isPowerOf2((unsigned)nx_) && Math::isPowerOf2((unsigned)ny_),
                    "Pyramid2d: Base dimensions must be powers of 2");

      nx = nx_;
      ny = ny_;
      num_levels = Math::floorLog2((uint32)std::min(nx, ny)) + 1;

      levels.resize(num_levels);

      if (smooth_base_level)
      {
        levels[0] = Array2d(nx, ny, 0);

        Real * dst = &levels[0].data[0];

        intx i = 0;
        for (int y = 0; y < ny; ++y)
          for (int x = 0; x < nx; ++x, ++i)
            PyramidInternal::smoothIncrement2d(dst, nx, ny, x, y, means[i].x(), means[i].y(), base_data[i]);
      }
      else
        levels[0] = Array2d(base_data, nx, ny);

      createPyramid(means);
    }

    /** Create the pyramid by recursively downsampling the base level. */
    void createPyramid(Vector2 * means);

    /** Resize a 2D array to half its size in each specified dimension. */
    void downsample(Array2d const & src, Array2d & dst, int num_dims_to_downsample, Vector2 * means) const;

    int num_levels;
    int nx, ny;
    Array<Array2d> levels;

    friend class PyramidMatch;

}; // class Pyramid2d

/** A 3D pyramid. */
class THEA_API Pyramid3d : public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(Pyramid3d)

    /**
     * Construct a pyramid from the base (highest-resolution) 3D array by recursive downsampling. The input array is assumed to
     * be in z-major order, i.e. z is the major index, y is the submajor index and x is the minor index.
     */
    template <typename T> Pyramid3d(T const * base_data, int nx_, int ny_, int nz_)
    { construct(base_data, nx_, ny_, nz_, nullptr, false); }

    /**
     * Construct a pyramid from the base (highest-resolution) 3D array by recursive downsampling with gaussian smoothing. The
     * input array is assumed to be in z-major order, i.e. z is the major index, y is the submajor index and x is the minor
     * index. The mean of each bin should be specified relative to the limits of the bin, in the range [0, 1] x [0, 1] x [0, 1].
     */
    template <typename T> Pyramid3d(T const * base_data, int nx_, int ny_, int nz_, Vector3 const * means,
                                    bool smooth_base_level = false)
    {
      Array<Vector3> means_copy(means, means + nx_ * ny_ * nz_);
      construct(base_data, nx_, ny_, nz_, &means_copy[0], smooth_base_level);
    }

    /** Load the pyramid from a binary input stream. */
    Pyramid3d(BinaryInputStream & input) { read(input); }

    /** Load the pyramid from a text input stream. */
    Pyramid3d(TextInputStream & input) { read(input); }

    /** Get the number of levels in the pyramid. */
    int numLevels() const { return num_levels; }

    /** Get the number of bins in the X direction in the base (input) level. */
    int baseSizeX() const { return nx; }

    /** Get the number of bins in the Y direction in the base (input) level. */
    int baseSizeY() const { return ny; }

    /** Get the number of bins in the Z direction in the base (input) level. */
    int baseSizeZ() const { return nz; }

    void read(BinaryInputStream & input, Codec const & codec = CodecAuto(), bool read_block_header = false);
    void write(BinaryOutputStream & output, Codec const & codec = CodecAuto(), bool write_block_header = false) const;
    void read(TextInputStream & input, Codec const & codec = CodecAuto());
    void write(TextOutputStream & output, Codec const & codec = CodecAuto()) const;

  protected:
    /** 3D array of scalars. */
    struct Array3d
    {
      Array<Real> data;
      int sx, sy, sz;

      Array3d() : sx(0), sy(0), sz(0) {}
      Array3d(int sx_, int sy_, int sz_) : data(sx_ * sy_ * sz_), sx(sx_), sy(sy_), sz(sz_) {}
      Array3d(int sx_, int sy_, int sz_, Real value)
      : data(sx_ * sy_ * sz_, value), sx(sx_), sy(sy_), sz(sz_) {}

      template <typename T>
      Array3d(T const * data_, int sx_, int sy_, int sz_) : data(data_, data_ + sx_ * sy_ * sz_), sx(sx_), sy(sy_), sz(sz_) {}

      Real operator()(int x, int y, int z) const { return data[(z * sy + y) * sx + x]; }
      Real & operator()(int x, int y, int z) { return data[(z * sy + y) * sx + x]; }
    };

    /** Default constructor. */
    Pyramid3d() {}

    /** Construct a pyramid from the base (highest-resolution) 3D array by recursive downsampling. */
    template <typename T>
    void construct(T const * base_data, int nx_, int ny_, int nz_, Vector3 * means, bool smooth_base_level)
    {
      alwaysAssertM(Math::isPowerOf2((unsigned)nx_) && Math::isPowerOf2((unsigned)ny_) && Math::isPowerOf2((unsigned)nz_),
                    "Pyramid3d: Base dimensions must be powers of 2");

      nx = nx_;
      ny = ny_;
      nz = nz_;
      num_levels = Math::floorLog2((uint32)std::min(nz, std::min(nx, ny))) + 1;

      levels.resize(num_levels);

      if (smooth_base_level)
      {
        levels[0] = Array3d(nx, ny, nz, 0);

        Real * dst = &levels[0].data[0];

        int i = 0;
        for (int z = 0; z < nz; ++z)
          for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x, ++i)
              PyramidInternal::smoothIncrement3d(dst, nx, ny, nz, x, y, z, means[i], base_data[i]);
      }
      else
        levels[0] = Array3d(base_data, nx, ny, nz);

      createPyramid(means);
    }

    /** Create the pyramid by recursively downsampling the base level. */
    void createPyramid(Vector3 * means);

    /** Resize a 3D array to half its size in each specified dimension. */
    void downsample(Array3d const & src, Array3d & dst, int num_dims_to_downsample, Vector3 * means) const;

    int num_levels;
    int nx, ny, nz;
    Array<Array3d> levels;

    friend class PyramidMatch;

}; // class Pyramid3d

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
    static Real similarity(Pyramid1d const & pyramid1, Pyramid1d const & pyramid2, Kernel kernel = Kernel::HIK, int levels = -1,
                           Real attenuation_factor = -1);

    /**
     * Compute the similarity of two pyramids with 2D supports. More similar pyramids will return a greater value (need not be
     * positive).
     */
    static Real similarity(Pyramid2d const & pyramid1, Pyramid2d const & pyramid2, Kernel kernel = Kernel::HIK, int levels = -1,
                           Real attenuation_factor = -1);

    /**
     * Compute the similarity of two pyramids with 3D supports. More similar pyramids will return a greater value (need not be
     * positive).
     */
    static Real similarity(Pyramid3d const & pyramid1, Pyramid3d const & pyramid2, Kernel kernel = Kernel::HIK, int levels = -1,
                           Real attenuation_factor = -1);

    /** Run some unit tests. Throws an error if a test fails. */
    static void test();

}; // class PyramidMatch

} // namespace Algorithms
} // namespace Thea

#endif
