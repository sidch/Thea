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

#include "PyramidMatch.hpp"

namespace Thea {
namespace Algorithms {

namespace PyramidInternal {

// Integrals of the one-dimensional gaussian with standard deviation 0.741 (chosen to give an integral of 0.5 within distance
// 0.5 of the mean). GAUSSIAN_INTEGRALS[i][0], [i][1] and [i][2] give the integrals over the ranges (-infty, 0], [0, 1] and
// [1, infty) respectively of the gaussian with mean (i + 0.5) / 10.
int const NUM_GAUSSIAN_INTEGRAL_BINS = 10;
Real const GAUSSIAN_INTEGRALS[NUM_GAUSSIAN_INTEGRAL_BINS][3] = { {0.473112f, 0.426885f, 0.100003f},
                                                                 {0.419823f, 0.454411f, 0.125766f},
                                                                 {0.367966f, 0.476202f, 0.155832f},
                                                                 {0.318412f, 0.491301f, 0.190287f},
                                                                 {0.271912f, 0.499026f, 0.229062f},
                                                                 {0.229062f, 0.499026f, 0.271912f},
                                                                 {0.190287f, 0.491301f, 0.318412f},
                                                                 {0.155832f, 0.476202f, 0.367966f},
                                                                 {0.125766f, 0.454411f, 0.419823f},
                                                                 {0.100003f, 0.426885f, 0.473112f} };

void
smoothIncrement1D(Real * dst, int dst_size, int index, Real mean, Real value)
{
  int sub_bin = Math::clamp((int)std::floor(NUM_GAUSSIAN_INTEGRAL_BINS * mean), 0, NUM_GAUSSIAN_INTEGRAL_BINS - 1);
  dst[std::max(index - 1, 0)           ]  +=  value * GAUSSIAN_INTEGRALS[sub_bin][0];  // (-infty, 0]
  dst[index                            ]  +=  value * GAUSSIAN_INTEGRALS[sub_bin][1];  // [0, 1]
  dst[std::min(index + 1, dst_size - 1)]  +=  value * GAUSSIAN_INTEGRALS[sub_bin][2];  // [1, infty)
}

void
smoothIncrement2D(Real * dst, int dst_sx, int dst_sy, int x, int y, Real mean_x, Real mean_y, Real value)
{
  int sub_bin_x = Math::clamp((int)std::floor(NUM_GAUSSIAN_INTEGRAL_BINS * mean_x), 0, NUM_GAUSSIAN_INTEGRAL_BINS - 1);
  int sub_bin_y = Math::clamp((int)std::floor(NUM_GAUSSIAN_INTEGRAL_BINS * mean_y), 0, NUM_GAUSSIAN_INTEGRAL_BINS - 1);

  int ix[3] = { std::max(x - 1, 0), x, std::min(x + 1, dst_sx - 1) };
  int iy[3] = { std::max(y - 1, 0), y, std::min(y + 1, dst_sy - 1) };

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      dst[iy[i] * dst_sx + ix[j]] += value * GAUSSIAN_INTEGRALS[sub_bin_x][j] * GAUSSIAN_INTEGRALS[sub_bin_y][i];
}

void
smoothIncrement3D(Real * dst, int dst_sx, int dst_sy, int dst_sz, int x, int y, int z, Vector3 const & mean, Real value)
{
  int sub_bin_x = Math::clamp((int)std::floor(NUM_GAUSSIAN_INTEGRAL_BINS * mean.x()), 0, NUM_GAUSSIAN_INTEGRAL_BINS - 1);
  int sub_bin_y = Math::clamp((int)std::floor(NUM_GAUSSIAN_INTEGRAL_BINS * mean.y()), 0, NUM_GAUSSIAN_INTEGRAL_BINS - 1);
  int sub_bin_z = Math::clamp((int)std::floor(NUM_GAUSSIAN_INTEGRAL_BINS * mean.z()), 0, NUM_GAUSSIAN_INTEGRAL_BINS - 1);

  int ix[3] = { std::max(x - 1, 0), x, std::min(x + 1, dst_sx - 1) };
  int iy[3] = { std::max(y - 1, 0), y, std::min(y + 1, dst_sy - 1) };
  int iz[3] = { std::max(z - 1, 0), z, std::min(z + 1, dst_sz - 1) };

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        dst[(iz[i] * dst_sy + iy[j]) * dst_sx + ix[k]] += value * GAUSSIAN_INTEGRALS[sub_bin_x][k]
                                                                * GAUSSIAN_INTEGRALS[sub_bin_y][j]
                                                                * GAUSSIAN_INTEGRALS[sub_bin_z][i];
}

inline Real toReal(Real            t) { return t;   }
inline Real toReal(Vector2 const & t) { return t.x(); }
inline Real toReal(Vector3 const & t) { return t.x(); }

inline Real    add1D(Real            t, Real inc) { return t + inc;                            }
inline Vector2 add1D(Vector2 const & t, Real inc) { return Vector2(t.x() + inc, t.y());        }
inline Vector3 add1D(Vector3 const & t, Real inc) { return Vector3(t.x() + inc, t.y(), t.z()); }

inline Vector2 add2D(Vector2 const & t, Vector2 const & inc) { return t + inc;                                          }
inline Vector3 add2D(Vector3 const & t, Vector2 const & inc) { return Vector3(t.x() + inc.x(), t.y() + inc.y(), t.z()); }

inline Real    halve1D(Real            t) { return 0.5f * t;                            }
inline Vector2 halve1D(Vector2 const & t) { return Vector2(0.5f * t.x(), t.y());        }
inline Vector3 halve1D(Vector3 const & t) { return Vector3(0.5f * t.x(), t.y(), t.z()); }

inline Vector2 halve2D(Vector2 const & t) { return 0.5f * t;                                   }
inline Vector3 halve2D(Vector3 const & t) { return Vector3(0.5f * t.x(), 0.5f * t.y(), t.z()); }

template <typename MeanT>
void
downsample1D(Real const * src, MeanT const * src_means, Real * dst, MeanT * dst_means, int dst_size)
{
  if (src_means)
  {
    Real w0, w1, sum_weights;

    for (int i = 0, j = 0; i < dst_size; ++i, j += 2)
    {
      w0 = src[j];
      w1 = src[j + 1];
      sum_weights = w0 + w1;

      dst_means[i] = Math::fuzzyEq(sum_weights, (Real)0) ? 0.5f * halve1D(src_means[j]
                                                                        + add1D(src_means[j + 1], 1))
                                                         : halve1D(w0 * src_means[j]
                                                                 + w1 * add1D(src_means[j + 1], 1)) / sum_weights;
      smoothIncrement1D(dst, dst_size, (int)i, toReal(dst_means[i]), sum_weights);
    }
  }
  else
  {
    for (int i = 0, j = 0; i < dst_size; ++i, j += 2)
      dst[i] = src[j] + src[j + 1];
  }
}

template <typename MeanT>
void
downsample2D(Real const * src, MeanT const * src_means, Real * dst, MeanT * dst_means, int dst_sx, int dst_sy)
{
  static Vector2 const offset01(0, 1);
  static Vector2 const offset10(1, 0);
  static Vector2 const offset11(1, 1);

  int src_sx = 2 * dst_sx;
  long i = 0;
  long j00 = 0, j10 = 1, j01 = src_sx, j11 = src_sx + 1;

  if (src_means)
  {
    Real w00, w10, w01, w11, sum_weights;

    for(int y = 0; y < dst_sy; ++y, j00 += src_sx, j10 += src_sx, j01 += src_sx, j11 += src_sx)
      for(int x = 0; x < dst_sx; ++x, ++i, j00 += 2, j10 += 2, j01 += 2, j11 += 2)
      {
        w00 = src[j00];
        w10 = src[j10];
        w01 = src[j01];
        w11 = src[j11];
        sum_weights = w00 + w10 + w01 + w11;

        dst_means[i] = Math::fuzzyEq(sum_weights, (Real)0) ? 0.25f * halve2D(      src_means[j00]
                                                                           + add2D(src_means[j10], offset10)
                                                                           + add2D(src_means[j01], offset01)
                                                                           + add2D(src_means[j11], offset11))
                                                           : halve2D(w00 *       src_means[j00]
                                                                   + w10 * add2D(src_means[j10], offset10)
                                                                   + w01 * add2D(src_means[j01], offset01)
                                                                   + w11 * add2D(src_means[j11], offset11)) / sum_weights;

        smoothIncrement2D(dst, dst_sx, dst_sy, x, y, dst_means[i].x(), dst_means[i].y(), sum_weights);
      }
  }
  else
  {
    for(int y = 0; y < dst_sy; ++y, j00 += src_sx, j10 += src_sx, j01 += src_sx, j11 += src_sx)
      for(int x = 0; x < dst_sx; ++x, ++i, j00 += 2, j10 += 2, j01 += 2, j11 += 2)
        dst[i] = src[j00] + src[j10] + src[j01] + src[j11];
  }
}

void
downsample3D(Real const * src, Vector3 const * src_means, Real * dst, Vector3 * dst_means, int dst_sx, int dst_sy, int dst_sz)
{
  static Vector3 const offset001(0, 0, 1);
  static Vector3 const offset010(0, 1, 0);
  static Vector3 const offset011(0, 1, 1);
  static Vector3 const offset100(1, 0, 0);
  static Vector3 const offset101(1, 0, 1);
  static Vector3 const offset110(1, 1, 0);
  static Vector3 const offset111(1, 1, 1);

  int src_sx = 2 * dst_sx;
  int src_sxy = src_sx * (2 * dst_sy);
  long i = 0;
  long j000 = 0,       j100 = 1,           j010 = src_sx,           j110 = src_sx + 1;
  long j001 = src_sxy, j101 = src_sxy + 1, j011 = src_sxy + src_sx, j111 = src_sxy + src_sx + 1;

  if (src_means)
  {
    Real w000, w100, w010, w110, w001, w101, w011, w111, sum_weights;

    for(long z = 0; z < dst_sz; ++z, j000 += src_sxy, j100 += src_sxy, j010 += src_sxy, j110 += src_sxy,
                                     j001 += src_sxy, j101 += src_sxy, j011 += src_sxy, j111 += src_sxy)
      for(long y = 0; y < dst_sy; ++y, j000 += src_sx, j100 += src_sx, j010 += src_sx, j110 += src_sx,
                                       j001 += src_sx, j101 += src_sx, j011 += src_sx, j111 += src_sx)
        for(long x = 0; x < dst_sx; ++x, ++i, j000 += 2, j100 += 2, j010 += 2, j110 += 2,
                                              j001 += 2, j101 += 2, j011 += 2, j111 += 2)
        {
          w000 = src[j000];
          w100 = src[j100];
          w010 = src[j010];
          w110 = src[j110];
          w001 = src[j001];
          w101 = src[j101];
          w011 = src[j011];
          w111 = src[j111];
          sum_weights = w000 + w100 + w010 + w110 + w001 + w101 + w011 + w111;

          dst_means[i] = Math::fuzzyEq(sum_weights, (Real)0) ? 0.125f * 0.5f * (src_means[j000]
                                                                              + src_means[j100] + offset100
                                                                              + src_means[j010] + offset010
                                                                              + src_means[j110] + offset110
                                                                              + src_means[j001] + offset001
                                                                              + src_means[j101] + offset101
                                                                              + src_means[j011] + offset011
                                                                              + src_means[j111] + offset111)
                                                             : (0.5f / sum_weights) * (w000 * src_means[j000]
                                                                                     + w100 * src_means[j100] + offset100
                                                                                     + w010 * src_means[j010] + offset010
                                                                                     + w110 * src_means[j110] + offset110
                                                                                     + w001 * src_means[j001] + offset001
                                                                                     + w101 * src_means[j101] + offset101
                                                                                     + w011 * src_means[j011] + offset011
                                                                                     + w111 * src_means[j111] + offset111);

          smoothIncrement3D(dst, dst_sx, dst_sy, dst_sz, x, y, z, dst_means[i], sum_weights);
        }
  }
  else
  {
    for(long z = 0; z < dst_sz; ++z, j000 += src_sxy, j100 += src_sxy, j010 += src_sxy, j110 += src_sxy,
                                     j001 += src_sxy, j101 += src_sxy, j011 += src_sxy, j111 += src_sxy)
      for(long y = 0; y < dst_sy; ++y, j000 += src_sx, j100 += src_sx, j010 += src_sx, j110 += src_sx,
                                       j001 += src_sx, j101 += src_sx, j011 += src_sx, j111 += src_sx)
        for(long x = 0; x < dst_sx; ++x, ++i, j000 += 2, j100 += 2, j010 += 2, j110 += 2,
                                              j001 += 2, j101 += 2, j011 += 2, j111 += 2)
          dst[i] = src[j000] + src[j100] + src[j010] + src[j110]
                 + src[j001] + src[j101] + src[j011] + src[j111];
  }
}

} // namespace PyramidInternal

void
Pyramid1D::createPyramid(Real * means)
{
  for(int i = 1; i < num_levels; ++i)
    downsample(levels[i - 1], levels[i], means);
}

void
Pyramid2D::createPyramid(Vector2 * means)
{
  for(int i = 1; i < num_levels; ++i)
    downsample(levels[i - 1], levels[i], 2, means);
}

void
Pyramid3D::createPyramid(Vector3 * means)
{
  for(int i = 1; i < num_levels; ++i)
    downsample(levels[i - 1], levels[i], 3, means);
}

void
Pyramid1D::downsample(TheaArray<Real> const & src, TheaArray<Real> & dst, Real * means) const
{
  array_size_t dst_size = src.size() / 2;

  if (means)
    dst.resize(dst_size, 0);
  else
    dst.resize(dst_size);

  PyramidInternal::downsample1D(&src[0], means, &dst[0], means, (int)dst_size);
}

void
Pyramid2D::downsample(Array2D const & src, Array2D & dst, int num_dims_to_downsample, Vector2 * means) const
{
  alwaysAssertM(num_dims_to_downsample > 0, "Pyramid2D: At least one dimension must reduce in size during downsampling");

  int dst_sx = src.sx / 2;
  int dst_sy = num_dims_to_downsample > 1 ? src.sy / 2 : src.sy;

  if (means)
    dst = Array2D(dst_sx, dst_sy, 0);
  else
    dst = Array2D(dst_sx, dst_sy);

  if (num_dims_to_downsample >= 2)
    PyramidInternal::downsample2D(&src.data[0], means, &dst.data[0], means, dst_sx, dst_sy);
  else
  {
    Real const * src_ptr = &src.data[0];
    Real * dst_ptr = &dst.data[0];

    if (means)
    {
      Vector2 * src_means_ptr = means;
      Vector2 * dst_means_ptr = means;

      for (int y = 0; y < dst_sy; ++y, src_ptr += src.sx, src_means_ptr += src.sx, dst_ptr += dst_sx, dst_means_ptr += dst_sx)
        PyramidInternal::downsample1D(src_ptr, src_means_ptr, dst_ptr, dst_means_ptr, dst_sx);
    }
    else
    {
      for (int y = 0; y < dst_sy; ++y, src_ptr += src.sx, dst_ptr += dst_sx)
        PyramidInternal::downsample1D(src_ptr, means, dst_ptr, means, dst_sx);
    }
  }
}

void
Pyramid3D::downsample(Array3D const & src, Array3D & dst, int num_dims_to_downsample, Vector3 * means) const
{
  alwaysAssertM(num_dims_to_downsample > 0, "Pyramid3D: At least one dimension must reduce in size during downsampling");

  int dst_sx = src.sx / 2;
  int dst_sy = num_dims_to_downsample > 1 ? src.sy / 2 : src.sy;
  int dst_sz = num_dims_to_downsample > 2 ? src.sz / 2 : src.sz;

  if (means)
    dst = Array3D(dst_sx, dst_sy, dst_sz, 0);
  else
    dst = Array3D(dst_sx, dst_sy, dst_sz);

  if (num_dims_to_downsample >= 3)
    PyramidInternal::downsample3D(&src.data[0], means, &dst.data[0], means, dst_sx, dst_sy, dst_sz);
  else if (num_dims_to_downsample == 2)
  {
    Real const * src_ptr = &src.data[0];
    Real * dst_ptr = &dst.data[0];
    long src_step = src.sx * src.sy;
    long dst_step = dst_sx * dst_sy;

    if (means)
    {
      Vector3 * src_means_ptr = means;
      Vector3 * dst_means_ptr = means;

      for (int z = 0; z < dst_sz; ++z, src_ptr += src_step, src_means_ptr += src_step, dst_ptr += dst_step,
                                       dst_means_ptr += dst_step)
        PyramidInternal::downsample2D(src_ptr, src_means_ptr, dst_ptr, dst_means_ptr, dst_sx, dst_sy);
    }
    else
    {
      for (int z = 0; z < dst_sz; ++z, src_ptr += src_step, dst_ptr += dst_step)
        PyramidInternal::downsample2D(src_ptr, means, dst_ptr, means, dst_sx, dst_sy);
    }
  }
  else
  {
    Real const * src_ptr = &src.data[0];
    Real * dst_ptr = &dst.data[0];

    if (means)
    {
      Vector3 * src_means_ptr = means;
      Vector3 * dst_means_ptr = means;

      for (int z = 0; z < dst_sz; ++z)
        for (int y = 0; y < dst_sy; ++y, src_ptr += src.sx, src_means_ptr += src.sx, dst_ptr += dst_sx, dst_means_ptr += dst_sx)
          PyramidInternal::downsample1D(src_ptr, src_means_ptr, dst_ptr, dst_means_ptr, dst_sx);
    }
    else
    {
      for (int z = 0; z < dst_sz; ++z)
        for (int y = 0; y < dst_sy; ++y, src_ptr += src.sx, dst_ptr += dst_sx)
          PyramidInternal::downsample1D(src_ptr, means, dst_ptr, means, dst_sx);
    }
  }
}

void
Pyramid1D::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  output.writeUInt32(1);                   // dimensions
  output.writeUInt32((uint32)num_levels);  // number of levels

  for (int i = 0; i < num_levels; ++i)
  {
    output.writeUInt32((uint32)levels[i].size());  // number of bins in this level

    TheaArray<Real> const & level_data = levels[i];
    for (array_size_t j = 0; j < level_data.size(); ++j)
      output.writeFloat64(level_data[j]);
  }
}

void
Pyramid1D::deserialize(BinaryInputStream & input, Codec const & codec)
{
  uint32 dims = input.readUInt32();
  if (dims != 1) throw Error("Pyramid1D: Number of dimensions must be 1");

  num_levels = input.readUInt32();
  levels.clear();
  if (num_levels > 0) levels.resize(num_levels);

  for (int i = 0; i < num_levels; ++i)
  {
    array_size_t level_size = (array_size_t)input.readUInt32();
    levels[i].resize(level_size);

    TheaArray<Real> & level_data = levels[i];
    for (array_size_t j = 0; j < level_data.size(); ++j)
      level_data[j] = input.readFloat64();
  }
}

void
Pyramid1D::serialize(TextOutputStream & output, Codec const & codec) const
{
  output.printf("1\n");               // dimensions
  output.printf("%d\n", num_levels);  // number of levels
  output.writeNewline();

  for (int i = 0; i < num_levels; ++i)
  {
    output.printf("%d\n", (int)levels[i].size());  // number of bins in this level

    TheaArray<Real> const & level_data = levels[i];
    for (array_size_t j = 0; j < level_data.size(); ++j)
      output.printf("%lf\n", level_data[j]);

    output.writeNewline();
  }
}

void
Pyramid1D::deserialize(TextInputStream & input, Codec const & codec)
{
  int dims = (int)input.readNumber();
  if (dims != 1) throw Error("Pyramid1D: Number of dimensions must be 1");

  num_levels = (int)input.readNumber();
  levels.clear();
  if (num_levels > 0) levels.resize(num_levels);

  for (int i = 0; i < num_levels; ++i)
  {
    array_size_t level_size = (array_size_t)input.readNumber();
    levels[i].resize(level_size);

    TheaArray<Real> & level_data = levels[i];
    for (array_size_t j = 0; j < level_data.size(); ++j)
      level_data[j] = input.readNumber();
  }
}

void
Pyramid2D::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  output.writeUInt32(2);                   // dimensions
  output.writeUInt32((uint32)num_levels);  // number of levels

  for (int i = 0; i < num_levels; ++i)
  {
    Array2D const & level = levels[i];
    output.writeUInt32((uint32)level.sx);
    output.writeUInt32((uint32)level.sy);

    TheaArray<Real> const & level_data = level.data;
    for (array_size_t j = 0; j < level_data.size(); ++j)
      output.writeFloat64(level_data[j]);
  }
}

void
Pyramid2D::deserialize(BinaryInputStream & input, Codec const & codec)
{
  uint32 dims = input.readUInt32();
  if (dims != 2) throw Error("Pyramid2D: Number of dimensions must be 2");

  num_levels = input.readUInt32();
  levels.clear();
  if (num_levels > 0) levels.resize(num_levels);

  nx = ny = 0;

  for (int i = 0; i < num_levels; ++i)
  {
    Array2D & level = levels[i];
    level.sx = (array_size_t)input.readUInt32();
    level.sy = (array_size_t)input.readUInt32();

    if (i == 0)
    {
      nx = level.sx;
      ny = level.sy;
    }

    TheaArray<Real> & level_data = level.data;
    level_data.resize(level.sx * level.sy);

    for (array_size_t j = 0; j < level_data.size(); ++j)
      level_data[j] = input.readFloat64();
  }
}

void
Pyramid2D::serialize(TextOutputStream & output, Codec const & codec) const
{
  output.printf("2\n");               // dimensions
  output.printf("%d\n", num_levels);  // number of levels
  output.writeNewline();

  for (int i = 0; i < num_levels; ++i)
  {
    Array2D const & level = levels[i];
    output.printf("%d\n", (int)level.sx);
    output.printf("%d\n", (int)level.sy);

    TheaArray<Real> const & level_data = level.data;
    for (array_size_t j = 0; j < level_data.size(); ++j)
      output.printf("%lf\n", level_data[j]);
  }
}

void
Pyramid2D::deserialize(TextInputStream & input, Codec const & codec)
{
  int dims = (int)input.readNumber();
  if (dims != 2) throw Error("Pyramid2D: Number of dimensions must be 2");

  num_levels = (int)input.readNumber();
  levels.clear();
  if (num_levels > 0) levels.resize(num_levels);

  nx = ny = 0;

  for (int i = 0; i < num_levels; ++i)
  {
    Array2D & level = levels[i];
    level.sx = (int)input.readNumber();
    level.sy = (int)input.readNumber();

    if (i == 0)
    {
      nx = level.sx;
      ny = level.sy;
    }

    TheaArray<Real> & level_data = level.data;
    level_data.resize(level.sx * level.sy);

    for (array_size_t j = 0; j < level_data.size(); ++j)
      level_data[j] = input.readNumber();
  }
}

void
Pyramid3D::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  output.writeUInt32(3);                   // dimensions
  output.writeUInt32((uint32)num_levels);  // number of levels

  for (int i = 0; i < num_levels; ++i)
  {
    Array3D const & level = levels[i];
    output.writeUInt32((uint32)level.sx);
    output.writeUInt32((uint32)level.sy);
    output.writeUInt32((uint32)level.sz);

    TheaArray<Real> const & level_data = level.data;
    for (array_size_t j = 0; j < level_data.size(); ++j)
      output.writeFloat64(level_data[j]);
  }
}

void
Pyramid3D::deserialize(BinaryInputStream & input, Codec const & codec)
{
  uint32 dims = input.readUInt32();
  if (dims != 3) throw Error("Pyramid3D: Number of dimensions must be 3");

  num_levels = input.readUInt32();
  levels.clear();
  if (num_levels > 0) levels.resize(num_levels);

  for (int i = 0; i < num_levels; ++i)
  {
    Array3D & level = levels[i];
    level.sx = (array_size_t)input.readUInt32();
    level.sy = (array_size_t)input.readUInt32();
    level.sz = (array_size_t)input.readUInt32();

    if (i == 0)
    {
      nx = level.sx;
      ny = level.sy;
      nz = level.sz;
    }

    TheaArray<Real> & level_data = level.data;
    level_data.resize(level.sx * level.sy * level.sz);

    for (array_size_t j = 0; j < level_data.size(); ++j)
      level_data[j] = input.readFloat64();
  }
}

void
Pyramid3D::serialize(TextOutputStream & output, Codec const & codec) const
{
  output.printf("3\n");               // dimensions
  output.printf("%d\n", num_levels);  // number of levels
  output.writeNewline();

  for (int i = 0; i < num_levels; ++i)
  {
    Array3D const & level = levels[i];
    output.printf("%d\n", (int)level.sx);
    output.printf("%d\n", (int)level.sy);
    output.printf("%d\n", (int)level.sz);

    TheaArray<Real> const & level_data = level.data;
    for (array_size_t j = 0; j < level_data.size(); ++j)
      output.printf("%lf\n", level_data[j]);
  }
}

void
Pyramid3D::deserialize(TextInputStream & input, Codec const & codec)
{
  int dims = (int)input.readNumber();
  if (dims != 3) throw Error("Pyramid3D: Number of dimensions must be 3");

  num_levels = (int)input.readNumber();
  levels.clear();
  if (num_levels > 0) levels.resize(num_levels);

  for (int i = 0; i < num_levels; ++i)
  {
    Array3D & level = levels[i];
    level.sx = (int)input.readNumber();
    level.sy = (int)input.readNumber();
    level.sz = (int)input.readNumber();

    if (i == 0)
    {
      nx = level.sx;
      ny = level.sy;
      nz = level.sz;
    }

    TheaArray<Real> & level_data = level.data;
    level_data.resize(level.sx * level.sy * level.sz);

    for (array_size_t j = 0; j < level_data.size(); ++j)
      level_data[j] = input.readNumber();
  }
}

namespace PyramidInternal {

inline Real
similarity(PyramidMatch::Kernel kernel, Real x, Real y)
{
  switch (kernel)
  {
    case PyramidMatch::Kernel::L1:           return -std::fabs(x - y);
    case PyramidMatch::Kernel::L2:           return -Math::square(x - y);
    case PyramidMatch::Kernel::CHI_SQUARED:  return -Math::square(x - y) / (x < 1 ? 1 : x);
    default: /* HIK */                       return  std::min(x, y);
  }
}

} // namespace PyramidInternal

Real
PyramidMatch::similarity(Pyramid1D const & pyramid1, Pyramid1D const & pyramid2, Kernel kernel, int levels,
                         Real attenuation_factor)
{
  alwaysAssertM(pyramid1.num_levels == pyramid2.num_levels, "PyramidMatch: Can't compare pyramids of different dimensions");

  int num_levels = (levels > 0 ? levels : pyramid1.num_levels);
  if (attenuation_factor <= 0) attenuation_factor = 0.5;

  Real ans = 0;
  Real scale = 1;
  for (int i = 0; i < num_levels; ++i)
  {
    TheaArray<Real> const & level_data1 = pyramid1.levels[i];
    TheaArray<Real> const & level_data2 = pyramid2.levels[i];

    debugAssertM(level_data1.size() == level_data2.size(), "PyramidMatch: Levels of comparable pyramids should have same size");

    Real level_ans = 0;
    for(array_size_t j = 0; j < level_data1.size(); ++ j)
      level_ans += PyramidInternal::similarity(kernel, level_data1[j], level_data2[j]);

    ans += (level_ans * scale);
    scale *= attenuation_factor;
  }

  return ans;
}

Real
PyramidMatch::similarity(Pyramid2D const & pyramid1, Pyramid2D const & pyramid2, Kernel kernel, int levels,
                         Real attenuation_factor)
{
  alwaysAssertM(pyramid1.nx == pyramid2.nx && pyramid1.ny == pyramid2.ny,
                "PyramidMatch: Can't compare pyramids of different dimensions");

  int num_levels = (levels > 0 ? levels : pyramid1.num_levels);
  if (attenuation_factor <= 0) attenuation_factor = 0.5;

  Real ans = 0;
  Real scale = 1;
  for (int i = 0; i < num_levels; ++i)
  {
    TheaArray<Real> const & level_data1 = pyramid1.levels[i].data;
    TheaArray<Real> const & level_data2 = pyramid2.levels[i].data;

    debugAssertM(level_data1.size() == level_data2.size(), "PyramidMatch: Levels of comparable pyramids should have same size");

    Real level_ans = 0;
    for(array_size_t j = 0; j < level_data1.size(); ++ j)
      level_ans += PyramidInternal::similarity(kernel, level_data1[j], level_data2[j]);

    ans += (level_ans * scale);
    scale *= attenuation_factor;
  }

  return ans;
}

Real
PyramidMatch::similarity(Pyramid3D const & pyramid1, Pyramid3D const & pyramid2, Kernel kernel, int levels,
                         Real attenuation_factor)
{
  alwaysAssertM(pyramid1.nx == pyramid2.nx && pyramid1.ny == pyramid2.ny && pyramid1.nz == pyramid2.nz,
                "PyramidMatch: Can't compare pyramids of different dimensions");

  int num_levels = (levels > 0 ? levels : pyramid1.num_levels);
  if (attenuation_factor <= 0) attenuation_factor = 0.5;

  Real ans = 0;
  Real scale = 1;
  for (int i = 0; i < num_levels; ++i)
  {
    TheaArray<Real> const & level_data1 = pyramid1.levels[i].data;
    TheaArray<Real> const & level_data2 = pyramid2.levels[i].data;

    debugAssertM(level_data1.size() == level_data2.size(), "PyramidMatch: Levels of comparable pyramids should have same size");

    Real level_ans = 0;
    for(array_size_t j = 0; j < level_data1.size(); ++ j)
      level_ans += PyramidInternal::similarity(kernel, level_data1[j], level_data2[j]);

    ans += (level_ans * scale);
    scale *= attenuation_factor;
  }

  return ans;
}

static void
unitTest1()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 1" << std::endl;

  TheaArray<int> x;
  x.push_back(1);

  Pyramid1D p1(x), p2(x);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 1;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 1 failed");
}

static void
unitTest2()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 2" << std::endl;

  TheaArray<int> x;
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);

  Pyramid1D p1(x), p2(x);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 7;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 2 failed");
}

static void
unitTest3()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 3" << std::endl;

  TheaArray<int> x, y;
  for(int i = 0; i < 8; ++i) {
    x.push_back(i);
    y.push_back(7 - i);
  }

  Pyramid1D p1(x), p2(y);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 24.5;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 3 failed");
}

static void
unitTest4()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 4" << std::endl;

  TheaArray<int> x;
  x.push_back(1);

  Pyramid2D p1(&x[0], 1, 1), p2(&x[0], 1, 1);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 1;
  alwaysAssertM(Math::fuzzyEq(result, expected), "PyramidMatch: Unit 4 test failed");
}

static void
unitTest5()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 5" << std::endl;

  TheaArray<int> x;
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);

  Pyramid2D p1(&x[0], 2, 2), p2(&x[0], 2, 2);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 6;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 5 failed");
}

static void
unitTest6()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 6" << std::endl;

  TheaArray<int> x, y;
  x.push_back(1);
  x.push_back(2);
  x.push_back(3);
  x.push_back(4);
  y.push_back(4);
  y.push_back(3);
  y.push_back(2);
  y.push_back(1);

  Pyramid2D p1(&x[0], 2, 2), p2(&y[0], 2, 2);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 11;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 6 failed");
}

static void
unitTest7()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 7" << std::endl;

  TheaArray<int> x;
  x.push_back(1);

  Pyramid3D p1(&x[0], 1, 1, 1), p2(&x[0], 1, 1, 1);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 1;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 7 failed");
}

static void
unitTest8()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 8" << std::endl;

  TheaArray<int> x;
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);
  x.push_back(1);

  Pyramid3D p1(&x[0], 2, 2, 2), p2(&x[0], 2, 2, 2);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 12;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 8 failed");
}

static void
unitTest9()
{
  std::cout << "=============================\n"
               "Pyramid Match: Unit test 9" << std::endl;

  TheaArray<int> x, y;
  for(int i = 0; i < 8; ++i) {
    x.push_back(i + 1);
    y.push_back(8 - i);
  }

  Pyramid3D p1(&x[0], 2, 2, 2), p2(&y[0], 2, 2, 2);

  Real result = PyramidMatch::similarity(p1, p2, PyramidMatch::Kernel::HIK, -1, 0.5f);
  std::cout << "Result = " << result << '\n' << std::endl;

  Real expected = 38;
  if (!Math::fuzzyEq(result, expected))
    throw Error("PyramidMatch: Unit test 9 failed");
}

void
PyramidMatch::test()
{
  unitTest1();
  unitTest2();
  unitTest3();
  unitTest4();
  unitTest5();
  unitTest6();
  unitTest7();
  unitTest8();
  unitTest9();
}

} // namespace Algorithms
} // namespace Thea
