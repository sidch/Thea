//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#include "SVD.hpp"
#include "../Math.hpp"
#include <cmath>

namespace Thea {
namespace Algorithms {
namespace SVDInternal {

// Code from G3D Matrix::svd(). Based on Dianne Cook's implementation, which is adapted from svdecomp.c in XLISP-STAT 2.1, which
// is code from Numerical Recipes adapted by Luke Tierney and David Betz. The Numerical Recipes code is adapted from Forsythe et
// al., who based their code on Golub and Reinsch's original implementation.

// Helper for svdCore
double
pythag(double a, double b)
{
  double at = std::fabs(a), bt = std::fabs(b), ct, result;

  if (at > bt)
  {
    ct = bt / at;
    result = at * std::sqrt(1.0 + Math::square(ct));
  }
  else if (bt > 0.0)
  {
    ct = at / bt;
    result = bt * std::sqrt(1.0 + Math::square(ct));
  }
  else
    result = 0.0;

  return result;
}

#define SVD_SIGN(a, b) ((b) >= 0.0 ? std::fabs(a) : -std::fabs(a))

template <typename R, typename S, typename T>
bool
svdCoreT(AddressableMatrix<R> & U, long rows, long cols, S * D, AddressableMatrix<T> & V)
{
  static int const MAX_ITERATIONS = 30;

  if (rows == 0 && cols == 0)
    return true;

  debugAssertM(rows >= cols, "SVD: Matrix must have more rows than columns");

  long i, j, jj, k, l = 0, nm = 0;
  int flag, its;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;

  // Temp row vector
  TheaArray<double> rv1((array_size_t)cols);

  // Householder reduction to bidiagonal form
  for (i = 0; i < cols; ++i)
  {
    // Left-hand reduction
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;

    if (i < rows)
    {

      for (k = i; k < rows; ++k)
      {
        scale += std::fabs((double)U.get(k, i));
      }

      if (scale)
      {
        for (k = i; k < rows; ++k)
        {
          double x = (double)U.get(k, i) / scale;
          U.set(k, i, x);
          s += (x * x);
        }

        f = (double)U.get(i, i);

        g = -SVD_SIGN(std::sqrt(s), f);
        h = f * g - s;
        U.set(i, i, (R)(f - g));

        if (i != cols - 1)
        {
          for (j = l; j < cols; j++)
          {
            for (s = 0.0, k = i; k < rows; ++k)
            {
              s += ((double)U.get(k, i) * (double)U.get(k, j));
            }

            f = s / h;
            for (k = i; k < rows; ++k)
            {
              U.getMutable(k, j) += (R)(f * (double)U.get(k, i));
            }
          }
        }
        for (k = i; k < rows; ++k)
        {
          U.getMutable(k, i) *= (R)scale;
        }
      }
    }
    D[i] = (S)(scale * g);

    // right-hand reduction
    g = s = scale = 0.0;
    if (i < rows && i != cols - 1)
    {
      for (k = l; k < cols; ++k)
      {
        scale += std::fabs((double)U.get(i, k));
      }

      if (scale)
      {
        for (k = l; k < cols; ++k)
        {
          double x = (double)U.get(i, k) / scale;
          U.set(i, k, (R)x);
          s += (x * x);
        }

        f = (double)U.get(i, l);
        g = -SVD_SIGN(std::sqrt(s), f);
        h = f * g - s;
        U.set(i, l, (R)(f - g));

        for (k = l; k < cols; ++k)
        {
          rv1[k] = (double)U.get(i, k) / h;
        }

        if (i != rows - 1)
        {
          for (j = l; j < rows; ++j)
          {
            for (s = 0.0, k = l; k < cols; ++k)
            {
              s += ((double)U.get(j, k) * (double)U.get(i, k));
            }

            for (k = l; k < cols; ++k)
            {
              U.getMutable(j, k) += (R)(s * rv1[k]);
            }
          }
        }

        for (k = l; k < cols; ++k)
        {
          U.getMutable(i, k) *= (R)scale;
        }
      }
    }

    anorm = std::max(anorm, std::fabs((double)D[i]) + std::fabs(rv1[i]));
  }

  // accumulate the right-hand transformation
  for (i = cols - 1; i >= 0; --i)
  {
    if (i < cols - 1)
    {
      if (g != 0.0)
      {
        for (j = l; j < cols; j++)
        {
          V.set(j, i, (T)(((double)U.get(i, j) / (double)U.get(i, l)) / g));
        }

        // double division to avoid underflow
        for (j = l; j < cols; ++j)
        {
          for (s = 0.0, k = l; k < cols; k++) {
            s += ((double)U.get(i, k) * (double)V.get(k, j));
          }

          for (k = l; k < cols; ++k)
          {
            V.getMutable(k, j) += (T)(s * (double)V.get(k, i));
          }
        }
      }

      for (j = l; j < cols; ++j)
      {
        V.set(i, j, 0.0);
        V.set(j, i, 0.0);
      }
    }

    V.set(i, i, 1.0);
    g = rv1[i];
    l = i;
  }

  // accumulate the left-hand transformation
  for (i = cols - 1; i >= 0; --i)
  {
    l = i + 1;
    g = (double)D[i];
    if (i < cols - 1)
    {
      for (j = l; j < cols; ++j)
      {
        U.set(i, j, 0.0);
      }
    }

    if (g)
    {
      g = 1.0 / g;
      if (i != cols - 1)
      {
        for (j = l; j < cols; ++j)
        {
          for (s = 0.0, k = l; k < rows; ++k)
          {
            s += ((double)U.get(k, i) * (double)U.get(k, j));
          }

          f = (s / (double)U.get(i, i)) * g;

          for (k = i; k < rows; ++k)
          {
            U.getMutable(k, j) += (R)(f * (double)U.get(k, i));
          }
        }
      }

      for (j = i; j < rows; ++j)
      {
        U.getMutable(j, i) *= (R)g;
      }
    }
    else
    {
      for (j = i; j < rows; ++j)
      {
        U.set(j, i, 0.0);
      }
    }
    U.getMutable(i, i)++;
  }

  // diagonalize the bidiagonal form
  for (k = cols - 1; k >= 0; --k)
  {
    // loop over singular values
    for (its = 0; its < MAX_ITERATIONS; ++its)
    {
      // loop over allowed iterations
      flag = 1;

      for (l = k; l >= 0; --l)
      {
        // test for splitting
        nm = l - 1;
        if (std::fabs(rv1[l]) + anorm == anorm)
        {
          flag = 0;
          break;
        }

        if (std::fabs((double)D[nm]) + anorm == anorm)
        {
          break;
        }
      }

      if (flag)
      {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; ++i)
        {
          f = s * rv1[i];
          if (std::fabs(f) + anorm != anorm)
          {
            g = (double)D[i];
            h = pythag(f, g);
            D[i] = (S)h;
            h = 1.0 / h;
            c = g * h;
            s = (- f * h);
            for (j = 0; j < rows; ++j)
            {
              y = (double)U.get(j, nm);
              z = (double)U.get(j, i);
              U.set(j, nm, (R)(y * c + z * s));
              U.set(j, i, (R)(z * c - y * s));
            }
          }
        }
      }

      z = (double)D[k];
      if (l == k)
      {
        // convergence
        if (z < 0.0)
        {
          // make singular value nonnegative
          D[k] = (S)(-z);

          for (j = 0; j < cols; ++j)
          {
            V.set(j, k, -V.get(j, k));
          }
        }
        break;
      }

      if (its >= MAX_ITERATIONS)
      {
        THEA_DEBUG << "SVD: Failed to converge";
        return false;
      }

      // shift from bottom 2 x 2 minor
      x = (double)D[l];
      nm = k - 1;
      y = (double)D[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = pythag(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SVD_SIGN(g, f))) - h)) / x;

      // next QR transformation
      c = s = 1.0;
      for (j = l; j <= nm; ++j)
      {
        i = j + 1;
        g = rv1[i];
        y = (double)D[i];
        h = s * g;
        g = c * g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;

        for (jj = 0; jj < cols; ++jj)
        {
          x = (double)V.get(jj, j);
          z = (double)V.get(jj, i);
          V.set(jj, j, (T)(x * c + z * s));
          V.set(jj, i, (T)(z * c - x * s));
        }
        z = pythag(f, h);
        D[j] = (S)z;
        if (z)
        {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = (c * g) + (s * y);
        x = (c * y) - (s * g);
        for (jj = 0; jj < rows; jj++)
        {
          y = (double)U.get(jj, j);
          z = (double)U.get(jj, i);
          U.set(jj, j, (R)(y * c + z * s));
          U.set(jj, i, (R)(z * c - y * s));
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      D[k] = (S)x;
    }
  }

  // So far V actually contains V^T, so we transpose V in-place
  for (i = 0; i < rows; ++i)
    for (j = i + 1; j < cols; ++j)
      std::swap(V.getMutable(i, j), V.getMutable(j, i));

  return true;
}

#undef SVD_SIGN

// Define non-template functions for all legal combinations of floating-point types
#define SVD_DEF_SVDCORE(R, S, T) \
  bool svdCore(AddressableMatrix<R> & U,  long rows, long cols, S * D, AddressableMatrix<T> & V) \
  { return svdCoreT(U, rows, cols, D, V); }

SVD_DEF_SVDCORE(float,  float,  float )
SVD_DEF_SVDCORE(float,  float,  double)
SVD_DEF_SVDCORE(float,  double, float )
SVD_DEF_SVDCORE(float,  double, double)
SVD_DEF_SVDCORE(double, float,  float )
SVD_DEF_SVDCORE(double, float,  double)
SVD_DEF_SVDCORE(double, double, float )
SVD_DEF_SVDCORE(double, double, double)

} // namespace SVDInternal

} // namespace Algorithms
} // namespace Thea
