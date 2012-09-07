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

#ifndef __Thea_Algorithms_LinearScalarFunction_hpp__
#define __Thea_Algorithms_LinearScalarFunction_hpp__

#include "AnalyticD2ScalarFunction.hpp"
#include "Array.hpp"
#include "FastCopy.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for a linear mapping from R^n to R. The function is computed as the inner product of a point in R^n with a
 * coefficient vector.
 */
class THEA_API LinearScalarFunction : public AnalyticD2ScalarFunction
{
  public:
    THEA_DEF_POINTER_TYPES(LinearScalarFunction, shared_ptr, weak_ptr)

    /**
     * Get the coefficients of the function mapping R^n to R. The function is computed by valueAt() as the inner product of the
     * point <b>p</b> with the coefficient vector <b>coeffs</b>.
     */
    virtual TheaArray<double> const & getCoefficients() const = 0;

    /**
     * {@inheritDoc}
     *
     * The function is computed as the inner product of the point <b>p</b> with the coefficient vector <b>coeffs</b> returned by
     * getCoefficients().
     */
    double valueAt(double const * p) const
    {
      TheaArray<double> const & coeffs = getCoefficients();

      double val = 0;
      for (array_size_t i = 0; i < coeffs.size(); ++i)
        val += (p[i] * coeffs[i]);

      return val;
    }

    // Gradient of a linear function is the coefficient vector

    void gradientAt(float const * p, float * result) const
    {
      TheaArray<double> const & coeffs = getCoefficients();
      for (array_size_t i = 0; i < coeffs.size(); ++i)
        result[i] = static_cast<float>(coeffs[i]);
    }

    void gradientAt(double const * p, double * result) const
    {
      TheaArray<double> const & coeffs = getCoefficients();
      if (!coeffs.empty()) fastCopy(&coeffs[0], &coeffs[0] + coeffs.size(), result);
    }

    // Hessian of a linear function is identically zero everywhere

    bool hasSparseHessian() const { return true; }

    void hessianAt(float const * p, AddressableMatrix<float> & result) const
    {
      checkHessianResultDims(result);
      result.makeZero();
    }

    void hessianAt(double const * p, AddressableMatrix<double> & result) const
    {
      checkHessianResultDims(result);
      result.makeZero();
    }

    void hessianAt(float const * p, CompressedRowMatrix<float> & result) const
    {
      checkHessianResultDims(result);
      result.makeZero();
    }

    void hessianAt(double const * p, CompressedRowMatrix<double> & result) const
    {
      checkHessianResultDims(result);
      result.makeZero();
    }

    void hessianAt(float const * p, CompressedColumnMatrix<float> & result) const
    {
      checkHessianResultDims(result);
      result.makeZero();
    }

    void hessianAt(double const * p, CompressedColumnMatrix<double> & result) const
    {
      checkHessianResultDims(result);
      result.makeZero();
    }

  private:
    /** Check dimensions of matrix supplied to store Hessian. */
    template <typename MatrixT> void checkHessianResultDims(MatrixT const & m) const
    {
      alwaysAssertM(m.numRows() == numDimensions() && m.numColumns() != numDimensions(),
                    "LinearScalarFunction: Hessian result matrix has wrong dimensions");
    }

}; // class LinearScalarFunction

} // namespace Algorithms
} // namespace Thea

#endif
