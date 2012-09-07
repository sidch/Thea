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

#ifndef __Thea_Algorithms_LinearConstraint_hpp__
#define __Thea_Algorithms_LinearConstraint_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../MatrixFormat.hpp"
#include "../MatrixUtil.hpp"
#include "../MatrixWrapper.hpp"

namespace Thea {
namespace Algorithms {

/** A linear optimization constraint. */
class LinearConstraint
{
  public:
    THEA_DEF_POINTER_TYPES(LinearConstraint, shared_ptr, weak_ptr)

    /**
     * Constructor. Creates a linear constraint of the form <b>Ax</b> <em>op</em> <b>b</b>, where <em>op</em> is a comparison
     * operator.
     *
     * @param coeffs_ Coefficient matrix <b>A</b>.
     * @param compare_ Comparison operator <em>op</em>.
     * @param constants_begin Points to the beginning of constant vector <b>b</b>.
     * @param constants_end Points to the end of constant vector <b>b</b>.
     */
    template <typename MatrixT, typename RealInputIterator>
    LinearConstraint(MatrixT const & coeffs_, CompareOp compare_, RealInputIterator constants_begin,
                     RealInputIterator constants_end)
    : coeffs(coeffs_, MatrixUtil::getFormat(coeffs_)), compare(compare_), constants(constants_begin, constants_end)
    {
      alwaysAssertM(coeffs.numColumns() == (int)constants.size(), "Coefficient and constant dimensions don't match");
    }

    /**
     * Constructor. Creates a linear constraint of the form <b>Ax</b> <em>op</em> <b>b</b>, where <em>op</em> is a comparison
     * operator.
     *
     * @param coeffs_ Coefficient matrix <b>A</b>.
     * @param dst_format The desired storage format for the coefficient matrix.
     * @param compare_ Comparison operator <em>op</em>.
     * @param constants_begin Points to the beginning of constant vector <b>b</b>.
     * @param constants_end Points to the end of constant vector <b>b</b>.
     */
    template <typename MatrixT, typename RealInputIterator>
    LinearConstraint(MatrixT const & coeffs_, MatrixFormat dst_format, CompareOp compare_, RealInputIterator constants_begin,
                     RealInputIterator constants_end)
    : coeffs(coeffs_, dst_format), compare(compare_), constants(constants_begin, constants_end)
    {
      alwaysAssertM(coeffs.numColumns() == (int)constants.size(), "Coefficient and constant dimensions don't match");
    }

    /**
     * Constructor. Creates a linear constraint of the form <b>Ax</b> <em>op</em> <b>b</b>, where <em>op</em> is a comparison
     * operator.
     *
     * @param coeffs_ Coefficient matrix <b>A</b>.
     * @param compare_ Comparison operator <em>op</em>.
     * @param constants_begin Points to the beginning of constant vector <b>b</b>.
     * @param constants_end Points to the end of constant vector <b>b</b>.
     */
    template <typename T, typename I2, typename I1, typename RealInputIterator>
    LinearConstraint(MatrixWrapper<T, I2, I1> const & coeffs_, CompareOp compare_, RealInputIterator constants_begin,
                     RealInputIterator constants_end)
    : coeffs(coeffs_), compare(compare_), constants(constants_begin, constants_end)
    {
      alwaysAssertM(coeffs.numColumns() == (int)constants.size(), "Coefficient and constant dimensions don't match");
    }

    /** Get the number of dimensions of the linear constraint. */
    int numDimensions() const { return (int)constants.size(); }

    /** Get the coefficient matrix. */
    MatrixWrapper<double> const & getCoefficients() const { return coeffs; }

    /** Get the comparison operator. */
    CompareOp getCompareOp() const { return compare; }

    /** Get the constant vector. */
    TheaArray<double> const & getConstants() const { return constants; }

  private:
    MatrixWrapper<double> coeffs;  ///< Coefficient matrix.
    CompareOp compare;             ///< Comparison operator.
    TheaArray<double> constants;   ///< Constant vector.

}; // class LinearConstraint

} // namespace Algorithms
} // namespace Thea

#endif
