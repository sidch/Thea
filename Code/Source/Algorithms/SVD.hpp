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

#ifndef __Thea_Algorithms_SVD_hpp__
#define __Thea_Algorithms_SVD_hpp__

#include "../Common.hpp"
#include "../AddressableMatrix.hpp"
#include "../Array.hpp"
#include "../Matrix.hpp"
#include "../ResizableMatrix.hpp"
#include "../TransposedMatrix.hpp"
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>
#include <algorithm>
#include <cmath>

namespace Thea {
namespace Algorithms {

/**
 * Singular value decomposition of matrices. A matrix A is decomposed into U * D * V, where U and V are orthogonal matrices and
 * D is a diagonal matrix.
 */
class THEA_API SVD
{
  public:
    /**
     * Compute the singular value decomposition of a matrix. An m x n matrix A is decomposed into U * D * V.
     * - If m >= n, U is an m x n matrix with orthogonal columns, D is an n x n diagonal matrix, and V is an n x n orthogonal
     *   matrix.
     * - If m < n, U is an m x m matrix with orthogonal columns, D is an m x m diagonal matrix, and V is an m x n matrix with
     *   orthogonal rows.
     *
     * @param a The input matrix A.
     * @param u Used to store the matrix U.
     * @param d Used to store the diagonal of the matrix D.
     * @param v Used to store the matrix V.
     *
     * @return True if the SVD was successfully computed, else false.
     *
     * @note The matrices \a a, \a u, \a d and \a v should all point to strictly different objects!
     * @note Full/partial specializations of the helper class SVDImpl define the decomposition for particular input/output
     *   types.
     *
     * @see SVDImpl
     */
    template <typename InputMatrixT, typename MatrixUT, typename T, typename MatrixVT>
    static bool compute(InputMatrixT const & a, MatrixUT & u, TheaArray<T> & d, MatrixVT & v);

    /**
     * Compute the Moore-Penrose pseudo-inverse of a (possibly rectangular) matrix A. If (A.transpose() * A) exists, the
     * pseudo-inverse is defined as (A.transpose() * A).inverse() * A.transpose(). However, even if the inverse does not exist,
     * the pseudo-inverse exists and is typically computed via a singular value decomposition.
     *
     * @param a The matrix whose pseudo-inverse is required.
     * @param result The computed pseudo-inverse is placed here.
     * @param tolerance Singular values less than this are discarded.
     *
     * @note The matrices \a a and \a result should point to strictly different objects!
     * @note Full/partial specializations of the helper class SVDPseudoInverseImpl define the decomposition for particular
     *   input/output types.
     *
     * @see SVDPseudoInverseImpl
     */
    template <typename InputMatrixT, typename OutputMatrixT>
    static bool pseudoInverse(InputMatrixT const & a, OutputMatrixT & result, double tolerance = -1);

}; // class SVD

namespace SVDInternal {

THEA_API bool svdCore(AddressableMatrix<float>  & U, long rows, long cols, float  * D, AddressableMatrix<float>  & V);
THEA_API bool svdCore(AddressableMatrix<float>  & U, long rows, long cols, float  * D, AddressableMatrix<double> & V);
THEA_API bool svdCore(AddressableMatrix<float>  & U, long rows, long cols, double * D, AddressableMatrix<float>  & V);
THEA_API bool svdCore(AddressableMatrix<float>  & U, long rows, long cols, double * D, AddressableMatrix<double> & V);
THEA_API bool svdCore(AddressableMatrix<double> & U, long rows, long cols, float  * D, AddressableMatrix<float>  & V);
THEA_API bool svdCore(AddressableMatrix<double> & U, long rows, long cols, float  * D, AddressableMatrix<double> & V);
THEA_API bool svdCore(AddressableMatrix<double> & U, long rows, long cols, double * D, AddressableMatrix<float>  & V);
THEA_API bool svdCore(AddressableMatrix<double> & U, long rows, long cols, double * D, AddressableMatrix<double> & V);

template <typename MU>
typename boost::disable_if< boost::is_base_of< ResizableMatrix<typename MU::Value>, MU > >::type
checkU(long m, long n, MU & u)
{
  if ((m >= n && !(u.numRows() == m && u.numColumns() == n))
   || (m <  n && !(u.numRows() == m && u.numColumns() == m)))
  {
    throw Error("SVD: Supplied matrix is of incorrect dimensions and cannot be resized");
  }
}

template <typename MU>
typename boost::enable_if< boost::is_base_of< ResizableMatrix<typename MU::Value>, MU > >::type
checkU(long m, long n, MU & u)
{
  if (m >= n && !(u.numRows() == m && u.numColumns() == n))
    u.resize(m, n);
  else if (m < n && !(u.numRows() == m && u.numColumns() == m))
    u.resize(m, m);
}

template <typename MV>
typename boost::disable_if< boost::is_base_of< ResizableMatrix<typename MV::Value>, MV > >::type
checkV(long m, long n, MV & v)
{
  if ((m >= n && !(v.numRows() == n && v.numColumns() == n))
   || (m <  n && !(v.numRows() == m && v.numColumns() == n)))
  {
    throw Error("SVD: Supplied matrix is of incorrect dimensions and cannot be resized");
  }
}

template <typename MV>
typename boost::enable_if< boost::is_base_of< ResizableMatrix<typename MV::Value>, MV > >::type
checkV(long m, long n, MV & v)
{
  if (m >= n && !(v.numRows() == n && v.numColumns() == n))
    v.resize(n, n);
  else if (m < n && !(v.numRows() == m && v.numColumns() == n))
    v.resize(m, n);
}

template <typename M>
typename boost::disable_if< boost::is_base_of< ResizableMatrix<typename M::Value>, M > >::type
checkDims(long m, long n, M & a)
{
  if (a.numRows() != m || a.numColumns() != n)
    throw Error("SVD: Supplied matrix is of incorrect dimensions and cannot be resized");
}

template <typename M>
typename boost::enable_if< boost::is_base_of< ResizableMatrix<typename M::Value>, M > >::type
checkDims(long m, long n, M & a)
{
  if (a.numRows() != m || a.numColumns() != n)
    a.resize(m, n);
}

} // namespace SVDInternal

/**
 * Helper class for SVD. Specializations of this class actually compute the singular value decomposition. This is required
 * because C++ does not allow partial specialization of function templates, and does unexpected things with specialized and
 * overloaded function templates (see http://www.gotw.ca/publications/mill17.htm). To extend SVD to other input/output types,
 * specialize this class (not SVD::compute()).
 *
 * @see SVD
 */
template <typename InputMatrixT, typename MatrixUT, typename T, typename MatrixVT, typename Enable = void>
class /* THEA_API */ SVDImpl
{
  public:
    /**
     * Compute the singular value decomposition of a matrix. An m x n matrix A is decomposed into U * D * V.
     * - If m >= n, U is an m x n matrix with orthogonal columns, D is an n x n diagonal matrix, and V is an n x n orthogonal
     *   matrix.
     * - If m < n, U is an m x m matrix with orthogonal columns, D is an m x m diagonal matrix, and V is an m x n matrix with
     *   orthogonal rows.
     *
     * @param a The input matrix A.
     * @param u Used to store the matrix U.
     * @param d Used to store the diagonal of the matrix D.
     * @param v Used to store the matrix V.
     *
     * @return True if the SVD was successfully computed, else false.
     *
     * @note The matrices \a a, \a u, \a d and \a v should all point to strictly different objects!
     */
    static bool compute(InputMatrixT const & a, MatrixUT & u, TheaArray<T> & d, MatrixVT & v)
    {
      // Check dimensions
      long m = a.numRows(), n = a.numColumns();
      SVDInternal::checkU(m, n, u);
      SVDInternal::checkV(m, n, v);

      bool status = false;
      if (m >= n)
      {
        a.copyTo(u);  // for in-place SVD
        d.resize((array_size_t)n);
        status = SVDInternal::svdCore(u, m, n, &d[0], v);
      }
      else
      {
        // Use SVD(A) = {U, D, V}, SVD(A^T) = {V^T, D^T, U^T}
        a.copyTo(v);
        d.resize((array_size_t)m);

        // Quick transpose: just wrap the matrix with a TransposedMatrix wrapper
        TransposedMatrix<typename MatrixVT::Value> vt(&v);
        TransposedMatrix<typename MatrixUT::Value> ut(&u);
        status = SVDInternal::svdCore(vt, n, m, &d[0], ut);

        // U, V wrap the same data as U^T, V^T, so now have the correct values (U^T)^T and (V^T)^T
      }

      return status;
    }

}; // class SVDImpl

/**
 * Helper class for SVD::pseudoInverse(). Specializations of this class actually compute the pseudo-inverse. This is required
 * because C++ does not allow partial specialization of function templates, and does unexpected things with specialized and
 * overloaded function templates (see http://www.gotw.ca/publications/mill17.htm). To extend SVD to other input/output types,
 * specialize this class (not SVD::pseudoInverse()).
 *
 * @see SVD
 */
template <typename InputMatrixT, typename OutputMatrixT, typename Enable = void>
class /* THEA_API */ SVDPseudoInverseImpl
{
  public:
    /**
     * Compute the Moore-Penrose pseudo-inverse of a (possibly rectangular) matrix A. If (A.transpose() * A) exists, the
     * pseudo-inverse is defined as (A.transpose() * A).inverse() * A.transpose(). However, even if the inverse does not exist,
     * the pseudo-inverse exists and is typically computed via a singular value decomposition.
     *
     * @param a The matrix whose pseudo-inverse is required.
     * @param result Singular values less than this are discarded.
     * @param tolerance Singular values less than this are discarded.
     */
    static bool pseudoInverse(InputMatrixT const & a, OutputMatrixT & result, double tolerance = -1);

}; // class SVDPseudoInverseImpl

// Specializations

template <typename InputMatrixT, typename OutputMatrixT>
class /* THEA_API */ SVDPseudoInverseImpl< InputMatrixT, OutputMatrixT,
                                           typename boost::enable_if<
                                                        boost::is_base_of< AddressableMatrix<typename OutputMatrixT::Value>,
                                                                           OutputMatrixT > >::type >
{
  public:
    static bool pseudoInverse(InputMatrixT const & a, OutputMatrixT & result, double tolerance = -1)
    {
      typedef typename OutputMatrixT::Value Value;
      Matrix<Value> u, v;
      TheaArray<Value> d;
      if (!SVD::compute(a, u, d, v))
        return false;

      if (tolerance < 0)
      {
        // Comment from G3D: TODO: Should be eps(d[0]), which is the largest diagonal
        tolerance = std::max(a.numRows(), a.numColumns()) * 1.0e-4;
      }

      // Invert D
      for (array_size_t i = 0; i < d.size(); ++i)
      {
        if (std::abs(d[i]) < tolerance)
          d[i] = static_cast<Value>(0);
        else
          d[i] = static_cast<Value>(1) / d[i];
      }

      SVDInternal::checkDims(a.numColumns(), a.numRows(), result);

      if (a.numRows() >= a.numColumns())
      {
        long vrows = v.numRows(), vcols = v.numColumns();
        long rrows = result.numRows(), rcols = result.numColumns();

        // Compute (V^T * D^{-1})^T = D^{-1} * V
        for (long r = 0; r < vrows; ++r)
          for (long c = 0; c < vcols; ++c)
            v(r, c) *= d[(array_size_t)r];

        // Compute (D^{-1} * V)^T * U^T = V^T * D^{-1} * U^T, without explicitly computing either transpose
        Value val;
        for (long i = 0; i < rrows; ++i)
          for (long k = 0; k < rcols; ++k)
          {
            val = v(0, i) * u(k, 0);  // invert usual order of indices

            for (long j = 1; j < vrows; ++j)
              val += v(j, i) * u(k, j);  // invert usual order of indices

            result.set(i, k, val);
          }
      }
      else
      {
        long urows = u.numRows(), ucols = u.numColumns();
        long rrows = result.numRows(), rcols = result.numColumns();

        // Compute (D^{-1} * U^T)^T = U * D^{-1}
        for (long c = 0; c < ucols; ++c)
          for (long r = 0; r < urows; ++r)
            u(r, c) *= d[(array_size_t)c];

        // Compute V^T * (U * D^{-1})^T = V^T * D^{-1} * U^T, without explicitly computing either transpose
        Value val;
        for (long i = 0; i < rrows; ++i)
          for (long k = 0; k < rcols; ++k)
          {
            val = v(0, i) * u(k, 0);  // invert usual order of indices

            for (long j = 1; j < ucols; ++j)
              val += v(j, i) * u(k, j);  // invert usual order of indices

            result.set(i, k, val);
          }
      }

      return true;
    }

}; // class SVDPseudoInverseImpl< InputMatrixT, AddressableMatrix >

template <typename InputMatrixT, typename MatrixUT, typename T, typename MatrixVT>
bool
SVD::compute(InputMatrixT const & a, MatrixUT & u, TheaArray<T> & d, MatrixVT & v)
{
  return SVDImpl<InputMatrixT, MatrixUT, T, MatrixVT>::compute(a, u, d, v);
}

template <typename InputMatrixT, typename OutputMatrixT>
bool
SVD::pseudoInverse(InputMatrixT const & a, OutputMatrixT & result, double tolerance)
{
  return SVDPseudoInverseImpl<InputMatrixT, OutputMatrixT>::pseudoInverse(a, result, tolerance);
}

} // namespace Algorithms
} // namespace Thea

#endif
