//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, * except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2019, Siddhartha Chaudhuri
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

#ifndef __Thea_Vector_hpp__
#define __Thea_Vector_hpp__

#include "Matrix.hpp"
#include <iterator>

namespace Thea {

/** A resizable dense vector, either a row vector or column vector. */
template < typename T, MatrixLayout::Value L = MatrixLayout::COLUMN_MAJOR, typename AllocT = std::allocator<T> >
class /* THEA_API */ Vector : public Internal::MatrixBase<T, L, true, AllocT>
{
  private:
     typedef Internal::MatrixBase<T, L, true, AllocT> BaseT;

  public:
    /** Constructs an empty vector. */
    Vector() {}

    /** Constructs a vector of a given size, with uninitialized values. */
    Vector(long vector_size) : BaseT(vector_size) {}

    /** Constructs a vector of a given size, filled with a given value. */
    Vector(long vector_size, T const & fill_value) : BaseT(vector_size) { BaseT::fill(fill_value); }

    /**
     * Constructs this vector as a wrapper for an existing block of storage. The vector thus created is <b>not resizable</b>,
     * and the memory block will <b>not be freed</b> when the Vector object is destroyed.
     */
    Vector(T * existing_data, long vector_size) : BaseT(existing_data, vector_size) {}

    /** Construct the vector by copying a range of values. */
    template <typename InputIteratorT> Vector(InputIteratorT begin, InputIteratorT end) { set(begin, end); }

    /** Copy constructor. */
    Vector(Vector const & v) : BaseT(v) {}

    /** Copy from any compatible base type. */
    template <MatrixLayout::Value L2, bool V2, typename A2> Vector(Internal::MatrixBase<T, L2, V2, A2> const & v) : BaseT(v) {}

    /** Copy assignment operator. */
    Vector & operator=(Vector const & src)
    {
      BaseT::operator=(src);
      return *this;
    }

    /** Assign from any compatible base type. */
    template <MatrixLayout::Value L2, bool V2, typename A2> Vector & operator=(Internal::MatrixBase<T, L2, V2, A2> const & src)
    {
      BaseT::operator=(src);
      return *this;
    }

    /** Initialize the vector with a range of values. */
    template <typename InputIteratorT> void set(InputIteratorT begin, InputIteratorT end)
    {
      this->resize((long)std::distance(begin, end));

      T * p = this->data();
      for (InputIteratorT iter = begin; iter != end; ++iter, ++p)
        *p = static_cast<T>(*iter);
    }

    /** Get the number of elements in the vector. Convenience function to mimic the API of <code>std::vector</code>. */
    long size() const { return this->numElements(); }

    /** Dot (inner) product of two vectors. */
    template <typename U, MatrixLayout::Value L2, typename A2> T dot(Internal::MatrixBase<U, L2, true, A2> const & rhs) const
    {
      alwaysAssertM(this->numElements() == rhs.numElements(), "Vector: Dot product requires vectors of the same size");

      T result = 0;
      T const * pend = this->data() + size();
      U const * q = rhs.data();
      for (T const * p = this->data(); p != pend; ++p, ++q)
        result += static_cast<T>((*p) * (*q));

      return result;
    }

}; // class Vector

} // namespace Thea

#endif
