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

#ifndef __Thea_TransposedMatrix_hpp__
#define __Thea_TransposedMatrix_hpp__

#include "AddressableMatrix.hpp"

namespace Thea {

/**
 * A matrix class that acts as the transpose of another matrix, without actually copying data. An object of this class holds a
 * pointer to another addressable matrix. Element (r, c) of this matrix corresponds to element (c, r) of the wrapped matrix. The
 * wrapped matrix must persist until this matrix is destroyed.
 */
template <typename T>
class /* THEA_API */ TransposedMatrix : public AddressableMatrix<T>
{
  public:
    THEA_DEF_POINTER_TYPES(TransposedMatrix, shared_ptr, weak_ptr)

    /** Constructor. Stores the passed pointer to a matrix, which should remain valid until this object is destroyed. */
    TransposedMatrix(AddressableMatrix<T> * m_) : m(m_)
    {
      alwaysAssertM(m_, "TransposedMatrix: The wrapped matrix must exist");
    }

    long numRows() const { return m->numColumns(); }
    long numColumns() const { return m->numRows(); }

    T const & get(long row, long col) const { return m->get(col, row); }
    T & getMutable(long row, long col) { return m->getMutable(col, row); }
    void set(long row, long col, T const & value) { m->set(col, row, value); }

  private:
    AddressableMatrix<T> * m;

}; // class TransposedMatrix

} // namespace Thea

#endif
