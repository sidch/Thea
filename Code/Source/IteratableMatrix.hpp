//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_IteratableMatrix_hpp__
#define __Thea_IteratableMatrix_hpp__

#include "BasicMatrix.hpp"

namespace Thea {

/**
 * Base class for a 2D matrix with an iterator that visits (at least) all the non-zero entries in the matrix. The iterator
 * is assumed to be read-only, and <b>cannot</b> in general be used to change the visited elements. Any class that inherits this
 * should have the following members.
 *
 * \code
 *  // Supporting types
 *  typedef std::pair<integer-type, integer-type> IndexPair;  ///< A (row, column) index pair.
 *  typedef std::pair<IndexPair, T> Entry;  ///< An entry in the matrix, mapping a (row, column) pair to a value.
 *
 *  // Read-only iterator
 *  typedef <type> ConstIterator;
 *
 *  // Access the beginning and end of the matrix
 *  ConstIterator begin() const;
 *  ConstIterator end() const;
 * \endcode
 *
 * <tt>ConstIterator</tt> should dereference to Entry.
 */
template <typename T>
class /* THEA_API */ IteratableMatrix : public virtual BasicMatrix<T>
{
  private:
    typedef BasicMatrix<T> BaseType;

  public:
    THEA_DEF_POINTER_TYPES(IteratableMatrix, shared_ptr, weak_ptr)

}; // class IteratableMatrix

} // namespace Thea

#endif
