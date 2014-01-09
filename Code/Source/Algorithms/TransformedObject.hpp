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

#ifndef __Thea_Algorithms_TransformedObject_hpp__
#define __Thea_Algorithms_TransformedObject_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Wraps pointers to an object and an associated transformation. The pointer values will be dereferenced as necessary by the
 * <code>get...()</code> functions. The original data must persist until this class is destroyed. One should always check
 * the return value of the <code>has...()</code> function before calling the corresponding <code>get...()</code> functions.
 */
template <typename ObjectT, typename TransformT>
class /* THEA_API */ TransformedObject
{
  public:
    typedef ObjectT Object;  ///< The type of wrapped object.
    typedef TransformT Transform;  ///< The type of wrapped transformation.

    /**
     * Constructor. Both pointers are stored directly, so the original objects must persist as long as the TransformedObject is
     * being used. This is the preferred constructor for efficiency.
     */
    explicit TransformedObject(Object const * obj_ = NULL, Transform const * trans_ = NULL) : obj(obj_), trans(trans_) {}

    /** Check if this class wraps an existing object or not. */
    bool hasObject() const { return obj != NULL; }

    /** Get the wrapped object. */
    Object const & getObject() const { return *obj; }

    /** Check if this class wraps an existing transformation or not. */
    bool hasTransform() const { return trans != NULL; }

    /** Get the wrapped transformation. */
    Transform const & getTransform() const { return *trans; }

  private:
    Object const * obj;
    Transform const * trans;

}; // class TransformedObject

/** Utility function to wrap an object-transformation pair. Allows the compiler to automatically infer template parameters. */
template <typename ObjectT, typename TransformT>
TransformedObject<ObjectT, TransformT>
makeTransformedObject(ObjectT const * obj, TransformT const * trans)
{
  return TransformedObject<ObjectT, TransformT>(obj, trans);
}

} // namespace Algorithms
} // namespace Thea

#endif
