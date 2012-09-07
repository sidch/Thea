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

#ifndef __Thea_Transformable_hpp__
#define __Thea_Transformable_hpp__

#include "Common.hpp"

namespace Thea {

/** An object with an associated transformation, stored separately. */
template <typename TransformT>
class /* THEA_API */ Transformable
{
  public:
    THEA_DEF_POINTER_TYPES(Transformable, shared_ptr, weak_ptr)

    typedef TransformT Transform;  ///< The type of transformation.

    /** Destructor. */
    virtual ~Transformable() = 0;

    /** Default constructor. */
    Transformable() : has_trans(false) {}

    /** Constructor. */
    Transformable(Transform const & trans_) : has_trans(true), trans(trans_) {}

    /** Copy constructor. */
    Transformable(Transformable const & src) : has_trans(src.has_trans), trans(src.trans) {}

    /** Check if a transform has been set. */
    bool hasTransform() const { return has_trans; }

    /** Get the transformation. */
    Transform const & getTransform() const { return trans; }

    /** Get the transformation. */
    Transform & getTransform() { return trans; }

    /** Set the transformation. */
    virtual void setTransform(Transform const & trans_) { has_trans = true; trans = trans_; }

    /** Clear any existing transform. */
    virtual void clearTransform() { has_trans = false; }

  private:
    bool has_trans;
    TransformT trans;

}; // Transformable

template <typename TransformT> inline Transformable<TransformT>::~Transformable() {}

} // namespace Thea

#endif
