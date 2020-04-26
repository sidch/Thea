//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2009
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
    THEA_DECL_SMART_POINTERS(Transformable)

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

    /** Get the transformation, if one has been set. Else, the return value is undefined. */
    Transform const & getTransform() const { return trans; }

    /** Get the transformation, if one has been set. Else, the return value is undefined. */
    Transform & getTransform() { return trans; }

    /** Set the transformation. */
    virtual void setTransform(Transform const & trans_) { has_trans = true; trans = trans_; }

    /** Clear any existing transform. */
    virtual void clearTransform() { has_trans = false; }

  private:
    bool has_trans;
    TransformT trans;

}; // Transformable

template <typename TransformT>
Transformable<TransformT>::~Transformable()
{
  // Pure virtual destructor should have a body
  // http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
}

} // namespace Thea

#endif
