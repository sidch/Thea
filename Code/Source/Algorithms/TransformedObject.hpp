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
    explicit TransformedObject(Object const * obj_ = nullptr, Transform const * trans_ = nullptr) : obj(obj_), trans(trans_) {}

    /** Check if this class wraps an existing object or not. */
    bool hasObject() const { return obj != nullptr; }

    /** Get the wrapped object. */
    Object const & getObject() const { return *obj; }

    /** Check if this class wraps an existing transformation or not. */
    bool hasTransform() const { return trans != nullptr; }

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
