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

#ifndef __Thea_AttributedObject_hpp__
#define __Thea_AttributedObject_hpp__

#include "Common.hpp"

namespace Thea {

/** An object with an attached attribute (which must be default-constructible). */
template <typename AttributeType>
class /* THEA_API */ AttributedObject
{
  public:
    typedef AttributeType Attribute;

    /** Default constructor. */
    AttributedObject() {}

    /** Constructor. */
    AttributedObject(Attribute const & attrib_) : attrib(attrib_) {}

    /** Get the attribute of the object. */
    Attribute const & attr() const { return attrib; }

    /** Get the attribute of the object. */
    Attribute & attr() { return attrib; }

    /** Set the attribute of the object. */
    void setAttr(Attribute const & attrib_) { attrib = attrib_; }

  private:
    Attribute attrib;

}; // class AttributedObject

/** A completely empty attribute. */
struct THEA_API NullAttribute {};

} // namespace Thea

#endif
