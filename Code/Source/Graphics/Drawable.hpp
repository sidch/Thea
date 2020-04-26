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

#ifndef __Thea_Graphics_Drawable_hpp__
#define __Thea_Graphics_Drawable_hpp__

#include "../Common.hpp"
#include "RenderOptions.hpp"
#include "RenderSystem.hpp"

namespace Thea {

/** %Graphics functionality. */
namespace Graphics {

/** Abstract base class for objects that can be displayed onscreen. */
class THEA_API Drawable
{
  public:
    THEA_DECL_SMART_POINTERS(Drawable)

    /** Destructor. */
    virtual ~Drawable() = 0;

    /** Draw the object using the specified options via the specified rendering system. */
    virtual void draw(RenderSystem & render_system,
                      AbstractRenderOptions const & options = RenderOptions::defaults()) const = 0;

}; // class Drawable

// Pure virtual destructor should have a body
// http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
inline Drawable::~Drawable() {}

} // namespace Graphics
} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Graphics::Drawable)

#endif
