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

#ifndef __Thea_Graphics_IDrawable_hpp__
#define __Thea_Graphics_IDrawable_hpp__

#include "../Common.hpp"
#include "RenderOptions.hpp"
#include "IRenderSystem.hpp"

namespace Thea {

/** %Graphics functionality. */
namespace Graphics {

/** Interface for objects that can be displayed onscreen. */
class THEA_API IDrawable
{
  public:
    THEA_DECL_SMART_POINTERS(IDrawable)

    /** Destructor. */
    virtual ~IDrawable() = 0;

    /** Draw the object using the specified options via the specified rendering system. */
    virtual void THEA_ICALL draw(IRenderSystem * render_system, IRenderOptions const * options = nullptr) const = 0;

}; // class IDrawable

// Pure virtual destructor should have a body
// http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
inline IDrawable::~IDrawable() {}

} // namespace Graphics
} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::Graphics::IDrawable)

#endif
