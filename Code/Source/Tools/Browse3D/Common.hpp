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
// First version: 2011
//
//============================================================================

#ifndef __Browse3D_Common_hpp__
#define __Browse3D_Common_hpp__

#include "../../Common.hpp"
#include "../../Colors.hpp"
#include <wx/event.h>

// Forward declarations
class wxPoint;
class wxRealPoint;
class wxKeyEvent;
class wxMouseEvent;
class wxSizeEvent;

/** Namespace for 3D file browser project. */
namespace Browse3D {

using namespace Thea;

// Put quotes around the result of a macro expansion.
#define BROWSE3D_STRINGIFY_(x) #x
#define BROWSE3D_STRINGIFY(x) BROWSE3D_STRINGIFY_(x)

/** A dummy event that is a default argument for callbacks. */
extern wxCommandEvent DUMMY_EVENT;

} // namespace Browse3D

#endif
