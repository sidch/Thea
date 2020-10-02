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

#ifndef __Thea_Graphics_MeshType_hpp__
#define __Thea_Graphics_MeshType_hpp__

#include "../Common.hpp"
#include "../Concept.hpp"

namespace Thea {
namespace Graphics {

/** Concept of a general-purpose mesh. */
THEA_HAS_TYPE(IsGeneralMesh, GENERAL_MESH_TAG)

/** Concept of a DCEL mesh. */
THEA_HAS_TYPE(IsDcelMesh, DCEL_MESH_TAG)

/** Concept of a display mesh. */
THEA_HAS_TYPE(IsDisplayMesh, DISPLAY_MESH_TAG)

} // namespace Graphics
} // namespace Thea

#endif
