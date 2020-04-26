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

#ifndef __Browse3D_MeshFwd_hpp__
#define __Browse3D_MeshFwd_hpp__

#include "Common.hpp"

namespace Thea {
namespace Graphics {

// Forward declarations
struct NullAttribute;

template <typename VA, typename EA, typename FA, template <typename T> class Alloc> class GeneralMeshVertex;
template <typename VA, typename EA, typename FA, template <typename T> class Alloc> class GeneralMeshFace;
template <typename MeshType> class MeshGroup;

} // namespace Graphics
} // namespace Thea

namespace Browse3D {

// Forward declarations
class VertexAttribute;
class FaceAttribute;
class Mesh;

/**< Handle to a mesh vertex. */
typedef Graphics::GeneralMeshVertex<VertexAttribute, Graphics::NullAttribute, FaceAttribute, std::allocator> MeshVertex;

/**< Handle to a mesh face. */
typedef Graphics::GeneralMeshFace<VertexAttribute, Graphics::NullAttribute, FaceAttribute, std::allocator> MeshFace;

typedef Graphics::MeshGroup<Mesh> MeshGroup;            ///< A hierarchical group of meshes.
typedef std::shared_ptr<Mesh> MeshPtr;                       ///< A shared pointer to a mesh.
typedef std::shared_ptr<Mesh const> MeshConstPtr;            ///< A shared pointer to an immutable mesh.
typedef std::shared_ptr<MeshGroup> MeshGroupPtr;             ///< A shared pointer to a group of meshes.
typedef std::shared_ptr<MeshGroup const> MeshGroupConstPtr;  ///< A shared pointer to an immutable group of meshes.

} // namespace Browse3D

#endif
