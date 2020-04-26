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
// First version: 2015
//
//============================================================================

#include "Segment.hpp"
#include "Mesh.hpp"

namespace Browse3D {


void
Segment::removeMesh(Mesh const * mesh, intx depth_promotion)
{
  if (!mesh)
    return;

  // Always remove this mesh
  meshes.erase(const_cast<Mesh *>(mesh));

  if (depth_promotion <= 0)
    return;

  // Remove every mesh in the segment that has a common ancestor with this mesh
  for (MeshSet::iterator mi = meshes.begin(); mi != meshes.end(); )
  {
    MeshGroup const * anc = (*mi)->getAncestor(depth_promotion);
    if (anc && mesh->hasAncestor(anc))
    {
      MeshSet::iterator to_remove = mi;
      ++mi;
      meshes.erase(to_remove);
    }
    else
      ++mi;
  }
}

namespace SegmentInternal {

struct MeshGroupAdder
{
  MeshGroupAdder(Segment & seg_) : seg(seg_) {}
  bool operator()(Mesh & mesh) { seg.addMesh(&mesh); return false; }

  Segment & seg;
};

struct MeshGroupRemover
{
  MeshGroupRemover(Segment & seg_) : seg(seg_) {}
  bool operator()(Mesh const & mesh) { seg.removeMesh(&mesh); return false; }

  Segment & seg;
};

} // namespace SegmentInternal

void
Segment::addMeshGroup(MeshGroup * mg)
{
  if (mg) mg->forEachMeshUntil(SegmentInternal::MeshGroupAdder(*this));
}

void
Segment::removeMeshGroup(MeshGroup const * mg)
{
  if (mg) mg->forEachMeshUntil(SegmentInternal::MeshGroupRemover(*this));
}

bool
Segment::hasMesh(Mesh const * mesh, intx depth_promotion) const
{
  if (!mesh)
    return false;

  if (depth_promotion <= 0)
    return meshes.find(const_cast<Mesh *>(mesh)) != meshes.end();

  // Check if any mesh in the segment has a common ancestor with the query mesh
  for (MeshSet::const_iterator mi = meshes.begin(); mi != meshes.end(); ++mi)
  {
    MeshGroup * anc = (*mi)->getAncestor(depth_promotion);
    if (anc && mesh->hasAncestor(anc))
      return true;
  }

  return false;
}

intx
Segment::minDepth() const
{
  intx min_depth = 0;
  for (MeshSet::const_iterator mi = meshes.begin(); mi != meshes.end(); ++mi)
  {
    MeshGroup const * p = (*mi)->getParent();
    if (!p)
      continue;

    intx d = p->getDepth() + 1;  // mesh is one level below parent
    if (min_depth <= 0 || d < min_depth)
      min_depth = d;
  }

  return min_depth;
}

} // namespace Browse3D
