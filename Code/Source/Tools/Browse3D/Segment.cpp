//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2015, Siddhartha Chaudhuri
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

#include "Segment.hpp"
#include "Mesh.hpp"

namespace Browse3D {


void
Segment::removeMesh(Mesh const * mesh, long depth_promotion)
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
  SegmentInternal::MeshGroupAdder adder(*this);
  if (mg) mg->forEachMeshUntil(&adder);
}

void
Segment::removeMeshGroup(MeshGroup const * mg)
{
  SegmentInternal::MeshGroupRemover remover(*this);
  if (mg) mg->forEachMeshUntil(&remover);
}

bool
Segment::hasMesh(Mesh const * mesh, long depth_promotion) const
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

long
Segment::minDepth() const
{
  long min_depth = 0;
  for (MeshSet::const_iterator mi = meshes.begin(); mi != meshes.end(); ++mi)
  {
    MeshGroup const * p = (*mi)->getParent();
    if (!p)
      continue;

    long d = p->getDepth() + 1;  // mesh is one level below parent
    if (min_depth <= 0 || d < min_depth)
      min_depth = d;
  }

  return min_depth;
}

} // namespace Browse3D
