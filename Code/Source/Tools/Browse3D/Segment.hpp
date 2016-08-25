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

#ifndef __Browse3D_Segment_hpp__
#define __Browse3D_Segment_hpp__

#include "Common.hpp"
#include "MeshFwd.hpp"
#include "../../UnorderedSet.hpp"

namespace Browse3D {

/** A labeled segment. */
class Segment
{
  public:
    typedef TheaUnorderedSet<Mesh *> MeshSet;  ///< A set of meshes

    /** Default constructor. */
    Segment() : label("AnonymousSegment") {}

    /** Create a segment with the given label. */
    Segment(std::string const & label_) : label(label_) {}

    /** Get the segment label. */
    std::string const & getLabel() const { return label; }

    /** Set the segment label. */
    void setLabel(std::string const & label_) { label = label_; }

    /** Get the number of meshes in the segment. */
    long numMeshes() const { return (long)meshes.size(); }

    /** Get the set of meshes. */
    MeshSet const & getMeshes() const { return meshes; }

    /** Add a mesh to the segment. */
    void addMesh(Mesh * mesh) { meshes.insert(mesh); }

    /** Remove a mesh from the segment. */
    void removeMesh(Mesh const * mesh, long depth_promotion = 0);

    /** Add a mesh group to the segment. */
    void addMeshGroup(MeshGroup * mg);

    /** Remove a mesh group from the segment. */
    void removeMeshGroup(MeshGroup const * mg);

    /** Check if the (possibly hierarchically expanded) segment contains a given mesh. */
    bool hasMesh(Mesh const * mesh, long depth_promotion = 0) const;

    /** Get the minimum depth (from the root) of a mesh in the segment. */
    long minDepth() const;

    /** Clear the segment. */
    void clear() { meshes.clear(); }

  private:
    MeshSet meshes;
    std::string label;

}; // class Segment

} // namespace Browse3D

#endif
