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
    typedef UnorderedSet<Mesh *> MeshSet;  ///< A set of meshes

    /** Default constructor. */
    Segment() : label("AnonymousSegment") {}

    /** Create a segment with the given label. */
    Segment(std::string const & label_) : label(label_) {}

    /** Get the segment label. */
    std::string const & getLabel() const { return label; }

    /** Set the segment label. */
    void setLabel(std::string const & label_) { label = label_; }

    /** Get the number of meshes in the segment. */
    intx numMeshes() const { return (intx)meshes.size(); }

    /** Get the set of meshes. */
    MeshSet const & getMeshes() const { return meshes; }

    /** Add a mesh to the segment. */
    void addMesh(Mesh * mesh) { meshes.insert(mesh); }

    /** Remove a mesh from the segment. */
    void removeMesh(Mesh const * mesh, intx depth_promotion = 0);

    /** Add a mesh group to the segment. */
    void addMeshGroup(MeshGroup * mg);

    /** Remove a mesh group from the segment. */
    void removeMeshGroup(MeshGroup const * mg);

    /** Check if the (possibly hierarchically expanded) segment contains a given mesh. */
    bool hasMesh(Mesh const * mesh, intx depth_promotion = 0) const;

    /** Get the minimum depth (from the root) of a mesh in the segment. */
    intx minDepth() const;

    /** Clear the segment. */
    void clear() { meshes.clear(); }

  private:
    MeshSet meshes;
    std::string label;

}; // class Segment

} // namespace Browse3D

#endif
