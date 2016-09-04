//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Graphics_DCELVertex_hpp__
#define __Thea_Graphics_DCELVertex_hpp__

#include "../Common.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../AttributedObject.hpp"
#include "DCELHalfedge.hpp"
#include "GraphicsAttributes.hpp"

namespace Thea {

namespace Graphics {

/**
 * Vertex of DCELMesh.
 *
 * Adapted from: DCELVertex class. Part of an example DCEL implementation.
 * - Webpage: http://www.holmes3d.net/graphics/dcel/
 * - Author: Ryan Holmes
 * - E-mail: ryan [at] holmes3d [dot] net
 * - Usage: Use freely. Please cite the website as the source if you use it substantially unchanged. Please leave this
 *   documentation in the code.
 */
template <typename VertexAttribute, typename HalfedgeAttribute, typename FaceAttribute>
class /* THEA_API */ DCELVertex
: public PositionAttribute<Vector3>,
  public NormalAttribute<Vector3>,
  public AttributedObject<VertexAttribute>
{
  public:
    typedef DCELHalfedge<VertexAttribute, HalfedgeAttribute, FaceAttribute> Halfedge;  ///< Halfedge of the mesh.

  private:
    typedef PositionAttribute<Vector3>  PositionBaseType;
    typedef NormalAttribute<Vector3>    NormalBaseType;

  public:
    /** Default constructor. */
    DCELVertex()
    : NormalBaseType(Vector3::zero()), leaving(NULL), has_precomputed_normal(false), normal_normalization_factor(0) {}

    /** Sets the vertex to have a location. */
    explicit DCELVertex(Vector3 const & p)
    : PositionBaseType(p), NormalBaseType(Vector3::zero()), leaving(NULL), has_precomputed_normal(false),
      normal_normalization_factor(0)
    {}

    /**
     * Sets the vertex to have a location and a precomputed normal. Adding faces will <b>not</b> change the vertex normal if
     * this constructor is used.
     */
    DCELVertex(Vector3 const & p, Vector3 const & n)
    : PositionBaseType(p), NormalBaseType(n), leaving(NULL), has_precomputed_normal(true), normal_normalization_factor(0)
    {}

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Halfedge const * getEdgeTo(DCELVertex const * v) const { return const_cast<DCELVertex *>(this)->getEdgeTo(v); }

    /** Get the edge from this vertex to another, if it exists, else return null. */
    Halfedge * getEdgeTo(DCELVertex const * v)
    {
      Halfedge * rval = NULL;

      if (leaving)
      {
        if (leaving->twin()->getOrigin() == v)
          rval = leaving;
        else
        {
          Halfedge * test = leaving->twin()->next();
          while (rval == NULL && test != leaving)
          {
            Halfedge * twin = test->twin();
            if (twin->getOrigin() == v)
              rval = test;
            else
              test = twin->next();
          }
        }
      }

      return rval;
    }

    /** Get a canonical halfedge leaving this vertex, or null if the vertex is isolated. */
    Halfedge const * getHalfedge() const { return leaving; }

    /** Get a canonical halfedge leaving this vertex, or null if the vertex is isolated. */
    Halfedge * getHalfedge() { return leaving; }

    /** Check if the vertex lies on a mesh boundary. */
    bool isBoundary() const
    {
      // Cycle round the halfedges emanating from this vertex, looking for one on the boundary
      if (leaving)
      {
        Halfedge * e = leaving;
        do
        {
          if (e->isBoundary())
            return true;

          e = e->nextAroundOrigin();

        } while (e != leaving);

        return false;
      }
      else
        return true;
    }

  private:
    template <typename V, typename H, typename F> friend class DCELMesh;

    /**
     * Add normal information from a new face at this vertex.
     *
     * @param n Unit (or weighted) normal of the new face.
     */
    void addFaceNormal(Vector3 const & n)
    {
      if (!has_precomputed_normal)
      {
        Vector3 sum_normals = normal_normalization_factor * getNormal() + n;
        normal_normalization_factor = sum_normals.length();
        setNormal(normal_normalization_factor < 1e-20 ? Vector3::zero() : sum_normals / normal_normalization_factor);
      }
    }

    /** Estimate the result of addFaceNormal() without actually updating any data. */
    Vector3 estimateUpdatedNormal(Vector3 const & new_normal) const
    {
      if (has_precomputed_normal)
        return getNormal();  // it won't change
      else
      {
        Vector3 sum_normals = normal_normalization_factor * getNormal() + new_normal;
        return sum_normals.unit();
      }
    }

    /** Get the packing_index of the vertex. */
    uint32 getPackingIndex() const { return packing_index; }

    /** Unmark the vertex. */
    void setPackingIndex(uint32 index_) { packing_index = index_; }

    Halfedge * leaving;
    bool has_precomputed_normal;
    float normal_normalization_factor;
    uint32 packing_index;

#if 0
    /** Check if one or more internal bits are set. */
    bool areInternalBitsSet(unsigned char mask) const { return ((internal_bits & mask) == mask); };

    /** Set one or more internal bits on or off. */
    void setInternalBits(unsigned char mask, bool value) const
    { internal_bits = (value ? (internal_bits | mask) : (internal_bits & ~mask)); };

    /** Set all internal bits to off. */
    void clearAllInternalBits() const { internal_bits = 0; };

    mutable unsigned char internal_bits;  // only for use by DCELMesh
#endif

}; // class DCELVertex

} // namespace Graphics

namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <typename VT, typename HT, typename FT>
class IsPointN<Graphics::DCELVertex<VT, HT, FT>, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <typename VT, typename HT, typename FT>
class PointTraitsN<Graphics::GeneralMeshVertex<VT, HT, FT>, 3>
{
  public:
    static Vector3 const & getPosition(Graphics::DCELVertex<VT, ET, FT, AT> const & t) { return t.getPosition(); }
};

} // namespace Algorithms

} // namespace Thea

#endif
