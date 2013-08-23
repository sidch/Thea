//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Graphics_SurfaceSample_hpp__
#define __Thea_Graphics_SurfaceSample_hpp__

#include "../Common.hpp"
#include "../Algorithms/KDTree3.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "RenderOptions.hpp"
#include "RenderSystem.hpp"
#include <Thea/BoundedSortedArray.hpp>
#include <boost/type_traits/has_trivial_assign.hpp>
#include <functional>

namespace Thea {
namespace Graphics {

class SurfaceSample;

namespace SurfaceSampleInternal {

class Neighbor
{
  private:
    SurfaceSample * sample;
    Real separation;

  public:
    Neighbor(SurfaceSample * sample_ = NULL, Real separation_ = 0) : sample(sample_), separation(separation_) {}

    bool operator<(Neighbor const & rhs) const
    {
      return separation < rhs.separation || (separation == rhs.separation && std::less<SurfaceSample *>()(sample, rhs.sample));
    }

}; // struct Neighbor

} // namespace SurfaceSampleInternal

} // namespace AE

namespace boost
{

template <> struct has_trivial_assign<AE::SurfaceSampleInternal::Neighbor> : public boost::true_type {};

} // namespace boost

namespace AE {

/** A point sample on the surface of a shape. */
class SurfaceSample
{
  public:
    static int const MAX_NEIGHBORS = 8;

    typedef SurfaceSampleInternal::Neighbor Neighbor;
    typedef BoundedSortedArrayN<MAX_NEIGHBORS, Neighbor> NeighborSet;

    SurfaceSample(PartWeakPtr part_ = PartWeakPtr(), Vector3 const & position_ = Vector3::zero(),
                  Vector3 const & normal_ = Vector3::zero(), Slot * src_slot_= NULL)
    : part(part_), position(position_), normal(normal_), src_slot(src_slot_), original(NULL), linked(NULL),
      nearest_vertex(-1), deformable(false), neighbors_found(false), data(NULL)
    {}

    SurfaceSample clone() const
    {
      SurfaceSample dst = *this;
      dst.original = this;
      return dst;
    }

    PartWeakPtr getPart() const { return part; }
    void setPart(PartWeakPtr part_) { part = part_; }

    Vector3 const & getPosition() const { return position; }
    void setPosition(Vector3 const & position_) { position = position_; }

    Vector3 const & getNormal() const { return normal; }
    void setNormal(Vector3 const & normal_) { normal = normal_; }

    Slot * getSourceSlot() const { return src_slot; }
    void setSourceSlot(Slot * src_slot_) { src_slot = src_slot_; }

    bool isClone() const { return (bool)original; }
    SurfaceSample const * getOriginal() const { return original; }

    SurfaceSample * getLinked() const { return linked; }
    void linkTo(SurfaceSample * sample) { linked = sample; }

    long getNearestVertex() const { return nearest_vertex; }
    bool isDeformable() const { return deformable; }

    NeighborSet const & getNeighbors() const { return neighbors; }
    NeighborSet & getNeighbors() { return neighbors; }

    bool neighborsFound() const { return neighbors_found; }
    bool & neighborsFound() { return neighbors_found; }

    int numNeighbors() const { return neighbors.size(); }
    Real maxNeighborSeparation() const { return neighbors.isEmpty() ? -1 : neighbors.last().separation; }
    void addNeighbor(SurfaceSample * neighbor);
    void addNeighbor(Neighbor neighbor) { neighbors.insert(neighbor); }
    void clearNeighbors() { neighbors.clear(); neighbors_found = false; }

    void setNearestVertex(long v) { nearest_vertex = v; }
    void setDeformable(bool value) { deformable = value; }

    // Temporary pointer to bookkeeping information, generated for example during parametrization.
    void * getData() const { return data; }
    void setData(void * data_) const { data = data_; }

    void transform(AffineTransform3 const & trans, Matrix3 const & trans_inverse_transpose)
    {
      position = trans * position;
      normal = (trans_inverse_transpose * normal).fastUnit();
    }

    template <typename InputIterator>
    static void drawSamples(InputIterator begin, InputIterator end, Graphics::RenderSystem & render_system)
    {
      render_system.pushShader();
      render_system.pushColorFlags();
      render_system.pushShapeFlags();

        render_system.setShader(NULL);
        render_system.setPointSize(4);

        render_system.beginPrimitive(Graphics::RenderSystem::Primitive::POINTS);
          for (InputIterator si = begin; si != end; ++si)
            drawSample(*si, render_system);
        render_system.endPrimitive();

      render_system.popShapeFlags();
      render_system.popColorFlags();
      render_system.popShader();
    }

  private:
    static void drawSample(SurfaceSample const & sample, Graphics::RenderSystem & render_system)
    {
      render_system.sendVertex(sample.position);
    }

    static void drawSample(SurfaceSample const * sample, Graphics::RenderSystem & render_system)
    {
      drawSample(*sample, render_system);
    }

    PartWeakPtr part;
    Vector3 position;
    Vector3 normal;
    Slot * src_slot;
    SurfaceSample const * original;
    SurfaceSample * linked;
    long nearest_vertex;
    bool deformable;
    bool neighbors_found;
    NeighborSet neighbors;

    mutable void * data;

}; // class SurfaceSample

} // namespace AE

namespace Thea {
namespace Algorithms {

template <>
class IsPointN<AE::SurfaceSample, 3>
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN<AE::SurfaceSample, 3>::getPosition(AE::SurfaceSample const & sample)
{
  return sample.getPosition();
}

} // namespace Algorithms
} // namespace Thea

namespace AE {

typedef TheaArray<SurfaceSample> SampleArray;  ///< Array of surface sample points.
typedef Algorithms::KDTree3<SurfaceSample *> SampleKDTree;  ///< KD-tree of surface sample points.

} // namespace AE

#endif
