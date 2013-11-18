//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_Algorithms_MeshFeatures_Curvature_hpp__
#define __Thea_Algorithms_MeshFeatures_Curvature_hpp__

#include "../../Common.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../KDTree3.hpp"
#include "../IntersectionTester.hpp"
#include "../MeshSampler.hpp"
#include "../MetricL2.hpp"
#include "../PointTraitsN.hpp"
#include "../../Math.hpp"
#include "../../Vector3.hpp"
#include <algorithm>

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {

namespace CurvatureInternal {

// A point sample and an associated normal
// TODO: Replace with SurfaceSample class if/when latter is completed?
class SurfaceSample
{
  public:
    SurfaceSample() {}
    SurfaceSample(Vector3 const & p, Vector3 const & n) : position(p), normal(n) {}

    Vector3 const & getPosition() const { return position; }
    Vector3 const & getNormal() const { return normal; }

  private:
    Vector3 position;
    Vector3 normal;
};

} // namespace CurvatureInternal
} // namespace MeshFeatures

template <>
class IsPointN<MeshFeatures::CurvatureInternal::SurfaceSample, 3>
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN<MeshFeatures::CurvatureInternal::SurfaceSample, 3>::getPosition(
  MeshFeatures::CurvatureInternal::SurfaceSample const & sample
)
{
  return sample.getPosition();
}

namespace MeshFeatures {

/** Computes the curvature at a given set of positions on a mesh. */
template < typename MeshT,
           typename ExternalSampleKDTreeT = KDTree3<CurvatureInternal::SurfaceSample> >
class Curvature
{
  public:
    typedef MeshT Mesh;  ///< The mesh class.
    typedef ExternalSampleKDTreeT ExternalSampleKDTree;  ///< A precomputed kd-tree on mesh samples.

  private:
    typedef CurvatureInternal::SurfaceSample SurfaceSample;  ///< A point plus a normal.
    typedef KDTree3<SurfaceSample> SampleKDTree;  ///< A kd-tree on mesh samples.

    static long const DEFAULT_NUM_SAMPLES = 50000;

    /** Get the smoothed normal at a point on a triangle with vertex normals. */
    template <typename TriangleT>
    Vector3 smoothNormal(TriangleT const & tri, Vector3 const & p)
    {
      Vector3 b = tri.barycentricCoordinates(p);

      Vector3 n0 = tri.getVertices().getVertexNormal(0);
      Vector3 n1 = tri.getVertices().getVertexNormal(1);
      Vector3 n2 = tri.getVertices().getVertexNormal(2);

      return b[0] * n0 + b[1] * n1 + b[2] * n2;
    }

    /** Compute sample points with normals on the mesh. */
    void computeSamples(MeshSampler<Mesh> & sampler, long num_samples, TheaArray<SurfaceSample> & samples)
    {
      typedef typename MeshSampler<Mesh>::Triangle Triangle;

      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      TheaArray<Vector3> positions;
      TheaArray<Triangle const *> tris;
      sampler.sampleEvenlyByArea(num_samples, positions, NULL, &tris);

      samples.clear();
      for (array_size_t i = 0; i < positions.size(); ++i)
        samples.push_back(SurfaceSample(positions[i], smoothNormal(*tris[i], positions[i])));
    }

    /** Get the normal at a sample point. */
    template <typename T> struct NormalTraits { static Vector3 getNormal(T const & t) { return t.getNormal(); } };

    /** Get the normal at a sample point referenced by a pointer. */
    template <typename T> struct NormalTraits<T *> { static Vector3 getNormal(T const * t) { return t->getNormal(); } };

  public:
    /**
     * Constructs the object to compute curvature at sample points on a given mesh. The mesh must persist as long as this object
     * does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     */
    Curvature(Mesh const & mesh, long num_samples = -1) : kdtree(new SampleKDTree), precomp_kdtree(NULL), scale(0)
    {
      MeshSampler<Mesh> sampler(mesh);
      TheaArray<SurfaceSample> samples;
      computeSamples(sampler, num_samples, samples);

      kdtree->init(samples.begin(), samples.end());
      scale = kdtree->getBounds().getExtent().length();
    }

    /**
     * Constructs the object to compute curvature at sample points on a given mesh group. The mesh group must persist as long as
     * this object does. Initializes internal data structures that do not need to be recomputed for successive calls to
     * compute().
     */
    Curvature(Graphics::MeshGroup<Mesh> const & mesh_group, long num_samples = -1)
    : kdtree(new SampleKDTree), precomp_kdtree(NULL), scale(0)
    {
      MeshSampler<Mesh> sampler(mesh_group);
      TheaArray<SurfaceSample> samples;
      computeSamples(sampler, num_samples, samples);

      kdtree->init(samples.begin(), samples.end());
      scale = kdtree->getBounds().getExtent().length();
    }

    /**
     * Constructs the object to compute curvature at sample points of a shape with a precomputed kd-tree. The kd-tree must
     * persist as long as this object does.
     */
    Curvature(ExternalSampleKDTree const * kdtree_) : kdtree(NULL), precomp_kdtree(kdtree_), scale(0)
    {
      alwaysAssertM(precomp_kdtree, "Curvature: Precomputed KD-tree cannot be null");
      scale = precomp_kdtree->getBounds().getExtent().length();
    }

    /** Destructor. */
    ~Curvature()
    {
      delete kdtree;
    }

    /**
     * Compute the <em>projected</em> curvatures at a given set of query points on the mesh. The projected curvature is an
     * approximation to the actual curvature, obtained by projecting sample points in the neighborhood of the query point onto
     * the normal.
     *
     * This version of the function explicitly computes the normal at each point -- the other version of the function should be
     * used if the normals are known in advance. The normals are computed from samples, so <b>may be quite inaccurate</b>
     * especially in thin areas.
     *
     * @note The returned curvatures are signed. Positive curvature surfaces curve <em>away</em> from the normal.
     */
    void computeProjectedCurvatures(TheaArray<Vector3> const & positions, TheaArray<Real> & curvatures, Real nbd_radius = -1)
         const
    {
      TheaArray<Vector3> const & normals(positions.size());
      for (array_size_t i = 0; i < positions.size(); ++i)
      {
        long nn_index = precomp_kdtree ? precomp_kdtree->template closestElement<MetricL2>(positions[i])
                                       : kdtree->template closestElement<MetricL2>(positions[i]);
        if (nn_index < 0)
        {
          THEA_WARNING << "Curvature: Query point cannot be mapped to mesh, all curvature values set to zero";
          curvatures.resize(positions.size());
          std::fill(curvatures.begin(), curvatures.end(), 0);
          return;
        }

        normals[i] = precomp_kdtree ? NormalTraits<typename ExternalSampleKDTree::Element>
                                              ::getNormal(precomp_kdtree->getElements()[(array_size_t)nn_index])
                                    : kdtree->getElements()[(array_size_t)nn_index].getNormal();
      }

      computeProjectedCurvatures(positions, normals, curvatures, nbd_radius);
    }

    /**
     * Compute the <em>projected</em> curvatures at a given set of query points with known normals on the mesh. The projected
     * curvature is an approximation to the actual curvature, obtained by projecting sample points in the neighborhood of the
     * query point onto the normal.
     *
     * @note The returned curvatures are signed. Positive curvature surfaces curve <em>away</em> from the normal.
     */
    void computeProjectedCurvatures(TheaArray<Vector3> const & positions, TheaArray<Vector3> const & normals,
                                    TheaArray<Real> & curvatures, Real nbd_radius = -1) const
    {
      alwaysAssertM(positions.size() == normals.size(), "Curvature: Number of sample positions and normals do not match");
      alwaysAssertM(kdtree || precomp_kdtree, "Curvature: KD-tree not initialized for mesh");

      if (nbd_radius <= 0)
        nbd_radius = 0.025f * scale;

      curvatures.resize(positions.size());

      for (array_size_t i = 0; i < curvatures.size(); ++i)
        curvatures[i] = computeProjectedCurvature(positions[i], normals[i], nbd_radius);
    }

  private:
    /** Called for each point in neighborhood. */
    struct ProjectedCurvatureFunctor
    {
      ProjectedCurvatureFunctor(Vector3 const & p, Vector3 const & n)
      : position(p), normal(n), num_offsets(0), sum_offsets(Vector3::zero())
      {}

      template <typename SampleT> bool operator()(long index, SampleT & t)
      {
        if (NormalTraits<SampleT>::getNormal(t).dot(normal) > -1.0e-05f)  // ignore points on hidden side
        {
          Vector3 offset = (PointTraitsN<SampleT, 3>::getPosition(t) - position).fastUnit();
          sum_offsets += offset;
          num_offsets++;
        }

        return false;
      }

      Real getCurvature() const
      {
        return num_offsets > 0 ? -sum_offsets.dot(normal) / num_offsets : 0;
      }

      Vector3 position, normal;
      long num_offsets;
      Vector3 sum_offsets;

    }; // struct ProjectedCurvatureFunctor

    /** Compute the projected curvature at a point. */
    Real computeProjectedCurvature(Vector3 const & p, Vector3 const & n, Real nbd_radius) const
    {
      ProjectedCurvatureFunctor func(p, n);
      Ball3 range(p, nbd_radius);

      if (precomp_kdtree)
        const_cast<ExternalSampleKDTree *>(precomp_kdtree)->template processRangeUntil<IntersectionTester>(range, &func);
      else
        kdtree->template processRangeUntil<IntersectionTester>(range, &func);

      return func.getCurvature();
    }

    SampleKDTree * kdtree;  ///< KD-tree on mesh samples.
    ExternalSampleKDTree const * precomp_kdtree;  ///< Precomputed KD-tree on mesh samples.
    Real scale;  ///< The normalization length.

}; // class Curvature

} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
