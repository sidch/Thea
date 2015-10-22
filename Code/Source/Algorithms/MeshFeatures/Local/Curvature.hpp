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

#ifndef __Thea_Algorithms_MeshFeatures_Local_Curvature_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_Curvature_hpp__

#include "../../../Common.hpp"
#include "../../../Graphics/MeshGroup.hpp"
#include "../../BestFitSphere3.hpp"
#include "../../IntersectionTester.hpp"
#include "../../KDTreeN.hpp"
#include "../../MeshSampler.hpp"
#include "../../MetricL2.hpp"
#include "../../PointCollectorN.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Vector3.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

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
} // namespace Local
} // namespace MeshFeatures

template <>
class IsPointN<MeshFeatures::Local::CurvatureInternal::SurfaceSample, 3>
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN<MeshFeatures::Local::CurvatureInternal::SurfaceSample, 3>::getPosition(
  MeshFeatures::Local::CurvatureInternal::SurfaceSample const & sample
)
{
  return sample.getPosition();
}

namespace MeshFeatures {
namespace Local {

/** Compute the curvature at a point on a shape. */
template < typename ExternalSampleKDTreeT = KDTreeN<CurvatureInternal::SurfaceSample, 3> >
class Curvature
{
  public:
    typedef ExternalSampleKDTreeT ExternalSampleKDTree;  ///< A precomputed kd-tree on mesh samples.

  private:
    typedef CurvatureInternal::SurfaceSample SurfaceSample;  ///< A point plus a normal.
    typedef KDTreeN<SurfaceSample, 3> SampleKDTree;  ///< A kd-tree on mesh samples.

    static long const DEFAULT_NUM_SAMPLES = 50000;  ///< Default number of points to sample from the shape.

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
    template <typename MeshT>
    void computeSamples(MeshSampler<MeshT> & sampler, long num_samples, TheaArray<SurfaceSample> & samples)
    {
      typedef typename MeshSampler<MeshT>::Triangle Triangle;

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
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    Curvature(MeshT const & mesh, long num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(new SampleKDTree), precomp_kdtree(NULL), scale(normalization_scale)
    {
      MeshSampler<MeshT> sampler(mesh);
      TheaArray<SurfaceSample> samples;
      computeSamples(sampler, num_samples, samples);

      sample_kdtree->init(samples.begin(), samples.end());

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh);
        scale = bsphere.getDiameter();
      }
    }

    /**
     * Constructs the object to compute curvature at sample points on a given mesh group. The mesh group must persist as long as
     * this object does. Initializes internal data structures that do not need to be recomputed for successive calls to
     * compute().
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    Curvature(Graphics::MeshGroup<MeshT> const & mesh_group, long num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(new SampleKDTree), precomp_kdtree(NULL), scale(normalization_scale)
    {
      MeshSampler<MeshT> sampler(mesh_group);
      TheaArray<SurfaceSample> samples;
      computeSamples(sampler, num_samples, samples);

      sample_kdtree->init(samples.begin(), samples.end());

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh_group);
        scale = bsphere.getDiameter();
      }
    }

    /**
     * Constructs the object to compute curvature of a shape with a precomputed kd-tree on points densely sampled from the
     * shape. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ A kd-tree on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    Curvature(ExternalSampleKDTree const * sample_kdtree_, Real normalization_scale = -1)
    : sample_kdtree(NULL), precomp_kdtree(sample_kdtree_), scale(normalization_scale)
    {
      alwaysAssertM(precomp_kdtree, "Curvature: Precomputed KD-tree cannot be null");

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addPoints(precomp_kdtree->getElements(),
                                                               precomp_kdtree->getElements() + precomp_kdtree->numElements());
        scale = bsphere.getDiameter();
      }
    }

    /** Destructor. */
    ~Curvature()
    {
      delete sample_kdtree;
    }

    /**
     * Compute the <em>projected</em> curvature at a query point on the mesh. The projected curvature is an approximation to the
     * actual curvature, obtained by projecting sample points in the neighborhood of the query point onto the normal.
     *
     * This version of the function explicitly computes the normal at the query point -- the other version of the function
     * should be used if the normal is known in advance. The normal is computed from samples, so <b>may be quite inaccurate</b>
     * especially in thin areas.
     *
     * @note The returned curvature is signed. Positive curvature surfaces curve <em>away</em> from the normal.
     */
    double computeProjectedCurvature(Vector3 const & position, Real nbd_radius = -1) const
    {
      long nn_index = precomp_kdtree ? precomp_kdtree->template closestElement<MetricL2>(position)
                                     : sample_kdtree->template closestElement<MetricL2>(position);
      if (nn_index < 0)
      {
        THEA_WARNING << "Curvature: Query point cannot be mapped to mesh, curvature value set to zero";
        return 0.0;
      }

      Vector3 normal = precomp_kdtree ? NormalTraits<typename ExternalSampleKDTree::Element>
                                            ::getNormal(precomp_kdtree->getElements()[(array_size_t)nn_index])
                                      : sample_kdtree->getElements()[(array_size_t)nn_index].getNormal();

      return computeProjectedCurvature(position, normal, nbd_radius);
    }

    /**
     * Compute the <em>projected</em> curvature at a query point with a known normal on the mesh. The projected curvature is an
     * approximation to the actual curvature, obtained by projecting sample points in the neighborhood of the query point onto
     * the normal.
     *
     * @note The returned curvature is signed. Positive curvature surfaces curve <em>away</em> from the normal.
     */
    double computeProjectedCurvature(Vector3 const & position, Vector3 const & normal, Real nbd_radius = -1) const
    {
      if (nbd_radius <= 0)
        nbd_radius = 0.05f * scale;

      ProjectedCurvatureFunctor func(position, normal);
      Ball3 range(position, nbd_radius);

      if (precomp_kdtree)
        const_cast<ExternalSampleKDTree *>(precomp_kdtree)->template processRangeUntil<IntersectionTester>(range, &func);
      else
        sample_kdtree->template processRangeUntil<IntersectionTester>(range, &func);

      return func.getCurvature();
    }

  private:
    /** Called for each point in the neighborhood. */
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

      double getCurvature() const
      {
        return num_offsets > 0 ? -sum_offsets.dot(normal) / num_offsets : 0;
      }

      Vector3 position, normal;
      long num_offsets;
      Vector3 sum_offsets;

    }; // struct ProjectedCurvatureFunctor

    SampleKDTree * sample_kdtree;  ///< KD-tree on mesh samples.
    ExternalSampleKDTree const * precomp_kdtree;  ///< Precomputed KD-tree on mesh samples.
    Real scale;  ///< The normalization length.

}; // class Curvature

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
