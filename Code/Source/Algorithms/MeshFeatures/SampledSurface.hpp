//============================================================================
//
// This file is part of the Thea project.
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

#ifndef __Thea_Algorithms_MeshFeatures_SampledSurface_hpp__
#define __Thea_Algorithms_MeshFeatures_SampledSurface_hpp__

#include "../../Common.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../BestFitSphere3.hpp"
#include "../KDTreeN.hpp"
#include "../MeshSampler.hpp"
#include "../PointCollectorN.hpp"
#include "../PointTraitsN.hpp"
#include "../../Vector3.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {

/** A point sample and an associated normal. */
class SurfaceSample
{
  public:
    /** Default constructor. */
    SurfaceSample() {}

    /** Construct a sample with a given position and normal. */
    SurfaceSample(Vector3 const & p, Vector3 const & n) : position(p), normal(n) {}

    /** Get the position of the sample. */
    Vector3 const & getPosition() const { return position; }

    /** Get the normal of the sample. */
    Vector3 const & getNormal() const { return normal; }

  private:
    Vector3 position;  ///< Sample position.
    Vector3 normal;    ///< Sample normal.
};

} // namespace MeshFeatures

template <>
class IsPointN<MeshFeatures::SurfaceSample, 3>
{
  public:
    static bool const value = true;
};

template <>
inline Vector3
PointTraitsN<MeshFeatures::SurfaceSample, 3>::getPosition(MeshFeatures::SurfaceSample const & sample)
{
  return sample.getPosition();
}

namespace MeshFeatures {

/** Get the normal at a sample point. */
template <typename T> struct NormalTraits { static Vector3 getNormal(T const & t) { return t.getNormal(); } };

// Get the normal at a sample point referenced by a pointer.
template <typename T> struct NormalTraits<T *> { static Vector3 getNormal(T const * t) { return t->getNormal(); } };

// Get the normal at a sample point represented as a Vector3.
template <> struct NormalTraits<Vector3> { static Vector3 getNormal(Vector3 const & t) { return Vector3::zero(); } };

/** A representation of a point-sampled surface, used for computing features. */
template <typename ExternalSampleKDTreeT>
class SampledSurface
{
  public:
    typedef ExternalSampleKDTreeT ExternalSampleKDTree;  ///< A precomputed kd-tree on mesh samples.

  protected:
    typedef KDTreeN<SurfaceSample, 3> InternalSampleKDTree;  ///< A kd-tree on mesh samples.

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

      static long const DEFAULT_NUM_SAMPLES = 5000;
      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      TheaArray<Vector3> positions;
      TheaArray<Triangle const *> tris;
      sampler.sampleEvenlyByArea(num_samples, positions, NULL, &tris);

      samples.clear();
      for (array_size_t i = 0; i < positions.size(); ++i)
        samples.push_back(SurfaceSample(positions[i], smoothNormal(*tris[i], positions[i])));
    }

    /**
     * Generates a point-sampled surface from a mesh.
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape. If zero, no samples are generated.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    SampledSurface(MeshT const & mesh, long num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(NULL), precomp_kdtree(NULL), scale(normalization_scale)
    {
      MeshSampler<MeshT> sampler(mesh);
      computeSamples(sampler, num_samples, samples);

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mesh);
        scale = bsphere.getDiameter();
      }
    }

    /**
     * Generates a point-sampled surface from a mesh group.
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape. If zero, no samples are generated.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    SampledSurface(Graphics::MeshGroup<MeshT> const & mesh_group, long num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(NULL), precomp_kdtree(NULL), scale(normalization_scale)
    {
      MeshSampler<MeshT> sampler(mesh_group);
      computeSamples(sampler, num_samples, samples);

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
    SampledSurface(ExternalSampleKDTree const * sample_kdtree_, Real normalization_scale = -1)
    : sample_kdtree(NULL), precomp_kdtree(sample_kdtree_), scale(normalization_scale)
    {
      alwaysAssertM(precomp_kdtree, "SampledSurface: Precomputed KD-tree cannot be null");

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addPoints(precomp_kdtree->getElements(),
                                                               precomp_kdtree->getElements() + precomp_kdtree->numElements());
        scale = bsphere.getDiameter();
      }
    }

    /** Get a non-const reference to the kd-tree on internally generated samples, or null if no such samples exist. */
    InternalSampleKDTree * getMutableInternalKDTree() const
    {
      if (!sample_kdtree && !samples.empty())
        sample_kdtree = new InternalSampleKDTree(samples.begin(), samples.end());

      return sample_kdtree;
    }

    /**
     * Get the kd-tree on internally generated samples, or null if no such samples exist. <b>Use with caution</b> -- you should
     * never actually change the tree even if you have a non-const handle to it!
     */
    InternalSampleKDTree const * getInternalKDTree() const
    {
      return const_cast<SampledSurface *>(this)->getMutableInternalKDTree();
    }

    /** Check if the shape has a kd-tree on externally generated samples. */
    bool hasExternalKDTree() const
    {
      return (bool)precomp_kdtree;
    }

    /**
     * Get a non-const reference to the kd-tree on externally generated samples, or null if no such samples exist. <b>Use with
     * caution</b> -- you should never actually change the tree even if you have a non-const handle to it!
     */
    ExternalSampleKDTree * getMutableExternalKDTree() const
    {
      return const_cast<ExternalSampleKDTree *>(precomp_kdtree);
    }

    /** Get the kd-tree on externally generated samples, or null if no such kd-tree exists. */
    ExternalSampleKDTree const * getExternalKDTree() const
    {
      return precomp_kdtree;
    }

  public:
    /** Destructor. */
    ~SampledSurface()
    {
      delete sample_kdtree;
    }

    /** Get the number of surface samples. */
    long numSamples() const
    {
      if (precomp_kdtree)
        return precomp_kdtree->numElements();
      else
        return (long)samples.size();
    }

    /** Get the position of the surface sample with index \a index. */
    Vector3 getSamplePosition(long index) const
    {
      debugAssertM(index >= 0 && index < numSamples(), format("SampledSurface: Sample index %ld out of bounds", index));

      if (precomp_kdtree)
      {
        typedef typename ExternalSampleKDTree::Element ExternalSample;
        return PointTraitsN<ExternalSample, 3>::getPosition(precomp_kdtree->getElements()[(array_size_t)index]);
      }
      else
        return samples[(array_size_t)index].getPosition();
    }

    /** Get the normal of the surface sample with index \a index. */
    Vector3 getSampleNormal(long index) const
    {
      debugAssertM(index >= 0 && index < numSamples(), format("SampledSurface: Sample index %ld out of bounds", index));

      if (precomp_kdtree)
      {
        typedef typename ExternalSampleKDTree::Element ExternalSample;
        ExternalSample const & e = precomp_kdtree->getElements()[(array_size_t)index];
        return NormalTraits<ExternalSample>::getNormal(e);
      }
      else
        return samples[(array_size_t)index].getNormal();
    }

    /** Get the normalization scale of the shape. */
    Real getNormalizationScale() const { return scale; }

  private:
    TheaArray<SurfaceSample> samples;  ///< Set of internally-generated surface samples.
    mutable InternalSampleKDTree * sample_kdtree;  ///< kd-tree on mesh samples.
    ExternalSampleKDTree const * precomp_kdtree;  ///< Precomputed kd-tree on mesh samples.
    Real scale;  ///< The normalization length.

}; // class SampledSurface

} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
