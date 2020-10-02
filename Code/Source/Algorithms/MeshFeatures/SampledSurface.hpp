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

#ifndef __Thea_Algorithms_MeshFeatures_SampledSurface_hpp__
#define __Thea_Algorithms_MeshFeatures_SampledSurface_hpp__

#include "../../Common.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../BestFitSphere3.hpp"
#include "../KdTreeN.hpp"
#include "../MeshSampler.hpp"
#include "../PointCollectorN.hpp"
#include "../PointTraitsN.hpp"
#include "../SampleGraph.hpp"
#include "../../MatVec.hpp"

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
template <> struct NormalTraits<Vector3> { static Vector3 getNormal(Vector3 const & t) { return Vector3::Zero(); } };

/** A representation of a point-sampled surface, used for computing features. */
template <typename ExternalSampleKdTreeT>
class SampledSurface
{
  public:
    typedef ExternalSampleKdTreeT ExternalSampleKdTree;  ///< A precomputed kd-tree on mesh samples.

  protected:
    typedef KdTreeN<SurfaceSample, 3> InternalSampleKdTree;  ///< A kd-tree on mesh samples.

  private:
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
    void computeSamples(MeshSampler<MeshT> & sampler, intx num_samples, Array<SurfaceSample> & samples)
    {
      typedef typename MeshSampler<MeshT>::Triangle Triangle;

      static intx const DEFAULT_NUM_SAMPLES = 5000;
      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      Array<Vector3> positions;
      Array<Triangle const *> tris;
      sampler.sampleEvenlyByArea(num_samples, positions, nullptr, &tris);

      samples.clear();
      for (size_t i = 0; i < positions.size(); ++i)
        samples.push_back(SurfaceSample(positions[i], smoothNormal(*tris[i], positions[i])));
    }

    /** Compute the scale of the shape as the diameter of the bounding sphere of the samples. */
    void computeScaleFromSamples()
    {
      BestFitSphere3 bsphere;
      PointCollectorN<BestFitSphere3, 3>(&bsphere).addPoints(samples.begin(), samples.end());
      scale = bsphere.getDiameter();
    }

  protected:
    /**
     * Initializes the object from a set of precomputed samples.
     *
     * @param num_samples The number of precomputed samples.
     * @param positions The positions of the samples.
     * @param normals The normals of the samples.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    SampledSurface(intx num_samples, Vector3 const * positions, Vector3 const * normals, Real normalization_scale = -1)
    : sample_kdtree(nullptr), precomp_kdtree(nullptr), owns_sample_graph(true), sample_graph(nullptr),
      scale(normalization_scale)
    {
      alwaysAssertM(num_samples >= 0,   "SampledSurface: Number of precomputed samples must be non-negative");
      alwaysAssertM(positions != nullptr,  "SampledSurface: Null array of sample positions");
      alwaysAssertM(normals != nullptr,    "SampledSurface: Null array of sample normals");

      samples.resize((size_t)num_samples);
      for (size_t i = 0; i < samples.size(); ++i)
        samples[i] = SurfaceSample(positions[i], normals[i]);

      if (scale <= 0)
        computeScaleFromSamples();
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
    SampledSurface(MeshT const & mesh, intx num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(nullptr), precomp_kdtree(nullptr), owns_sample_graph(true), sample_graph(nullptr),
      scale(normalization_scale)
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
    SampledSurface(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(nullptr), precomp_kdtree(nullptr), owns_sample_graph(true), sample_graph(nullptr),
      scale(normalization_scale)
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
     * Constructs the surface from a precomputed kd-tree on points sampled from the surface. The kd-tree must persist as long as
     * this object does.
     *
     * @param sample_kdtree_ A kd-tree on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    SampledSurface(ExternalSampleKdTree const * sample_kdtree_, Real normalization_scale = -1)
    : sample_kdtree(nullptr), precomp_kdtree(sample_kdtree_), owns_sample_graph(true), sample_graph(nullptr),
      scale(normalization_scale)
    {
      alwaysAssertM(precomp_kdtree, "SampledSurface: Precomputed KD-tree cannot be null");

      // Cache the external samples for quick access to positions and normals
      typedef typename ExternalSampleKdTree::Element ExternalSample;
      ExternalSample const * ext_samples = precomp_kdtree->getElements();
      intx num_samples = precomp_kdtree->numElements();

      samples.resize((size_t)num_samples);
      for (size_t i = 0; i < samples.size(); ++i)
      {
        samples[i] = SurfaceSample(PointTraitsN<ExternalSample, 3>::getPosition(ext_samples[i]),
                                   NormalTraits<ExternalSample>::getNormal(ext_samples[i]));
      }

      if (scale <= 0)
        computeScaleFromSamples();
    }

    /**
     * Constructs the surface from a precomputed adjacency graph on points sampled from the surface. The graph must persist as
     * long as this object does.
     *
     * @param sample_graph_ A graph on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    SampledSurface(SampleGraph const * sample_graph_, Real normalization_scale = -1)
    : sample_kdtree(nullptr), precomp_kdtree(nullptr), owns_sample_graph(false),
      sample_graph(const_cast<SampleGraph *>(sample_graph_)), scale(normalization_scale)
    {
      alwaysAssertM(sample_graph, "SampledSurface: Sample graph cannot be null");

      // Cache the external samples for quick access to positions and normals
      SampleGraph::SampleArray const & ext_samples = sample_graph->getSamples();
      samples.resize(ext_samples.size());
      for (size_t i = 0; i < samples.size(); ++i)
      {
        samples[i] = SurfaceSample(ext_samples[i].getPosition(), ext_samples[i].getNormal());
      }

      if (scale <= 0)
        computeScaleFromSamples();
    }

    /** Get the array of samples. */
    Array<SurfaceSample> const & getSamples() const { return samples; }

    /** Get a particular sample. */
    SurfaceSample const & getSample(intx index) const
    {
      debugAssertM(index >= 0 && index < (intx)samples.size(), "SampledSurface: Sample index out of bounds");
      return samples[(size_t)index];
    }

    /** Get a non-const reference to the kd-tree on internally generated samples, or null if no such samples exist. */
    InternalSampleKdTree * getMutableInternalKdTree() const
    {
      if (!sample_kdtree && !samples.empty())
        sample_kdtree = new InternalSampleKdTree(samples.begin(), samples.end());

      return sample_kdtree;
    }

    /**
     * Get the kd-tree on internally generated samples, or null if no such samples exist. <b>Use with caution</b> -- you should
     * never actually change the tree even if you have a non-const handle to it!
     */
    InternalSampleKdTree const * getInternalKdTree() const
    {
      return const_cast<SampledSurface *>(this)->getMutableInternalKdTree();
    }

    /** Check if the shape has a kd-tree on externally generated samples. */
    bool hasExternalKdTree() const
    {
      return (bool)precomp_kdtree;
    }

    /**
     * Get a non-const reference to the kd-tree on externally generated samples, or null if no such samples exist. <b>Use with
     * caution</b> -- you should never actually change the tree even if you have a non-const handle to it!
     */
    ExternalSampleKdTree * getMutableExternalKdTree() const
    {
      return const_cast<ExternalSampleKdTree *>(precomp_kdtree);
    }

    /** Get the kd-tree on externally generated samples, or null if no such kd-tree exists. */
    ExternalSampleKdTree const * getExternalKdTree() const
    {
      return precomp_kdtree;
    }

    /** Get an adjacency graph on surface samples, creating it from scratch if it was not specified in the constructor. */
    SampleGraph const * getSampleGraph(SampleGraph::Options const & options = SampleGraph::Options::defaults()) const
    {
      if (sample_graph)
        return sample_graph;

      size_t n = (size_t)numSamples();
      Array<Vector3> positions(n);
      Array<Vector3> normals(n);

      for (size_t i = 0; i < n; ++i)
      {
        positions[i] = getSamplePosition((intx)i);
        normals[i] = getSampleNormal((intx)i);
      }

      sample_graph = new SampleGraph(options);

      if (n > 0)
        sample_graph->setSamples((intx)n, &positions[0], &normals[0]);

      sample_graph->init();

      return sample_graph;
    }

  public:
    /** Destructor. */
    ~SampledSurface()
    {
      delete sample_kdtree;
      if (owns_sample_graph) delete sample_graph;
    }

    /** Get the number of surface samples. */
    intx numSamples() const
    {
      return (intx)samples.size();
    }

    /** Get the position of the surface sample with index \a index. */
    Vector3 getSamplePosition(intx index) const
    {
      debugAssertM(index >= 0 && index < numSamples(), format("SampledSurface: Sample index %ld out of bounds", index));
      return samples[(size_t)index].getPosition();
    }

    /** Get the normal of the surface sample with index \a index. */
    Vector3 getSampleNormal(intx index) const
    {
      debugAssertM(index >= 0 && index < numSamples(), format("SampledSurface: Sample index %ld out of bounds", index));
      return samples[(size_t)index].getNormal();
    }

    /** Get the normalization scale of the shape. */
    Real getNormalizationScale() const { return scale; }

  private:
    Array<SurfaceSample> samples;  ///< Set of internally-generated surface samples.
    mutable InternalSampleKdTree * sample_kdtree;  ///< kd-tree on surface samples.
    ExternalSampleKdTree const * precomp_kdtree;  ///< Precomputed kd-tree on surface samples.
    bool owns_sample_graph;  ///< Was the sample graph precomputed?
    mutable SampleGraph * sample_graph;  ///< Precomputed or internally computed graph on surface samples.
    Real scale;  ///< The normalization length.

}; // class SampledSurface

} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
