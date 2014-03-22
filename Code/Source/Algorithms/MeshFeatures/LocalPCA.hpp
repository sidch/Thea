//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2014, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_MeshFeatures_LocalPCA_hpp__
#define __Thea_Algorithms_MeshFeatures_LocalPCA_hpp__

#include "../../Common.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../BestFitSphere3.hpp"
#include "../IntersectionTester.hpp"
#include "../KDTreeN.hpp"
#include "../MeshSampler.hpp"
#include "../PCA_N.hpp"
#include "../PointCollectorN.hpp"
#include "../../Vector3.hpp"

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {

/**
 * Compute local PCA features at a point on a mesh. The local PCA features are eigenvalues of the distribution of samples in the
 * neighborhood of the point, sorted in decreasing order. Optionally, the eigenvectors (principal components) of the
 * distribution may also be returned.
 */
template < typename MeshT,
           typename ExternalSampleKDTreeT = KDTreeN<Vector3, 3> >
class LocalPCA
{
  public:
    typedef MeshT Mesh;  ///< The mesh class.
    typedef ExternalSampleKDTreeT ExternalSampleKDTree;  ///< A precomputed kd-tree on mesh samples.

  private:
    typedef KDTreeN<Vector3, 3> SampleKDTree;  ///< A kd-tree on mesh samples.

    static long const DEFAULT_NUM_SAMPLES = 50000;  ///< Default number of points to sample from the shape.

    /** Compute sample points with normals on the mesh. */
    void computeSamples(MeshSampler<Mesh> & sampler, long num_samples, TheaArray<Vector3> & samples)
    {
      if (num_samples < 0) num_samples = DEFAULT_NUM_SAMPLES;

      samples.clear();
      sampler.sampleEvenlyByArea(num_samples, samples);
    }

  public:
    /**
     * Constructs the object to compute PCA features at sample points on a given mesh. The mesh must persist as long as this
     * object does. Initializes internal data structures that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    LocalPCA(Mesh const & mesh, long num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(new SampleKDTree), precomp_kdtree(NULL), scale(normalization_scale)
    {
      MeshSampler<Mesh> sampler(mesh);
      TheaArray<Vector3> samples;
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
     * Constructs the object to compute PCA features at sample points on a given mesh group. The mesh group must persist as long
     * as this object does. Initializes internal data structures that do not need to be recomputed for successive calls to
     * compute().
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    LocalPCA(Graphics::MeshGroup<Mesh> const & mesh_group, long num_samples = -1, Real normalization_scale = -1)
    : sample_kdtree(new SampleKDTree), precomp_kdtree(NULL), scale(normalization_scale)
    {
      MeshSampler<Mesh> sampler(mesh_group);
      TheaArray<Vector3> samples;
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
     * Constructs the object to compute PCA features of a shape with a precomputed kd-tree on points densely sampled from the
     * shape. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ A kd-tree on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    LocalPCA(ExternalSampleKDTree const * sample_kdtree_, Real normalization_scale = -1)
    : sample_kdtree(NULL), precomp_kdtree(sample_kdtree_), scale(normalization_scale)
    {
      alwaysAssertM(precomp_kdtree, "LocalPCA: Precomputed KD-tree cannot be null");

      if (scale <= 0)
      {
        BestFitSphere3 bsphere;
        PointCollectorN<BestFitSphere3, 3>(&bsphere).addPoints(precomp_kdtree->getElements(),
                                                               precomp_kdtree->getElements() + precomp_kdtree->numElements());
        scale = bsphere.getDiameter();
      }
    }

    /** Destructor. */
    ~LocalPCA()
    {
      delete sample_kdtree;
    }

    /**
     * Compute the PCA features at a query point on the mesh.
     *
     * @param position Point at which to compute features.
     * @param eigenvectors If non-null, used to return eigenvectors of samples in the neighborhood, sorted in order of
     *   decreasing eigenvalue. Must be pre-allocated to (at least) 3 elements.
     *
     * @return Eigenvalues of samples in the neighborhood, sorted in decreasing order.
     */
    Vector3 compute(Vector3 const & position, Vector3 * eigenvectors = NULL, Real nbd_radius = -1) const
    {
      if (nbd_radius <= 0)
        nbd_radius = 0.05f * scale;

      Ball3 range(position, nbd_radius);
      func.reset();

      if (precomp_kdtree)
        const_cast<ExternalSampleKDTree *>(precomp_kdtree)->template processRangeUntil<IntersectionTester>(range, &func);
      else
        sample_kdtree->template processRangeUntil<IntersectionTester>(range, &func);

      return func.getPCAFeatures(eigenvectors);
    }

  private:
    /** Aggregates points in the neighborhood and computes PCA features. */
    struct LocalPCAFunctor
    {
      void reset() { nbd_pts.clear(); }

      template <typename SampleT> bool operator()(long index, SampleT & t)
      {
        nbd_pts.push_back(PointTraitsN<SampleT, 3>::getPosition(t));
        return false;
      }

      Vector3 getPCAFeatures(Vector3 * eigenvectors) const
      {
        Real eval[3];
        Vector3 evec[3];
        PCA_N<Vector3, 3>::compute(nbd_pts.begin(), nbd_pts.end(), eval, evec);  // returns ordered by decreasing eigenvalue

        if (eigenvectors)
        {
          eigenvectors[0] = evec[0];
          eigenvectors[1] = evec[1];
          eigenvectors[2] = evec[2];
        }

        return Vector3(eval[0], eval[1], eval[2]);
      }

      TheaArray<Vector3> nbd_pts;

    }; // struct LocalPCAFunctor

    SampleKDTree * sample_kdtree;  ///< KD-tree on mesh samples.
    ExternalSampleKDTree const * precomp_kdtree;  ///< Precomputed KD-tree on mesh samples.
    Real scale;  ///< The normalization length.

    mutable LocalPCAFunctor func;

}; // class LocalPCA

} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
