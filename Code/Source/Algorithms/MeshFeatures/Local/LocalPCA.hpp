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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_MeshFeatures_Local_LocalPCA_hpp__
#define __Thea_Algorithms_MeshFeatures_Local_LocalPCA_hpp__

#include "../../../Common.hpp"
#include "../SampledSurface.hpp"
#include "../../IntersectionTester.hpp"
#include "../../PcaN.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Noncopyable.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {
namespace MeshFeatures {
namespace Local {

/**
 * Compute local PCA features at a point on a mesh. The local PCA features are eigenvalues of the distribution of samples in the
 * neighborhood of the point, sorted in decreasing order. Optionally, the eigenvectors (principal components) of the
 * distribution may also be returned.
 */
template < typename ExternalSampleKdTreeT = KdTreeN<Vector3, 3> >
class LocalPCA : public SampledSurface<ExternalSampleKdTreeT>
{
  private:
    typedef SampledSurface<ExternalSampleKdTreeT> BaseT;  ///< Base class.
    static intx const DEFAULT_NUM_SAMPLES = 100000;  ///< Default number of points to sample from the shape.

  public:
    /**
     * Constructs the object to compute PCA features at sample points on a given mesh. Initializes internal data structures that
     * do not need to be recomputed for successive calls to compute().
     *
     * @param mesh The mesh representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    LocalPCA(MeshT const & mesh, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute PCA features at sample points on a given mesh group. Initializes internal data
     * structures that do not need to be recomputed for successive calls to compute().
     *
     * @param mesh_group The mesh group representing the shape.
     * @param num_samples The number of samples to compute on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    template <typename MeshT>
    LocalPCA(Graphics::MeshGroup<MeshT> const & mesh_group, intx num_samples = -1, Real normalization_scale = -1)
    : BaseT(mesh_group, (num_samples < 0 ? DEFAULT_NUM_SAMPLES : num_samples), normalization_scale)
    {}

    /**
     * Constructs the object to compute PCA features of a shape with a precomputed kd-tree on points densely sampled from the
     * shape. The kd-tree must persist as long as this object does.
     *
     * @param sample_kdtree_ A kd-tree on a dense set of samples on the shape.
     * @param normalization_scale The scale of the shape, used to define neighborhood sizes. If <= 0, the bounding sphere
     *   diameter will be used.
     */
    LocalPCA(ExternalSampleKdTreeT const * sample_kdtree_, Real normalization_scale = -1)
    : BaseT(sample_kdtree_, normalization_scale)
    {}

    /**
     * Compute the PCA features at a query point on the mesh.
     *
     * @param position Point at which to compute features.
     * @param eigenvectors If non-null, used to return eigenvectors of samples in the neighborhood, sorted in order of
     *   decreasing eigenvalue. Must be pre-allocated to (at least) 3 elements.
     * @param nbd_radius The size of the local neighborhood for which the features are computed, specified as a multiple of the
     *   shape scale. A negative argument selects a default size.
     *
     * @return Eigenvalues of samples in the neighborhood, sorted in decreasing order.
     */
    Vector3 compute(Vector3 const & position, Vector3 * eigenvectors = nullptr, Real nbd_radius = -1) const
    {
      if (nbd_radius <= 0)
        nbd_radius = 0.1f;

      nbd_radius *= this->getNormalizationScale();

      Ball3 range(position, nbd_radius);
      func.reset();

      if (this->hasExternalKdTree())
        this->getMutableExternalKdTree()->template processRangeUntil<IntersectionTester>(range, std::ref(func));
      else
        this->getMutableInternalKdTree()->template processRangeUntil<IntersectionTester>(range, std::ref(func));

      return func.getPCAFeatures(eigenvectors);
    }

  private:
    /** Aggregates points in the neighborhood and computes PCA features. */
    struct LocalPCAFunctor : public Noncopyable
    {
      void reset() { nbd_pts.clear(); }

      template <typename SampleT> bool operator()(intx index, SampleT & t)
      {
        nbd_pts.push_back(PointTraitsN<SampleT, 3>::getPosition(t));
        return false;
      }

      Vector3 getPCAFeatures(Vector3 * eigenvectors) const
      {
        Real eval[3];
        Vector3 evec[3];
        PcaN<Vector3, 3>::compute(nbd_pts.begin(), nbd_pts.end(), eval, evec);  // returns ordered by decreasing eigenvalue

        if (eigenvectors)
        {
          eigenvectors[0] = evec[0];
          eigenvectors[1] = evec[1];
          eigenvectors[2] = evec[2];
        }

        return Vector3(eval[0], eval[1], eval[2]);
      }

      Array<Vector3> nbd_pts;

    }; // struct LocalPCAFunctor

    mutable LocalPCAFunctor func;

}; // class LocalPCA

} // namespace Local
} // namespace MeshFeatures
} // namespace Algorithms
} // namespace Thea

#endif
