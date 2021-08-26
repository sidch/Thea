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

#ifndef __Thea_Algorithms_SurfaceFeatures_Local_LocalPca_hpp__
#define __Thea_Algorithms_SurfaceFeatures_Local_LocalPca_hpp__

#include "../../PointCloud3.hpp"
#include "../../PointTraitsN.hpp"
#include "../../../Noncopyable.hpp"

namespace Thea {
namespace Algorithms {
namespace SurfaceFeatures {
namespace Local {

/**
 * Compute local PCA features at a point on a surface. The local PCA features are eigenvalues of the distribution of samples in
 * the neighborhood of the point, sorted in decreasing order. Optionally, the eigenvectors (principal components) of the
 * distribution may also be returned.
 */
class LocalPca
{
  public:
    /**
     * Constructs the object to compute PCA features at points on a given surface. The sampled surface must persist as long as
     * this object does.
     */
    LocalPca(PointCloud3 const * surf_);

    /** Get the underlying point-sampled surface. */
    PointCloud3 const * getSurface() const { return surf; }

    /**
     * Compute the PCA features at a query point on the surface.
     *
     * @param position Point at which to compute features.
     * @param eigenvectors If non-null, used to return eigenvectors of samples in the neighborhood, sorted in order of
     *   decreasing eigenvalue. Must be pre-allocated to (at least) 3 elements.
     * @param nbd_radius The size of the local neighborhood for which the features are computed, specified as a multiple of the
     *   shape scale. A negative argument selects a default size.
     *
     * @return Eigenvalues of samples in the neighborhood, sorted in decreasing order.
     */
    Vector3 compute(Vector3 const & position, Vector3 * eigenvectors = nullptr, Real nbd_radius = -1) const;

  private:
    /** Aggregates points in the neighborhood and computes PCA features. */
    class LocalPcaFunctor : private Noncopyable
    {
      public:
        void reset() { nbd_pts.clear(); /* no reallocation by spec of std::vector */ }

        template <typename SampleT> bool operator()(intx index, SampleT & t)
        {
          nbd_pts.push_back(PointTraitsN<SampleT, 3>::getPosition(t));
          return false;
        }

        Vector3 getPcaFeatures(Vector3 * eigenvectors) const;

      private:
        Array<Vector3> nbd_pts;

    }; // struct LocalPcaFunctor

    PointCloud3 const * surf;   ///< The point-sampled surface.
    mutable LocalPcaFunctor func;  ///< Functor that is constructed just once to reduce array reallocation costs.

}; // class LocalPca

} // namespace Local
} // namespace SurfaceFeatures
} // namespace Algorithms
} // namespace Thea

#endif
