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
// First version: 2011
//
//============================================================================

#ifndef __Thea_HyperplaneN_hpp__
#define __Thea_HyperplaneN_hpp__

#include "Common.hpp"
#include "Array.hpp"
#include "Math.hpp"
#include "MatVec.hpp"
#include "RayIntersectableN.hpp"

namespace Thea {

// Forward declarations
template <int N, typename T> class HyperplaneN;

namespace Internal {

/**
 * <b>[Internal]</b> Base class for hyperplanes ((N - 1)-flats) in N-dimensional space, where N is any <b>positive</b>
 * (non-zero) integer and T is a field.
 *
 * @note This class is <b>INTERNAL</b>! Don't use it directly.
 */
template <int N, typename T>
class /* THEA_DLL_LOCAL */ HyperplaneNBase : public RayIntersectableN<N, T>
{
  public:
    typedef HyperplaneN<N, T>  HyperplaneT;  ///< N-dimensional hyperplane.
    typedef Vector<N, T>       VectorT;      ///< N-dimensional vector.

    THEA_DECL_SMART_POINTERS(HyperplaneT)

    /** Default constructor. */
    HyperplaneNBase() : normal(VectorT::Zero()), dist(0) {}

    /**
     * Construct a hyperplane from its signed distance from the origin, and its normal vector (need not be a unit vector). The
     * \a normalize argument suppresses rescaling of the normal to unit length if set to false.
     */
    static HyperplaneT fromDistanceAndNormal(T dist_, VectorT const & normal_, bool normalize = true)
    {
      HyperplaneT hyperplane;
      hyperplane.normal  =  normal_; if (normalize) hyperplane.normal.normalize();
      hyperplane.dist    =  dist_;

      return hyperplane;
    }

    /**
     * Construct a hyperplane from a point on it, and its normal vector (need not be a unit vector). The \a normalize argument
     * suppresses rescaling of the normal to unit length if set to false.
     */
    static HyperplaneT fromPointAndNormal(VectorT const & point_, VectorT const & normal_, bool normalize = true)
    {
      HyperplaneT hyperplane;
      hyperplane.normal  =  normal_; if (normalize) hyperplane.normal.normalize();
      hyperplane.dist    =  hyperplane.normal.dot(point_);

      return hyperplane;
    }

    /** Construct a hyperplane from N points on it. */
    static HyperplaneT fromNPoints(Array<VectorT> const & points)
    {
      alwaysAssertM(points.size() >= N,
                    format("HyperplaneN: Too few points specified (provided %ld points, hyperplane requires %ld points)",
                           (intx)points.size(), N));

      Matrix<N, N, T> a;
      for (intx i = 0; i < N; ++i)
        a.setRow(i, points[i]);

      try
      {
        a = a.inverse();
      }
      catch (...)
      {
        throw Error("HyperplaneN: Points are degenerate");
      }

      HyperplaneT hyperplane;
      hyperplane.normal  =  -(a * VectorT::Ones()).normalized();
      hyperplane.dist    =  hyperplane.normal.dot(points[0]);
      return hyperplane;
    }

    /** Cast the hyperplane to a different scalar type. */
    template <typename U> HyperplaneN<N, U> cast() const
    {
      return HyperplaneN<N, U>::fromDistanceAndNormal(static_cast<U>(dist), normal.template cast<U>(), /* normalize = */ false);
    }

    /** Get a point on the hyperplane. */
    VectorT getPoint() const { return dist * normal; }

    /** Get the unit normal vector of the hyperplane. */
    VectorT const & getNormal() const { return normal; }

    /**
     * Get the coefficients {a_i} of the hyperplane equation a_0 * x_0 + a_1 * x_1 + ... + a_{N - 1} * x_{N - 1} + a_N * 1 = 0.
     */
    Vector<N + 1, T> getEquation() const
    {
      Vector<N + 1, T> coeffs;

      for (intx i = 0; i < N; ++i)
        coeffs[i] = normal[i];

      coeffs[N] = -dist;

      return coeffs;
    }

    /** Flip the hyperplane so that the normal points the other way. */
    void flip()
    {
      normal = -normal;
      dist = -dist;
    }

    /** Get the (unsigned) distance of a given point from the hyperplane. */
    T distance(VectorT const & p) const
    {
      return std::abs(signedDistance(p));
    }

    /**
     * Get the signed distance of a given point from the hyperplane. This is positive if the point is on the side of the
     * hyperplane containing the normal, and negative if the point is on the other side.
     */
    T signedDistance(VectorT const & p) const
    {
      return p.dot(normal) - dist;
    }

    /** Get the square of the distance of the hyperplane from a given point. */
    T squaredDistance(VectorT const & p) const
    {
      return Math::square(signedDistance(p));
    }

    /** Get the point on the hyperplane closest to a given point. */
    VectorT closestPoint(VectorT const & p) const
    {
      return p - signedDistance(p) * normal;
    }

    /**
     * Check if the positive half space (the side of the hyperplane containing the normal) contains a given point. Returns true
     * if the point lies on the hyperplane.
     */
    bool positiveHalfSpaceContains(VectorT const & p) const
    {
      return signedDistance(p) >= 0;
    }

    /**
     * Check if the negative half space (the side of the hyperplane <b>not</b> containing the normal) contains a given point.
     * Returns true if the point lies on the hyperplane.
     */
    bool negativeHalfSpaceContains(VectorT const & p) const
    {
      return signedDistance(p) <= 0;
    }

    /** Reflect a point in the hyperplane. */
    VectorT reflect(VectorT const & p) const
    {
      return p - 2 * signedDistance(p) * normal;
    }

    bool rayIntersects(RayN<N, T> const & ray, T max_time = -1) const
    {
      return rayIntersectionTime(ray, max_time) >= 0;
    }

    T rayIntersectionTime(RayN<N, T> const & ray, T max_time = -1) const
    {
      T numer = dist - normal.dot(ray.getOrigin());
      T denom = normal.dot(ray.getDirection());
      if (std::abs(denom) < Math::eps(numer, denom))
        return -1;
      else
      {
        T t = numer / denom;
        return (max_time < 0 || t <= max_time) ? t : -1;
      }
    }

    RayIntersectionN<N, T> rayIntersection(RayN<N, T> const & ray, T max_time = -1) const
    {
      T t = rayIntersectionTime(ray, max_time);
      if (t >= 0)
        return RayIntersectionN<N, T>(t, &normal);
      else
        return RayIntersectionN<N, T>(-1);
    }

    /** Get a textual description of the hyperplane. */
    std::string toString() const
    {
      std::ostringstream oss;
      oss << "[N: " << Thea::toString(normal) << ", D: " << dist << ']';
      return oss.str();
    }

  protected:
    VectorT  normal;  ///< The unit normal vector of the hyperplane.
    T        dist;    ///< The signed distance of the hyperplane from the origin.

}; // class HyperplaneNBase

} // namespace Internal

/** A hyperplane ((N - 1)-flat) in N-dimensional space, where N is any <b>positive</b> (non-zero) integer and T is a field. */
template <int N, typename T = Real>
class /* THEA_API */ HyperplaneN : public Internal::HyperplaneNBase<N, T>
{
  public:
    /** Default constructor. */
    HyperplaneN() {}

}; // class HyperplaneN

/** Pipe a textual representation of a hyperplane to a <code>std::ostream</code>. */
template <int N, typename T>
std::ostream &
operator<<(std::ostream & os, HyperplaneN<N, T> const & plane)
{
  return os << plane.toString();
}

} // namespace Thea

#include "Hyperplane3.hpp"

#endif
