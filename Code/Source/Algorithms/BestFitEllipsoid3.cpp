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
// First version: 2009
//
//============================================================================

#include "BestFitEllipsoid3.hpp"

#if THEA_ENABLE_CGAL
// sprintf in boost::lexical_cast triggers a deprecation warning
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
#    include <CGAL/Cartesian.h>
#    include <CGAL/MP_Float.h>
#    include <CGAL/Approximate_min_ellipsoid_d.h>
#    include <CGAL/Approximate_min_ellipsoid_d_traits_3.h>
#  pragma clang diagnostic pop
#  include <algorithm>
#endif

namespace Thea {
namespace Algorithms {

BestFitEllipsoid3::BestFitEllipsoid3(Real eps_)
: eps(eps_), center(Vector3::Zero()), updated(true)
{
  axis[0] = axis[1] = axis[2] = Vector3::Zero();
}

void
BestFitEllipsoid3::addPoint(Vector3 const & point)
{
  points.push_back(point);
  updated = false;
}

void
BestFitEllipsoid3::clear()
{
  points.clear();
  updated = false;
}

void
BestFitEllipsoid3::getAxes(Vector3 & axis0, Vector3 & axis1, Vector3 & axis2) const
{
  update();

  axis0 = axis[0];
  axis1 = axis[1];
  axis2 = axis[2];
}

Vector3 const &
BestFitEllipsoid3::getCenter() const
{
  update();
  return center;
}

Box3 const &
BestFitEllipsoid3::getOrientedBoundingBox() const
{
  update();
  return obb;
}

void
BestFitEllipsoid3::update() const
{
  if (updated)
    return;

  if (points.empty())
    center = axis[0] = axis[1] = axis[2] = Vector3::Zero();
  else if (points.size() == 1)
  {
    center = points[0];
    axis[0] = axis[1] = axis[2] = Vector3::Zero();
  }
  else
  {
#if THEA_ENABLE_CGAL

    typedef CGAL::Cartesian<double>                                 Kernel;
    typedef CGAL::MP_Float                                          ET;
    typedef CGAL::Approximate_min_ellipsoid_d_traits_3<Kernel, ET>  Traits;
    typedef Traits::Point                                           Point;
    typedef CGAL::Approximate_min_ellipsoid_d<Traits>               AME;

    Array<Point> p(points.size());
    for (size_t i = 0; i < points.size(); ++i)
    {
      Vector3 const & point = points[i];
      p[i] = Point(point.x(), point.y(), point.z());
    }

    Traits traits;
    AME ame(eps, p.begin(), p.end(), traits);

    if (!ame.is_full_dimensional())
      axis[0] = axis[1] = axis[2] = Vector3::Zero();  // FIXME
    else
    {
      int coord = 0;
      for (AME::Center_coordinate_iterator c_it = ame.center_cartesian_begin(); coord < 3; ++c_it, ++coord)
        center[coord] = *c_it;

      AME::Axes_lengths_iterator l_it = ame.axes_lengths_begin();
      Real lengths[3];
      for (int i = 0; i < 3; ++i, ++l_it)
      {
        int coord = 0;
        for (AME::Axes_direction_coordinate_iterator d_it = ame.axis_direction_cartesian_begin(i); coord < 3; ++d_it, ++coord)
          axis[i][coord] = (Real)(*d_it);

        lengths[i] = (Real)(*l_it);
      }

      if (lengths[0] < lengths[2])  // CGAL bug? Lengths returned in ascending order although spec says otherwise
      {
        std::swap(lengths[0], lengths[2]);
        std::swap(axis[0], axis[2]);
      }

      // Compute oriented bounding box before we scale the unit axes by the dimensions of the ellipsoid
      CoordinateFrame3 cframe(RigidTransform3::_fromAffine(
                                  AffineTransform3((Matrix3() << axis[0], axis[1], axis[2]).finished())));
      Vector3 min, max;
      for (size_t i = 0; i < points.size(); ++i)
      {
        Vector3 op = cframe.pointToObjectSpace(points[i]);

        if (i == 0)
          min = max = op;
        else
        {
          min = op.cwiseMin(min);
          max = op.cwiseMax(max);
        }
      }

      Vector3 obb_center = cframe.pointToWorldSpace((Real)0.5 * (min + max));
      cframe.setTranslation(obb_center);

      Vector3 half_extent = (Real)0.5 * (max - min);
      obb = Box3(AxisAlignedBox3(-half_extent, half_extent), cframe);

      // Now scale the axes to the dimensions of the ellipsoid
      axis[0] *= lengths[0];
      axis[1] *= lengths[1];
      axis[2] *= lengths[2];

      // Make sure the axes, in order, form a right-handed system
      if ((axis[0].cross(axis[1])).dot(axis[2]) < 0)
        axis[2] = -axis[2];
    }

#else  // THEA_ENABLE_CGAL

    throw Error("BestFitEllipsoid3: Requires CGAL for 2 or more points");

#endif // THEA_ENABLE_CGAL
  }

  updated = true;
}

} // namespace Algorithms
} // namespace Thea
