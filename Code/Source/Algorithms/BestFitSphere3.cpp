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

#include "BestFitSphere3.hpp"
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_spheres_d_traits_3.h>
#include <algorithm>

namespace Thea {
namespace Algorithms {

BestFitSphere3::BestFitSphere3()
: ball(Vector3::zero(), 0), updated(true)
{
}

void
BestFitSphere3::addPoint(Vector3 const & point)
{
  points.push_back(point);
  updated = false;
}

void
BestFitSphere3::clear()
{
  points.clear();
  updated = false;
}

Real
BestFitSphere3::getRadius() const
{
  update();
  return ball.getRadius();
}

Real
BestFitSphere3::getDiameter() const
{
  update();
  return ball.getDiameter();
}

Vector3 const &
BestFitSphere3::getCenter() const
{
  update();
  return ball.getCenter();
}

Ball3 const &
BestFitSphere3::getBall() const
{
  update();
  return ball;
}

void
BestFitSphere3::update() const
{
  if (updated)
    return;

  if (points.empty())
  {
    ball = Ball3(Vector3::zero(), 0);
  }
  else if (points.size() == 1)
  {
    ball = Ball3(points[0], 0);
  }
  else
  {
    typedef CGAL::Cartesian<double>                              Kernel;
    typedef double                                               FT;
    typedef CGAL::Min_sphere_of_spheres_d_traits_3<Kernel, FT>   Traits;
    typedef Traits::Point                                        Point;
    typedef Traits::Sphere                                       Sphere;
    typedef CGAL::Min_sphere_of_spheres_d<Traits>                MSS;

    TheaArray<Sphere> p(points.size());
    for (array_size_t i = 0; i < points.size(); ++i)
    {
      Vector3 const & point = points[i];
      p[i] = Sphere(Point(point.x(), point.y(), point.z()), 0);
    }

    Traits traits;
    MSS mss(p.begin(), p.end(), traits);

    Real r = (Real)mss.radius();
    MSS::Cartesian_const_iterator ci = mss.center_cartesian_begin();
    Vector3 c;
    c[0] = (Real)(*ci);
    c[1] = (Real)(*(++ci));
    c[2] = (Real)(*(++ci));

    ball = Ball3(c, r);
  }

  updated = true;
}

} // namespace Thea
} // namespace Algorithms
