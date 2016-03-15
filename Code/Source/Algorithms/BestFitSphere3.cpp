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
#include "CentroidN.hpp"
#include "SVD.hpp"
#include "../Math.hpp"
#include "../Matrix.hpp"
#include "../TransposedMatrix.hpp"
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
    Vector3 centroid = CentroidN<Vector3, 3>::compute(points.begin(), points.end());
    double x = centroid[0];
    double y = centroid[1];
    double z = centroid[2];

    double r = 0;
    for (array_size_t i = 0; i < points.size(); ++i)
      r += (points[i] - centroid).length();

    r /= points.size();

    static double const THRESHOLD = 1.0e-10;
    double g_new = 100.0;
    double g_old = 1.0;

    Matrix<double> J((long)points.size(), 4);
    Matrix<double> D((long)points.size(), 1);
    Matrix<double> Inv(4, (long)points.size());

    while (std::fabs(g_new - g_old) > THRESHOLD)
    {
      g_old = g_new;

      for (array_size_t i = 0; i < points.size(); ++i)
      {
        double pX = points[i][0] - x;
        double pY = points[i][1] - y;
        double pZ = points[i][2] - z;
        double ri = std::sqrt(pX * pX + pY * pY + pZ * pZ);

        D((long)i, 0) = r - ri;
        J((long)i, 0) = -pX / ri;
        J((long)i, 1) = -pY / ri;
        J((long)i, 2) = -pZ / ri;
        J((long)i, 3) = -1;
      }

      SVD::pseudoInverse(J, Inv);  // 4 x n
      Matrix<double> X = Inv * D;  // 4 x 1

      x += X(0, 0);
      y += X(1, 0);
      z += X(2, 0);
      r += X(3, 0);

      TransposedMatrix<double> Jt(&J);  // 4 x n
      g_new = 0.0;
      for (int i = 0; i < 4; ++i)
      {
        for (long j = 0; j < Jt.numColumns(); ++j)
          g_new -= Jt.get(i, j) * D(j, 0);
      }
    }

    ball = Ball3(Vector3(x, y, z), (Real)r);
  }

  updated = true;
}

} // namespace Thea
} // namespace Algorithms
