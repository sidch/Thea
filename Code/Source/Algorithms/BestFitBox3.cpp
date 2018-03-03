//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2009, Siddhartha Chaudhuri/Stanford University
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

#include "BestFitBox3.hpp"
#include "CentroidN.hpp"
#include "LinearLeastSquares3.hpp"
#include "../Quat.hpp"

namespace Thea {
namespace Algorithms {

namespace BestFitBox3Internal {

//
// The algorithm (brute-force search of rotations in best-fit plane) is a reimplementation of John Ratcliff's code. The original
// license is:
//
// Copyright (c) 2009 by John W. Ratcliff mailto:jratcliffscarab@gmail.com
//
// Portions of this source has been released with the PhysXViewer application, as well as
// Rocket, CreateDynamics, ODF, and as a number of sample code snippets.
//
// If you find this code useful or you are feeling particularily generous I would
// ask that you please go to http://www.amillionpixels.us and make a donation
// to Troy DeMolay.
//
// DeMolay is a youth group for young men between the ages of 12 and 21.
// It teaches strong moral principles, as well as leadership skills and
// public speaking.  The donations page uses the 'pay for pixels' paradigm
// where, in this case, a pixel is only a single penny.  Donations can be
// made for as small as $4 or as high as a $100 block.  Each person who donates
// will get a link to their own site as well as acknowledgement on the
// donations blog located here http://www.amillionpixels.blogspot.com/
//
// If you wish to contact me you can use the following methods:
//
// Skype ID: jratcliff63367
// Yahoo: jratcliff63367
// AOL: jratcliff1961
// email: jratcliffscarab@gmail.com
//
//
// The MIT license:
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is furnished
// to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

static Quat
rotationArc(Vector3 const & v0, Vector3 const & v1)
{
  Vector3 cross = v0.cross(v1);
  Real d = v0.dot(v1);
  Real s = std::sqrt((1 + d) * 2);
  Real recip = 1.0f / s;

  Vector3 rc = recip * cross;
  return Quat(rc, 0.5f * s);
}

// Convert a plane equation to a rigid transform.
static void
planeToCFrame(Plane3 const & plane, Vector3 const & centroid, CoordinateFrame3 & cframe)
{
  Vector3 n;
  Real d;
  plane.getEquation(n.x(), n.y(), n.z(), d);

  Matrix3 rot;
  if (n.y() < -0.999)
  {
    rot = Matrix3(1,  0,  0,
                  0, -1,  0,
                  0,  0, -1);
  }
  else
  {
    Quat quat = rotationArc(Vector3::unitY(), n);
    rot = quat.toRotationMatrix();
  }

  cframe = CoordinateFrame3::_fromAffine(AffineTransform3(rot, centroid));
}

Matrix3
basisMatrix(Vector3 const & u, Vector3 const & v, Vector3 const & w)
{
  return Matrix3(u.x(), v.x(), w.x(),
                 u.y(), v.y(), w.y(),
                 u.z(), v.z(), w.z());
}

struct OBB
{
  Vector3 lo, hi;
  CoordinateFrame3 cframe;
  double volume;

  OBB() {}

  OBB(Vector3 const & lo_, Vector3 const & hi_, CoordinateFrame3 const & cframe_, double volume_)
  : lo(lo_), hi(hi_), cframe(cframe_), volume(volume_) {}
};

// Compute the OBB for a set of points relative to a transform matrix and see if it is smaller than the current best value. If
// so, return it.
static void
computeOBB(TheaArray<Vector3> const & points, CoordinateFrame3 const & cframe, OBB & current_best_result, bool overwrite)
{
  if (points.empty())
    return;

  Vector3 lo, hi;
  lo = hi = cframe.pointToObjectSpace(points[0]);

  Vector3 p;
  for (array_size_t i = 1; i < points.size(); ++i)
  {
    p = cframe.pointToObjectSpace(points[i]);
    lo = lo.min(p);
    hi = hi.max(p);
  }

  Vector3 e = hi - lo;
  double volume = e.x() * e.y() * e.z();
  if (overwrite || volume < current_best_result.volume)
  {
    Vector3 c  = 0.5f * (lo + hi);
    Vector3 he = 0.5f * e;
    current_best_result = OBB(c - he, c + he, cframe, volume);
  }
}

static void
computeBestFitOBB(TheaArray<Vector3> const & points, Box3 & result, bool has_up, Vector3 const & up)
{
  if (points.empty())
  {
    result = Box3();
    return;
  }
  else if (points.size() == 1)
  {
    result = Box3(AxisAlignedBox3(points[0]));
    return;
  }

  Vector3 centroid;
  CoordinateFrame3 cframe;

  if (has_up)
  {
    centroid = CentroidN<Vector3, 3>::compute(points.begin(), points.end());
    Vector3 u, v;
    up.createOrthonormalBasis(u, v);
    cframe = CoordinateFrame3::_fromAffine(AffineTransform3(basisMatrix(u, v, up), centroid));
  }
  else
  {
    Plane3 plane;
    LinearLeastSquares3<Vector3>::fitPlane(points.begin(), points.end(), plane, &centroid);
    planeToCFrame(plane, centroid, cframe);
  }

  OBB best_obb;
  computeOBB(points, cframe, best_obb, true);

  Matrix3 rot = Matrix3::rotationAxisAngle((has_up ? up : Vector3::unitY()), Math::degreesToRadians(10));
  for (int a = 10;  a < 180; a += 10)
  {
    cframe._setRotation(cframe.getRotation() * rot);
    computeOBB(points, cframe, best_obb, false);
  }

  result = Box3(AxisAlignedBox3(best_obb.lo, best_obb.hi), best_obb.cframe);
}

} // namespace BestFitBox3Internal

BestFitBox3::BestFitBox3()
: has_up(false), updated(true)
{}

void
BestFitBox3::addPoint(Vector3 const & point)
{
  points.push_back(point);
  updated = false;
}

void
BestFitBox3::clear()
{
  points.clear();
  updated = false;
}

void
BestFitBox3::releaseMemoryWithoutUpdate()
{
  points.clear();
}

Box3 const &
BestFitBox3::getBox() const
{
  update();
  return box;
}

void
BestFitBox3::update() const
{
  if (updated)
    return;

  BestFitBox3Internal::computeBestFitOBB(points, box, has_up, up);
  updated = true;
}

} // namespace Algorithms
} // namespace Thea
