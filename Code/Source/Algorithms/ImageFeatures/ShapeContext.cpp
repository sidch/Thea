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

#include "ShapeContext.hpp"
#include "../../Array.hpp"
#include "../../AxisAlignedBoxN.hpp"
#include "../../HyperplaneN.hpp"
#include "../../Image.hpp"
#include "../../Math.hpp"
#include "../../Vector2.hpp"

namespace Thea {
namespace Algorithms {
namespace ImageFeatures {

namespace ShapeContextInternal {

typedef AxisAlignedBoxN<2, Real> Rect;
typedef HyperplaneN<2, Real> Halfspace;

struct Sector
{
  Vector2 center;
  Real min_radius, max_radius, min_angle, max_angle;
  Vector2 corner[4];  // small arc min angle to max angle, large arc min angle to max angle
  Halfspace min_space, max_space;  // positive halfspace contains sector

  // min_angle, max_ang must be in range [0, 2 * pi]
  Sector(Vector2 const & center_, Real min_radius_, Real max_radius_, Real min_angle_, Real max_angle_)
  : center(center_), min_radius(min_radius_), max_radius(max_radius_), min_angle(min_angle_), max_angle(max_angle_)
  {
    Real cos_min = Math::fastCos(min_angle);
    Real sin_min = Math::fastSin(min_angle);

    Real cos_max = Math::fastCos(max_angle);
    Real sin_max = Math::fastSin(max_angle);

    corner[0] = center + Vector2(min_radius * cos_min, min_radius * sin_min);
    corner[1] = center + Vector2(min_radius * cos_max, min_radius * sin_max);

    corner[2] = center + Vector2(max_radius * cos_min, max_radius * sin_min);
    corner[3] = center + Vector2(max_radius * cos_max, max_radius * sin_max);

    Vector2 delta_min = corner[2] - center;
    min_space = Halfspace::fromPointAndNormal(center, Vector2(-delta_min.y(), delta_min.x()));

    Vector2 delta_max = corner[3] - center;
    max_space = Halfspace::fromPointAndNormal(center, Vector2(delta_max.y(), -delta_max.x()));
  }

  bool contains(Vector2 const & p) const
  {
    Vector2 delta = p - center;
    Real sqdist = delta.squaredLength();
    if (sqdist < min_radius * min_radius || sqdist > max_radius * max_radius)
      return false;

    Real ang = Math::fastArcTan2(delta.y(), delta.x());
    if (ang < 0) ang += Math::twoPi();

    return (ang >= min_angle && ang <= max_angle);
  }

  bool contains(Rect const & rect) const
  {
    if (rect.squaredDistance(center) < min_radius * min_radius)
      return false;

    return (contains(rect.getCorner(0))
         && contains(rect.getCorner(1))
         && contains(rect.getCorner(2))
         && contains(rect.getCorner(3)));
  }

  // A conservative intersection test -- can return true even when the shapes don't intersect
  bool maybeIntersects(Rect const & rect) const
  {
    // Any corner inside rectangle ==> intersect
    if (rect.contains(corner[0])
     || rect.contains(corner[1])
     || rect.contains(corner[2])
     || rect.contains(corner[3]))
      return true;

    // Rectangle outside ring ==> no intersect
    if (rect.squaredDistance(center) > max_radius * max_radius
     || rect.squaredMaxDistance(center) < min_radius * min_radius)
      return false;

    // Any rectangle corner inside sector ==> intersect
    if (contains(rect.getCorner(0))
     || contains(rect.getCorner(1))
     || contains(rect.getCorner(2))
     || contains(rect.getCorner(3)))
      return true;

    // Test against min_angle separating ray
    if (min_space.negativeHalfSpaceContains(rect.getCorner(0))
     && min_space.negativeHalfSpaceContains(rect.getCorner(1))
     && min_space.negativeHalfSpaceContains(rect.getCorner(2))
     && min_space.negativeHalfSpaceContains(rect.getCorner(3)))
      return false;

    // Test against max_angle separating ray
    if (max_space.negativeHalfSpaceContains(rect.getCorner(0))
     && max_space.negativeHalfSpaceContains(rect.getCorner(1))
     && max_space.negativeHalfSpaceContains(rect.getCorner(2))
     && max_space.negativeHalfSpaceContains(rect.getCorner(3)))
      return false;

    return true;
  }

}; // struct Sector

struct QuadTreeNode
{
  Rect rect;
  Real weight;
  TheaArray<array_size_t> children;

  void getBounds(int xmin, int ymin, int xmax, int ymax) const
  {
    xmin = (int)rect.getLow().x();
    ymin = (int)rect.getLow().y();

    xmax = (int)rect.getHigh().x();
    ymax = (int)rect.getHigh().y();
  }

}; // struct QuadTreeNode

struct QuadTree
{
  static int const MIN_NODE_WIDTH = 10;

  TheaArray<QuadTreeNode> nodes;
  Image owned_image;
  Image const * image;

  void clear()
  {
    nodes.clear();
  }

  bool init(Image const & image_)
  {
    clear();

    image = &image_;
    if (image_.getType() != Image::Type::LUMINANCE_8U)
    {
      if (!image_.convert(Image::Type::LUMINANCE_8U, owned_image))
        return false;

      image = &owned_image;
    }

    int w = image->getWidth();
    int h = image->getHeight();

    nodes.push_back(QuadTreeNode());
    nodes[0].rect.set(Vector2(0, 0), Vector2(w - 1, h - 1));
    buildTree(&nodes[0]);
  }

  void buildTree(QuadTreeNode * node)
  {
    int xmin, ymin, xmax, ymax;
    node->getBounds(xmin, ymin, xmax, ymax);

    int num_xsteps = (xmax - xmin >= MIN_NODE_WIDTH ? 2 : 1);
    int num_ysteps = (ymax - ymin >= MIN_NODE_WIDTH ? 2 : 1);

    int xmid = (num_xsteps > 1 ? (xmin + xmax) / 2 : xmax);
    int ymid = (num_ysteps > 1 ? (ymin + ymax) / 2 : ymax);

    for (int i = 0; i < num_xsteps; ++i)
    {
      int cell_xmin = (i == 0 ? xmin : xmid + 1);
      int cell_xmax = (i == 0 ? xmid : xmax);

      for (int j = 0; j < num_ysteps; ++j)
      {
        int cell_ymin = (i == 0 ? ymin : ymid + 1);
        int cell_ymax = (i == 0 ? ymid : ymax);

        QuadTreeNode child;
        child.rect.set(Vector2(cell_xmin, cell_ymin), Vector2(cell_xmax, cell_ymax));
        initWeight(&child);  // do this brute force here so we can aggressively prune recursion when weight --> 0

        if (child.weight > 1e-10)
        {
          nodes.push_back(child);
          node->children.push_back(nodes.size() - 1);
        }
      }
    }

    for (array_size_t i = 0; i < node->children.size(); ++i)
      buildTree(&nodes[node->children[i]]);
  }

  void initWeight(QuadTreeNode * node)
  {
    node->weight = 0;

    int xmin, ymin, xmax, ymax;
    node->getBounds(xmin, ymin, xmax, ymax);

    for (int i = ymin; i <= ymax; ++i)
    {
      uint8 const * scanline = (uint8 const *)image->getScanLine(i);

      for (int j = xmin; j <= xmax; ++j)
        node->weight += scanline[j];
    }

    node->weight /= 255.0;  // normalize each pixel's weight to [0, 1]
  }

  Real rangeWeight(Sector const & range)
  {
    return rangeWeight(&nodes[0], range);
  }

  Real rangeWeight(QuadTreeNode const * node, Sector const & range)
  {
    if (range.contains(node->rect))
      return node->weight;
    else if (!range.maybeIntersects(node->rect))
      return 0;

    Real weight = 0;

    if (node->children.empty())
    {
      int xmin, ymin, xmax, ymax;
      node->getBounds(xmin, ymin, xmax, ymax);

      for (int i = ymin; i <= ymax; ++i)
      {
        uint8 const * scanline = (uint8 const *)image->getScanLine(i);

        for (int j = xmin; j <= xmax; ++j)
        {
          if (range.contains(Vector2(i, j)))
            weight += scanline[j];
        }
      }

      weight /= 255.0;  // normalize each pixel's weight to [0, 1]
    }
    else
    {
      for (array_size_t i = 0; i < node->children.size(); ++i)
        weight += rangeWeight(&nodes[node->children[i]], range);
    }

    return weight;
  }

}; // QuadTree

} // namespace ShapeContextInternal

bool
ShapeContext::compute(Image const & image, long num_radial_bins, long num_polar_bins, bool ignore_empty_pixels, Real max_radius)
const
{
  return false;
}

} // namespace ImageFeatures
} // namespace Algorithms
} // namespace Thea
