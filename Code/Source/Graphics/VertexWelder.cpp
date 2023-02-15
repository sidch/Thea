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

#include "VertexWelder.hpp"
#include "../Hash.hpp"
#include "../UnorderedMap.hpp"
#include <cmath>

namespace Thea {
namespace Graphics {

class VertexWelderImpl
{
  private:
    struct WeldableVertex
    {
      WeldableVertex(void * v, Vector3 const & p) : vertex(v), position(p) {}

      void * vertex;
      Vector3 position;
    };

    struct Long3
    {
      intx x, y, z;

      Long3(intx x_, intx y_, intx z_) : x(x_), y(y_), z(z_) {}

      bool operator==(Long3 const & other) const { return x == other.x && y == other.y && z == other.z; }
    };

    struct Long3Hasher
    {
      size_t operator()(Long3 const & a) const
      {
        size_t s = hashValue(a.x);
        hashCombine(s, a.y);
        hashCombine(s, a.z);
        return s;
      }
    };

    typedef UnorderedMultiMap<Long3, WeldableVertex, Long3Hasher> VertexMap;

    Long3 toGrid(Vector3 const & pos) const
    {
      return Long3(static_cast<intx>(std::floor(pos.x() / (double)weld_radius)),
                   static_cast<intx>(std::floor(pos.y() / (double)weld_radius)),
                   static_cast<intx>(std::floor(pos.z() / (double)weld_radius)));
    }

    Real weld_radius;
    Real squared_weld_radius;
    VertexMap vertex_map;

  public:
    VertexWelderImpl(Real weld_radius_) : weld_radius(weld_radius_), squared_weld_radius(weld_radius_ * weld_radius_) {}

    void addVertex(void * vertex, Vector3 const & position)
    {
      WeldableVertex value(vertex, position);
      Long3 key = toGrid(position);
      vertex_map.insert(VertexMap::value_type(key, value));
    }

    void * getVertex(Vector3 const & position) const
    {
      typedef std::pair<VertexMap::const_iterator, VertexMap::const_iterator> IteratorRange;
      Long3 base_key = toGrid(position);

      // Every vertex within the welding distance of the given point must be inside a neighbor of the grid cell containing the
      // point
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz)
          {
            Long3 key(base_key.x + dx, base_key.y + dy, base_key.z + dz);
            IteratorRange nbrs = vertex_map.equal_range(key);
            for (VertexMap::const_iterator ni = nbrs.first; ni != nbrs.second; ++ni)
              if ((ni->second.position - position).squaredNorm() <= squared_weld_radius)
                return ni->second.vertex;
          }

      return nullptr;
    }

}; // class VertexWelderImpl

VertexWelder::VertexWelder(Real weld_radius)
: impl(new VertexWelderImpl(weld_radius))
{
  alwaysAssertM(weld_radius >= 0, "VertexWelder: Weld radius cannot be negative");
}

VertexWelder::~VertexWelder()
{
  delete impl;
}

void
VertexWelder::addVertex(void * vertex, Vector3 const & position)
{
  impl->addVertex(vertex, position);
}

void *
VertexWelder::getVertex(Vector3 const & position) const
{
  return impl->getVertex(position);
}

} // namespace Graphics
} // namespace Thea
