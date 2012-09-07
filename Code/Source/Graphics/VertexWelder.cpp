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

#include "VertexWelder.hpp"
#include "../UnorderedMap.hpp"
#include <boost/functional/hash.hpp>
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
      long x, y, z;

      Long3(long x_, long y_, long z_) : x(x_), y(y_), z(z_) {}

      bool operator==(Long3 const & other) const { return x == other.x && y == other.y && z == other.z; }
    };

    struct Long3Hasher
    {
      std::size_t operator()(Long3 const & a) const
      {
        std::size_t s = boost::hash_value(a.x);
        boost::hash_combine(s, a.y);
        boost::hash_combine(s, a.z);
        return s;
      }
    };

    typedef TheaUnorderedMultiMap<Long3, WeldableVertex, Long3Hasher> VertexMap;

    Long3 toGrid(Vector3 const & pos) const
    {
      return Long3(static_cast<long>(std::floor(pos.x() / (double)weld_radius)),
                   static_cast<long>(std::floor(pos.y() / (double)weld_radius)),
                   static_cast<long>(std::floor(pos.z() / (double)weld_radius)));
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
              if ((ni->second.position - position).squaredLength() <= squared_weld_radius)
                return ni->second.vertex;
          }

      return NULL;
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
