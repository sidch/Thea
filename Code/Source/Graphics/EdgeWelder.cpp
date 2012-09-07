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

#include "EdgeWelder.hpp"
#include "../UnorderedMap.hpp"
#include <boost/functional/hash.hpp>
#include <cmath>
#include <utility>

namespace Thea {
namespace Graphics {

class EdgeWelderImpl
{
  private:
    struct WeldableEdge
    {
      WeldableEdge(void * edge_, Vector3 const & e0_, Vector3 const & e1_) : edge(edge_), e0(e0_), e1(e1_) {}

      void * edge;
      Vector3 e0, e1;
    };

    struct Long3
    {
      long x, y, z;

      Long3(long x_, long y_, long z_) : x(x_), y(y_), z(z_) {}

      bool operator==(Long3 const & other) const { return x == other.x && y == other.y && z == other.z; }
    };

    typedef std::pair<Long3, Long3> Long3Pair;

    struct Long3PairHasher
    {
      std::size_t operator()(Long3Pair const & a) const
      {
        std::size_t s = boost::hash_value(a.first.x);
        boost::hash_combine(s, a.first.y);
        boost::hash_combine(s, a.first.z);
        boost::hash_combine(s, a.second.x);
        boost::hash_combine(s, a.second.y);
        boost::hash_combine(s, a.second.z);
        return s;
      }
    };

    typedef TheaUnorderedMultiMap<Long3Pair, WeldableEdge, Long3PairHasher> EdgeMap;

    Long3 toGrid(Vector3 const & pos) const
    {
      return Long3(static_cast<long>(std::floor(pos.x() / (double)weld_radius)),
                   static_cast<long>(std::floor(pos.y() / (double)weld_radius)),
                   static_cast<long>(std::floor(pos.z() / (double)weld_radius)));
    }

    Real weld_radius;
    Real squared_weld_radius;
    EdgeMap edge_map;

  public:
    EdgeWelderImpl(Real weld_radius_) : weld_radius(weld_radius_), squared_weld_radius(weld_radius_ * weld_radius_) {}

    void addEdge(void * edge, Vector3 const & e0, Vector3 const & e1)
    {
      WeldableEdge value(edge, e0, e1);
      Long3Pair key(toGrid(e0), toGrid(e1));
      edge_map.insert(EdgeMap::value_type(key, value));
    }

    void * getUndirectedEdge(Vector3 const & e0, Vector3 const & e1) const
    {
      void * edge = getDirectedEdge(e0, e1);
      if (edge)
        return edge;
      else
        return getDirectedEdge(e1, e0);
    }

    void * getDirectedEdge(Vector3 const & e0, Vector3 const & e1) const
    {
      typedef std::pair<EdgeMap::const_iterator, EdgeMap::const_iterator> IteratorRange;
      Long3Pair base_key(toGrid(e0), toGrid(e1));

      // Every edge within welding distance of the given edge must be inside a neighbor of the grid cell containing the edge
      for (int dx0 = -1; dx0 <= 1; ++dx0)
        for (int dy0 = -1; dy0 <= 1; ++dy0)
          for (int dz0 = -1; dz0 <= 1; ++dz0)
            for (int dx1 = -1; dx1 <= 1; ++dx1)
              for (int dy1 = -1; dy1 <= 1; ++dy1)
                for (int dz1 = -1; dz1 <= 1; ++dz1)
                {
                  Long3Pair key(Long3(base_key.first.x + dx0, base_key.first.y + dy0, base_key.first.z + dz0),
                                Long3(base_key.second.x + dx1, base_key.second.y + dy1, base_key.second.z + dz1));
                  IteratorRange nbrs = edge_map.equal_range(key);
                  for (EdgeMap::const_iterator ni = nbrs.first; ni != nbrs.second; ++ni)
                  {
                    if ((ni->second.e0 - e0).squaredLength() <= squared_weld_radius
                     && (ni->second.e1 - e1).squaredLength() <= squared_weld_radius)
                      return ni->second.edge;
                  }
                }

      return NULL;
    }

}; // class EdgeWelderImpl

EdgeWelder::EdgeWelder(Real weld_radius)
: impl(new EdgeWelderImpl(weld_radius))
{
  alwaysAssertM(weld_radius >= 0, "EdgeWelder: Weld radius cannot be negative");
}

EdgeWelder::~EdgeWelder()
{
  delete impl;
}

void
EdgeWelder::addEdge(void * edge, Vector3 const & e0, Vector3 const & e1)
{
  impl->addEdge(edge, e0, e1);
}

void *
EdgeWelder::getDirectedEdge(Vector3 const & e0, Vector3 const & e1) const
{
  return impl->getDirectedEdge(e0, e1);
}

void *
EdgeWelder::getUndirectedEdge(Vector3 const & e0, Vector3 const & e1) const
{
  return impl->getUndirectedEdge(e0, e1);
}

} // namespace Graphics
} // namespace Thea
