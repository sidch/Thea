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
// First version: 2016
//
//============================================================================

#ifndef __Thea_UnionFind_hpp__
#define __Thea_UnionFind_hpp__

#include "Common.hpp"
#include "UnorderedMap.hpp"

namespace Thea {

/**
 * Union-find data structure. Original code by Kartik Kukreja, https://github.com/kartikkukreja/ . The data structure will save
 * a copy of each object to support getObjectId(), so T should be easily copyable (a small/POD class, or a pointer).
 */
template <typename T = intx>
class UnionFind
{
  public:
    /**
     * Create a union-find data structure for \a n objects, which will be assumed to have sequential IDs from 0 to \a n - 1. If
     * you use this constructor, note that getObjectId() will return the expected results (operating as an identity function)
     * only if <code>T = intx</code>. Else, avoid calling getObjectId().
     */
    UnionFind(intx n)
    : has_objects(false)
    {
      init(n);
    }

    /**
     * Create a union-find data structure for objects in the range [\a begin, \a end). The objects will be assigned sequential
     * IDs 0, 1, 2, ..., n - 1 in the same order as the input. Note that the data structure will save a copy of each object to
     * support getObjectId(), so T should be easily copyable (a small/POD class, or a pointer).
     */
    template <typename InputIterator> UnionFind(InputIterator begin, InputIterator end)
    : has_objects(true)
    {
      intx n = 0;
      for (InputIterator iter = begin; iter != end; ++iter, ++n)
        objects[*iter] = n;

      init(n);
    }

    /** Destructor. */
    ~UnionFind()
    {
      delete [] id;
      delete [] sz;
    }

    /**
     * Get the ID of an object, or a negative value if it is not present in the data structure. Note that if you used the
     * UnionFind(intx) constructor, this function will give the expected results (operating as an identity function) only if
     * <code>T = intx</code>.
     */
    intx getObjectId(T const & obj) const
    {
      typename ObjectMap::const_iterator existing = objects.find(obj);
      return (existing == objects.end() ? -1 : existing->second);
    }

    /** Return the ID of the representative element of the set corresponding to the object with ID \a p. */
    intx find(intx p) const
    {
      intx root = p;

      while (root != id[root])
        root = id[root];

      while (p != root)
      {
        intx newp = id[p];
        id[p] = root;
        p = newp;
      }

      return root;
    }

    /** Replace sets containing \a x and \a y with their union. */
    void merge(intx x, intx y)
    {
      intx i = find(x);
      intx j = find(y);

      if (i == j)
        return;

      // Make smaller root point to larger one
      if (sz[i] < sz[j])
      {
        id[i] = j;
        sz[j] += sz[i];
      }
      else
      {
        id[j] = i;
        sz[i] += sz[j];
      }

      count--;
    }

    /** Are objects \a x and \a y in the same set? */
    bool sameSet(intx x, intx y) const
    {
      return find(x) == find(y);
    }

    /** Get the number of disjoint sets. */
    intx numSets() const
    {
      return count;
    }

    /** Get the size of the set containing object \a p. */
    intx sizeOfSet(intx p) const
    {
      return sz[find(p)];
    }

  private:
    /** Initialize buffers for a set of \a n objects. */
    void init(intx n)
    {
      count = n;
      id = new intx[n];
      sz = new intx[n];

      for (intx i = 0; i < n; ++i)
      {
        id[i] = i;
        sz[i] = 1;
      }
    }

    typedef UnorderedMap<T, intx> ObjectMap;  ///< Map from objects to their IDs.

    bool has_objects;
    ObjectMap objects;
    intx * id;
    intx * sz;
    intx count;

}; // class UnionFind

// Specialization for T = intx
template <>
inline intx
UnionFind<intx>::getObjectId(intx const & obj) const
{
  if (has_objects)
  {
    ObjectMap::const_iterator existing = objects.find(obj);
    return (existing == objects.end() ? -1 : existing->second);
  }
  else
    return obj;
}

} // namespace Thea

#endif
