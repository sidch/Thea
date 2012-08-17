//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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

#ifndef __Thea_Algorithms_KDTree3_hpp__
#define __Thea_Algorithms_KDTree3_hpp__

#include "../Common.hpp"
#include "../AffineTransform3.hpp"
#include "../Array.hpp"
#include "../AttributedObject.hpp"
#include "../Math.hpp"
#include "../Noncopyable.hpp"
#include "../Transformable.hpp"
#include "BoundedObjectTraits3.hpp"
#include "Filter.hpp"
#include "ProximityQueryStructure3.hpp"
#include "RangeQueryStructure3.hpp"
#include "RayQueryStructure3.hpp"
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>

namespace Thea {
namespace Algorithms {

/**
 * A kd-tree for a set of objects in 3-space. The requirements on the optional template parameters are:
 *
 * - <code>NodeAttributeT</code> must be default-constructible.
 * - <code>BoundedObjectTraitsT</code> must support obtaining bounding volumes for objects of type T.
 *
 * An affine transformation may be applied to the kd-tree. The tree does <em>not</em> need to be recomputed after the
 * transformation, though all operations may be somewhat slower (the precise overhead depends on how difficult it is to compute
 * distances, intersections etc after a transform). Normally, this means that the elements of type T should be affine
 * transformable via Transformer::transform().
 */
template < typename T,
           typename NodeAttributeT = NullAttribute,
           typename BoundedObjectTraitsT = BoundedObjectTraits3<T> >
class /* THEA_API */ KDTree3
: public RangeQueryStructure3<T>,
  public ProximityQueryStructure3,
  public RayQueryStructure3,
  public Transformable<AffineTransform3>,
  private Noncopyable
{
  private:
    typedef RangeQueryStructure3<T>          RangeQueryBaseT;
    typedef ProximityQueryStructure3         ProximityQueryBaseT;
    typedef RayQueryStructure3               RayQueryBaseT;
    typedef Transformable<AffineTransform3>  TransformableBaseT;

  public:
    typedef array_size_t ElementIndex;  ///< Index of an element in the kd-tree.
    using ProximityQueryBaseT::NeighborPair;

  private:
    /**
     * A memory pool for fast allocation of arrays. A pointer to memory allocated from the pool remains valid, and the pool
     * never reduces in size, until the pool is cleared, reinitialized, or destroyed.
     */
    template <typename U> class MemoryPool : private Noncopyable
    {
      private:
        /** A single buffer in the memory pool. */
        class Buffer : private Noncopyable
        {
          public:
            /** Constructor. */
            Buffer(array_size_t capacity_ = 0) : data(NULL), capacity(capacity_), current_end(0)
            {
              alwaysAssertM(capacity_ > 0, "KDTree3: Memory pool buffer capacity must be positive");
              if (capacity > 0)
              {
                data = new U[capacity];
                // THEA_CONSOLE << "Allocated data block " << data << " for buffer " << this;
              }
            }

            /** Destructor. */
            ~Buffer()
            {
              // THEA_CONSOLE << "Deleting data " << data << " of buffer " << this;
              delete [] data;
            }

            /** Clears the buffer without deallocating buffer memory. */
            void reset()
            {
              current_end = 0;
            }

            /**
             * Allocate a block of elements and return a pointer to the first allocated element, or null if the allocation
             * exceeded buffer capacity.
             */
            U * alloc(array_size_t num_elems)
            {
              // THEA_CONSOLE << "KDTree3: Allocating " << num_elems << " elements from buffer of capacity " << capacity;

              if (current_end + num_elems > capacity)
                return NULL;

              U * ret = &data[current_end];
              current_end += num_elems;
              return ret;
            }

            /** Deallocate a block of elements and return the number of elements successfully deallocated. */
            array_size_t free(array_size_t num_elems)
            {
              if (num_elems > current_end)
              {
                array_size_t num_freed = current_end;
                current_end = 0;
                return num_freed;
              }

              current_end -= num_elems;
              return num_elems;
            }

            /** Element access. */
            U const & operator[](array_size_t i) const { return data[i]; }

            /** Element access. */
            U & operator[](array_size_t i) { return data[i]; }

          private:
            U * data;
            array_size_t capacity;
            array_size_t current_end;

        }; // Buffer

        TheaArray<Buffer *> buffers;
        array_size_t buffer_capacity;
        int current_buffer;

        /** Get the current buffer. */
        Buffer & getCurrentBuffer() { return *buffers[(array_size_t)current_buffer]; }

        /** Get the next available buffer, creating it if necessary and making it the current one. */
        Buffer & getNextBuffer()
        {
          int next_buffer = current_buffer < 0 ? 0 : current_buffer + 1;
          if ((array_size_t)next_buffer >= buffers.size())
          {
            buffers.push_back(new Buffer(buffer_capacity));

            // THEA_CONSOLE << "KDTree3: Added buffer to memory pool " << this << ", current_buffer = " << current_buffer
            //              << ", next_buffer = " << next_buffer;
          }

          current_buffer = next_buffer;
          return *buffers[(array_size_t)current_buffer];
        }

      public:
        /** Constructor. */
        MemoryPool() : buffer_capacity(0), current_buffer(-1)
        {
          // THEA_CONSOLE << "KDTree3: Creating memory pool " << this;
        }

        /** Destructor. */
        ~MemoryPool()
        {
          clear(true);
          // THEA_CONSOLE << "KDTree3: Destroyed memory pool " << this;
        }

        /** Initialize the memory pool to hold buffers of a given capacity. Previous data in the pool is deallocated. */
        void init(array_size_t buffer_capacity_)
        {
          // THEA_CONSOLE << "KDTree3: Initializing memory pool " << this << " with buffer capacity " << buffer_capacity_
          //              << " elements";

          clear(true);

          buffer_capacity = buffer_capacity_;
        }

        /** Get the maximum number of elements of type T that a single buffer can hold. */
        array_size_t getBufferCapacity() const
        {
          return buffer_capacity;
        }

        /** Reset the memory pool, optionally deallocating and removing all buffers. */
        void clear(bool deallocate_all_memory = true)
        {
          // THEA_CONSOLE << "KDTree3: Clearing memory pool " << this;

          if (deallocate_all_memory)
          {
            for (array_size_t i = 0; i < buffers.size(); ++i)
              delete buffers[i];

            buffers.clear();
          }
          else
          {
            for (array_size_t i = 0; i < buffers.size(); ++i)
              buffers[i]->reset();
          }

          current_buffer = -1;
        }

        /** Allocate a block of elements from the pool and return a pointer to the first allocated element. */
        U * alloc(array_size_t num_elems)
        {
          alwaysAssertM(num_elems <= buffer_capacity, "KDTree3: A single memory pool allocation cannot exceed buffer capacity");

          if (current_buffer >= 0)
          {
            U * ret = getCurrentBuffer().alloc(num_elems);
            if (!ret)
              return getNextBuffer().alloc(num_elems);
            else
              return ret;
          }

          return getNextBuffer().alloc(num_elems);
        }

        /** Deallocate a block of elements from the pool. */
        void free(array_size_t num_elems)
        {
          long n = (long)num_elems;
          while (n > 0 && current_buffer >= 0)
          {
            array_size_t num_freed = getCurrentBuffer().free(num_elems);
            n -= (long)num_freed;

            if (n > 0)
              current_buffer--;
          }
        }

    }; // class MemoryPool

    /** A functor to add results of a range query to an array. */
    class RangeQueryFunctor
    {
      public:
        RangeQueryFunctor(TheaArray<T> & result_) : result(result_) {}
        bool operator()(long index, T & t) { result.push_back(t); return false; }

      private:
        TheaArray<T> & result;
    };

    /** A functor to add the indices of results of a range query to an array. */
    class RangeQueryIndicesFunctor
    {
      public:
        RangeQueryIndicesFunctor(TheaArray<long> & result_) : result(result_) {}
        bool operator()(long index, T & t) { result.push_back(index); return false; }

      private:
        TheaArray<long> & result;
    };

  public:
    THEA_DEF_POINTER_TYPES(KDTree3, shared_ptr, weak_ptr)

    typedef T                     Element;              ///< The type of elements in the kd-tree.
    typedef T                     value_type;           ///< The type of elements in the kd-tree (STL convention).
    typedef NodeAttributeT        NodeAttribute;        ///< Attributes attached to nodes.
    typedef BoundedObjectTraitsT  BoundedObjectTraits;  ///< Gives bounding volumes for elements.

    /** A node of the kd-tree. Only immutable objects of this class should be exposed by the external kd-tree interface. */
    class Node : public AttributedObject<NodeAttributeT>
    {
      private:
        int depth;
        AxisAlignedBox3 bounds;
        array_size_t num_elems;
        ElementIndex * elems;
        Node * lo;
        Node * hi;

        friend class KDTree3;

        void init(int depth_)
        {
          depth = depth_;
          bounds = AxisAlignedBox3();
          num_elems = 0;
          elems = NULL;
          lo = hi = NULL;
        }

      public:
        /** Iterator over immutable element indices. Dereferences to an array index. */
        typedef ElementIndex const * ElementIndexConstIterator;

        /** Constructor. */
        Node(int depth_ = 0) : depth(depth_), lo(NULL), hi(NULL) {}

        /** Get the depth of the node in the tree (the root is at depth 0). */
        int getDepth() const { return depth; }

        /** Get the bounding box of the node. */
        AxisAlignedBox3 const & getBounds() const { return bounds; }

        /**
         * Get the number of element indices stored at this node. This is <b>not</b> the number of elements within the node's
         * bounding box: in memory-saving mode, indices of all such elements are only held at the leaves of the subtree rooted
         * at this node.
         */
        long numElementIndices() const { return (long)num_elems; }

        /** Get an iterator to the first element index stored at the node. */
        ElementIndexConstIterator elementIndicesBegin() const { return elems; }

        /** Get an iterator to one past the last element index stored at the node. */
        ElementIndexConstIterator elementIndicesEnd() const { return elems + num_elems; }

        /** Get the child corresponding to the lower half of the range. */
        Node const * getLowChild() const { return lo; }

        /** Get the child containing the upper half of the range. */
        Node const * getHighChild() const { return hi; }

        /** Check if the node is a leaf (both children are null) are not. */
        bool isLeaf() const { return !(lo || hi); }

    }; // Node

  private:
    typedef MemoryPool<Node> NodePool;  ///< A pool for quickly allocating kd-tree nodes.
    typedef MemoryPool<ElementIndex> IndexPool;  ///< A pool for quickly allocating element indices.
    typedef TheaArray<T> ElementArray;  ///< An array of elements.

  public:
    /** Default constructor. */
    KDTree3() : root(NULL), num_elems(0), num_nodes(0), max_depth(0), max_elems_in_leaf(0) {}

    /**
     * Construct from a list of elements. InputIterator must dereference to type T.
     *
     * @param begin Points to the first element to be added.
     * @param end Points to one position beyond the last element to be added.
     * @param max_depth_ Maximum depth of the tree. The root is at depth zero. Use a negative argument to auto-select a suitable
     *   value.
     * @param max_elems_in_leaf_ Maximum number of elements in a leaf (unless the depth exceeds the maximum). Use a negative
     *   argument to auto-select a suitable value.
     * @param save_memory If true, element references at inner nodes of the tree are deleted to save memory. This could slow
     *   down range searches since every positive result will only be obtained at the leaves.
     */
    template <typename InputIterator>
    KDTree3(InputIterator begin, InputIterator end, int max_depth_ = -1, int max_elems_in_leaf_ = -1, bool save_memory = false)
    : root(NULL), num_elems(0), num_nodes(0), max_depth(0), max_elems_in_leaf(0), valid_bounds(true)
    {
      init(begin, end, max_elems_in_leaf_, max_depth_, save_memory, false /* no previous data to deallocate */);
    }

    /**
     * Construct from a list of elements. InputIterator must dereference to type T. Any previous data is discarded. If any
     * filters are active at this time, only those input elements that pass the filters will be retained in the tree.
     *
     * @param begin Points to the first element to be added.
     * @param end Points to one position beyond the last element to be added.
     * @param max_depth_ Maximum depth of the tree. The root is at depth zero. Use a negative argument to auto-select a suitable
     *   value.
     * @param max_elems_in_leaf_ Maximum number of elements in a leaf (unless the depth exceeds the maximum). Use a negative
     *   argument to auto-select a suitable value.
     * @param save_memory If true, element references at inner nodes of the tree are deleted to save memory. This could slow
     *   down range searches since every positive result will only be obtained at the leaves.
     * @param deallocate_previous_memory If true, all previous data held in internal memory pools is explicitly deallocated.
     *   Else, all such space is reused and overwritten when possible. If \a save_memory is true, or some filters are active,
     *   this flag may not be quite as effective since it's more likely that some space will be allocated/deallocated. Note that
     *   if this flag is set to false, the space used internally by the kd-tree will not decrease except in some special
     *   implementation-specific cases.
     */
    template <typename InputIterator>
    void init(InputIterator begin, InputIterator end, int max_depth_ = -1, int max_elems_in_leaf_ = -1,
              bool save_memory = false, bool deallocate_previous_memory = true)
    {
      clear(deallocate_previous_memory);

      if (deallocate_previous_memory)
      {
        for (InputIterator ii = begin; ii != end; ++ii, ++num_elems)
          if (elementPassesFilters(*ii))
            elems.push_back(*ii);
      }
      else
      {
        array_size_t max_new_elems = (array_size_t)std::distance(begin, end);
        bool resized = false;
        if (max_new_elems > elems.size())
        {
          if (filters.empty())
            elems.resize((array_size_t)std::ceil(1.2 * max_new_elems));  // add a little more space to avoid future reallocs
          else
            elems.clear();  // we don't know how many elements will pass the filter

          resized = true;
        }

        if (filters.empty())
        {
          std::copy(begin, end, elems.begin());
          num_elems = max_new_elems;
        }
        else
        {
          if (resized)
          {
            for (InputIterator ii = begin; ii != end; ++ii)
              if (elementPassesFilters(*ii))
              {
                elems.push_back(*ii);
                ++num_elems;
              }
          }
          else
          {
            typename ElementArray::iterator ei = elems.begin();
            for (InputIterator ii = begin; ii != end; ++ii)
              if (elementPassesFilters(*ii))
              {
                *(ei++) = *ii;
                ++num_elems;
              }
          }
        }
      }

      if (num_elems <= 0)
        return;

      static int const DEFAULT_MAX_ELEMS_IN_LEAF = 10;
      max_elems_in_leaf = max_elems_in_leaf_ < 0 ? DEFAULT_MAX_ELEMS_IN_LEAF : max_elems_in_leaf_;

      // The fraction of elements held by the larger node at each split is 0.5
      static double const SPLIT_FRACTION = 0.5;
      int est_depth = Math::kdtreeDepth(num_elems, max_elems_in_leaf, SPLIT_FRACTION);
      max_depth = max_depth_;
      if (max_depth < 0)
        max_depth = est_depth;
      else if (max_depth < est_depth)
        est_depth = max_depth;

      // THEA_CONSOLE << "KDTree3: max_depth = " << max_depth << ", est_depth = " << est_depth;

      // Each index is stored at most once at each level
      array_size_t BUFFER_SAFETY_MARGIN = 10;
      array_size_t index_buffer_capacity = num_elems + BUFFER_SAFETY_MARGIN;
      if (!save_memory)
        index_buffer_capacity *= (array_size_t)(1 + est_depth);  // reserve space for all levels at once

      if (deallocate_previous_memory || index_buffer_capacity > 1.3 * index_pool.getBufferCapacity())
      {
        // THEA_CONSOLE << "KDTree3: Resizing index pool: old buffer capacity = " << index_pool.getBufferCapacity()
        //              << ", new buffer capacity = " << index_buffer_capacity;
        index_pool.init(index_buffer_capacity);
      }

      // Assume a complete, balanced binary tree upto the estimated depth to guess the number of leaf nodes
      array_size_t node_buffer_capacity = (array_size_t)(1 << est_depth) + BUFFER_SAFETY_MARGIN;
      if (deallocate_previous_memory || node_buffer_capacity > 1.3 * node_pool.getBufferCapacity())
      {
        // THEA_CONSOLE << "KDTree3: Resizing node pool: old buffer capacity = " << node_pool.getBufferCapacity()
        //              << ", new buffer capacity = " << node_buffer_capacity;
        node_pool.init(node_buffer_capacity);
      }

      // Create the root node
      root = node_pool.alloc(1);
      root->init(0);
      num_nodes = 1;

      // THEA_CONSOLE << "Allocated root " << root << " from mempool " << &node_pool;

      root->num_elems = num_elems;
      root->elems = index_pool.alloc(root->num_elems);

      AxisAlignedBox3 elem_bounds;
      for (array_size_t i = 0; i < (array_size_t)num_elems; ++i)
      {
        root->elems[i] = i;

        BoundedObjectTraitsT::getBounds(elems[i], elem_bounds);
        root->bounds.merge(elem_bounds);
      }

      // Expand the bounding box slightly to handle numerical error
      root->bounds.scaleCentered(BOUNDS_EXPANSION_FACTOR);

      if (save_memory)
      {
        // Estimate the maximum number of indices that will need to be held in the scratch pool at any time during depth-first
        // traversal with earliest-possible deallocation. This is
        //
        //      #elements * (sum of series 1 + 1 + SPLIT_FRACTION + SPLIT_FRACTION^2 + ... + SPLIT_FRACTION^(est_depth - 1))
        //  <=  #elements * (1 + 1 / (1 - SPLIT_FRACTION))
        //
        array_size_t est_max_path_indices = (array_size_t)(num_elems * (1 + 1 / (1 - SPLIT_FRACTION)));
        // THEA_CONSOLE << "KDTree3: Estimated maximum number of indices on a single path = " << est_max_path_indices;

        // Create a temporary pool for scratch storage
        IndexPool tmp_index_pool;
        tmp_index_pool.init(est_max_path_indices + BUFFER_SAFETY_MARGIN);

        createTree(root, true, &tmp_index_pool, &index_pool);
      }
      else
        createTree(root, false, &index_pool, NULL);

      invalidateBounds();
    }

    /** Destructor. */
    ~KDTree3() { clear(true); }

    void setTransform(Transform const & trans_)
    {
      TransformableBaseT::setTransform(trans_);
      transform_inverse_transpose = trans_.getLinear().inverse().transpose();
      invalidateBounds();
    }

    /**
     * Clear the tree. If \a deallocate_all_memory is false, memory allocated in pools is held to be reused if possible by the
     * next init() operation.
     */
    virtual void clear(bool deallocate_all_memory = true)
    {
      num_elems = 0;
      if (deallocate_all_memory)
        elems.clear();

      node_pool.clear(deallocate_all_memory);
      index_pool.clear(deallocate_all_memory);

      root = NULL;

      invalidateBounds();
    }

    /** Check if the tree is empty. */
    bool isEmpty() const { return num_elems <= 0; }

    /** Get the number of elements in the tree. The elements themselves can be obtained with getElements(). */
    long numElements() const { return num_elems; }

    /** Get a pointer to an array of the elements in the tree. The number of elements can be obtained with numElements(). */
    T const * getElements() const { return &elems[0]; }

    /**
     * Get the node corresponding to the root of the kd-tree. This function is provided so that users can implement their own
     * tree traversal procedures without being restricted by the interface of RangeQueryStructure3.
     *
     * This function cannot be used to change the structure of the tree, or any value in it (unless <code>const_cast</code> is
     * used, which is not recommended).
     *
     * @note An empty tree has a null root.
     */
    Node const * getRoot() const { return root; }

    /** Get the number of nodes in the tree. */
    long numNodes() const { return num_nodes; }

    AxisAlignedBox3 const & getBounds() const
    {
      updateBounds();
      return bounds;
    }

    /**
     * Push an element filter onto the filter stack. Elements in the tree that are not passed by all filters currently on the
     * stack are ignored for all operations, including init().
     *
     * The filter must persist until it is popped off. Must be matched with popFilter().
     *
     * @see popFilter()
     */
    void pushFilter(Filter<T> * filter)
    {
      filters.push_back(filter);
    }

    /**
     * Pops the last pushed element filter off the filter stack. Must be matched with a preceding pushFilter().
     *
     * @see pushFilter()
     */
    void popFilter()
    {
      filters.pop_back();
    }

    template <typename MetricT, typename QueryT> double distance(QueryT const & query, double dist_bound = -1) const
    {
      double result = -1;
      if (closestElement<MetricT>(query, dist_bound, &result) >= 0)
        return result;
      else
        return -1;
    }

    template <typename MetricT, typename QueryT>
    long closestElement(QueryT const & query, double dist_bound = -1, double * dist = NULL, Vector3 * closest_point = NULL)
    const
    {
      NeighborPair pair = closestPair<MetricT>(query, dist_bound, (bool)closest_point);

      if (pair.isValid())
      {
        if (dist) *dist = MetricT::invertMonotoneApprox(pair.getMonotoneApproxDistance());
        if (closest_point) *closest_point = pair.getTargetPoint();
      }

      return pair.getTargetIndex();
    }

    // BoundedObjectTraits3<QueryT> must be defined.
    template <typename MetricT, typename QueryT>
    NeighborPair closestPair(QueryT const & query, double dist_bound = -1, bool get_closest_points = false) const
    {
      if (!root) return NeighborPair(-1);

      AxisAlignedBox3 query_bounds;
      BoundedObjectTraits3<QueryT>::getBounds(query, query_bounds);

      // Early pruning if the entire structure is too far away from the query
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = MetricT::monotoneApproxDistance(getBoundsWorldSpace(*root), query_bounds);
        if (lower_bound > mon_approx_dist_bound)
          return NeighborPair(-1);
      }

      NeighborPair pair(-1, -1, mon_approx_dist_bound);
      closestPair<MetricT>(root, query, query_bounds, pair, get_closest_points);

      return pair;
    }

    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSet>
    long kClosestPairs(QueryT const & query, BoundedNeighborPairSet & k_closest_pairs, double dist_bound = -1,
                       bool get_closest_points = false, bool clear_set = true, long use_as_query_index_and_swap = -1) const
    {
      if (clear_set) k_closest_pairs.clear();

      if (!root) return 0;

      AxisAlignedBox3 query_bounds;
      BoundedObjectTraits3<QueryT>::getBounds(query, query_bounds);

      // Early pruning if the entire structure is too far away from the query
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = MetricT::monotoneApproxDistance(getBoundsWorldSpace(*root), query_bounds);
        if (lower_bound > mon_approx_dist_bound)
          return 0;

        if (!k_closest_pairs.isInsertable(NeighborPair(0, 0, lower_bound)))
          return 0;
      }

      kClosestPairs<MetricT>(root, query, query_bounds, k_closest_pairs, dist_bound, get_closest_points,
                             use_as_query_index_and_swap);

      return k_closest_pairs.size();
    }

    template <typename IntersectionTesterT, typename RangeT>
    void rangeQuery(RangeT const & range, TheaArray<T> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root)
      {
        RangeQueryFunctor functor(result);
        const_cast<KDTree3 *>(this)->processRangeUntil<IntersectionTesterT>(root, range, &functor);
      }
    }

    template <typename IntersectionTesterT, typename RangeT>
    void rangeQueryIndices(RangeT const & range, TheaArray<long> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root)
      {
        RangeQueryIndicesFunctor functor(result);
        const_cast<KDTree3 *>(this)->processRangeUntil<IntersectionTesterT>(root, range, &functor);
      }
    }

    /**
     * Apply a functor to all elements in a range, until the functor returns true. See the base class documentation
     * (RangeQueryStructure3::processRangeUntil()) for more information.
     *
     * The RangeT class should support intersection queries with AxisAlignedBox3 and containment queries with Vector3 and
     * AxisAlignedBox3.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    bool processRangeUntil(RangeT const & range, FunctorT * functor)
    {
      if (!functor)  // no-op
        return false;

      return root ? processRangeUntil<IntersectionTesterT>(root, range, functor) : false;
    }

    template <typename RayIntersectionTesterT> bool rayIntersects(Ray3 const & ray, Real max_time = -1) const
    {
      return rayIntersectionTime<RayIntersectionTesterT>(ray, max_time) >= 0;
    }

    template <typename RayIntersectionTesterT> Real rayIntersectionTime(Ray3 const & ray, Real max_time = -1) const
    {
      if (root)
      {
        if (hasTransform())
        {
          Ray3 tr_ray = toObjectSpace(ray);
          if (root->bounds.rayIntersects(tr_ray, max_time))
            return rayIntersectionTime<RayIntersectionTesterT>(root, tr_ray, max_time);
        }
        else
        {
          if (root->bounds.rayIntersects(ray, max_time))
            return rayIntersectionTime<RayIntersectionTesterT>(root, ray, max_time);
        }
      }

      return -1;
    }

    template <typename RayIntersectionTesterT>
    RayStructureIntersection3 rayStructureIntersection(Ray3 const & ray, Real max_time = -1) const
    {
      if (root)
      {
        if (hasTransform())
        {
          Ray3 tr_ray = toObjectSpace(ray);
          if (root->bounds.rayIntersects(tr_ray, max_time))
          {
            RayStructureIntersection3 isec = rayStructureIntersection<RayIntersectionTesterT>(root, tr_ray, max_time);
            if (isec.isValid() && isec.hasNormal())
              isec.setNormal(normalToWorldSpace(isec.getNormal()));

            return isec;
          }
        }
        else
        {
          if (root->bounds.rayIntersects(ray, max_time))
            return rayStructureIntersection<RayIntersectionTesterT>(root, ray, max_time);
        }
      }

      return RayStructureIntersection3(-1);
    }

  private:
    /** Comparator for sorting elements along an axis. */
    struct ObjectLess
    {
      int coord;
      KDTree3 const * tree;

      /** Constructor. Axis 0 = X, 1 = Y, 2 = Z. */
      ObjectLess(int coord_, KDTree3 const * tree_) : coord(coord_), tree(tree_) {}

      /** Less-than operator, along the specified axis. */
      bool operator()(ElementIndex a, ElementIndex b)
      {
        // Compare object min coords
        return BoundedObjectTraitsT::getLow(tree->elems[a], coord)
             < BoundedObjectTraitsT::getLow(tree->elems[b], coord);
      }
    };

    // Allow the comparator unrestricted access to the kd-tree.
    friend struct ObjectLess;

    /** A stack of element filters. */
    typedef TheaArray<Filter<T> *> FilterStack;

    void moveIndicesToLeafPool(Node * leaf, IndexPool * main_index_pool, IndexPool * leaf_index_pool)
    {
      if (leaf)
      {
        ElementIndex * leaf_indices_start = leaf_index_pool->alloc(leaf->num_elems);
        std::memcpy(leaf_indices_start, leaf->elems, leaf->num_elems * sizeof(ElementIndex));
        leaf->elems = leaf_indices_start;
        main_index_pool->free(leaf->num_elems);
      }
    }

    /** Recursively construct the tree. */
    void createTree(Node * start, bool save_memory, IndexPool * main_index_pool, IndexPool * leaf_index_pool)
    {
      // Assume the start node is fully constructed at this stage.
      //
      // If we are in memory-saving mode, then we assume the node's indices are last in the main index pool (but not yet in the
      // leaf pool, since we don't yet know if this node will turn out to be a leaf). In this case, after the function finishes,
      // the node's indices will be deallocated from the main index pool (and possibly moved to the leaf pool).

      if (!start || start->depth >= max_depth || (long)start->num_elems <= max_elems_in_leaf)
      {
        if (save_memory)
          moveIndicesToLeafPool(start, main_index_pool, leaf_index_pool);

        return;
      }

      // Find a splitting plane
#define THEA_KDTREE3_SPLIT_LONGEST
#ifdef THEA_KDTREE3_SPLIT_LONGEST
      Vector3 ext = start->bounds.getExtent();  // split longest dimension
      int coord = (ext.x() > ext.y() ? (ext.x() > ext.z() ? 0 : 2) : (ext.y() > ext.z() ? 1 : 2));
#else
      int coord = (int)(start->depth % 3);
#endif

      // Split elements into lower and upper halves
      array_size_t mid = start->num_elems / 2;
      std::nth_element(start->elems, start->elems + mid, start->elems + start->num_elems, ObjectLess(coord, this));

      // Create child nodes
      start->lo = node_pool.alloc(1);
      start->lo->init(start->depth + 1);
      num_nodes++;

      start->hi = node_pool.alloc(1);
      start->hi->init(start->depth + 1);
      num_nodes++;

      // THEA_CONSOLE << "num_nodes = " << num_nodes;

      // Allocate element arrays for the children
      start->lo->elems = main_index_pool->alloc(start->num_elems - mid);
      start->hi->elems = main_index_pool->alloc(mid);

      // Add first half of array (elems less than median) to low child
      AxisAlignedBox3 elem_bounds;
      bool lo_first = true;
      for (ElementIndex i = 0; i < start->num_elems - mid; ++i)
      {
        ElementIndex index = start->elems[i];
        BoundedObjectTraitsT::getBounds(elems[index], elem_bounds);

        start->lo->elems[start->lo->num_elems++] = index;

        if (lo_first)
        {
          start->lo->bounds = elem_bounds;
          lo_first = false;
        }
        else
          start->lo->bounds.merge(elem_bounds);
      }

      // Add second half of array (elems greater than median) to high child
      bool hi_first = true;
      for (ElementIndex i = start->num_elems - mid; i < start->num_elems; ++i)
      {
        ElementIndex index = start->elems[i];
        BoundedObjectTraitsT::getBounds(elems[index], elem_bounds);

        start->hi->elems[start->hi->num_elems++] = index;

        if (hi_first)
        {
          start->hi->bounds = elem_bounds;
          hi_first = false;
        }
        else
          start->hi->bounds.merge(elem_bounds);
      }

      // Expand the bounding boxes slightly to handle numerical error
      start->lo->bounds.scaleCentered(BOUNDS_EXPANSION_FACTOR);
      start->hi->bounds.scaleCentered(BOUNDS_EXPANSION_FACTOR);

      // Recurse on the high child first, since its indices are at the end of the main index pool and can be freed first if
      // necessary
      createTree(start->hi, save_memory, main_index_pool, leaf_index_pool);

      // Recurse on the low child next, if we are in memory-saving mode its indices are now the last valid entries in the main
      // index pool
      createTree(start->lo, save_memory, main_index_pool, leaf_index_pool);

      // If we are in memory-saving mode, deallocate the indices stored at this node, which are currently the last entries in
      // the main index pool
      if (save_memory)
      {
        main_index_pool->free(start->num_elems);
        start->num_elems = 0;
        start->elems = NULL;
      }
    }

    /** Mark that the bounding box requires an update. */
    void invalidateBounds()
    {
      valid_bounds = false;
    }

    /** Recompute the bounding box if it has been invalidated. */
    void updateBounds() const
    {
      if (valid_bounds) return;

      if (root)
        bounds = hasTransform() ? root->bounds.transformAndBound(getTransform()) : root->bounds;
      else
        bounds = AxisAlignedBox3();

      valid_bounds = true;
    }

    /** Check if an element passes all filters currently on the stack. */
    bool elementPassesFilters(T const & elem) const
    {
      if (filters.empty()) return true;  // early exit

      for (typename FilterStack::const_iterator fi = filters.begin(); fi != filters.end(); ++fi)
        if (!(*fi)->allows(elem))
          return false;

      return true;
    }

    /** Get a bounding box for a node, in world space. */
    AxisAlignedBox3 getBoundsWorldSpace(Node const & node) const
    {
      return hasTransform() ? node.bounds.transformAndBound(getTransform()) : node.bounds;
    }

    /**
     * Recursively look for the closest pair of points between two elements. Only pairs separated by less than the current
     * minimum distance will be considered.
     */
    template <typename MetricT, typename QueryT>
    void closestPair(Node const * start, QueryT const & query, AxisAlignedBox3 const & query_bounds, NeighborPair & pair,
                     bool get_closest_points) const
    {
      if (!start->lo)  // leaf
        closestPairLeaf<MetricT>(start, query, pair, get_closest_points);
      else  // not leaf
      {
        // Figure out which child is closer (optimize for point queries?)
        Node const * n[2] = { start->lo, start->hi };
        double mad[2] = { MetricT::monotoneApproxDistance(getBoundsWorldSpace(*n[0]), query_bounds),
                          MetricT::monotoneApproxDistance(getBoundsWorldSpace(*n[1]), query_bounds) };

        if (mad[1] > mad[0])
        {
          std::swap(n[0], n[1]);
          std::swap(mad[0], mad[1]);
        }

        for (int i = 0; i < 2; ++i)
          if (pair.getMonotoneApproxDistance() < 0 || mad[i] <= pair.getMonotoneApproxDistance())
            closestPair<MetricT>(n[i], query, query_bounds, pair, get_closest_points);
      }
    }

    /**
     * Search the elements in a leaf node for the one closest to another element, when the latter is a proximity query
     * structure.
     */
    template <typename MetricT, typename QueryT>
    void closestPairLeaf(
      Node const * leaf,
      QueryT const & query,
      NeighborPair & pair,
      bool get_closest_points,
      typename boost::enable_if<boost::is_base_of<ProximityQueryBaseT, QueryT>, void>::type * dummy = NULL) const
    {
      for (array_size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        NeighborPair swapped;
        if (hasTransform())
          swapped = query.closestPair<MetricT>(makeTransformedObject(&elem, &getTransform()), pair.getMonotoneApproxDistance(),
                                               get_closest_points);
        else
          swapped = query.closestPair<MetricT>(elem, pair.getMonotoneApproxDistance(), get_closest_points);

        if (swapped.isValid())
        {
          pair = swapped.swapped();
          pair.setTargetIndex((long)index);
        }
      }
    }

    /**
     * Search the elements in a leaf node for the one closest to another element, when the latter is NOT a proximity query
     * structure.
     */
    template <typename MetricT, typename QueryT>
    void closestPairLeaf(
      Node const * leaf,
      QueryT const & query,
      NeighborPair & pair,
      bool get_closest_points,
      typename boost::disable_if<boost::is_base_of<ProximityQueryBaseT, QueryT>, void>::type * dummy = NULL) const
    {
      Vector3 qp, tp;
      double mad;

      for (array_size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        if (hasTransform())
          mad = MetricT::closestPoints(makeTransformedObject(&elem, &getTransform()), query, tp, qp);
        else
          mad = MetricT::closestPoints(elem, query, tp, qp);

        if (pair.getMonotoneApproxDistance() < 0 || mad <= pair.getMonotoneApproxDistance())
          pair = NeighborPair(0, (long)index, mad, qp, tp);
      }
    }

    /**
     * Recursively look for the k closest elements to a query object. Only elements at less than the specified maximum distance
     * will be considered.
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSet>
    void kClosestPairs(Node const * start, QueryT const & query, AxisAlignedBox3 const & query_bounds,
                       BoundedNeighborPairSet & k_closest_pairs, double dist_bound, bool get_closest_points,
                       long use_as_query_index_and_swap)
    const
    {
      if (!start->lo)  // leaf
        kClosestPairsLeaf<MetricT>(start, query, k_closest_pairs, dist_bound, get_closest_points, use_as_query_index_and_swap);
      else  // not leaf
      {
        // Figure out which child is closer (optimize for point queries?)
        Node const * n[2] = { start->lo, start->hi };
        double d[2] = { MetricT::monotoneApproxDistance(getBoundsWorldSpace(*n[0]), query_bounds),
                        MetricT::monotoneApproxDistance(getBoundsWorldSpace(*n[1]), query_bounds) };

        if (d[1] > d[0])
        {
          std::swap(n[0], n[1]);
          std::swap(d[0], d[1]);
        }

        double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);

        for (int i = 0; i < 2; ++i)
          if ((mon_approx_dist_bound < 0 || d[i] <= mon_approx_dist_bound)
            && k_closest_pairs.isInsertable(NeighborPair(0, 0, d[i])))
          {
            kClosestPairs<MetricT>(n[i], query, query_bounds, k_closest_pairs, dist_bound, get_closest_points,
                                   use_as_query_index_and_swap);
          }
      }
    }

    /**
     * Search the elements in a leaf node for the k nearest neighbors of an object, when the latter is a proximity query
     * structure.
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSet>
    void kClosestPairsLeaf(
      Node const * leaf,
      QueryT const & query,
      BoundedNeighborPairSet & k_closest_pairs,
      double dist_bound,
      bool get_closest_points,
      long use_as_query_index_and_swap,
      typename boost::enable_if<boost::is_base_of<ProximityQueryBaseT, QueryT>, void>::type * dummy = NULL) const
    {
      for (array_size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        if (hasTransform())
          query.kClosestPairs<MetricT>(makeTransformedObject(&elem, &getTransform()), k_closest_pairs, dist_bound,
                                       get_closest_points, false, (long)index);
        else
          query.kClosestPairs<MetricT>(elem, k_closest_pairs, dist_bound, get_closest_points, false, (long)index);
      }
    }

    /**
     * Search the elements in a leaf node for the one closest to another element, when the latter is NOT a proximity query
     * structure.
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSet>
    void kClosestPairsLeaf(
      Node const * leaf,
      QueryT const & query,
      BoundedNeighborPairSet & k_closest_pairs,
      double dist_bound,
      bool get_closest_points,
      long use_as_query_index_and_swap,
      typename boost::disable_if<boost::is_base_of<ProximityQueryBaseT, QueryT>, void>::type * dummy = NULL) const
    {
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);

      std::equal_to<NeighborPair> eq_comp;
      Vector3 qp, tp;
      double mad;

      for (array_size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        // Check if the element is already in the set of neighbors or not
        NeighborPair pair = (use_as_query_index_and_swap >= 0 ? NeighborPair((long)index, use_as_query_index_and_swap)
                                                              : NeighborPair(0, (long)index));
        if (k_closest_pairs.contains(pair, eq_comp))  // already found
          continue;

        if (hasTransform())
          mad = MetricT::closestPoints(makeTransformedObject(&elem, &getTransform()), query, tp, qp);
        else
          mad = MetricT::closestPoints(elem, query, tp, qp);

        if (mon_approx_dist_bound < 0 || mad <= mon_approx_dist_bound)
        {
          pair.setMonotoneApproxDistance(mad);

          if (get_closest_points)
          {
            if (use_as_query_index_and_swap >= 0)
            {
              pair.setQueryPoint(tp);
              pair.setTargetPoint(qp);
            }
            else
            {
              pair.setQueryPoint(qp);
              pair.setTargetPoint(tp);
            }
          }

          k_closest_pairs.insert(pair);
        }
      }
    }

    /**
     * Apply a functor to all elements of a subtree within a range, stopping when the functor returns true on any point. The
     * RangeT class should support containment queries with AxisAlignedBox3.
     *
     * @return True if the functor evaluated to true on any point in the range, else false.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    bool processRangeUntil(Node const * start, RangeT const & range, FunctorT * functor)
    {
      // Early exit if the range and node are disjoint
      AxisAlignedBox3 tr_start_bounds = getBoundsWorldSpace(*start);
      if (!IntersectionTesterT::intersects(range, tr_start_bounds))
        return false;

      // If the entire node is contained in the range AND there are element references at this node (so it's either a leaf or we
      // have not saved memory by flushing references at internal nodes), then process all these elems.
      if (start->num_elems > 0 && range.contains(tr_start_bounds))
      {
        for (array_size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          if ((*functor)(static_cast<long>(index), elem))
            return true;
        }
      }
      else if (!start->lo)  // leaf
      {
        for (array_size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          bool intersects = hasTransform()
                          ? IntersectionTesterT::intersects(makeTransformedObject(&elem, &getTransform()), range)
                          : IntersectionTesterT::intersects(elem, range);
          if (intersects)
            if ((*functor)(static_cast<long>(index), elem))
              return true;
        }
      }
      else  // not leaf
      {
        if (processRangeUntil<IntersectionTesterT>(start->lo, range, functor)) return true;
        if (processRangeUntil<IntersectionTesterT>(start->hi, range, functor)) return true;
      }

      return false;
    }

    /** Transform a ray to local/object space. */
    Ray3 toObjectSpace(Ray3 const & ray) const
    {
      AffineTransform3 inv_trans = getTransform().inverse();  // TODO: cache this
      return Ray3(inv_trans * ray.getOrigin(), inv_trans.getLinear() * ray.getDirection());
    }

    /** Transform a normal to world space. */
    Vector3 normalToWorldSpace(Vector3 const & n) const
    {
      return (transform_inverse_transpose * n).fastUnit();  // the slower unit(), perhaps?
    }

    /**
     * Check if the ray intersection time \a new_time represents a closer, or equally close, valid hit than the previous best
     * time \a old_time.
     */
    static bool improvedRayTime(Real new_time, Real old_time)
    {
      return (new_time >= 0 && (old_time < 0 || new_time <= old_time));
    }

    /** Get the time taken for a ray to hit the nearest object in a node, in the forward direction. */
    template <typename RayIntersectionTesterT>
    Real rayIntersectionTime(Node const * start, Ray3 const & ray, Real max_time) const
    {
      if (!start->lo)  // leaf
      {
        Real best_time = max_time;
        bool found = false;
        for (array_size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element const & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          Real time = RayIntersectionTesterT::rayIntersectionTime(ray, elem, best_time);
          if (improvedRayTime(time, best_time))
          {
            best_time = time;
            found = true;
          }
        }

        return found ? best_time : -1;
      }
      else  // not leaf
      {
        // Figure out which child will be hit first
        Node const * n[2] = { start->lo, start->hi };
        Real t[2] = { n[0]->bounds.rayIntersectionTime(ray, max_time),
                      n[1]->bounds.rayIntersectionTime(ray, max_time) };

        if (t[0] < 0 && t[1] < 0)
          return -1;

        if (improvedRayTime(t[1], t[0]))
        {
          std::swap(n[0], n[1]);
          std::swap(t[0], t[1]);
        }

        Real best_time = max_time;
        bool found = false;
        for (int i = 0; i < 2; ++i)
        {
          if (improvedRayTime(t[i], best_time))
          {
            Real time = rayIntersectionTime<RayIntersectionTesterT>(n[i], ray, best_time);
            if (improvedRayTime(time, best_time))
            {
              best_time = time;
              found = true;
            }
          }
        }

        return found ? best_time : -1;
      }
    }

    /** Get the nearest intersection of a ray with a node in the forward direction. */
    template <typename RayIntersectionTesterT>
    RayStructureIntersection3 rayStructureIntersection(Node const * start, Ray3 const & ray, Real max_time) const
    {
      if (!start->lo)  // leaf
      {
        RayStructureIntersection3 best_isec(max_time);
        bool found = false;
        for (array_size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element const & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          RayIntersection3 isec = RayIntersectionTesterT::rayIntersection(ray, elem, best_isec.getTime());
          if (improvedRayTime(isec.getTime(), best_isec.getTime()))
          {
            best_isec = RayStructureIntersection3(isec, (long)index);
            found = true;
          }
        }

        return found ? best_isec : RayStructureIntersection3(-1);
      }
      else  // not leaf
      {
        // Figure out which child will be hit first
        Node const * n[2] = { start->lo, start->hi };
        Real t[2] = { n[0]->bounds.rayIntersectionTime(ray, max_time),
                      n[1]->bounds.rayIntersectionTime(ray, max_time) };

        if (t[0] < 0 && t[1] < 0)
          return -1;

        if (improvedRayTime(t[1], t[0]))
        {
          std::swap(n[0], n[1]);
          std::swap(t[0], t[1]);
        }

        RayStructureIntersection3 best_isec(max_time);
        bool found = false;
        for (int i = 0; i < 2; ++i)
        {
          if (improvedRayTime(t[i], best_isec.getTime()))
          {
            RayStructureIntersection3 isec = rayStructureIntersection<RayIntersectionTesterT>(n[i], ray, best_isec.getTime());
            if (improvedRayTime(isec.getTime(), best_isec.getTime()))
            {
              best_isec = isec;
              found = true;
            }
          }
        }

        return found ? best_isec : RayStructureIntersection3(-1);
      }
    }

    Node * root;

    long num_elems;  // elems.size() doesn't tell us how many elements there are, it's just the capacity of the elems array
    ElementArray elems;  // elems.size() is *not* the number of elements in the tree!!!

    long num_nodes;
    NodePool node_pool;

    IndexPool index_pool;

    int max_depth;
    int max_elems_in_leaf;

    Matrix3 transform_inverse_transpose;

    FilterStack filters;

    mutable bool valid_bounds;
    mutable AxisAlignedBox3 bounds;

    static Real const BOUNDS_EXPANSION_FACTOR;

}; // class KDTree3

// Static variables
template <typename T, typename N, typename B> Real const KDTree3<T, N, B>::BOUNDS_EXPANSION_FACTOR = 1.05f;

} // namespace Algorithms
} // namespace Thea

#endif
