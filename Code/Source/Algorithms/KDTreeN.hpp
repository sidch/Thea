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

#ifndef __Thea_Algorithms_KDTreeN_hpp__
#define __Thea_Algorithms_KDTreeN_hpp__

#include "../Common.hpp"
#include "../AffineTransformN.hpp"
#include "../Array.hpp"
#include "../AttributedObject.hpp"
#include "../Math.hpp"
#include "../Noncopyable.hpp"
#include "../Random.hpp"
#include "../Transformable.hpp"
#include "BoundedObjectTraitsN.hpp"
#include "Filter.hpp"
#include "ProximityQueryStructureN.hpp"
#include "RangeQueryStructure.hpp"
#include "RayQueryStructureN.hpp"
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>

namespace Thea {
namespace Algorithms {

namespace KDTreeNInternal {

/** A point sample drawn from a kd-tree element, used for accelerating nearest neighbor queries. */
template <long N, typename ScalarT>
struct ElementSample
{
  VectorN<N, ScalarT> position;
  void const * element;  // this cannot be pointer-to-T, else we get recursive instantiation overflow

  /** Default constructor. */
  ElementSample() {}

  /** Initialize a sample at point \a p drawn from an element \a e. */
  ElementSample(VectorN<N, ScalarT> const & p, void const * e) : position(p), element(e) {}

}; // struct ElementSample

/**
 * A filter that passes or rejects a sample point, depending on whether a base filter passes or rejects the point's parent
 * element.
 */
template <typename T, long N, typename ScalarT>
class SampleFilter : public Filter< ElementSample<N, ScalarT> >
{
  public:
    /** Constructor. */
    SampleFilter(Filter<T> * base_filter_) : base_filter(base_filter_) {}

    bool allows(ElementSample<N, ScalarT> const & sample) const
    {
      // Assume null pointer case won't happen by construction
      return base_filter->allows(*static_cast<T const *>(sample.element));
    }

  private:
    Filter<T> * base_filter;  ///< The underlying element filter.

}; // struct SampleFilter

} // namespace KDTreeNInternal

template <long N, typename ScalarT>
class IsPointN< KDTreeNInternal::ElementSample<N, ScalarT>, N >
{
  public:
    static bool const value = true;
};

template <long N, typename ScalarT>
class PointTraitsN< KDTreeNInternal::ElementSample<N, ScalarT>, N, ScalarT >
{
  public:
    typedef VectorN<N, ScalarT> VectorT;
    static VectorT getPosition(KDTreeNInternal::ElementSample<N, ScalarT> const & sample) { return sample.position; }
};

/**
 * A kd-tree for a set of objects in N-space. The requirements on the optional template parameters are:
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
           long N,
           typename ScalarT = Real,
           typename NodeAttributeT = NullAttribute,
           typename BoundedObjectTraitsT = BoundedObjectTraitsN<T, N, ScalarT> >
class /* THEA_API */ KDTreeN
: public RangeQueryStructure<T>,
  public ProximityQueryStructureN<N, ScalarT>,
  public RayQueryStructureN<N, ScalarT>,
  public Transformable< AffineTransformN<N, ScalarT> >,
  private Noncopyable
{
  private:
    typedef RangeQueryStructure<T>                         RangeQueryBaseT;
    typedef ProximityQueryStructureN<N, ScalarT>           ProximityQueryBaseT;
    typedef RayQueryStructureN<N, ScalarT>                 RayQueryBaseT;
    typedef Transformable< AffineTransformN<N, ScalarT> >  TransformableBaseT;

  public:
    typedef array_size_t ElementIndex;  ///< Index of an element in the kd-tree.
    typedef typename ProximityQueryBaseT::NeighborPair NeighborPair;  ///< A pair of neighboring elements.
    typedef typename TransformableBaseT::Transform Transform;  ///< Transform applied to the kd-tree.

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
              alwaysAssertM(capacity_ > 0, "KDTreeN: Memory pool buffer capacity must be positive");
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
              // THEA_CONSOLE << "KDTreeN: Allocating " << num_elems << " elements from buffer of capacity " << capacity;

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
        long current_buffer;

        /** Get the current buffer. */
        Buffer & getCurrentBuffer() { return *buffers[(array_size_t)current_buffer]; }

        /** Get the next available buffer, creating it if necessary and making it the current one. */
        Buffer & getNextBuffer()
        {
          long next_buffer = current_buffer < 0 ? 0 : current_buffer + 1;
          if ((array_size_t)next_buffer >= buffers.size())
          {
            buffers.push_back(new Buffer(buffer_capacity));

            // THEA_CONSOLE << "KDTreeN: Added buffer to memory pool " << this << ", current_buffer = " << current_buffer
            //              << ", next_buffer = " << next_buffer;
          }

          current_buffer = next_buffer;
          return *buffers[(array_size_t)current_buffer];
        }

      public:
        /** Constructor. */
        MemoryPool() : buffer_capacity(0), current_buffer(-1)
        {
          // THEA_CONSOLE << "KDTreeN: Creating memory pool " << this;
        }

        /** Destructor. */
        ~MemoryPool()
        {
          clear(true);
          // THEA_CONSOLE << "KDTreeN: Destroyed memory pool " << this;
        }

        /** Initialize the memory pool to hold buffers of a given capacity. Previous data in the pool is deallocated. */
        void init(array_size_t buffer_capacity_)
        {
          // THEA_CONSOLE << "KDTreeN: Initializing memory pool " << this << " with buffer capacity " << buffer_capacity_
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
          // THEA_CONSOLE << "KDTreeN: Clearing memory pool " << this;

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
          alwaysAssertM(num_elems <= buffer_capacity, "KDTreeN: A single memory pool allocation cannot exceed buffer capacity");

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
    THEA_DEF_POINTER_TYPES(KDTreeN, shared_ptr, weak_ptr)

    typedef T                     Element;              ///< Type of elements in the kd-tree.
    typedef T                     value_type;           ///< Type of elements in the kd-tree (STL convention).
    typedef NodeAttributeT        NodeAttribute;        ///< Attributes attached to nodes.
    typedef BoundedObjectTraitsT  BoundedObjectTraits;  ///< Gives bounding volumes for elements.

    typedef typename ProximityQueryBaseT::VectorT              VectorT;                    ///< Vector in N-space.
    typedef AxisAlignedBoxN<N, ScalarT>                        AxisAlignedBoxT;            ///< Axis-aligned box in N-space.
    typedef typename RayQueryBaseT::RayT                       RayT;                       ///< Ray in N-space.
    typedef typename RayQueryBaseT::RayStructureIntersectionT  RayStructureIntersectionT;  /**< Ray intersection structure in
                                                                                                N-space. */

    typedef KDTreeNInternal::ElementSample<N, ScalarT> ElementSample;  /**< A point sample drawn from a kd-tree element, used
                                                                            for accelerating nearest neighbor queries. */
    typedef KDTreeN<ElementSample, N, ScalarT> NearestNeighborAccelerationStructure;  /**< Structure to speed up nearest
                                                                                           neighbor queries. */

    /** A node of the kd-tree. Only immutable objects of this class should be exposed by the external kd-tree interface. */
    class Node : public AttributedObject<NodeAttributeT>
    {
      private:
        long depth;
        AxisAlignedBoxT bounds;
        array_size_t num_elems;
        ElementIndex * elems;
        Node * lo;
        Node * hi;

        friend class KDTreeN;

        void init(long depth_)
        {
          depth = depth_;
          bounds = AxisAlignedBoxT();
          num_elems = 0;
          elems = NULL;
          lo = hi = NULL;
        }

      public:
        /** Iterator over immutable element indices. Dereferences to an array index. */
        typedef ElementIndex const * ElementIndexConstIterator;

        /** Constructor. */
        Node(long depth_ = 0) : depth(depth_), lo(NULL), hi(NULL) {}

        /** Get the depth of the node in the tree (the root is at depth 0). */
        long getDepth() const { return depth; }

        /** Get the bounding box of the node. */
        AxisAlignedBoxT const & getBounds() const { return bounds; }

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
    typedef KDTreeNInternal::SampleFilter<T, N, ScalarT> SampleFilter;  ///< Filter for samples, wrapping a filter for elements.

  public:
    /** Default constructor. */
    KDTreeN()
    : root(NULL), num_elems(0), num_nodes(0), max_depth(0), max_elems_in_leaf(0), accelerate_nn_queries(false),
      valid_acceleration_structure(false), acceleration_structure(NULL), valid_bounds(true)
    {}

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
    KDTreeN(InputIterator begin, InputIterator end, long max_depth_ = -1, long max_elems_in_leaf_ = -1,
            bool save_memory = false)
    : root(NULL), num_elems(0), num_nodes(0), max_depth(0), max_elems_in_leaf(0), accelerate_nn_queries(false),
      valid_acceleration_structure(false), acceleration_structure(NULL), valid_bounds(true)
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
    void init(InputIterator begin, InputIterator end, long max_depth_ = -1, long max_elems_in_leaf_ = -1,
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
          num_elems = (long)max_new_elems;
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

      static long const DEFAULT_MAX_ELEMS_IN_LEAF = 10;
      max_elems_in_leaf = max_elems_in_leaf_ < 0 ? DEFAULT_MAX_ELEMS_IN_LEAF : max_elems_in_leaf_;

      // The fraction of elements held by the larger node at each split is 0.5
      static double const SPLIT_FRACTION = 0.5;
      long est_depth = Math::binaryTreeDepth(num_elems, max_elems_in_leaf, SPLIT_FRACTION);
      max_depth = max_depth_;
      if (max_depth < 0)
        max_depth = est_depth;
      else if (max_depth < est_depth)
        est_depth = max_depth;

      // THEA_CONSOLE << "KDTreeN: max_depth = " << max_depth << ", est_depth = " << est_depth;

      // Each index is stored at most once at each level
      array_size_t BUFFER_SAFETY_MARGIN = 10;
      array_size_t index_buffer_capacity = num_elems + BUFFER_SAFETY_MARGIN;
      if (!save_memory)
        index_buffer_capacity *= (array_size_t)(1 + est_depth);  // reserve space for all levels at once

      if (deallocate_previous_memory || index_buffer_capacity > 1.3 * index_pool.getBufferCapacity())
      {
        // THEA_CONSOLE << "KDTreeN: Resizing index pool: old buffer capacity = " << index_pool.getBufferCapacity()
        //              << ", new buffer capacity = " << index_buffer_capacity;
        index_pool.init(index_buffer_capacity);
      }

      // Assume a complete, balanced binary tree upto the estimated depth to guess the number of leaf nodes
      array_size_t node_buffer_capacity = (array_size_t)(1 << est_depth) + BUFFER_SAFETY_MARGIN;
      if (deallocate_previous_memory || node_buffer_capacity > 1.3 * node_pool.getBufferCapacity())
      {
        // THEA_CONSOLE << "KDTreeN: Resizing node pool: old buffer capacity = " << node_pool.getBufferCapacity()
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

      AxisAlignedBoxT elem_bounds;
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
        // THEA_CONSOLE << "KDTreeN: Estimated maximum number of indices on a single path = " << est_max_path_indices;

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
    ~KDTreeN() { clear(true); }

    /**
     * Enable acceleration of nearest neighbor queries with an auxiliary structure on a sparse set of points.
     *
     * @see disableNearestNeighborAcceleration()
     */
    void enableNearestNeighborAcceleration(long num_acceleration_samples_ = -1)
    {
      accelerate_nn_queries = true;

      if (valid_acceleration_structure && num_acceleration_samples != num_acceleration_samples_)
        valid_acceleration_structure = false;

      num_acceleration_samples = num_acceleration_samples_;
    }

    /** Disable acceleration of nearest neighbor queries with an auxiliary structure on a sparse set of points off. */
    void disableNearestNeighborAcceleration(bool deallocate_memory = true)
    {
      accelerate_nn_queries = false;
      clearAccelerationStructure(deallocate_memory);
    }

    /** Check if nearest neighbor queries are accelerated by an auxiliary structure. */
    bool hasNearestNeighborAcceleration() const
    {
      return accelerate_nn_queries;
    }

    /** Get the auxiliary structure to accelerate nearest neighbor queries, if available. */
    template <typename MetricT> NearestNeighborAccelerationStructure const * getNearestNeighborAccelerationStructure() const
    {
      if (hasNearestNeighborAcceleration())
      {
        buildAccelerationStructure<MetricT>();
        return acceleration_structure;
      }
      else
        return NULL;
    }

    void setTransform(Transform const & trans_)
    {
      TransformableBaseT::setTransform(trans_);
      transform_inverse_transpose = trans_.getLinear().inverse().transpose();
      invalidateBounds();

      if (valid_acceleration_structure)
        acceleration_structure->setTransform(trans_);
    }

    void clearTransform()
    {
      TransformableBaseT::clearTransform();
      invalidateBounds();

      if (valid_acceleration_structure)
        acceleration_structure->clearTransform();
    }

    /**
     * Clear the tree. If \a deallocate_all_memory is false, memory allocated in pools is held to be reused if possible by the
     * next init() operation.
     */
    virtual void clear(bool deallocate_all_memory = true)
    {
      clearAccelerationStructure(deallocate_all_memory);

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

    /** Get a bounding box for all the objects in the tree. */
    AxisAlignedBoxT const & getBounds() const
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
      alwaysAssertM(filter, "KDTreeN: Filter must be non-null");

      filters.push_back(filter);

      if (accelerate_nn_queries && valid_acceleration_structure)
      {
        sample_filters.push_back(SampleFilter(filter));
        acceleration_structure->pushFilter(&sample_filters.back());
      }
    }

    /**
     * Pops the last pushed element filter off the filter stack. Must be matched with a preceding pushFilter().
     *
     * @see pushFilter()
     */
    void popFilter()
    {
      filters.pop_back();

      if (accelerate_nn_queries && valid_acceleration_structure)
      {
        acceleration_structure->popFilter();
        sample_filters.pop_back();
      }
    }

    /** Get the minimum distance between this structure and a query object. */
    template <typename MetricT, typename QueryT> double distance(QueryT const & query, double dist_bound = -1) const
    {
      double result = -1;
      if (closestElement<MetricT>(query, dist_bound, &result) >= 0)
        return result;
      else
        return -1;
    }

    /**
     * Get the closest element in this structure to a query object, within a specified distance bound.
     *
     * @param query Query object.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param dist The distance to the query object is placed here. Ignored if null.
     * @param closest_point The coordinates of the closest point are placed here. Ignored if null.
     *
     * @return A non-negative handle to the closest element, if one was found, else a negative number.
     */
    template <typename MetricT, typename QueryT>
    long closestElement(QueryT const & query, double dist_bound = -1, double * dist = NULL, VectorT * closest_point = NULL)
    const
    {
      NeighborPair pair = closestPair<MetricT>(query, dist_bound, closest_point != NULL);

      if (pair.isValid())
      {
        if (dist) *dist = MetricT::invertMonotoneApprox(pair.getMonotoneApproxDistance());
        if (closest_point) *closest_point = pair.getTargetPoint();
      }

      return pair.getTargetIndex();
    }

    /**
     * Get the closest pair of elements between this structure and another structure, whose separation is less than a specified
     * upper bound.
     *
     * @param query Query object. BoundedObjectTraitsN<QueryT, N, ScalarT> must be defined.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param get_closest_points If true, the coordinates of the closest pair of points on the respective elements is computed
     *   and stored in the returned structure.
     *
     * @return Non-negative handles to the closest pair of elements in their respective objects, if such a pair was found. Else
     *   returns a pair of negative numbers.
     */
    template <typename MetricT, typename QueryT>
    NeighborPair closestPair(QueryT const & query, double dist_bound = -1, bool get_closest_points = false) const
    {
      if (!root) return NeighborPair(-1);

      AxisAlignedBoxT query_bounds;
      BoundedObjectTraitsN<QueryT, N, ScalarT>::getBounds(query, query_bounds);

      // Early pruning if the entire structure is too far away from the query
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = MetricT::template monotoneApproxDistance<N, ScalarT>(getBoundsWorldSpace(*root), query_bounds);
        if (lower_bound > mon_approx_dist_bound)
          return NeighborPair(-1);
      }

      // If acceleration is enabled, set an upper limit to the distance to the nearest object
      double accel_bound = accelerationBound<MetricT>(query, dist_bound);
      if (accel_bound >= 0)
      {
        double fudge = 0.001 * getBoundsWorldSpace(*root).getExtent().fastLength();
        mon_approx_dist_bound = MetricT::computeMonotoneApprox(accel_bound + fudge);
      }

      NeighborPair pair(-1, -1, mon_approx_dist_bound);
      closestPair<MetricT>(root, query, query_bounds, pair, get_closest_points);

      return pair;
    }

    /**
     * Get the k elements closest to a query object. The returned elements are placed in a set of bounded size (k). The template
     * type BoundedNeighborPairSetT should typically be BoundedSortedArray<NeighborPair> or BoundedSortedArrayN<k, NeighborPair>
     * if only a few neighbors are requested. BoundedObjectTraitsN<QueryT, N, ScalarT> must be defined.
     *
     * @param query Query object. BoundedObjectTraitsN<QueryT, N, ScalarT> must be defined.
     * @param k_closest_pairs The k (or fewer) nearest neighbors are placed here.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param get_closest_points If true, the coordinates of the closest pair of points on each pair of neighboring elements is
     *   computed and stored in the returned pairs.
     * @param clear_set If true (default), this function discards prior data in \a k_closest_pairs. This is chiefly for internal
     *   use and the default value of true should normally be left as is.
     * @param use_as_query_index_and_swap If non-negative, the supplied index is used as the index of the query object (instead
     *   of the default 0), following which query and target indices/points are swapped in the returned pairs of neighbors. This
     *   is chiefly for internal use and the default value of -1 should normally be left as is.
     *
     * @return The number of neighbors found (i.e. the size of \a k_closest_pairs).
     *
     * @note k-closest pairs <b>cannot</b> be accelerated by the auxiliary structure created by
     *   enableNearestNeighborAcceleration().
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSetT>
    long kClosestPairs(QueryT const & query, BoundedNeighborPairSetT & k_closest_pairs, double dist_bound = -1,
                       bool get_closest_points = false, bool clear_set = true, long use_as_query_index_and_swap = -1) const
    {
      if (clear_set) k_closest_pairs.clear();

      if (!root) return 0;

      AxisAlignedBoxT query_bounds;
      BoundedObjectTraitsN<QueryT, N, ScalarT>::getBounds(query, query_bounds);

      // Early pruning if the entire structure is too far away from the query
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = MetricT::template monotoneApproxDistance<N, ScalarT>(getBoundsWorldSpace(*root), query_bounds);
        if (lower_bound > mon_approx_dist_bound)
          return 0;

        if (!k_closest_pairs.isInsertable(NeighborPair(0, 0, lower_bound)))
          return 0;
      }

      kClosestPairs<MetricT>(root, query, query_bounds, k_closest_pairs, dist_bound, get_closest_points,
                             use_as_query_index_and_swap);

      return k_closest_pairs.size();
    }

    /**
     * Get all objects intersecting a range.
     *
     * @param range The range to search in.
     * @param result The objects intersecting the range are stored here.
     * @param discard_prior_results If true, the contents of \a results are cleared before the range query proceeds. If false,
     *   the previous results are retained and new objects are appended to the array (this is useful for range queries over a
     *   union of simpler ranges).
     */
    template <typename IntersectionTesterT, typename RangeT>
    void rangeQuery(RangeT const & range, TheaArray<T> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root)
      {
        RangeQueryFunctor functor(result);
        const_cast<KDTreeN *>(this)->processRangeUntil<IntersectionTesterT>(root, range, &functor);
      }
    }

    /**
     * Get the indices of all objects intersecting a range.
     *
     * @param range The range to search in.
     * @param result The indices of objects intersecting the range are stored here.
     * @param discard_prior_results If true, the contents of \a results are cleared before the range query proceeds. If false,
     *   the previous results are retained and indices of new objects are appended to the array (this is useful for range
     *   queries over a union of simpler ranges).
     */
    template <typename IntersectionTesterT, typename RangeT>
    void rangeQueryIndices(RangeT const & range, TheaArray<long> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root)
      {
        RangeQueryIndicesFunctor functor(result);
        const_cast<KDTreeN *>(this)->processRangeUntil<IntersectionTesterT>(root, range, &functor);
      }
    }

    /**
     * Apply a functor to all objects in a range, until the functor returns true. The functor should provide the member function
     * (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(long index, T & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object).
     *
     * The RangeT class should support intersection queries with AxisAlignedBoxT and containment queries with VectorT and
     * AxisAlignedBoxT.
     *
     * @return True if the functor evaluated to true on any object in the range (and hence stopped immediately after processing
     *   this object), else false.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    bool processRangeUntil(RangeT const & range, FunctorT * functor)
    {
      if (!functor)  // no-op
        return false;

      return root ? processRangeUntil<IntersectionTesterT>(root, range, functor) : false;
    }

    template <typename RayIntersectionTesterT> bool rayIntersects(RayT const & ray, Real max_time = -1) const
    {
      return rayIntersectionTime<RayIntersectionTesterT>(ray, max_time) >= 0;
    }

    template <typename RayIntersectionTesterT> Real rayIntersectionTime(RayT const & ray, Real max_time = -1) const
    {
      if (root)
      {
        if (TransformableBaseT::hasTransform())
        {
          RayT tr_ray = toObjectSpace(ray);
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
    RayStructureIntersectionT rayStructureIntersection(RayT const & ray, Real max_time = -1) const
    {
      if (root)
      {
        if (TransformableBaseT::hasTransform())
        {
          RayT tr_ray = toObjectSpace(ray);
          if (root->bounds.rayIntersects(tr_ray, max_time))
          {
            RayStructureIntersectionT isec = rayStructureIntersection<RayIntersectionTesterT>(root, tr_ray, max_time);
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

      return RayStructureIntersectionT(-1);
    }

  private:
    /** Comparator for sorting elements along an axis. */
    struct ObjectLess
    {
      long coord;
      KDTreeN const * tree;

      /** Constructor. Axis 0 = X, 1 = Y, 2 = Z. */
      ObjectLess(long coord_, KDTreeN const * tree_) : coord(coord_), tree(tree_) {}

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

    typedef TheaArray<Filter<T> *> FilterStack;  ///< A stack of element filters.
    typedef TheaArray<SampleFilter> SampleFilterStack;  ///< A stack of point sample filters.

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
#define THEA_KDTREEN_SPLIT_LONGEST
#ifdef THEA_KDTREEN_SPLIT_LONGEST
      long coord = start->bounds.getExtent().maxAxis();  // split longest dimension
#else
      long coord = (long)(start->depth % N);  // cycle between dimensions
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
      AxisAlignedBoxT elem_bounds;
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
      {
        bounds = TransformableBaseT::hasTransform()
               ? root->bounds.transformAndBound(TransformableBaseT::getTransform())
               : root->bounds;
      }
      else
      {
        bounds = AxisAlignedBoxT();
      }

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
    AxisAlignedBoxT getBoundsWorldSpace(Node const & node) const
    {
      return TransformableBaseT::hasTransform()
           ? node.bounds.transformAndBound(TransformableBaseT::getTransform())
           : node.bounds;
    }

    /**
     * Recursively look for the closest pair of points between two elements. Only pairs separated by less than the current
     * minimum distance will be considered.
     */
    template <typename MetricT, typename QueryT>
    void closestPair(Node const * start, QueryT const & query, AxisAlignedBoxT const & query_bounds, NeighborPair & pair,
                     bool get_closest_points) const
    {
      if (!start->lo)  // leaf
        closestPairLeaf<MetricT>(start, query, pair, get_closest_points);
      else  // not leaf
      {
        // Figure out which child is closer (optimize for point queries?)
        Node const * n[2] = { start->lo, start->hi };
        double mad[2] = { MetricT::template monotoneApproxDistance<N, ScalarT>(getBoundsWorldSpace(*n[0]), query_bounds),
                          MetricT::template monotoneApproxDistance<N, ScalarT>(getBoundsWorldSpace(*n[1]), query_bounds) };

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
        if (TransformableBaseT::hasTransform())
          swapped = query.template closestPair<MetricT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                        pair.getMonotoneApproxDistance(), get_closest_points);
        else
          swapped = query.template closestPair<MetricT>(elem, pair.getMonotoneApproxDistance(), get_closest_points);

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
      VectorT qp, tp;
      double mad;

      for (array_size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        if (TransformableBaseT::hasTransform())
          mad = MetricT::template closestPoints<N, ScalarT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                            query, tp, qp);
        else
          mad = MetricT::template closestPoints<N, ScalarT>(elem, query, tp, qp);

        if (pair.getMonotoneApproxDistance() < 0 || mad <= pair.getMonotoneApproxDistance())
          pair = NeighborPair(0, (long)index, mad, qp, tp);
      }
    }

    /**
     * Recursively look for the k closest elements to a query object. Only elements at less than the specified maximum distance
     * will be considered.
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSet>
    void kClosestPairs(Node const * start, QueryT const & query, AxisAlignedBoxT const & query_bounds,
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
        double d[2] = { MetricT::template monotoneApproxDistance<N, ScalarT>(getBoundsWorldSpace(*n[0]), query_bounds),
                        MetricT::template monotoneApproxDistance<N, ScalarT>(getBoundsWorldSpace(*n[1]), query_bounds) };

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
     * structure of compatible type.
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

        if (TransformableBaseT::hasTransform())
          query.template kClosestPairs<MetricT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                                      k_closest_pairs, dist_bound, get_closest_points, false,
                                                                      (long)index);
        else
          query.template kClosestPairs<MetricT>(elem, k_closest_pairs, dist_bound, get_closest_points, false, (long)index);
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
      VectorT qp, tp;
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

        if (TransformableBaseT::hasTransform())
          mad = MetricT::template closestPoints<N, ScalarT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                            query, tp, qp);
        else
          mad = MetricT::template closestPoints<N, ScalarT>(elem, query, tp, qp);

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
     * RangeT class should support containment queries with AxisAlignedBoxT.
     *
     * @return True if the functor evaluated to true on any point in the range, else false.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    bool processRangeUntil(Node const * start, RangeT const & range, FunctorT * functor)
    {
      // Early exit if the range and node are disjoint
      AxisAlignedBoxT tr_start_bounds = getBoundsWorldSpace(*start);
      if (!IntersectionTesterT::template intersects<N, ScalarT>(range, tr_start_bounds))
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

          bool intersects = TransformableBaseT::hasTransform()
                          ? IntersectionTesterT::template intersects<N, ScalarT>(
                                makeTransformedObject(&elem, &TransformableBaseT::getTransform()), range)
                          : IntersectionTesterT::template intersects<N, ScalarT>(elem, range);
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
    RayT toObjectSpace(RayT const & ray) const
    {
      AffineTransform3 inv_trans = TransformableBaseT::getTransform().inverse();  // TODO: cache this
      return RayT(inv_trans * ray.getOrigin(), inv_trans.getLinear() * ray.getDirection());
    }

    /** Transform a normal to world space. */
    VectorT normalToWorldSpace(VectorT const & n) const
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
    Real rayIntersectionTime(Node const * start, RayT const & ray, Real max_time) const
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

          Real time = RayIntersectionTesterT::template rayIntersectionTime<N, ScalarT>(ray, elem, best_time);
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
    RayStructureIntersectionT rayStructureIntersection(Node const * start, RayT const & ray, Real max_time) const
    {
      if (!start->lo)  // leaf
      {
        RayStructureIntersectionT best_isec(max_time);
        bool found = false;
        for (array_size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element const & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          RayIntersection3 isec = RayIntersectionTesterT::template rayIntersection<N, ScalarT>(ray, elem, best_isec.getTime());
          if (improvedRayTime(isec.getTime(), best_isec.getTime()))
          {
            best_isec = RayStructureIntersectionT(isec, (long)index);
            found = true;
          }
        }

        return found ? best_isec : RayStructureIntersectionT(-1);
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

        RayStructureIntersectionT best_isec(max_time);
        bool found = false;
        for (int i = 0; i < 2; ++i)
        {
          if (improvedRayTime(t[i], best_isec.getTime()))
          {
            RayStructureIntersectionT isec = rayStructureIntersection<RayIntersectionTesterT>(n[i], ray, best_isec.getTime());
            if (improvedRayTime(isec.getTime(), best_isec.getTime()))
            {
              best_isec = isec;
              found = true;
            }
          }
        }

        return found ? best_isec : RayStructureIntersectionT(-1);
      }
    }

    /** Build a structure to accelerate nearest neighbor queries. */
    template <typename MetricT> void buildAccelerationStructure() const
    {
      if (valid_acceleration_structure) return;

      delete acceleration_structure; acceleration_structure = NULL;

      static int const DEFAULT_NUM_ACCELERATION_SAMPLES = 250;
      TheaArray<ElementSample> acceleration_samples(num_acceleration_samples <= 0 ? DEFAULT_NUM_ACCELERATION_SAMPLES
                                                                                  : num_acceleration_samples);
      VectorT src_cp, dst_cp;
      for (array_size_t i = 0; i < acceleration_samples.size(); ++i)
      {
        long elem_index = Random::common().integer(0, (int32)num_elems - 1);
        VectorT p = BoundedObjectTraitsT::getCenter(elems[elem_index]);

        // Snap point to element, else it's not a valid NN proxy
        MetricT::template closestPoints<N, ScalarT>(p, elems[elem_index], src_cp, dst_cp);

        acceleration_samples[i] = ElementSample(dst_cp, &elems[elem_index]);
      }

      acceleration_structure = new NearestNeighborAccelerationStructure;
      acceleration_structure->disableNearestNeighborAcceleration();
      acceleration_structure->init(acceleration_samples.begin(), acceleration_samples.end());

      if (this->hasTransform())
        acceleration_structure->setTransform(this->getTransform());

      sample_filters.clear();
      for (array_size_t i = 0; i < filters.size(); ++i)
      {
        sample_filters.push_back(SampleFilter(filters[i]));
        acceleration_structure->pushFilter(&sample_filters.back());  // careful, make sure order is not reversed!
      }

      valid_acceleration_structure = true;
    }

    /** Destroy any existing structure to accelerate nearest neighbor queries. */
    void clearAccelerationStructure(bool deallocate_memory = true) const
    {
      valid_acceleration_structure = false;

      if (deallocate_memory)
      {
        delete acceleration_structure;
        acceleration_structure = NULL;

        sample_filters.clear();
      }
    }

    /** Get an upper bound on the distance to a query object, using the acceleration structure if it exists. */
    template <typename MetricT, typename QueryT>
    double accelerationBound(QueryT const & query, double dist_bound) const
    {
      NearestNeighborAccelerationStructure const * accel = getNearestNeighborAccelerationStructure<MetricT>();
      return accel ? accel->template distance<MetricT>(query, dist_bound) : -1;
    }

    /**
     * Get an upper bound on the distance to a query kd-tree, using the acceleration structures of both the query and of this
     * object if they exist.
     */
    template <typename MetricT, typename E, typename S, typename A, typename B>
    double accelerationBound(KDTreeN<E, N, S, A, B> const & query, double dist_bound) const
    {
      NearestNeighborAccelerationStructure const * accel = getNearestNeighborAccelerationStructure<MetricT>();
      if (accel)
      {
        if (query.hasNearestNeighborAcceleration())
        {
          typename KDTreeN<E, N, S, A, B>::NearestNeighborAccelerationStructure const * query_accel
              = query.template getNearestNeighborAccelerationStructure<MetricT>();

          if (query_accel)
            return accel->template distance<MetricT>(*query_accel, dist_bound);
        }

        return accel->template distance<MetricT>(query, dist_bound);
      }
      else
      {
        if (query.hasNearestNeighborAcceleration())
        {
          typename KDTreeN<E, N, S, A, B>::NearestNeighborAccelerationStructure const * query_accel
              = query.template getNearestNeighborAccelerationStructure<MetricT>();

          if (query_accel)
            return distance<MetricT>(*query_accel, dist_bound);
        }

        return -1;
      }
    }

    Node * root;

    long num_elems;  // elems.size() doesn't tell us how many elements there are, it's just the capacity of the elems array
    ElementArray elems;  // elems.size() is *not* the number of elements in the tree!!!

    long num_nodes;
    NodePool node_pool;

    IndexPool index_pool;

    long max_depth;
    long max_elems_in_leaf;

    MatrixMN<3, 3, ScalarT> transform_inverse_transpose;

    FilterStack filters;

    bool accelerate_nn_queries;
    long num_acceleration_samples;
    mutable bool valid_acceleration_structure;
    mutable NearestNeighborAccelerationStructure * acceleration_structure;
    mutable SampleFilterStack sample_filters;

    mutable bool valid_bounds;
    mutable AxisAlignedBoxT bounds;

    static Real const BOUNDS_EXPANSION_FACTOR;

}; // class KDTreeN

// Static variables
template <typename T, long N, typename S, typename A, typename B>
Real const KDTreeN<T, N, S, A, B>::BOUNDS_EXPANSION_FACTOR = 1.05f;

} // namespace Algorithms
} // namespace Thea

#endif
