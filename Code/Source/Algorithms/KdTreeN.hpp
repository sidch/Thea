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

#ifndef __Thea_Algorithms_KdTreeN_hpp__
#define __Thea_Algorithms_KdTreeN_hpp__

#include "../Common.hpp"
#include "../AffineTransformN.hpp"
#include "../Array.hpp"
#include "../AttributedObject.hpp"
#include "../Math.hpp"
#include "../Noncopyable.hpp"
#include "../Random.hpp"
#include "../Transformable.hpp"
#include "BoundedTraitsN.hpp"
#include "Filter.hpp"
#include "MetricL2.hpp"
#include "ProximityQueryStructureN.hpp"
#include "RangeQueryStructure.hpp"
#include "RayQueryStructureN.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <type_traits>

namespace Thea {
namespace Algorithms {

namespace KdTreeNInternal {

/** A point sample drawn from a kd-tree element, used for accelerating nearest neighbor queries. */
template <int N, typename ScalarT>
struct ElementSample
{
  Vector<N, ScalarT> position;
  void const * element;  // this cannot be pointer-to-T, else we get recursive instantiation overflow

  /** Default constructor. */
  ElementSample() {}

  /** Initialize a sample at point \a p drawn from an element \a e. */
  ElementSample(Vector<N, ScalarT> const & p, void const * e) : position(p), element(e) {}

}; // struct ElementSample

/**
 * A filter that passes or rejects a sample point, depending on whether a base filter passes or rejects the point's parent
 * element.
 */
template <typename T, int N, typename ScalarT>
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

} // namespace KdTreeNInternal

template <int N, typename ScalarT>
class IsPointN< KdTreeNInternal::ElementSample<N, ScalarT>, N >
{
  public:
    static bool const value = true;
};

template <int N, typename ScalarT>
class PointTraitsN< KdTreeNInternal::ElementSample<N, ScalarT>, N, ScalarT >
{
  public:
    typedef Vector<N, ScalarT> VectorT;
    static VectorT getPosition(KdTreeNInternal::ElementSample<N, ScalarT> const & sample) { return sample.position; }
};

/**
 * A kd-tree for a set of bounded objects in N-space. IsBoundedNN<T, N> must evaluate to true, and BoundedTraitsN<T, N>
 * appropriately defined to compute bounding boxes of T objects. The optional template parameter <code>NodeAttributeT</code>
 * must be default-constructible.
 *
 * An affine transformation may be applied to the kd-tree. The tree does <em>not</em> need to be recomputed after the
 * transformation, though all operations may be somewhat slower (the precise overhead depends on how difficult it is to compute
 * distances, intersections etc after a transform). Normally, this means that the elements of type T should be affine
 * transformable via Transformer::transform().
 */
template < typename T,
           int N,
           typename ScalarT = Real,
           typename NodeAttributeT = NullAttribute >
class /* THEA_API */ KdTreeN
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
    typedef BoundedTraitsN<T, N, ScalarT>                  BoundedTraitsT;

  public:
    typedef size_t ElementIndex;  ///< Index of an element in the kd-tree.
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
            Buffer(size_t capacity_ = 0) : data(nullptr), capacity(capacity_), current_end(0)
            {
              alwaysAssertM(capacity_ > 0, "KdTreeN: Memory pool buffer capacity must be positive");
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
            U * alloc(size_t num_elems)
            {
              // THEA_CONSOLE << "KdTreeN: Allocating " << num_elems << " elements from buffer of capacity " << capacity;

              if (current_end + num_elems > capacity)
                return nullptr;

              U * ret = &data[current_end];
              current_end += num_elems;
              return ret;
            }

            /** Deallocate a block of elements and return the number of elements successfully deallocated. */
            size_t free(size_t num_elems)
            {
              if (num_elems > current_end)
              {
                size_t num_freed = current_end;
                current_end = 0;
                return num_freed;
              }

              current_end -= num_elems;
              return num_elems;
            }

            /** Element access. */
            U const & operator[](size_t i) const { return data[i]; }

            /** Element access. */
            U & operator[](size_t i) { return data[i]; }

          private:
            U * data;
            size_t capacity;
            size_t current_end;

        }; // Buffer

        Array<Buffer *> buffers;
        size_t buffer_capacity;
        intx current_buffer;

        /** Get the current buffer. */
        Buffer & getCurrentBuffer() { return *buffers[(size_t)current_buffer]; }

        /** Get the next available buffer, creating it if necessary and making it the current one. */
        Buffer & getNextBuffer()
        {
          intx next_buffer = current_buffer < 0 ? 0 : current_buffer + 1;
          if ((size_t)next_buffer >= buffers.size())
          {
            buffers.push_back(new Buffer(buffer_capacity));

            // THEA_CONSOLE << "KdTreeN: Added buffer to memory pool " << this << ", current_buffer = " << current_buffer
            //              << ", next_buffer = " << next_buffer;
          }

          current_buffer = next_buffer;
          return *buffers[(size_t)current_buffer];
        }

      public:
        /** Constructor. */
        MemoryPool() : buffer_capacity(0), current_buffer(-1)
        {
          // THEA_CONSOLE << "KdTreeN: Creating memory pool " << this;
        }

        /** Destructor. */
        ~MemoryPool()
        {
          clear(true);
          // THEA_CONSOLE << "KdTreeN: Destroyed memory pool " << this;
        }

        /** Initialize the memory pool to hold buffers of a given capacity. Previous data in the pool is deallocated. */
        void init(size_t buffer_capacity_)
        {
          // THEA_CONSOLE << "KdTreeN: Initializing memory pool " << this << " with buffer capacity " << buffer_capacity_
          //              << " elements";

          clear(true);

          buffer_capacity = buffer_capacity_;
        }

        /** Get the maximum number of elements of type T that a single buffer can hold. */
        size_t getBufferCapacity() const
        {
          return buffer_capacity;
        }

        /** Reset the memory pool, optionally deallocating and removing all buffers. */
        void clear(bool deallocate_all_memory = true)
        {
          // THEA_CONSOLE << "KdTreeN: Clearing memory pool " << this;

          if (deallocate_all_memory)
          {
            for (size_t i = 0; i < buffers.size(); ++i)
              delete buffers[i];

            buffers.clear();
          }
          else
          {
            for (size_t i = 0; i < buffers.size(); ++i)
              buffers[i]->reset();
          }

          current_buffer = -1;
        }

        /** Allocate a block of elements from the pool and return a pointer to the first allocated element. */
        U * alloc(size_t num_elems)
        {
          alwaysAssertM(num_elems <= buffer_capacity, "KdTreeN: A single memory pool allocation cannot exceed buffer capacity");

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
        void free(size_t num_elems)
        {
          intx n = (intx)num_elems;
          while (n > 0 && current_buffer >= 0)
          {
            size_t num_freed = getCurrentBuffer().free(num_elems);
            n -= (intx)num_freed;

            if (n > 0)
              current_buffer--;
          }
        }

    }; // class MemoryPool

    /** A functor to add results of a range query to an array. */
    class RangeQueryFunctor
    {
      public:
        RangeQueryFunctor(Array<T> & result_) : result(result_) {}
        bool operator()(intx index, T & t) { result.push_back(t); return false; }

      private:
        Array<T> & result;
    };

    /** A functor to add the indices of results of a range query to an array. */
    class RangeQueryIndicesFunctor
    {
      public:
        RangeQueryIndicesFunctor(Array<intx> & result_) : result(result_) {}
        bool operator()(intx index, T & t) { result.push_back(index); return false; }

      private:
        Array<intx> & result;
    };

  public:
    THEA_DECL_SMART_POINTERS(KdTreeN)

    typedef T                                    Element;        ///< Type of elements in the kd-tree.
    typedef T                                    value_type;     ///< Type of elements in the kd-tree (STL convention).
    typedef NodeAttributeT                       NodeAttribute;  ///< Attributes attached to nodes.

    typedef typename ProximityQueryBaseT::VectorT              VectorT;                    ///< Vector in N-space.
    typedef AxisAlignedBoxN<N, ScalarT>                        AxisAlignedBoxT;            ///< Axis-aligned box in N-space.
    typedef typename RayQueryBaseT::RayT                       RayT;                       ///< Ray in N-space.
    typedef typename RayQueryBaseT::RayStructureIntersectionT  RayStructureIntersectionT;  /**< Ray intersection structure in
                                                                                                N-space. */

    typedef KdTreeNInternal::ElementSample<N, ScalarT> ElementSample;  /**< A point sample drawn from a kd-tree element, used
                                                                            for accelerating nearest neighbor queries. */
    typedef KdTreeN<ElementSample, N, ScalarT> NearestNeighborAccelerationStructure;  /**< Structure to speed up nearest
                                                                                           neighbor queries. */

    /** A node of the kd-tree. Only immutable objects of this class should be exposed by the external kd-tree interface. */
    class Node : public AttributedObject<NodeAttributeT>
    {
      private:
        intx depth;
        AxisAlignedBoxT bounds;
        size_t num_elems;
        ElementIndex * elems;
        Node * lo;
        Node * hi;

        friend class KdTreeN;

        void init(intx depth_)
        {
          depth = depth_;
          bounds = AxisAlignedBoxT();
          num_elems = 0;
          elems = nullptr;
          lo = hi = nullptr;
        }

      public:
        /** Iterator over immutable element indices. Dereferences to an array index. */
        typedef ElementIndex const * ElementIndexConstIterator;

        /** Constructor. */
        Node(intx depth_ = 0) : depth(depth_), lo(nullptr), hi(nullptr) {}

        /** Get the depth of the node in the tree (the root is at depth 0). */
        intx getDepth() const { return depth; }

        /** Get the bounding box of the node. */
        AxisAlignedBoxT const & getBounds() const { return bounds; }

        /**
         * Get the number of element indices stored at this node. This is <b>not</b> the number of elements within the node's
         * bounding box: in memory-saving mode, indices of all such elements are only held at the leaves of the subtree rooted
         * at this node.
         */
        intx numElementIndices() const { return (intx)num_elems; }

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

  protected:
    typedef MemoryPool<Node> NodePool;  ///< A pool for quickly allocating kd-tree nodes.
    typedef MemoryPool<ElementIndex> IndexPool;  ///< A pool for quickly allocating element indices.

  private:
    typedef Array<T> ElementArray;  ///< An array of elements.
    typedef KdTreeNInternal::SampleFilter<T, N, ScalarT> SampleFilter;  ///< Filter for samples, wrapping a filter for elements.

  public:
    /** Default constructor. */
    KdTreeN()
    : root(nullptr), num_elems(0), num_nodes(0), max_depth(0), max_elems_per_leaf(0),
      transform_inverse(AffineTransformN<N, ScalarT>::identity()),
      transform_inverse_transpose(Matrix<N, N, ScalarT>::Identity()),
      accelerate_nn_queries(false), valid_acceleration_structure(false), acceleration_structure(nullptr), valid_bounds(true)
    {}

    /**
     * Construct from a list of elements. InputIterator must dereference to type T.
     *
     * @param begin Points to the first element to be added.
     * @param end Points to one position beyond the last element to be added.
     * @param max_depth_ Maximum depth of the tree. The root is at depth zero. Use a negative argument to auto-select a suitable
     *   value.
     * @param max_elems_per_leaf_ Maximum number of elements in a leaf (unless the depth exceeds the maximum). Use a negative
     *   argument to auto-select a suitable value.
     * @param save_memory If true, element references at inner nodes of the tree are deleted to save memory. This could slow
     *   down range searches since every positive result will only be obtained at the leaves.
     */
    template <typename InputIterator>
    KdTreeN(InputIterator begin, InputIterator end, intx max_depth_ = -1, intx max_elems_per_leaf_ = -1,
            bool save_memory = false)
    : root(nullptr), num_elems(0), num_nodes(0), max_depth(0), max_elems_per_leaf(0),
      transform_inverse(AffineTransformN<N, ScalarT>::identity()),
      transform_inverse_transpose(Matrix<N, N, ScalarT>::Identity()),
      accelerate_nn_queries(false), valid_acceleration_structure(false), acceleration_structure(nullptr), valid_bounds(true)
    {
      init(begin, end, max_elems_per_leaf_, max_depth_, save_memory, false /* no previous data to deallocate */);
    }

    /**
     * Construct from a list of elements. InputIterator must dereference to type T. Any previous data is discarded. If any
     * filters are active at this time, only those input elements that pass the filters will be retained in the tree.
     *
     * @param begin Points to the first element to be added.
     * @param end Points to one position beyond the last element to be added.
     * @param max_depth_ Maximum depth of the tree. The root is at depth zero. Use a negative argument to auto-select a suitable
     *   value.
     * @param max_elems_per_leaf_ Maximum number of elements in a leaf (unless the depth exceeds the maximum). Use a negative
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
    void init(InputIterator begin, InputIterator end, intx max_depth_ = -1, intx max_elems_per_leaf_ = -1,
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
        size_t max_new_elems = (size_t)std::distance(begin, end);
        bool resized = false;
        if (max_new_elems > elems.size())
        {
          if (filters.empty())
            elems.resize((size_t)std::ceil(1.2 * max_new_elems));  // add a little more space to avoid future reallocs
          else
            elems.clear();  // we don't know how many elements will pass the filter

          resized = true;
        }

        if (filters.empty())
        {
          std::copy(begin, end, elems.begin());
          num_elems = (intx)max_new_elems;
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

      static intx const DEFAULT_MAX_ELEMS_IN_LEAF = 10;
      max_elems_per_leaf = max_elems_per_leaf_ < 0 ? DEFAULT_MAX_ELEMS_IN_LEAF : max_elems_per_leaf_;

      // The fraction of elements held by the larger node at each split is 0.5
      static double const SPLIT_FRACTION = 0.5;
      intx est_depth = Math::binaryTreeDepth(num_elems, max_elems_per_leaf, SPLIT_FRACTION);
      max_depth = max_depth_;
      if (max_depth < 0)
        max_depth = est_depth;
      else if (max_depth < est_depth)
        est_depth = max_depth;

      // THEA_CONSOLE << "KdTreeN: max_depth = " << max_depth << ", est_depth = " << est_depth;

      // Each index is stored at most once at each level
      size_t BUFFER_SAFETY_MARGIN = 10;
      size_t index_buffer_capacity = num_elems + BUFFER_SAFETY_MARGIN;
      if (!save_memory)
        index_buffer_capacity *= (size_t)(1 + est_depth);  // reserve space for all levels at once

      if (deallocate_previous_memory || index_buffer_capacity > 1.3 * index_pool.getBufferCapacity())
      {
        // THEA_CONSOLE << "KdTreeN: Resizing index pool: old buffer capacity = " << index_pool.getBufferCapacity()
        //              << ", new buffer capacity = " << index_buffer_capacity;
        index_pool.init(index_buffer_capacity);
      }

      // Assume a complete, balanced binary tree upto the estimated depth to guess the number of leaf nodes
      size_t node_buffer_capacity = (size_t)(1 << est_depth) + BUFFER_SAFETY_MARGIN;
      if (deallocate_previous_memory || node_buffer_capacity > 1.3 * node_pool.getBufferCapacity())
      {
        // THEA_CONSOLE << "KdTreeN: Resizing node pool: old buffer capacity = " << node_pool.getBufferCapacity()
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
      for (size_t i = 0; i < (size_t)num_elems; ++i)
      {
        root->elems[i] = i;

        BoundedTraitsT::getBounds(elems[i], elem_bounds);
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
        size_t est_max_path_indices = (size_t)(num_elems * (1 + 1 / (1 - SPLIT_FRACTION)));
        // THEA_CONSOLE << "KdTreeN: Estimated maximum number of indices on a single path = " << est_max_path_indices;

        // Create a temporary pool for scratch storage
        IndexPool tmp_index_pool;
        tmp_index_pool.init(est_max_path_indices + BUFFER_SAFETY_MARGIN);

        createTree(root, true, &tmp_index_pool, &index_pool);
      }
      else
        createTree(root, false, &index_pool, nullptr);

      invalidateBounds();
    }

    /**
     * Recompute the bounding box of every node in the tree (or a specific subtree) from the element data at the leaves. The
     * tree topology is unchanged. This function is useful for quickly updating the tree when relative element sizes and
     * element-to-element proximity do not change much (e.g. under character articulations). Of course, this also requires that
     * the leaf elements (objects of type <tt>T</tt>) be pointers or other handles to external data that can be changed outside
     * this object.
     *
     * @param start The root of the subtree to update. If null, the entire tree is updated.
     */
    void updateNodeBounds(Node * start = nullptr)
    {
      if (!start)
      {
        if (!root) return;  // nothing to process
        start = root;
      }

      updateNodeBoundsRecursive(start);

      if (accelerate_nn_queries && valid_acceleration_structure && acceleration_structure)
      {
        // Would be nicer to touch only the smallest subtree of points from updated elements, but this is too complicated for
        // now, so we'll just resample all the points

        // Ok to cast away the const, this is a privately held object and we have full control over it
        samplePointsFromElements(acceleration_structure->numElements(),
                                 const_cast<ElementSample *>(acceleration_structure->getElements()));
        acceleration_structure->updateNodeBounds();
      }
    }

    /** Destructor. */
    ~KdTreeN() { clear(true); }

    /**
     * Enable acceleration of nearest neighbor queries with an auxiliary structure on a sparse set of points.
     *
     * @param num_acceleration_samples_ Number of sampled points in the acceleration structure. If this number is zero, no
     *   structure will be created. If it is negative, the number will be automatically determined as a fraction of the number
     *   of elements. In the latter case, if the number is very small, no structure will be created.
     *
     * @see disableNearestNeighborAcceleration()
     */
    void enableNearestNeighborAcceleration(intx num_acceleration_samples_ = -1)
    {
      if (num_acceleration_samples_ < 0)
      {
        static intx const MIN_ACCEL_SAMPLES = 10;

        num_acceleration_samples_ = std::min(250L, (intx)std::ceil(0.1 * numElements()));
        accelerate_nn_queries = (num_acceleration_samples_ >= MIN_ACCEL_SAMPLES);
      }
      else
        accelerate_nn_queries = (num_acceleration_samples_ > 0);

      if (accelerate_nn_queries)
      {
        if (valid_acceleration_structure && num_acceleration_samples != num_acceleration_samples_)
          valid_acceleration_structure = false;

        num_acceleration_samples = num_acceleration_samples_;
      }
    }

    /**
     * Disable acceleration of nearest neighbor queries with an auxiliary structure on a sparse set of points.
     *
     * @see enableNearestNeighborAcceleration()
     */
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
        return nullptr;
    }

    void setTransform(Transform const & trans_)
    {
      TransformableBaseT::setTransform(trans_);
      transform_inverse = trans_.inverse();
      transform_inverse_transpose = transform_inverse.getLinear().transpose();
      invalidateBounds();

      if (valid_acceleration_structure)
        acceleration_structure->setTransform(trans_);
    }

    void clearTransform()
    {
      TransformableBaseT::clearTransform();
      transform_inverse = AffineTransformN<N, ScalarT>::identity();
      transform_inverse_transpose = Matrix<N, N, ScalarT>::Identity();
      invalidateBounds();

      if (valid_acceleration_structure)
        acceleration_structure->clearTransform();
    }

    /** Get the inverse of the transform applied to the kd-tree. */
    AffineTransformN<N, ScalarT> const & getTransformInverse() const { return transform_inverse; }

    /** Get the transpose of the inverse of the linear part of the transform applied to the kd-tree. */
    Matrix<N, N, ScalarT> const & getLinearTransformInverseTranspose() const { return transform_inverse_transpose; }

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

      root = nullptr;

      invalidateBounds();
    }

    /** Check if the tree is empty. */
    bool isEmpty() const { return num_elems <= 0; }

    /** Get the number of elements in the tree. The elements themselves can be obtained with getElements(). */
    intx numElements() const { return num_elems; }

    /** Get a pointer to an array of the elements in the tree. The number of elements can be obtained with numElements(). */
    T const * getElements() const { return &elems[0]; }

    /**
     * Get the node corresponding to the root of the kd-tree. This function is provided so that users can implement their own
     * tree traversal procedures.
     *
     * This function cannot be used to change the structure of the tree, or any value in it (unless <code>const_cast</code> is
     * used, which is not recommended).
     *
     * @note An empty tree has a null root.
     */
    Node const * getRoot() const { return root; }

    /** Get the number of nodes in the tree. */
    intx numNodes() const { return num_nodes; }

    /** Get the maximum subdivision depth (number of levels not counting the root) of the kd-tree. */
    intx maxDepth() const { return max_depth; }

    /** Get the maximum number of elements in each leaf of the kd-tree. */
    intx maxElementsPerLeaf() const { return max_elems_per_leaf; }

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
      alwaysAssertM(filter, "KdTreeN: Filter must be non-null");

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
    intx closestElement(QueryT const & query, double dist_bound = -1, double * dist = nullptr,
                        VectorT * closest_point = nullptr) const
    {
      NeighborPair pair = closestPair<MetricT>(query, dist_bound, closest_point != nullptr);

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
     * @param query Query object. BoundedTraitsN<QueryT, N, ScalarT> must be defined.
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

      // Early pruning if the entire structure is too far away from the query
      AxisAlignedBoxT query_bounds;
      getObjectBounds(query, query_bounds);
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = monotonePruningDistance<MetricT>(root, query, query_bounds);
        if (lower_bound >= 0 && lower_bound > mon_approx_dist_bound)
          return NeighborPair(-1);
      }

      // If acceleration is enabled, set an upper limit to the distance to the nearest object
      double accel_bound = accelerationBound<MetricT>(query, dist_bound);
      if (accel_bound >= 0)
      {
        double fudge = 0.001 * getBoundsWorldSpace(*root).getExtent().norm();
        mon_approx_dist_bound = MetricT::computeMonotoneApprox(accel_bound + fudge);
      }

      NeighborPair pair(-1, -1, mon_approx_dist_bound);
      closestPair<MetricT>(root, query, query_bounds, pair, get_closest_points);

      return pair;
    }

    /**
     * Get the k elements closest to a query object. The returned elements are placed in a set of bounded size (k). The template
     * type BoundedNeighborPairSetT should typically be BoundedSortedArray<NeighborPair> or BoundedSortedArrayN<k, NeighborPair>
     * if only a few neighbors are requested. BoundedTraitsN<QueryT, N, ScalarT> must be defined.
     *
     * @param query Query object. BoundedTraitsN<QueryT, N, ScalarT> must be defined.
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
    intx kClosestPairs(QueryT const & query, BoundedNeighborPairSetT & k_closest_pairs, double dist_bound = -1,
                       bool get_closest_points = false, bool clear_set = true, intx use_as_query_index_and_swap = -1) const
    {
      if (clear_set) k_closest_pairs.clear();

      if (!root) return 0;

      // Early pruning if the entire structure is too far away from the query
      AxisAlignedBoxT query_bounds;
      getObjectBounds(query, query_bounds);
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = monotonePruningDistance<MetricT>(root, query, query_bounds);
        if (lower_bound >= 0)
        {
          if (lower_bound > mon_approx_dist_bound)
            return 0;

          if (!k_closest_pairs.isInsertable(NeighborPair(0, 0, lower_bound)))
            return 0;
        }
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
    void rangeQuery(RangeT const & range, Array<T> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root) const_cast<KdTreeN *>(this)->processRangeUntil<IntersectionTesterT>(range, RangeQueryFunctor(result));
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
    void rangeQueryIndices(RangeT const & range, Array<intx> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root) const_cast<KdTreeN *>(this)->processRangeUntil<IntersectionTesterT>(range, RangeQueryIndicesFunctor(result));
    }

    /**
     * Apply a functor to all objects in a range, until the functor returns true. The functor should provide the member function
     * (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T const & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * The RangeT class should support intersection queries with AxisAlignedBoxT and containment queries with VectorT and
     * AxisAlignedBoxT.
     *
     * @return The index of the first object in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this object), else a negative value.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    intx processRangeUntil(RangeT const & range, FunctorT functor) const
    {
      return root ? const_cast<KdTreeN *>(this)->processRangeUntil<IntersectionTesterT, T const>(root, range, functor) : -1;
    }

    /**
     * Apply a functor to all objects in a range, until the functor returns true. The functor should provide the member function
     * (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T [const] & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * The RangeT class should support intersection queries with AxisAlignedBoxT and containment queries with VectorT and
     * AxisAlignedBoxT.
     *
     * @return The index of the first object in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this object), else a negative value.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    intx processRangeUntil(RangeT const & range, FunctorT functor)
    {
      return root ? processRangeUntil<IntersectionTesterT, T>(root, range, functor) : -1;
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
      intx coord;
      KdTreeN const * tree;

      /** Constructor. Axis 0 = X, 1 = Y, 2 = Z. */
      ObjectLess(intx coord_, KdTreeN const * tree_) : coord(coord_), tree(tree_) {}

      /** Less-than operator, along the specified axis. */
      bool operator()(ElementIndex a, ElementIndex b)
      {
        // Compare object min coords
        return BoundedTraitsT::getLow(tree->elems[a], coord)
             < BoundedTraitsT::getLow(tree->elems[b], coord);
      }
    };

    // Allow the comparator unrestricted access to the kd-tree.
    friend struct ObjectLess;

    typedef Array<Filter<T> *> FilterStack;  ///< A stack of element filters.
    typedef Array<SampleFilter> SampleFilterStack;  ///< A stack of point sample filters.

  protected:
    /** Get the main memory pool for nodes. */
    NodePool & getNodePool() { return node_pool; }

    /** Get the main memory pool for indices. */
    IndexPool & getIndexPool() { return index_pool; }

  private:
    /** Move the set of element indices for a leaf from the end of the main index pool to an index pool for leaves. */
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

  protected:
    /**
     * Recursively construct a (sub-)tree. If \a save_memory is true, element indices for this subtree are assumed to be in a
     * block at the end of \a main_index_pool. They are subsequently moved to arrays associated with the leaves, allocated in
     * \a leaf_index_pool. They are deleted from \a main_index_pool when this function exits, thus deallocating index arrays
     * held by internal nodes. If \a save_memory is false, each node maintains its own list of all elements in its subtree.
     *
     * @note To use this function with a node accessed via getRoot(), you'll have to const_cast it to a non-const type first.
     */
    void createTree(Node * start, bool save_memory, IndexPool * main_index_pool, IndexPool * leaf_index_pool)
    {
      // Assume the start node is fully constructed at this stage.
      //
      // If we are in memory-saving mode, then we assume the node's indices are last in the main index pool (but not yet in the
      // leaf pool, since we don't yet know if this node will turn out to be a leaf). In this case, after the function finishes,
      // the node's indices will be deallocated from the main index pool (and possibly moved to the leaf pool).

      if (!start || start->depth >= max_depth || (intx)start->num_elems <= max_elems_per_leaf)
      {
        if (save_memory)
          moveIndicesToLeafPool(start, main_index_pool, leaf_index_pool);

        return;
      }

      // Find a splitting plane
#define THEA_KDTREEN_SPLIT_LONGEST
#ifdef THEA_KDTREEN_SPLIT_LONGEST
      intx coord = Math::maxAxis(start->bounds.getExtent());  // split longest dimension
#else
      intx coord = (intx)(start->depth % N);  // cycle between dimensions
#endif

      // Split elements into lower and upper halves
      size_t mid = start->num_elems / 2;
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
        BoundedTraitsT::getBounds(elems[index], elem_bounds);

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
        BoundedTraitsT::getBounds(elems[index], elem_bounds);

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
        start->elems = nullptr;
      }
    }

    /**
     * Recompute the bounding box of a leaf node from its elements. Can be overriden in derived classes for additional speed if
     * possible.
     *
     * @param leaf A pointer to a leaf node. This is guaranteed to be non-null, and an actual leaf without children.
     */
    virtual void updateLeafBounds(Node * leaf)
    {
      leaf->bounds.setNull();

      AxisAlignedBoxT elem_bounds;
      for (ElementIndex i = 0; i < leaf->num_elems; ++i)
      {
        BoundedTraitsT::getBounds(elems[leaf->elems[i]], elem_bounds);
        leaf->bounds.merge(elem_bounds);
      }

      // Expand the bounding box slightly to handle numerical error
      leaf->bounds.scaleCentered(BOUNDS_EXPANSION_FACTOR);
    }

  private:
    /** Recursively applied helper function for updateNodeBounds(). */
    void updateNodeBoundsRecursive(Node * start)
    {
      if (!start->lo)  // leaf
        updateLeafBounds(start);
      else  // recurse
      {
        updateNodeBoundsRecursive(start->lo);
        updateNodeBoundsRecursive(start->hi);
      }
    }

  protected:
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

    /** Get the bounding box for an object, if it is bounded. */
    template < typename U, typename std::enable_if< IsBoundedN<U, N>::value, int >::type = 0 >
    static void getObjectBounds(U const & u, AxisAlignedBoxT & bounds)
    {
      BoundedTraitsN<U, N, ScalarT>::getBounds(u, bounds);
    }

    /** Returns a null bounding box for unbounded objects. */
    template < typename U, typename std::enable_if< !IsBoundedN<U, N>::value, int >::type = 0 >
    static void getObjectBounds(U const & u, AxisAlignedBoxT & bounds)
    {
      bounds.setNull();
    }

    /**
     * Get a lower bound on (the monotone approximation to) the distance between the bounding box of a kd-tree node and a
     * bounded query object, or a negative value if no such lower bound can be calculated.
     */
    template < typename MetricT, typename QueryT, typename std::enable_if< IsBoundedN<QueryT, N>::value, int >::type = 0 >
    double monotonePruningDistance(Node const * node, QueryT const & query, AxisAlignedBoxT query_bounds) const
    {
      if (node && !query_bounds.isNull())
        return MetricT::template monotoneApproxDistance<N, ScalarT>(query_bounds, getBoundsWorldSpace(*node));
      else
        return -1;
    }

    /**
     * Get a lower bound on (the monotone approximation to) the distance between the bounding box of a kd-tree node and an
     * unbounded query object, or a negative value if no such lower bound can be calculated. \a query_bounds is ignored.
     */
    template < typename MetricT, typename QueryT, typename std::enable_if< !IsBoundedN<QueryT, N>::value, int >::type = 0 >
    double monotonePruningDistance(Node const * node, QueryT const & query, AxisAlignedBoxT query_bounds) const
    {
      // Assume the following specialization exists
      return node ? MetricT::template monotoneApproxDistance<N, ScalarT>(query, getBoundsWorldSpace(*node)) : -1;
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
     * minimum distance (as stored in \a pair) will be considered. If \a get_closest_points is true, the positions of the
     * closest pair of points will be stored in \a pair, not just the distance between them.
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
        double mad[2] = { monotonePruningDistance<MetricT>(n[0], query, query_bounds),
                          monotonePruningDistance<MetricT>(n[1], query, query_bounds) };

        // The smaller non-negative value should be first
        if (mad[1] >= 0 && (mad[0] < 0 || mad[0] > mad[1]))
        {
          std::swap(n[0], n[1]);
          std::swap(mad[0], mad[1]);
        }

        for (int i = 0; i < 2; ++i)
          if (pair.getMonotoneApproxDistance() < 0 || mad[i] <= pair.getMonotoneApproxDistance())
            closestPair<MetricT>(n[i], query, query_bounds, pair, get_closest_points);
      }
    }

  private:
    /**
     * Search the elements in a leaf node for the one closest to another element, when the latter is a proximity query
     * structure.
     */
    template < typename MetricT, typename QueryT,
               typename std::enable_if< std::is_base_of<ProximityQueryBaseT, QueryT>::value, int >::type = 0 >
    void closestPairLeaf(
      Node const * leaf,
      QueryT const & query,
      NeighborPair & pair,
      bool get_closest_points) const
    {
      for (size_t i = 0; i < leaf->num_elems; ++i)
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
          pair.setTargetIndex((intx)index);
        }
      }
    }

    /**
     * Search the elements in a leaf node for the one closest to another element, when the latter is NOT a proximity query
     * structure.
     */
    template < typename MetricT, typename QueryT,
               typename std::enable_if< !std::is_base_of<ProximityQueryBaseT, QueryT>::value, int >::type = 0 >
    void closestPairLeaf(
      Node const * leaf,
      QueryT const & query,
      NeighborPair & pair,
      bool get_closest_points) const
    {
      VectorT qp = VectorT::Zero(), tp = VectorT::Zero();  // initialize to squash uninitialized variable warning
      double mad;

      for (size_t i = 0; i < leaf->num_elems; ++i)
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
          pair = NeighborPair(0, (intx)index, mad, qp, tp);
      }
    }

  protected:
    /**
     * Recursively look for the k closest elements to a query object. Only elements at less than the specified maximum distance
     * \a dist_bound will be considered. If \a get_closest_points is true, the positions of the closest pair of points will be
     * stored with each pair, not just the distance between them.
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSet>
    void kClosestPairs(Node const * start, QueryT const & query, AxisAlignedBoxT const & query_bounds,
                       BoundedNeighborPairSet & k_closest_pairs, double dist_bound, bool get_closest_points,
                       intx use_as_query_index_and_swap) const
    {
      if (!start->lo)  // leaf
        kClosestPairsLeaf<MetricT>(start, query, k_closest_pairs, dist_bound, get_closest_points, use_as_query_index_and_swap);
      else  // not leaf
      {
        // Figure out which child is closer (optimize for point queries?)
        Node const * n[2] = { start->lo, start->hi };
        double mad[2] = { monotonePruningDistance<MetricT>(n[0], query, query_bounds),
                          monotonePruningDistance<MetricT>(n[1], query, query_bounds) };

        // The smaller non-negative value should be first
        if (mad[1] >= 0 && (mad[0] < 0 || mad[0] > mad[1]))
        {
          std::swap(n[0], n[1]);
          std::swap(mad[0], mad[1]);
        }

        double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
        for (int i = 0; i < 2; ++i)
        {
          if ((mon_approx_dist_bound < 0 || mad[i] <= mon_approx_dist_bound)
            && k_closest_pairs.isInsertable(NeighborPair(0, 0, mad[i])))
          {
            kClosestPairs<MetricT>(n[i], query, query_bounds, k_closest_pairs, dist_bound, get_closest_points,
                                   use_as_query_index_and_swap);
          }
        }
      }
    }

  private:
    /**
     * Search the elements in a leaf node for the k nearest neighbors of an object, when the latter is a proximity query
     * structure of compatible type.
     */
    template < typename MetricT, typename QueryT, typename BoundedNeighborPairSet,
               typename std::enable_if< std::is_base_of<ProximityQueryBaseT, QueryT>::value, int >::type = 0 >
    void kClosestPairsLeaf(
      Node const * leaf,
      QueryT const & query,
      BoundedNeighborPairSet & k_closest_pairs,
      double dist_bound,
      bool get_closest_points,
      intx use_as_query_index_and_swap) const
    {
      for (size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        if (TransformableBaseT::hasTransform())
          query.template kClosestPairs<MetricT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                                      k_closest_pairs, dist_bound, get_closest_points, false,
                                                                      (intx)index);
        else
          query.template kClosestPairs<MetricT>(elem, k_closest_pairs, dist_bound, get_closest_points, false, (intx)index);
      }
    }

    /**
     * Search the elements in a leaf node for the one closest to another element, when the latter is NOT a proximity query
     * structure.
     */
    template < typename MetricT, typename QueryT, typename BoundedNeighborPairSet,
               typename std::enable_if< !std::is_base_of<ProximityQueryBaseT, QueryT>::value, int >::type = 0 >
    void kClosestPairsLeaf(
      Node const * leaf,
      QueryT const & query,
      BoundedNeighborPairSet & k_closest_pairs,
      double dist_bound,
      bool get_closest_points,
      intx use_as_query_index_and_swap) const
    {
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);

      std::equal_to<NeighborPair> eq_comp;
      VectorT qp, tp;
      double mad;

      for (size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        // Check if the element is already in the set of neighbors or not
        NeighborPair pair = (use_as_query_index_and_swap >= 0 ? NeighborPair((intx)index, use_as_query_index_and_swap)
                                                              : NeighborPair(0, (intx)index));
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

  protected:
    /**
     * Apply a functor to all elements of a subtree within a range, until the functor returns true. The functor should provide
     * the member function (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T & t)
     * \endcode
     * and will be passed the index of each object contained in the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object). To pass a functor by reference, wrap it in <tt>std::ref</tt>
     *
     * The RangeT class should support intersection queries with AxisAlignedBoxT and containment queries with VectorT and
     * AxisAlignedBoxT.
     *
     * @return The index of the first object in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this object), else a negative value.
     */
    template <typename IntersectionTesterT, typename FunctorArgT, typename RangeT, typename FunctorT>
    intx processRangeUntil(Node const * start, RangeT const & range, FunctorT functor)
    {
      // Early exit if the range and node are disjoint
      AxisAlignedBoxT tr_start_bounds = getBoundsWorldSpace(*start);
      if (!IntersectionTesterT::template intersects<N, ScalarT>(range, tr_start_bounds))
        return -1;

      // If the entire node is contained in the range AND there are element references at this node (so it's either a leaf or we
      // have not saved memory by flushing references at internal nodes), then process all these elems.
      if (start->num_elems > 0 && range.contains(tr_start_bounds))
      {
        for (size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          if (functor(static_cast<intx>(index), static_cast<FunctorArgT &>(elem)))
            return static_cast<intx>(index);
        }
      }
      else if (!start->lo)  // leaf
      {
        for (size_t i = 0; i < start->num_elems; ++i)
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
            if (functor(static_cast<intx>(index), static_cast<FunctorArgT &>(elem)))
              return static_cast<intx>(index);
        }
      }
      else  // not leaf
      {
        intx index = processRangeUntil<IntersectionTesterT, FunctorArgT>(start->lo, range, functor);
        if (index >= 0) return index;
        return processRangeUntil<IntersectionTesterT, FunctorArgT>(start->hi, range, functor);
      }

      return -1;
    }

  private:
    /** Transform a ray to local/object space. */
    RayT toObjectSpace(RayT const & ray) const
    {
      return RayT(transform_inverse * ray.getOrigin(), transform_inverse.getLinear() * ray.getDirection());
    }

    /** Transform a normal to world space. */
    VectorT normalToWorldSpace(VectorT const & n) const
    {
      return (transform_inverse_transpose * n).normalized();
    }

    /**
     * Check if the ray intersection time \a new_time represents a closer, or equally close, valid hit than the previous best
     * time \a old_time.
     */
    static bool improvedRayTime(Real new_time, Real old_time)
    {
      return (new_time >= 0 && (old_time < 0 || new_time <= old_time));
    }

  protected:
    /** Get the time taken for a ray to hit the nearest object in a node, in the forward direction. */
    template <typename RayIntersectionTesterT>
    Real rayIntersectionTime(Node const * start, RayT const & ray, Real max_time) const
    {
      if (!start->lo)  // leaf
      {
        Real best_time = max_time;
        bool found = false;
        for (size_t i = 0; i < start->num_elems; ++i)
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
        for (size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element const & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          RayIntersectionN<N, ScalarT> isec = RayIntersectionTesterT::template rayIntersection<N, ScalarT>(ray, elem,
                                                                                                           best_isec.getTime());
          if (improvedRayTime(isec.getTime(), best_isec.getTime()))
          {
            best_isec = RayStructureIntersectionT(isec, (intx)index);
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

  protected:
    /**
     * Given a preallocated array of ElementSample objects, each pointing to a valid element, assign each object a point
     * (randomly or deterministically) sampled from the corresponding element. The default implementation finds the point
     * closest to the bounding box center of each element: subclasses can override this with faster sampling methods.
     *
     * @note This version finds nearest neighbors on elements using the L2 distance. Hence MetricL2 must support NN queries with
     *   the element type T.
     */
    virtual void samplePointsFromElements(intx num_samples, ElementSample * samples) const
    {
      VectorT src_cp;
      for (intx i = 0; i < num_samples; ++i)
      {
        T const * elem = static_cast<T const *>(samples[i].element);
        VectorT p = BoundedTraitsT::getCenter(*elem);
        MetricL2::closestPoints<N, ScalarT>(p, *elem, src_cp, samples[i].position);  // snap point to element
      }
    }

  private:
    /** Build a structure to accelerate nearest neighbor queries. */
    template <typename MetricT> void buildAccelerationStructure() const
    {
      if (valid_acceleration_structure) return;

      delete acceleration_structure; acceleration_structure = nullptr;

      static int const DEFAULT_NUM_ACCELERATION_SAMPLES = 250;
      Array<ElementSample> acceleration_samples(num_acceleration_samples <= 0 ? DEFAULT_NUM_ACCELERATION_SAMPLES
                                                                              : num_acceleration_samples);
      for (size_t i = 0; i < acceleration_samples.size(); ++i)
        acceleration_samples[i].element = &elems[Random::common().integer(0, (int32)num_elems - 1)];

      samplePointsFromElements((intx)acceleration_samples.size(), acceleration_samples.data());

      acceleration_structure = new NearestNeighborAccelerationStructure;
      acceleration_structure->disableNearestNeighborAcceleration();
      acceleration_structure->init(acceleration_samples.begin(), acceleration_samples.end());

      if (this->hasTransform())
        acceleration_structure->setTransform(this->getTransform());

      sample_filters.clear();
      for (size_t i = 0; i < filters.size(); ++i)
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
        acceleration_structure = nullptr;

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
    double accelerationBound(KdTreeN<E, N, S, A> const & query, double dist_bound) const
    {
      NearestNeighborAccelerationStructure const * accel = getNearestNeighborAccelerationStructure<MetricT>();
      if (accel)
      {
        if (query.hasNearestNeighborAcceleration())
        {
          typename KdTreeN<E, N, S, A>::NearestNeighborAccelerationStructure const * query_accel
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
          typename KdTreeN<E, N, S, A>::NearestNeighborAccelerationStructure const * query_accel
              = query.template getNearestNeighborAccelerationStructure<MetricT>();

          if (query_accel)
            return distance<MetricT>(*query_accel, dist_bound);
        }

        return -1;
      }
    }

  private:
    Node * root;

    intx num_elems;  // elems.size() doesn't tell us how many elements there are, it's just the capacity of the elems array
    ElementArray elems;  // elems.size() is *not* the number of elements in the tree!!!

    intx num_nodes;
    NodePool node_pool;

    IndexPool index_pool;

    intx max_depth;
    intx max_elems_per_leaf;

    AffineTransformN<N, ScalarT> transform_inverse;
    Matrix<N, N, ScalarT> transform_inverse_transpose;

    FilterStack filters;

    bool accelerate_nn_queries;
    intx num_acceleration_samples;
    mutable bool valid_acceleration_structure;
    mutable NearestNeighborAccelerationStructure * acceleration_structure;
    mutable SampleFilterStack sample_filters;

    mutable bool valid_bounds;
    mutable AxisAlignedBoxT bounds;

    static Real const BOUNDS_EXPANSION_FACTOR;

}; // class KdTreeN

// Static variables
template <typename T, int N, typename S, typename A>
Real const KdTreeN<T, N, S, A>::BOUNDS_EXPANSION_FACTOR = 1.05f;

// Mark the kd-tree and its (public) descendants as bounded objects. The default BoundedTraitsN implementation is good enough.
template <typename T, int N>
class IsBoundedN< T, N, typename std::enable_if< std::is_base_of< KdTreeN< typename T::Element, N,
                                                                           typename T::VectorT::value_type,
                                                                           typename T::NodeAttribute >,
                                                                  T >::value >::type >
{
  public:
    static bool const value = true;
};

} // namespace Algorithms
} // namespace Thea

#endif
