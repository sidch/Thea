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

#ifndef __Thea_Algorithms_BvhN_hpp__
#define __Thea_Algorithms_BvhN_hpp__

#include "../Common.hpp"
#include "../AffineTransformN.hpp"
#include "../Array.hpp"
#include "../AttributedObject.hpp"
#include "../BoundedSortedArrayN.hpp"
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
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <type_traits>

namespace Thea {
namespace Algorithms {

// Forward declaration
template <typename T, int N, typename S, typename A, int D, typename V> class BvhN;

namespace BvhNInternal {

// A point sample drawn from a BVH element, used for accelerating nearest neighbor queries.
template <typename T, int N, typename ScalarT>
struct ElementSample
{
  Vector<N, ScalarT> position;
  T const * element;

  /** Default constructor. */
  ElementSample() {}

  /** Initialize a sample at point \a p drawn from an element \a e. */
  ElementSample(Vector<N, ScalarT> const & p, T const * e) : position(p), element(e) {}

}; // struct ElementSample

// A filter that passes or rejects a sample point, depending on whether a base filter passes or rejects the point's parent
// element.
template <typename T, int N, typename ScalarT>
class SampleFilter : public Filter< ElementSample<T, N, ScalarT> >
{
  public:
    /** Constructor. */
    SampleFilter(Filter<T> * base_filter_) : base_filter(base_filter_) {}

    bool allows(ElementSample<T, N, ScalarT> const & sample) const
    {
      // Assume null pointer case won't happen by construction
      return base_filter->allows(*sample.element);
    }

  private:
    Filter<T> * base_filter;  ///< The underlying element filter.

}; // struct SampleFilter

// A dummy class with the interface of BvhN, whose functions are all no-ops. Used as the acceleration structure for acceleration
// structures, to avoid recursive instantiation overflow.
template <typename T, int N, typename ScalarT>
struct DummyAccelerationStructure
{
  intx numElements() const { return 0; }
  ElementSample<T, N, ScalarT> * getElements() const { return nullptr; }
  void updateNodeBounds() {}
  template <typename TransformT> void setTransform(TransformT const & trans_) {}
  void clearTransform() {}
  template <typename FilterT> void pushFilter(FilterT * filter) {}
  void popFilter() {}
  void disableNearestNeighborAcceleration() {}
  template <typename SampleIteratorT> void init(SampleIteratorT begin, SampleIteratorT end) {}

  template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
  double
  distance(QueryT const & query, double dist_bound = -1, CompatibilityFunctorT compatibility = CompatibilityFunctorT()) const
  {
    return Math::inf<double>();
  }

}; // DummyAccelerationStructure

// Defines the acceleration structure for a BVH.
template <typename T, int N, typename ScalarT, typename NodeAttributeT, int MaxDegree, typename BoundingVolumeT>
struct AccelerationTraits
{
  /** Structure to speed up nearest neighbor queries. */
  typedef BvhN< ElementSample<T, N, ScalarT>, N, ScalarT, NodeAttributeT, 2, AxisAlignedBoxN<N, ScalarT> >
          NearestNeighborAccelerationStructure;

}; // class AccelerationTraits

// Defines a dummy acceleration structure for a BVH that is itself an acceleration structure.
template <typename T, int N, typename ScalarT, typename NodeAttributeT, int MaxDegree, typename BoundingVolumeT >
struct AccelerationTraits< ElementSample<T, N, ScalarT>, N, ScalarT, NodeAttributeT, MaxDegree, BoundingVolumeT >
{
  typedef DummyAccelerationStructure< ElementSample<T, N, ScalarT>, N, ScalarT > NearestNeighborAccelerationStructure;

}; // class AccelerationTraits< ElementSample<T, N, ScalarT>, N, ScalarT >

// A utility class that exchanges the arguments to a compatibility functor.
template <typename CompatibilityFunctorT> class SwappedCompatibility
{
  public:
    SwappedCompatibility(CompatibilityFunctorT compatibility_) : base_compatibility(compatibility_) {}
    template <typename U, typename V> bool operator()(U const & u, V const & v) const { return base_compatibility(v, u); }

    CompatibilityFunctorT base_compatibility;  // The base compatibility functor.

}; // class SwappedCompatibility

// Utility function to construct a SwappedCompatibility object.
template <typename CompatibilityFunctorT>
SwappedCompatibility<CompatibilityFunctorT>
makeSwappedCompatibility(CompatibilityFunctorT compatibility)
{
  return SwappedCompatibility<CompatibilityFunctorT>(compatibility);
}

// Version of makeSwappedCompatibility that handles the nested swap case (swap(swap(c)) == c).
template <typename CompatibilityFunctorT>
CompatibilityFunctorT
makeSwappedCompatibility(SwappedCompatibility<CompatibilityFunctorT> swapped)
{
  return swapped.base_compatibility;
}

// A utility class that evaluates the compatibility of samples in terms of the compatibility of their parent elements.
template <typename CompatibilityFunctorT> class SampleCompatibility
{
  public:
    SampleCompatibility(CompatibilityFunctorT compatibility_) : compatibility(compatibility_) {}

    template <typename U, typename V> bool operator()(U const & u, V const & v) const
    {
      return compatibility(getObject(u), getObject(v));
    }

  private:
    template <typename T> static T const & getObject(T const & t) { return t; }

    template <typename T, int N, typename ScalarT>
    static T const & getObject(ElementSample<T, N, ScalarT> const & t)
    {
      return *t.element;
    }

    CompatibilityFunctorT compatibility;

}; // class SampleCompatibility

// Check if a distance value \a lhs is less than another distance value \a rhs. All negative distances are considered equal, and
// greater than any positive value.
template <typename U, typename V> bool distanceLessThan(U lhs, V rhs)
{
  return (lhs >= 0 && (rhs < 0 || lhs < rhs));
}

} // namespace BvhNInternal

template <typename T, int N, typename ScalarT>
class IsPointN< BvhNInternal::ElementSample<T, N, ScalarT>, N > { public: static bool const value = true; };

template <typename T, int N, typename ScalarT>
class PointTraitsN< BvhNInternal::ElementSample<T, N, ScalarT>, N, ScalarT >
{
  public:
    typedef Vector<N, ScalarT> VectorT;
    static VectorT getPosition(BvhNInternal::ElementSample<T, N, ScalarT> const & sample) { return sample.position; }
};

// Define a DummyAccelerationStructure as a logical point, to enable MetricL2 to handle it and thereby allow compilation to
// proceed. These functions will never actually be called.
template <typename T, int N, typename ScalarT>
class IsPointN< BvhNInternal::DummyAccelerationStructure<T, N, ScalarT>, N > { public: static bool const value = true; };

template <typename T, int N, typename ScalarT>
class PointTraitsN< BvhNInternal::DummyAccelerationStructure<T, N, ScalarT>, N, ScalarT >
{
  public:
    typedef Vector<N, ScalarT> VectorT;

    static VectorT getPosition(BvhNInternal::DummyAccelerationStructure<T, N, ScalarT> const & dummy)
    {
      (void)dummy;
      return VectorT::Constant(Math::inf<ScalarT>());
    }
};

/**
 * A bounding volume hierarchy (BVH) for a set of bounded objects in N-space. IsBoundedN<T, N> must evaluate to true, and
 * BoundedTraitsN<T, N> appropriately defined to compute bounding volumes of T objects. The optional template parameter
 * <code>NodeAttributeT</code> must be default-constructible.
 *
 * BvhN supports fast nearest neighbor (including BVH-to-BVH), range, and ray intersection queries. An optional auxiliary
 * structure based on point samples can be automatically computed to accelerate nearest neighbor queries
 * (see enableNearestNeighborAcceleration()).
 *
 * This class is designed to handle many different types of spatial hierarchies, including kd-trees (in any dimension),
 * quadtrees (in 2 dimensions), octrees (in 3 dimension) etc. The appropriate type of hierarchy is selected (or auto-selected)
 * via the <tt>method</tt> parameter of the init() function. (Currently only kd-tree construction is supported.) Note that the
 * MaxDegree parameter, which specifies the maximum number of children a tree node can have, must be at least 2 for kd-trees, 4
 * for quadtrees, and 8 for octrees.
 *
 * The default (and currently the only supported) bounding volume is an axis-aligned bounding box. The class is designed to
 * support other bounding volumes (e.g. spheres) in the future, but it currently makes some AAB-specific assumptions so
 * instantiating the class with (say) BallN will not yet work.
 *
 * An affine transformation may be applied to the BVH. The tree does <em>not</em> need to be recomputed after the
 * transformation, though all operations may be somewhat slower (the precise overhead depends on how difficult it is to compute
 * distances, intersections etc after a transform, and how much a transformation enlarges each bounding volume).
 * Normally, this means that the elements of type T should be affine-transformable via Transformer::transform().
 *
 * If the elements in the tree are modified in such a way that siblings stay relatively close together (e.g. an articulating
 * human), then the bounding volumes can be rapidly updated, without changing the tree topology, by calling updateNodeBounds().
 * This is much faster than recomputing the tree from scratch.
 *
 * @todo Support bounding volumes other than AxisAlignedBoxN. This should be fairly straightforward but needs a bunch of small
 *   fixes here and there.
 */
template < typename T,
           int N,
           typename ScalarT = Real,
           typename NodeAttributeT = NullAttribute,
           int MaxDegree = 2,
           typename BoundingVolumeT = AxisAlignedBoxN<N, ScalarT> >
class /* THEA_API */ BvhN
: public RangeQueryStructure<T>,
  public ProximityQueryStructureN<N, ScalarT>,
  public RayQueryStructureN<N, ScalarT>,
  public Transformable< AffineTransformN<N, ScalarT> >,
  private Noncopyable
{
  private:
    static_assert(N >= 1, "BvhN: Dimensionality must be at least 1");
    static_assert(MaxDegree >= 2, "BvhN: Maximum degree (number of children) must be at least 2");

    typedef RangeQueryStructure<T>                         RangeQueryBaseT;
    typedef ProximityQueryStructureN<N, ScalarT>           ProximityQueryBaseT;
    typedef RayQueryStructureN<N, ScalarT>                 RayQueryBaseT;
    typedef Transformable< AffineTransformN<N, ScalarT> >  TransformableBaseT;
    typedef BoundedTraitsN<T, N, ScalarT>                  BoundedTraitsT;

    typedef BvhNInternal::AccelerationTraits<T, N, ScalarT, NodeAttributeT, MaxDegree, BoundingVolumeT>  AccelerationTraitsT;

  public:
    typedef size_t ElementIndex;  ///< Index of an element in the BVH.
    typedef typename ProximityQueryBaseT::NeighborPair NeighborPair;  ///< A pair of neighboring elements.
    typedef typename TransformableBaseT::Transform Transform;  ///< Transform applied to the BVH.

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
              alwaysAssertM(capacity_ > 0, "BvhN: Memory pool buffer capacity must be positive");
              if (capacity > 0)
                data = new U[capacity];
            }

            /** Destructor. */
            ~Buffer() { delete [] data; }

            /** Clears the buffer without deallocating buffer memory. */
            void reset() { current_end = 0; }

            /**
             * Allocate a block of elements and return a pointer to the first allocated element, or null if the allocation
             * exceeded buffer capacity.
             */
            U * alloc(size_t num_elems)
            {
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
            buffers.push_back(new Buffer(buffer_capacity));

          current_buffer = next_buffer;
          return *buffers[(size_t)current_buffer];
        }

      public:
        /** Constructor. */
        MemoryPool() : buffer_capacity(0), current_buffer(-1) {}

        /** Destructor. */
        ~MemoryPool() { clear(true); }

        /** Initialize the memory pool to hold buffers of a given capacity. Previous data in the pool is deallocated. */
        void init(size_t buffer_capacity_)
        {
          clear(true);

          buffer_capacity = buffer_capacity_;
        }

        /** Get the maximum number of elements of type T that a single buffer can hold. */
        size_t getBufferCapacity() const { return buffer_capacity; }

        /** Reset the memory pool, optionally deallocating and removing all buffers. */
        void clear(bool deallocate_all_memory = true)
        {
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
          alwaysAssertM(num_elems <= buffer_capacity, "BvhN: A single memory pool allocation cannot exceed buffer capacity");

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

    /** A functor to get the element closest to a query. */
    template <typename MetricT>
    class ClosestPairFunctor
    {
      public:
        ClosestPairFunctor(double dist_bound, NeighborPair & result_)
        : result(result_)
        {
          if (!result_.isValid())  // not pre-initialized
            result_ = NeighborPair(-1, -1, (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1));
        }

        bool allows(NeighborPair const & pair) const
        {
          return true;  // do distance-based pruning in operator() -- not sure if we also want to do that here or not
        }

        bool operator()(NeighborPair const & pair)
        {
          auto best_mad = result.getMonotoneApproxDistance();
          if (best_mad < 0 || pair.getMonotoneApproxDistance() < best_mad)
            result = pair;

          return false;
        }

      private:
        NeighborPair & result;

    }; // class ClosestPairFunctor

    /** A functor to get the element closest to a query. */
    template <typename MetricT, typename BoundedNeighborPairSetT>
    class KClosestPairsFunctor
    {
      public:
        KClosestPairsFunctor(double dist_bound_, BoundedNeighborPairSetT & result_)
        : mad_bound(dist_bound_ >= 0 ? MetricT::computeMonotoneApprox(dist_bound_) : -1), result(result_)
        {}

        bool allows(NeighborPair const & pair) const
        {
          return (!pair.isValid()  // invalid indices implies distance-based check, which we will ignore here
               || !result.contains(pair, std::equal_to<NeighborPair>()));  // avoid duplicates
        }

        bool operator()(NeighborPair const & pair)
        {
          if ((mad_bound < 0 || pair.getMonotoneApproxDistance() < mad_bound) && result.isInsertable(pair))
            result.insert(pair);

          return false;
        }

      private:
        double mad_bound;
        BoundedNeighborPairSetT & result;

    }; // class KClosestPairsFunctor

    /** A functor to add results of a range query to an array. */
    class RangeQueryFunctor
    {
      public:
        RangeQueryFunctor(Array<T> & result_) : result(result_) {}
        bool operator()(intx index, T & t) { result.push_back(t); return false; }

      private:
        Array<T> & result;

    }; // class RangeQueryFunctor

    /** A functor to add the indices of results of a range query to an array. */
    class RangeQueryIndicesFunctor
    {
      public:
        RangeQueryIndicesFunctor(Array<intx> & result_) : result(result_) {}
        bool operator()(intx index, T & t) { result.push_back(index); return false; }

      private:
        Array<intx> & result;

    }; // class RangeQueryIndicesFunctor

  public:
    THEA_DECL_SMART_POINTERS(BvhN)

    typedef T                Element;         ///< Type of elements in the BVH.
    typedef T                value_type;      ///< Type of elements in the BVH (STL convention).
    typedef NodeAttributeT   NodeAttribute;   ///< Attributes attached to nodes.
    typedef BoundingVolumeT  BoundingVolume;  ///< Bounding volume for elements in N-space.

    typedef typename ProximityQueryBaseT::VectorT              VectorT;                    ///< Vector in N-space.
    typedef typename RayQueryBaseT::RayT                       RayT;                       ///< Ray in N-space.
    typedef typename RayQueryBaseT::RayStructureIntersectionT  RayStructureIntersectionT;  /**< Ray intersection structure in
                                                                                                N-space. */

    /** A point sample drawn from a BVH element, used for accelerating nearest neighbor queries. */
    typedef BvhNInternal::ElementSample<T, N, ScalarT> ElementSample;

    /** Structure to speed up nearest neighbor queries. */
    typedef typename AccelerationTraitsT::NearestNeighborAccelerationStructure NearestNeighborAccelerationStructure;

    /** A node of the BVH. Only immutable objects of this class should be exposed by the external BVH interface. */
    class Node : public AttributedObject<NodeAttributeT>
    {
      private:
        typedef std::array<Node *, MaxDegree> ChildArray;

        intx depth;
        BoundingVolume bounds;
        size_t num_elems;
        ElementIndex * elems;
        ChildArray children;

        template <typename E, int M, typename S, typename A, int D, typename V> friend class BvhN;

        void init(intx depth_)
        {
          depth = depth_;
          bounds = BoundingVolume();
          num_elems = 0;
          elems = nullptr;
          children.fill(nullptr);
        }

      public:
        /** Iterator over immutable element indices. Dereferences to an array index. */
        typedef ElementIndex const * ElementIndexConstIterator;

        /** Constructor. */
        Node(intx depth_ = 0) : depth(depth_) { children.fill(nullptr); }

        /** Get the depth of the node in the tree (the root is at depth 0). */
        intx getDepth() const { return depth; }

        /** Get the bounding volume of the node. */
        BoundingVolume const & getBounds() const { return bounds; }

        /**
         * Get the number of element indices stored at this node. This is <b>not</b> the number of elements within the node's
         * bounding volume: in memory-saving mode, indices of all such elements are only held at the leaves of the subtree
         * rooted at this node.
         */
        intx numElementIndices() const { return (intx)num_elems; }

        /** Get an iterator to the first element index stored at the node. */
        ElementIndexConstIterator elementIndicesBegin() const { return elems; }

        /** Get an iterator to one past the last element index stored at the node. */
        ElementIndexConstIterator elementIndicesEnd() const { return elems + num_elems; }

        /**
         * Get the maximum possible number of children of a node. It is possible no node in a tree actually has this many
         * children.
         */
        static constexpr int maxDegree() { return MaxDegree; }

        /**
         * Get the number of non-null children of the node. This involves a pass over the array of children and is not just a
         * value lookup. Also note that the non-null children may not be contiguously stored in the array.
         */
        int numChildren() const { return MaxDegree - (int)std::count(children.begin(), children.end(), nullptr); }

        /** Get the child with a given index. The returned handle may be null if no child exists with this index. */
        Node const * getChild(intx i) const { return children[(size_t)i]; }

        /** Check if the node is a leaf (all children are null) are not. */
        bool isLeaf() const
        {
          for (auto c : children)
            if (c) { return false; }

          return true;
        }

    }; // class Node

    /** The method used for tree construction (enum class). */
    struct Method
    {
      /** Supported values. */
      enum Value
      {
        AUTO,      ///< Auto-select a suitable method based on dimensionality, maximum degree etc.
        KDTREE,    ///< Construct a kd-tree (split node with a single axis-aligned plane in each subdivision step).
        QUADTREE,  ///< Construct a quadtree (split node into 4 quadrants in each subdivision step). Only for <tt>N == 2</tt>.
        OCTREE,    ///< Construct an octree (split node into 8 octants in each subdivision step). Only for <tt>N == 3</tt>.
      };

      THEA_ENUM_CLASS_BODY(Method)

      THEA_ENUM_CLASS_STRINGS_BEGIN(Method)
        THEA_ENUM_CLASS_STRING(AUTO,      "auto")
        THEA_ENUM_CLASS_STRING(KDTREE,    "kd-tree")
        THEA_ENUM_CLASS_STRING(QUADTREE,  "quadtree")
        THEA_ENUM_CLASS_STRING(OCTREE,    "octree")
      THEA_ENUM_CLASS_STRINGS_END(Method)
    };

  protected:
    typedef MemoryPool<Node> NodePool;  ///< A pool for quickly allocating BVH nodes.
    typedef MemoryPool<ElementIndex> IndexPool;  ///< A pool for quickly allocating element indices.

  private:
    typedef Array<T> ElementArray;  ///< An array of elements.
    typedef BvhNInternal::SampleFilter<T, N, ScalarT> SampleFilter;  ///< Filter for samples, wrapping a filter for elements.

  public:
    /** Default constructor. */
    BvhN()
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
     * @param method The tree construction method. Pass Method::AUTO to auto-select a suitable method.
     * @param max_depth_ Maximum depth of the tree. The root is at depth zero. Use a negative argument to auto-select a suitable
     *   value.
     * @param max_elems_per_leaf_ Maximum number of elements in a leaf (unless the depth exceeds the maximum). Use a negative
     *   argument to auto-select a suitable value.
     * @param save_memory If true, element references at inner nodes of the tree are deleted to save memory. This could slow
     *   down range searches since every positive result will only be obtained at the leaves.
     */
    template <typename InputIterator>
    BvhN(InputIterator begin, InputIterator end, Method method = Method::AUTO, intx max_depth_ = -1,
         intx max_elems_per_leaf_ = -1, bool save_memory = false)
    : root(nullptr), num_elems(0), num_nodes(0), max_depth(0), max_elems_per_leaf(0),
      transform_inverse(AffineTransformN<N, ScalarT>::identity()),
      transform_inverse_transpose(Matrix<N, N, ScalarT>::Identity()),
      accelerate_nn_queries(false), valid_acceleration_structure(false), acceleration_structure(nullptr), valid_bounds(true)
    {
      init(begin, end, method, max_elems_per_leaf_, max_depth_, save_memory, /* no previous data to deallocate, hence */ false);
    }

    /**
     * Construct from a list of elements. InputIterator must dereference to type T. Any previous data is discarded. If any
     * filters are active at this time, only those input elements that pass the filters will be retained in the tree.
     *
     * @param begin Points to the first element to be added.
     * @param end Points to one position beyond the last element to be added.
     * @param method The tree construction method. Pass Method::AUTO to auto-select a suitable method.
     * @param max_depth_ Maximum depth of the tree. The root is at depth zero. Use a negative argument to auto-select a suitable
     *   value.
     * @param max_elems_per_leaf_ Maximum number of elements in a leaf (unless the depth exceeds the maximum). Use a negative
     *   argument to auto-select a suitable value.
     * @param save_memory If true, element references at inner nodes of the tree are deleted to save memory. This could slow
     *   down range searches since every positive result will only be obtained at the leaves.
     * @param deallocate_previous_memory If true, all previous data held in internal memory pools is explicitly deallocated.
     *   Else, all such space is reused and overwritten when possible. If \a save_memory is true, or some filters are active,
     *   this flag may not be quite as effective since it's more likely that some space will be allocated/deallocated. Note that
     *   if this flag is set to false, the space used internally by the BVH will not decrease except in some special
     *   implementation-specific cases.
     */
    template <typename InputIterator>
    void init(InputIterator begin, InputIterator end, Method method = Method::AUTO, intx max_depth_ = -1,
              intx max_elems_per_leaf_ = -1, bool save_memory = false, bool deallocate_previous_memory = true)
    {
      if (method == Method::AUTO)
      {
        if (MaxDegree >= 4 && N == 2)
          method = Method::QUADTREE;
        else if (MaxDegree >= 8 && N == 3)
          method = Method::OCTREE;
        else if (MaxDegree >= 2)
          method = Method::KDTREE;
        else
          throw Error("BvhN: Could not infer tree construction method");
      }

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

      switch (method)
      {
        case Method::KDTREE:
          initKdTree(begin, end, max_depth_, max_elems_per_leaf_, save_memory, deallocate_previous_memory); break;

        default: throw Error("BvhN: Unsupported BVH construction method");
      }
    }

    /**
     * Update cached properties of all elements in the BVH, or in one of its subtrees. In this base class, this function
     * does nothing, but see MeshBvh for a non-trivial override.
     */
    virtual void updateElements(Node const * start = nullptr) {}

    /**
     * Recompute the bounding volume of every node in the tree (or a specific subtree) from the element data at the leaves. The
     * tree topology is unchanged. This function is useful for quickly updating the tree when relative element sizes and
     * element-to-element proximity do not change much (e.g. under character articulations). Of course, this also requires that
     * the leaf elements (objects of type <tt>T</tt>) be pointers or other handles to external data that can be changed outside
     * this object.
     *
     * @param start The root of the subtree to update. If null, the entire tree is updated.
     *
     * @note Since this function is meaningful only if the elements have changed, you may need to call updateElements() before
     *   calling this.
     */
    virtual void updateNodeBounds(Node * start = nullptr)
    {
      if (!start)
      {
        if (!root) return;  // nothing to process
        start = root;
      }

      (void)updateNodeBoundsRecursive(start);

      if (accelerate_nn_queries && valid_acceleration_structure && acceleration_structure)
      {
        // Would be nicer to touch only the smallest subtree of points from updated elements, but this is too complicated for
        // now, so we'll just resample all the points

        samplePointsFromElements(acceleration_structure->numElements(), acceleration_structure->getElements());
        acceleration_structure->updateNodeBounds();
      }
    }

    /** Destructor. */
    ~BvhN() { clear(true); }

    /**
     * Get the maximum possible number of children of a node. It is possible no node in a tree actually has this many children.
     */
    static constexpr int maxDegree() { return MaxDegree; }

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
      transform_inverse.setIdentity();
      transform_inverse_transpose.setIdentity();
      invalidateBounds();

      if (valid_acceleration_structure)
        acceleration_structure->clearTransform();
    }

    /** Get the inverse of the transform applied to the BVH. */
    AffineTransformN<N, ScalarT> const & getTransformInv() const { return transform_inverse; }

    /** Get the transpose of the inverse of the linear part of the transform applied to the BVH. */
    Matrix<N, N, ScalarT> const & getLinearTransformInvTr() const { return transform_inverse_transpose; }

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
    bool empty() const { return num_elems <= 0; }

    /** Get the number of elements in the tree. The elements themselves can be obtained with getElements(). */
    intx numElements() const { return num_elems; }

    /** Get a pointer to an array of the elements in the tree. The number of elements can be obtained with numElements(). */
    T const * getElements() const { return &elems[0]; }

    /**
     * Get a pointer to an array of the elements in the tree. The number of elements can be obtained with numElements().
     *
     * @warning This non-const version of the function can be used to change the elements of the tree, which can invalidate the
     *  tree and necessitate either a call to updateNodeBounds() or a complete tree re-initialization. It should be used with
     *  extreme caution.
     */
    T * getElements() { return &elems[0]; }

  public:
    /**
     * Get the node corresponding to the root of the BVH. This function is provided so that users can implement their own
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

    /** Get the maximum subdivision depth (number of levels not counting the root) of the BVH. */
    intx maxDepth() const { return max_depth; }

    /** Get the maximum number of elements in each leaf of the BVH. */
    intx maxElementsPerLeaf() const { return max_elems_per_leaf; }

    /** Get a bounding volume for all the objects in the tree. */
    BoundingVolume const & getBounds() const
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
      alwaysAssertM(filter, "BvhN: Filter must be non-null");

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

    /**
     * Get the minimum distance between this structure and a query object, or a negative number if no such distance can be
     * computed, for instance if the structure is empty.
     */
    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    double distance(QueryT const & query, double dist_bound = -1, CompatibilityFunctorT compatibility = CompatibilityFunctorT())
           const
    {
      double result = -1;
      if (closestElement<MetricT>(query, dist_bound, compatibility, &result) >= 0)
        return result;
      else
        return -1;
    }

    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    intx closestElement(QueryT const & query, double dist_bound = -1,
                        CompatibilityFunctorT compatibility = CompatibilityFunctorT(),
                        double * dist = nullptr, VectorT * closest_point = nullptr) const
    {
      NeighborPair pair = closestPair<MetricT>(query, dist_bound, compatibility, closest_point != nullptr);

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
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param get_closest_points If true, the coordinates of the closest pair of points on the respective elements is computed
     *   and stored in the returned structure.
     *
     * @return Non-negative handles to the closest pair of elements in their respective objects, if such a pair was found. Else
     *   returns a pair of negative numbers.
     */
    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    NeighborPair closestPair(QueryT const & query, double dist_bound = -1,
                             CompatibilityFunctorT compatibility = CompatibilityFunctorT(),
                             bool get_closest_points = false) const
    {
      NeighborPair result(-1, -1);
      if (!root) { return result; }

      // Early pruning if the entire structure is too far away from the query
      BoundingVolume query_bounds;
      getObjectBounds(query, query_bounds);
      double mon_approx_dist_bound = (dist_bound >= 0 ? MetricT::computeMonotoneApprox(dist_bound) : -1);
      if (mon_approx_dist_bound >= 0)
      {
        double lower_bound = monotonePruningDistance<MetricT>(root, query, query_bounds);
        if (lower_bound >= 0 && lower_bound > mon_approx_dist_bound)
          return result;
      }

      // If acceleration is enabled, set an upper limit to the distance to the nearest object
      double accel_bound = accelerationBound<MetricT>(query, dist_bound, compatibility);
      if (accel_bound >= 0)
      {
        double fudge = std::max((double)(0.001 * getBoundsWorldSpace(*root).getExtent().norm()),  // FIXME: AABB-specific
                                (double)(100 * Math::eps<ScalarT>()));
        dist_bound = accel_bound + fudge;
      }

      ClosestPairFunctor<MetricT> functor(dist_bound, result);
      processNeighbors<MetricT>(root, query, query_bounds, functor, compatibility, get_closest_points,
                                /* query_index = */ 0, /* swap_query_and_target = */ false);

      return result;
    }

    /**
     * Get the k elements closest to a query object. The returned elements are placed in a set of bounded size (<tt>k</tt>). The
     * template type BoundedNeighborPairSetT should typically be BoundedSortedArray<NeighborPair> or
     * BoundedSortedArrayN<k, NeighborPair> if only a few neighbors are requested. The set is <b>not</b> cleared at the start.
     *
     * @note The auxiliary nearest-neighbor acceleration structure, even if present (enableNearestNeighborAcceleration()), is
     *   <b>NOT</b> used for pruning by this function.
     *
     * @param query Query object. BoundedTraitsN<QueryT, N, ScalarT> must be defined.
     * @param k_closest_pairs The <tt>k</tt> (or fewer) nearest neighbors are placed here.
     * @param dist_bound Upper bound on the distance between any pair of points considered. Ignored if negative.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param get_closest_points If true, the coordinates of the closest pair of points on each pair of neighboring elements is
     *   computed and stored in the returned pairs.
     *
     * @return The number of neighbors found, i.e. the size of \a k_closest_pairs, which may be less than its capacity
     *   <tt>k</tt>.
     */
    template <typename MetricT, typename QueryT, typename BoundedNeighborPairSetT,
              typename CompatibilityFunctorT = UniversalCompatibility>
    intx kClosestPairs(QueryT const & query, BoundedNeighborPairSetT & k_closest_pairs, double dist_bound = -1,
                       CompatibilityFunctorT compatibility = CompatibilityFunctorT(), bool get_closest_points = false) const
    {
      KClosestPairsFunctor<MetricT, BoundedNeighborPairSetT> functor(dist_bound, k_closest_pairs);
      processNeighbors<MetricT>(query, functor, compatibility, get_closest_points);

      return (intx)k_closest_pairs.size();
    }

    /**
     * Apply a functor to all elements neighboring a query object, until the functor returns true. The functor should provide
     * the following two member functions:
     * \code
     * bool allows(NeighborPair const &) const;
     * bool operator()(NeighborPair const &)
     * \endcode
     * The first function will be used to check if a potential pair of neighbors will be accepted by the functor or not. Note
     * that this pair may be <i>incomplete</i>, e.g. one or both source/target element indices, or the element separation, could
     * be negative. Any such invalid fields will be ignored. This feature is used for early pruning of elements and bounding
     * volumes.
     *
     * The second function will be passed information about every valid pair of elements that passes the <tt>allows()</tt>
     * function. If the query is itself a proximity query structure, the corresponding field of the NeighborPair will point to
     * the relevant element in the structure. Otherwise, the corresponding field will always contain the index 0 and refer to
     * the entire query object. This function will always receive <i>complete</i> NeighborPair objects (no negative indices or
     * distances), though the closest point positions will be uninitialized if \a get_closest_points is false.
     *
     * This is a generic function that can be used to mimic the behavior of closestPair() or kClosestPairs(), though those
     * functions can have individual optimizations (e.g. closestPair() uses the auxiliary acceleration structure).
     *
     * @note Since two different leaf volumes can contain the same overlapping element, the same pair of elements may be found
     *   twice and passed to <tt>FunctorT::operator()</tt>. <tt>FunctorT::allows()</tt> should check for repetitions if it wants
     *   to avoid this.
     * @note The auxiliary nearest-neighbor acceleration structure, even if present (enableNearestNeighborAcceleration()), is
     *   <b>NOT</b> used for pruning by this function.
     *
     * @param query Query object. BoundedTraitsN<QueryT, N, ScalarT> must be defined.
     * @param functor The functor that will be called for each pair of neighbors.
     * @param compatibility A functor that checks if two elements being considered as near neighbors are compatible with each
     *   other or not. It should have a <tt>bool operator()(q, e) const</tt> function that returns true if the two objects
     *   <tt>q</tt> and <tt>e</tt> are compatible with each other, where <tt>q</tt> is (i) an element of <tt>QueryT</tt> if the
     *   latter is a ProximityQueryStructureN, or (ii) \a query otherwise; and <tt>e</tt> is an element of this structure.
     * @param get_closest_points If true, the coordinates of the closest pair of points on each pair of neighboring elements is
     *   computed and passed to the functor.
     * @param query_index The supplied index is passed to the functor as the index of the query object. This is chiefly used for
     *   internal processing and the default value of 0 should normally be left as is.
     * @param swap_query_and_target If true, the indices and other properties of neighboring query and target objects are
     *   swapped when they are passed to the functor. This is chiefly used for internal processing and the default value of
     *   false should normally be left as is.
     *
     * If the functor returns true on any object, the search will terminate immediately (this is useful for searching for a
     * particular pair). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     */
    template <typename MetricT, typename QueryT, typename FunctorT, typename CompatibilityFunctorT = UniversalCompatibility>
    void processNeighbors(QueryT const & query, FunctorT functor, CompatibilityFunctorT compatibility = CompatibilityFunctorT(),
                          bool get_closest_points = false, intx query_index = 0, bool swap_query_and_target = false) const
    {
      if (!root) { return; }

      // Early pruning if the entire structure is too far away from the query
      BoundingVolume query_bounds;
      getObjectBounds(query, query_bounds);
      if (!processNeighborsAllows(functor, NeighborPair(-1, -1, monotonePruningDistance<MetricT>(root, query, query_bounds))))
        return;

      processNeighbors<MetricT>(root, query, query_bounds, functor, compatibility, get_closest_points,
                                query_index, swap_query_and_target);
    }

    template <typename IntersectionTesterT, typename RangeT>
    void rangeQuery(RangeT const & range, Array<T> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root) { const_cast<BvhN *>(this)->processRange<IntersectionTesterT>(range, RangeQueryFunctor(result)); }
    }

    template <typename IntersectionTesterT, typename RangeT>
    void rangeQueryIndices(RangeT const & range, Array<intx> & result, bool discard_prior_results = true) const
    {
      if (discard_prior_results) result.clear();
      if (root) { const_cast<BvhN *>(this)->processRange<IntersectionTesterT>(range, RangeQueryIndicesFunctor(result)); }
    }

    /**
     * Apply a functor to all elements intersecting a range, until the functor returns true. The functor should provide the
     * member function (or be a function pointer with the equivalent signature):
     * \code
     * bool operator()(intx index, T const & t)
     * \endcode
     * and will be passed the index of each element intersecting in the range as well as a handle to the element itself. If the
     * functor returns true on any element, the search will terminate immediately (this is useful for searching for a particular
     * element). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * The RangeT class should support intersection queries with BoundingVolume and T, and containment queries with
     * BoundingVolume.
     *
     * @return The index of the first element in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this element), else a negative value.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    intx processRange(RangeT const & range, FunctorT functor) const
    {
      return root ? const_cast<BvhN *>(this)->processRange<IntersectionTesterT, T const>(root, range, functor) : -1;
    }

    /**
     * Apply a functor to all elements intersecting a range, until the functor returns true. The functor should provide the
     * member function (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T [const] & t)
     * \endcode
     * and will be passed the index of each element intersecting the range as well as a handle to the element itself. If the
     * functor returns true on any element, the search will terminate immediately (this is useful for searching for a particular
     * element). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * The RangeT class should support intersection queries with BoundingVolume and T, and containment queries with
     * BoundingVolume.
     *
     * @return The index of the first element in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this element), else a negative value.
     */
    template <typename IntersectionTesterT, typename RangeT, typename FunctorT>
    intx processRange(RangeT const & range, FunctorT functor)
    {
      return root ? processRange<IntersectionTesterT, T>(root, range, functor) : -1;
    }

    template <typename RayIntersectionTesterT> bool rayIntersects(RayT const & ray, ScalarT max_time = -1) const
    {
      return rayIntersectionTime<RayIntersectionTesterT>(ray, max_time) >= 0;
    }

    template <typename RayIntersectionTesterT> ScalarT rayIntersectionTime(RayT const & ray, ScalarT max_time = -1) const
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
    RayStructureIntersectionT rayStructureIntersection(RayT const & ray, ScalarT max_time = -1) const
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
      BvhN const * tree;

      /** Constructor. Axis 0 = X, 1 = Y, 2 = Z. */
      ObjectLess(intx coord_, BvhN const * tree_) : coord(coord_), tree(tree_) {}

      /** Less-than operator, along the specified axis. */
      bool operator()(ElementIndex a, ElementIndex b)
      {
        // Compare object min coords
        return BoundedTraitsT::getLow(tree->elems[a], coord)
             < BoundedTraitsT::getLow(tree->elems[b], coord);
      }
    };

    // Allow the comparator unrestricted access to the BVH.
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
    /** Initialize a kd-tree. Assumes that the global list of elements <tt>elems</tt> has already been initialized. */
    template <typename InputIterator>
    void initKdTree(InputIterator begin, InputIterator end, intx max_depth_ = -1, intx max_elems_per_leaf_ = -1,
                    bool save_memory = false, bool deallocate_previous_memory = true)
    {
      static intx const DEFAULT_MAX_ELEMS_IN_LEAF = 5;
      max_elems_per_leaf = max_elems_per_leaf_ < 0 ? DEFAULT_MAX_ELEMS_IN_LEAF : max_elems_per_leaf_;

      // The fraction of elements held by the larger node at each split is 0.5
      static double const SPLIT_FRACTION = 0.5;
      intx est_depth = Math::binaryTreeDepth(num_elems, max_elems_per_leaf, SPLIT_FRACTION);
      max_depth = max_depth_;
      if (max_depth < 0)
        max_depth = est_depth;
      else if (max_depth < est_depth)
        est_depth = max_depth;

      // Each index is stored at most once at each level
      size_t BUFFER_SAFETY_MARGIN = 10;
      size_t index_buffer_capacity = num_elems + BUFFER_SAFETY_MARGIN;
      if (!save_memory)
        index_buffer_capacity *= (size_t)(1 + est_depth);  // reserve space for all levels at once

      if (deallocate_previous_memory || index_buffer_capacity > 1.3 * index_pool.getBufferCapacity())
        index_pool.init(index_buffer_capacity);

      // Assume a complete, balanced binary tree upto the estimated depth to guess the number of leaf nodes
      size_t node_buffer_capacity = (size_t)(1 << est_depth) + BUFFER_SAFETY_MARGIN;
      if (deallocate_previous_memory || node_buffer_capacity > 1.3 * node_pool.getBufferCapacity())
        node_pool.init(node_buffer_capacity);

      // Create the root node
      root = node_pool.alloc(1);
      root->init(0);
      num_nodes = 1;

      root->num_elems = num_elems;
      root->elems = index_pool.alloc(root->num_elems);

      BoundingVolume elem_bounds;
      for (size_t i = 0; i < (size_t)num_elems; ++i)
      {
        root->elems[i] = i;

        BoundedTraitsT::getBounds(elems[i], elem_bounds);
        root->bounds.merge(elem_bounds);
      }

      // Expand the bounding volume slightly to handle numerical error
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

        // Create a temporary pool for scratch storage
        IndexPool tmp_index_pool;
        tmp_index_pool.init(est_max_path_indices + BUFFER_SAFETY_MARGIN);

        createKdTree(root, true, &tmp_index_pool, &index_pool);
      }
      else
        createKdTree(root, false, &index_pool, nullptr);

      invalidateBounds();
    }

    /**
     * Recursively construct a (sub-)kd-tree. If \a save_memory is true, element indices for this subtree are assumed to be in a
     * block at the end of \a main_index_pool. They are subsequently moved to arrays associated with the leaves, allocated in
     * \a leaf_index_pool. They are deleted from \a main_index_pool when this function exits, thus deallocating index arrays
     * held by internal nodes. If \a save_memory is false, each node maintains its own list of all elements in its subtree.
     *
     * @note To use this function with a node accessed via getRoot(), you'll have to const_cast it to a non-const type first.
     */
    void createKdTree(Node * start, bool save_memory, IndexPool * main_index_pool, IndexPool * leaf_index_pool)
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
#define THEA_BVHN_KDTREE_SPLIT_LONGEST
#ifdef THEA_BVHN_KDTREE_SPLIT_LONGEST
      intx coord = Math::maxAxis(start->bounds.getExtent());  // split longest dimension (FIXME: AABB-specific)
#else
      intx coord = (intx)(start->depth % N);  // cycle between dimensions
#endif

      // Split elements into lower and upper halves
      size_t mid = start->num_elems / 2;
      std::nth_element(start->elems, start->elems + mid, start->elems + start->num_elems, ObjectLess(coord, this));

      // Create child nodes
      for (size_t i = 0; i < 2; ++i)
      {
        Node * child = start->children[i] = node_pool.alloc(1);
        child->init(start->depth + 1);
        num_nodes++;

        ElementIndex elems_begin, elems_end;
        if (i == 0)
        { elems_begin = 0; elems_end = start->num_elems - mid; }
        else
        { elems_begin = start->num_elems - mid; elems_end = start->num_elems; }

        child->elems = main_index_pool->alloc(elems_end - elems_begin);

        // Add the appropriate half of the elements to the child node
        BoundingVolume elem_bounds;
        bool first = true;
        for (ElementIndex i = elems_begin; i < elems_end; ++i)
        {
          ElementIndex index = start->elems[i];
          BoundedTraitsT::getBounds(elems[index], elem_bounds);

          child->elems[child->num_elems++] = index;

          if (first)
          {
            child->bounds = elem_bounds;
            first = false;
          }
          else
            child->bounds.merge(elem_bounds);
        }

        // Expand the bounding volumes slightly to handle numerical error
        child->bounds.scaleCentered(BOUNDS_EXPANSION_FACTOR);
      }

      // Recurse on the high child first, since its indices are at the end of the main index pool and can be freed first if
      // necessary
      createKdTree(start->children[1], save_memory, main_index_pool, leaf_index_pool);

      // Recurse on the low child next, if we are in memory-saving mode its indices are now the last valid entries in the main
      // index pool
      createKdTree(start->children[0], save_memory, main_index_pool, leaf_index_pool);

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
     * Check if an element passes all filters currently on the stack.
     *
     * @see pushFilter(), popFilter()
     */
    bool elementPassesFilters(T const & elem) const
    {
      if (filters.empty()) return true;  // early exit

      for (typename FilterStack::const_iterator fi = filters.begin(); fi != filters.end(); ++fi)
        if (!(*fi)->allows(elem))
          return false;

      return true;
    }

    /**
     * Recompute the bounding volume of a leaf node from its elements. Can be overriden in derived classes for additional speed
     * if possible.
     *
     * @param leaf A pointer to a leaf node. This is guaranteed to be non-null, and an actual leaf without children.
     *
     * @return The updated bounding volume of \a leaf, <b>without</b> any padding. The bounding volume actually stored in
     *   \a leaf <b>will</b> be padded.
     */
    virtual BoundingVolume updateLeafBounds(Node * leaf)
    {
      BoundingVolume unpadded_bounds, elem_bounds;
      for (ElementIndex i = 0; i < leaf->num_elems; ++i)
      {
        BoundedTraitsT::getBounds(elems[leaf->elems[i]], elem_bounds);
        unpadded_bounds.merge(elem_bounds);
      }

      // Expand the bounding volume slightly to handle numerical error
      leaf->bounds = unpadded_bounds.scaleCenteredCopy(BOUNDS_EXPANSION_FACTOR);

      return unpadded_bounds;
    }

  private:
    /**
     * Recursively applied helper function for updateNodeBounds().
     *
     * @return The updated bounding volume of \a start, <b>without</b> any padding. The bounding volume actually stored in
     *   \a start <b>will</b> be padded.
     */
    BoundingVolume updateNodeBoundsRecursive(Node * start)
    {
      if (start->isLeaf())
        return updateLeafBounds(start);
      else  // recurse
      {
        BoundingVolume unpadded_bounds;
        for (auto c : start->children)
          unpadded_bounds.merge(updateNodeBoundsRecursive(c));

        start->bounds = unpadded_bounds.scaleCenteredCopy(BOUNDS_EXPANSION_FACTOR);

        return unpadded_bounds;
      }
    }

  protected:
    /** Mark that the bounding volume requires an update. */
    void invalidateBounds()
    {
      valid_bounds = false;
    }

    /** Recompute the bounding volume if it has been invalidated. */
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
        bounds = BoundingVolume();
      }

      valid_bounds = true;
    }

    /** Get the bounding volume for an object, if it is bounded. */
    template < typename U, typename std::enable_if< IsBoundedN<U, N>::value, int >::type = 0 >
    static void getObjectBounds(U const & u, BoundingVolume & bounds)
    {
      BoundedTraitsN<U, N, ScalarT>::getBounds(u, bounds);
    }

    /** Returns a null bounding volume for unbounded objects. */
    template < typename U, typename std::enable_if< !IsBoundedN<U, N>::value, int >::type = 0 >
    static void getObjectBounds(U const & u, BoundingVolume & bounds)
    {
      bounds.setNull();
    }

    /**
     * Get a bounding volume for a node of the tree, in world space. For the bounding volume in object coordinates,
     * Node::getBounds() can be used.
     */
    BoundingVolume getBoundsWorldSpace(Node const & node) const
    {
      return TransformableBaseT::hasTransform()
           ? node.bounds.transformAndBound(TransformableBaseT::getTransform())
           : node.bounds;
    }

    /** Get a bounding volume for an element, in world space. */
    BoundingVolume getBoundsWorldSpace(T const & t) const
    {
      BoundingVolume tb; getObjectBounds(t, tb);
      return TransformableBaseT::hasTransform() ? tb.transformAndBound(TransformableBaseT::getTransform()) : tb;
    }

    /**
     * Get a lower bound on (the monotone approximation to) the distance between the bounding volume of a BVH node and a
     * bounded query object, or a negative value if no such lower bound can be calculated.
     */
    template < typename MetricT, typename QueryT, typename std::enable_if< IsBoundedN<QueryT, N>::value, int >::type = 0 >
    double monotonePruningDistance(Node const * node, QueryT const & query, BoundingVolume query_bounds) const
    {
      if (node && !query_bounds.isNull())
        return MetricT::template monotoneApproxDistance<N, ScalarT>(query_bounds, getBoundsWorldSpace(*node));
      else
        return -1;
    }

    /**
     * Get a lower bound on (the monotone approximation to) the distance between the bounding volume of a BVH node and an
     * unbounded query object, or a negative value if no such lower bound can be calculated. \a query_bounds is ignored.
     */
    template < typename MetricT, typename QueryT, typename std::enable_if< !IsBoundedN<QueryT, N>::value, int >::type = 0 >
    double monotonePruningDistance(Node const * node, QueryT const & query, BoundingVolume query_bounds) const
    {
      // Assume the following specialization exists
      return node ? MetricT::template monotoneApproxDistance<N, ScalarT>(query, getBoundsWorldSpace(*node)) : -1;
    }

  private:
    /** Wrap a node and a distance to it. */
    struct NodeDistance
    {
      Node * node;      ///< The wrapped node.
      double distance;  ///< The distance to the node.

      /** Default constructor. */
      NodeDistance() : node(nullptr), distance(-1) {}

      /** Initializing constructor. */
      NodeDistance(Node * node_, double distance_) : node(node_), distance(distance_) {}

      /** Less-than operator. All negative distances are considered equal, and greater than any positive value. */
      bool operator<(NodeDistance const & rhs) const { return BvhNInternal::distanceLessThan(distance, rhs.distance); }

    }; // struct NodeDistance

    /** Sorted array of distances to the children of a tree node. */
    typedef BoundedSortedArrayN<MaxDegree, NodeDistance> ChildDistanceArray;

  protected:
    /** Check whether the type Q is a compatible BvhN. This default implementation sets the <tt>value</tt> field to false. */
    template <typename Q, typename Enable = void>
    struct IsCompatibleBvhN
    {
      static bool const value = false;
    };

    /**
     * Check whether the type Q is a compatible BvhN. This specialization is used when Q is actually a BvhN, and it sets the
     * <tt>value</tt> field to true.
     */
    template <typename Q>
    struct IsCompatibleBvhN< Q, typename std::enable_if< std::is_base_of< BvhN< typename Q::Element, N, ScalarT,
                                                                                typename Q::NodeAttribute,
                                                                                Q::maxDegree(),
                                                                                typename Q::BoundingVolume >,
                                                                          Q >::value >::type >
    {
      static bool const value = true;
    };

    /**
     * Recursively call \a functor for elements neighboring a query object in the subtree rooted at \a start, when the query
     * object is NOT itself a compatible BvhN. If \a get_closest_points is true, the positions of the closest pair of points
     * will be stored with each pair, not just the distance between them.
     */
    template < typename MetricT, typename QueryT, typename FunctorT, typename CompatibilityFunctorT = UniversalCompatibility,
               typename std::enable_if< !IsCompatibleBvhN<QueryT>::value, int >::type = 0 >
    void processNeighbors(Node const * start, QueryT const & query, BoundingVolume const & query_bounds, FunctorT functor,
                          CompatibilityFunctorT compatibility, bool get_closest_points,
                          intx query_index, bool swap_query_and_target) const
    {
      if (start->isLeaf())
      {
        processNeighborsLeaf<MetricT>(start, query, functor, compatibility, get_closest_points,
                                      query_index, swap_query_and_target);
      }
      else  // not leaf
      {
        // Sort the children by increasing distance to their bounding volumes (optimize for point queries?)
        ChildDistanceArray c;
        for (auto n : start->children)
          if (n) { c.insert(NodeDistance(n, monotonePruningDistance<MetricT>(n, query, query_bounds))); }

        for (size_t i = 0; i < c.size(); ++i)
        {
          if (processNeighborsAllows(functor, NeighborPair(-1, -1, c[i].distance)))
          {
            processNeighbors<MetricT>(c[i].node, query, query_bounds, functor, compatibility, get_closest_points,
                                      query_index, swap_query_and_target);
          }
        }
      }
    }

    /**
     * Recursively call \a functor for elements neighboring a query object in the subtree rooted at \a start, when the query
     * object IS a compatible BvhN. If \a get_closest_points is true, the positions of the closest pair of points will be stored
     * with each pair, not just the distance between them.
     */
    template < typename MetricT, typename QueryT, typename FunctorT, typename CompatibilityFunctorT = UniversalCompatibility,
               typename std::enable_if< IsCompatibleBvhN<QueryT>::value, int >::type = 0 >
    void processNeighbors(Node const * start, QueryT const & query, BoundingVolume const & query_bounds, FunctorT functor,
                          CompatibilityFunctorT compatibility, bool get_closest_points,
                          intx query_index, bool swap_query_and_target) const
    {
      processNeighborsBvh<MetricT>(start, query, query.getRoot(), functor, compatibility, get_closest_points,
                                   swap_query_and_target);
    }

  private:
    /**
     * Recursively call \a functor for neighboring elements between subtrees rooted at \a this_node (of this BvhN) and
     * \a query_node (of a different but compatible BvhN). If \a get_closest_points is true, the positions of the closest pair
     * of points will be stored with each pair, not just the distance between them.
     */
    template <typename MetricT, typename QueryBvhT, typename QueryNodeT, typename FunctorT, typename CompatibilityFunctorT>
    void processNeighborsBvh(Node const * this_node, QueryBvhT const & query, QueryNodeT const * query_node, FunctorT functor,
                             CompatibilityFunctorT compatibility, bool get_closest_points, bool swap_query_and_target) const
    {
      if (!this_node || !query_node) { return; }

      if (this_node->isLeaf() && query_node->isLeaf())  // both are leaves
      {
        // Loop over the elements in query_node and compare each one to the bounding box of this_node for pruning
        auto const * qelems = query.getElements();
        BoundingVolume qe_bounds_compatible;  // QueryBvhT::BoundingVolume may differ from BoundingVolume
        for (auto qei = query_node->elementIndicesBegin(); qei != query_node->elementIndicesEnd(); ++qei)
        {
          auto const & qe = qelems[*qei];
          if (!query.elementPassesFilters(qe))
            continue;

          auto qe_bounds = query.getBoundsWorldSpace(qe);
          getObjectBounds(qe_bounds, qe_bounds_compatible);

          auto mad = monotonePruningDistance<MetricT>(this_node, qe, qe_bounds_compatible);
          if (processNeighborsAllows(functor, NeighborPair(-1, -1, mad)))
            processNeighborsLeaf<MetricT>(this_node, qe, functor, compatibility, get_closest_points, *qei,
                                          swap_query_and_target);
        }
      }
      else  // at most one is a leaf
      {
        bool recurse_into_this = false;
        BoundingVolume query_bounds_compatible;  // QueryBvhT::BoundingVolume may differ from BoundingVolume
        typename QueryBvhT::BoundingVolume this_bounds_compatible;

        if (!this_node->isLeaf() && !query_node->isLeaf())  // neither are leaves
        {
          // Find the bigger node and recurse into it
          auto this_bounds = getBoundsWorldSpace(*this_node);
          auto query_bounds = query.getBoundsWorldSpace(*query_node);
          getObjectBounds(query_bounds, query_bounds_compatible);

          // FIXME: AABB-specific
          recurse_into_this = (this_bounds.getExtent().squaredNorm() > query_bounds_compatible.getExtent().squaredNorm());

          if (!recurse_into_this)
            query.getObjectBounds(this_bounds, this_bounds_compatible);
        }
        else  // exactly one is a leaf
        {
          recurse_into_this = !this_node->isLeaf();

          if (recurse_into_this)
            getObjectBounds(query.getBoundsWorldSpace(*query_node), query_bounds_compatible);
          else
            query.getObjectBounds(getBoundsWorldSpace(*this_node), this_bounds_compatible);
        }

        if (recurse_into_this)
        {
          processNeighborsBvhRecurse<MetricT>(this_node, query, query_node, query_bounds_compatible, functor, compatibility,
                                              get_closest_points, swap_query_and_target);
        }
        else
        {
          auto swapped_compatibility = BvhNInternal::makeSwappedCompatibility(compatibility);

          query.template processNeighborsBvhRecurse<MetricT>(query_node, *this, this_node, this_bounds_compatible, functor,
                                                             swapped_compatibility, get_closest_points, !swap_query_and_target);
        }
      }
    }

    /** Helper function for processNeighborsBvh that does one step of recursion into the children of \a this_node. */
    template <typename MetricT, typename QueryBvhT, typename QueryNodeT, typename FunctorT, typename CompatibilityFunctorT>
    void processNeighborsBvhRecurse(Node const * this_node, QueryBvhT const & query, QueryNodeT const * query_node,
                                    BoundingVolume const & query_node_bounds, FunctorT functor,
                                    CompatibilityFunctorT compatibility, bool get_closest_points,
                                    bool swap_query_and_target) const
    {
      ChildDistanceArray c;
      for (auto n : this_node->children)
        if (n)
        {
          // query_node_bounds twice is intentional
          c.insert(NodeDistance(n, monotonePruningDistance<MetricT>(n, query_node_bounds, query_node_bounds)));
        }

      for (size_t i = 0; i < c.size(); ++i)
      {
        if (processNeighborsAllows(functor, NeighborPair(-1, -1, c[i].distance)))
        {
          processNeighborsBvh<MetricT>(c[i].node, query, query_node, functor, compatibility, get_closest_points,
                                       swap_query_and_target);
        }
      }
    }

    /** Find elements in a leaf node neighboring another object, when the latter IS a proximity query structure. */
    template < typename MetricT, typename QueryT, typename FunctorT, typename CompatibilityFunctorT,
               typename std::enable_if< std::is_base_of<ProximityQueryBaseT, QueryT>::value, int >::type = 0 >
    void processNeighborsLeaf(Node const * leaf, QueryT const & query, FunctorT functor, CompatibilityFunctorT compatibility,
                              bool get_closest_points, intx query_index, bool swap_query_and_target) const
    {
      auto swapped_compatibility = BvhNInternal::makeSwappedCompatibility(compatibility);

      for (size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        if (TransformableBaseT::hasTransform())
          query.template processNeighbors<MetricT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                   functor, swapped_compatibility, get_closest_points, (intx)index, true);
        else
          query.template processNeighbors<MetricT>(elem, functor, swapped_compatibility, get_closest_points, (intx)index, true);
      }
    }

    /** Find elements in a leaf node neighboring another object, when the latter is NOT a proximity query structure. */
    template < typename MetricT, typename QueryT, typename FunctorT, typename CompatibilityFunctorT,
               typename std::enable_if< !std::is_base_of<ProximityQueryBaseT, QueryT>::value, int >::type = 0 >
    void processNeighborsLeaf(Node const * leaf, QueryT const & query, FunctorT functor, CompatibilityFunctorT compatibility,
                              bool get_closest_points, intx query_index, bool swap_query_and_target) const
    {
      VectorT qp, tp;
      double mad;
      for (size_t i = 0; i < leaf->num_elems; ++i)
      {
        ElementIndex index = leaf->elems[i];
        Element const & elem = elems[index];

        if (!elementPassesFilters(elem))
          continue;

        // Check if the element is already in the set of neighbors or not
        NeighborPair pair = (swap_query_and_target ? NeighborPair((intx)index, query_index)
                                                   : NeighborPair(query_index, (intx)index));
        if (!processNeighborsAllows(functor, pair))  // same pair may be found twice, so functor may want to ignore repetitions
          continue;

        if (!compatibility(query, elem))
          continue;

        if (TransformableBaseT::hasTransform())
        {
          mad = MetricT::template closestPoints<N, ScalarT>(makeTransformedObject(&elem, &TransformableBaseT::getTransform()),
                                                            query, tp, qp);
        }
        else
          mad = MetricT::template closestPoints<N, ScalarT>(elem, query, tp, qp);

        if (processNeighborsAllows(functor, NeighborPair(-1, -1, mad)))  // no need to check indices again
        {
          pair.setMonotoneApproxDistance(mad);

          if (get_closest_points)
          {
            if (swap_query_and_target)
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

          functor(pair);
        }
      }
    }

    /** Helper function to check if a functor passed to processNeighbors allows a pair. */
    template <typename FunctorT>
    static bool processNeighborsAllows(FunctorT functor, NeighborPair pair)
    { return functor.allows(pair); }

    // Specialization when the functor has been wrapped with std::ref
    template <typename FunctorT>
    static bool processNeighborsAllows(std::reference_wrapper<FunctorT> functor, NeighborPair pair)
    { return functor.get().allows(pair); }

  protected:
    /**
     * Apply a functor to all elements of a subtree intersecting a range, until the functor returns true. The functor should
     * provide the member function (or be a function pointer with the equivalent signature)
     * \code
     * bool operator()(intx index, T & t)
     * \endcode
     * and will be passed the index of each object intersecting the range as well as a handle to the object itself. If the
     * functor returns true on any object, the search will terminate immediately (this is useful for searching for a particular
     * object). To pass a functor by reference, wrap it in <tt>std::ref</tt>.
     *
     * The RangeT class should support intersection queries with BoundingVolume and T, and containment queries with
     * BoundingVolume.
     *
     * @return The index of the first object in the range for which the functor evaluated to true (the search stopped
     *   immediately after processing this object), else a negative value.
     */
    template <typename IntersectionTesterT, typename FunctorArgT, typename RangeT, typename FunctorT>
    intx processRange(Node const * start, RangeT const & range, FunctorT functor)
    {
      // Early exit if the range and node are disjoint
      BoundingVolume tr_start_bounds = getBoundsWorldSpace(*start);
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
      else if (start->isLeaf())
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
        for (auto c : start->children)
        {
          intx index = processRange<IntersectionTesterT, FunctorArgT>(c, range, functor);
          if (index >= 0) return index;
        }
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

  protected:
    /** Get the time taken for a ray to hit the nearest object in a node, in the forward direction. */
    template <typename RayIntersectionTesterT>
    ScalarT rayIntersectionTime(Node const * start, RayT const & ray, ScalarT max_time) const
    {
      if (start->isLeaf())
      {
        ScalarT best_time = max_time;
        bool found = false;
        for (size_t i = 0; i < start->num_elems; ++i)
        {
          ElementIndex index = start->elems[i];
          Element const & elem = elems[index];

          if (!elementPassesFilters(elem))
            continue;

          ScalarT time = RayIntersectionTesterT::template rayIntersectionTime<N, ScalarT>(ray, elem, best_time);
          if (BvhNInternal::distanceLessThan(time, best_time))
          {
            best_time = time;
            found = true;
          }
        }

        return found ? best_time : -1;
      }
      else  // not leaf
      {
        // Sort the children by increasing hit times to their bounding volumes
        ChildDistanceArray c;
        for (auto n : start->children)
          if (n) { c.insert(NodeDistance(n, n->bounds.rayIntersectionTime(ray, max_time))); }

        ScalarT best_time = max_time;
        bool found = false;
        for (size_t i = 0; i < c.size(); ++i)
        {
          if (BvhNInternal::distanceLessThan(c[i].distance, best_time))
          {
            ScalarT time = rayIntersectionTime<RayIntersectionTesterT>(c[i].node, ray, best_time);
            if (BvhNInternal::distanceLessThan(time, best_time))
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
    RayStructureIntersectionT rayStructureIntersection(Node const * start, RayT const & ray, ScalarT max_time) const
    {
      if (start->isLeaf())
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
          if (BvhNInternal::distanceLessThan(isec.getTime(), best_isec.getTime()))
          {
            best_isec = RayStructureIntersectionT(isec, (intx)index);
            found = true;
          }
        }

        return found ? best_isec : RayStructureIntersectionT(-1);
      }
      else  // not leaf
      {
        // Sort the children by increasing hit times to their bounding volumes
        ChildDistanceArray c;
        for (auto n : start->children)
          if (n) { c.insert(NodeDistance(n, n->bounds.rayIntersectionTime(ray, max_time))); }

        RayStructureIntersectionT best_isec(max_time);
        bool found = false;
        for (size_t i = 0; i < c.size(); ++i)
        {
          if (BvhNInternal::distanceLessThan(c[i].distance, best_isec.getTime()))
          {
            RayStructureIntersectionT isec = rayStructureIntersection<RayIntersectionTesterT>(c[i].node, ray,
                                                                                              best_isec.getTime());
            if (BvhNInternal::distanceLessThan(isec.getTime(), best_isec.getTime()))
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
     * closest to the center of the bounding volume each element: subclasses can override this with faster sampling methods.
     *
     * @note This version finds nearest neighbors on elements using the L2 distance. Hence MetricL2 must support NN queries with
     *   the element type T.
     */
    virtual void samplePointsFromElements(intx num_samples, ElementSample * samples) const
    {
      VectorT src_cp;
      for (intx i = 0; i < num_samples; ++i)
      {
        T const * elem = samples[i].element;
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

    /**
     * Get an upper bound on the distance to a query object, using the acceleration structure if it exists. If a non-negative
     * distance upper bound has been specified, the function just returns that value instead of using the acceleration
     * structure.
     */
    template <typename MetricT, typename QueryT, typename CompatibilityFunctorT = UniversalCompatibility>
    double accelerationBound(QueryT const & query, double dist_bound,
                             CompatibilityFunctorT compatibility = CompatibilityFunctorT()) const
    {
      if (dist_bound >= 0)
        return dist_bound;

      NearestNeighborAccelerationStructure const * accel = getNearestNeighborAccelerationStructure<MetricT>();
      return accel
             ? accel->template distance<MetricT>(query, dist_bound,
                                                 BvhNInternal::SampleCompatibility<CompatibilityFunctorT>(compatibility))
             : -1;
    }

    /**
     * Get an upper bound on the distance to a query BVH, using the acceleration structures of both the query and of this
     * object if they exist.
     */
    template <typename MetricT, typename E, typename S, typename A, int D, typename V,
              typename CompatibilityFunctorT = UniversalCompatibility>
    double accelerationBound(BvhN<E, N, S, A, D, V> const & query, double dist_bound,
                             CompatibilityFunctorT compatibility = CompatibilityFunctorT()) const
    {
      NearestNeighborAccelerationStructure const * accel = getNearestNeighborAccelerationStructure<MetricT>();
      if (accel)
      {
        if (query.hasNearestNeighborAcceleration())
        {
          typename BvhN<E, N, S, A, D, V>::NearestNeighborAccelerationStructure const * query_accel
              = query.template getNearestNeighborAccelerationStructure<MetricT>();

          if (query_accel)
            return accel->template distance<MetricT>(*query_accel, dist_bound,
                                                     BvhNInternal::SampleCompatibility<CompatibilityFunctorT>(compatibility));
        }

        return accel->template distance<MetricT>(query, dist_bound,
                                                 BvhNInternal::SampleCompatibility<CompatibilityFunctorT>(compatibility));
      }
      else
      {
        if (query.hasNearestNeighborAcceleration())
        {
          typename BvhN<E, N, S, A, D, V>::NearestNeighborAccelerationStructure const * query_accel
              = query.template getNearestNeighborAccelerationStructure<MetricT>();

          if (query_accel)
            return distance<MetricT>(*query_accel, dist_bound, compatibility);
        }

        return -1;
      }
    }

  private:
    template <typename E, int M, typename S, typename A, int D, typename V> friend class BvhN;

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
    mutable BoundingVolume bounds;

  protected:
    /** A factor by which node bounding volumes are upscaled to add a little padding for protection against numerical error. */
    static ScalarT const BOUNDS_EXPANSION_FACTOR;

}; // class BvhN

// Static variables
template <typename T, int N, typename S, typename A, int D, typename V>
S const BvhN<T, N, S, A, D, V>::BOUNDS_EXPANSION_FACTOR = static_cast<S>(1.05);

// Mark the BVH and its (public) descendants as bounded objects. The default BoundedTraitsN implementation is good enough.
template <typename T, int N>
class IsBoundedN< T, N, typename std::enable_if< std::is_base_of< BvhN< typename T::Element, N,
                                                                        typename T::VectorT::value_type,
                                                                        typename T::NodeAttribute,
                                                                        T::maxDegree(),
                                                                        typename T::BoundingVolume >,
                                                                  T >::value >::type >
{
  public:
    static bool const value = true;
};

} // namespace Algorithms
} // namespace Thea

#endif
