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
// First version: 2014
//
//============================================================================

#ifndef __Thea_Algorithms_SamplePoint3_hpp__
#define __Thea_Algorithms_SamplePoint3_hpp__

#include "../Common.hpp"
#include "NormalTraitsN.hpp"
#include "PointTraitsN.hpp"
#include "../BoundedSortedArray.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Algorithms {

// Forward declarations
class SamplePoint3;

namespace SamplePoint3Internal {

/** A sample neighboring another sample. */
class NeighboringSample
{
  public:
    /** Constructor. */
    NeighboringSample(SamplePoint3 * sample_ = nullptr, Real separation_ = 0) : sample(sample_), separation(separation_) {}

    /** Get a reference to this sample. */
    SamplePoint3 * getSample() const { return sample; }

    /** Set the reference to the neighboring sample. */
    void setSample(SamplePoint3 * sample_) { sample = sample_; }

    /** Get the separation of this sample from the source sample. */
    Real getSeparation() const { return separation; }

    /** Set the separation of this sample from the source sample. */
    void setSeparation(Real separation_) { separation = separation_; }

    /** Compare neighbors by their separation from the source. */
    bool operator<(NeighboringSample const & rhs) const
    {
      return separation < rhs.separation || (separation == rhs.separation && std::less<SamplePoint3 *>()(sample, rhs.sample));
    }

  private:
    SamplePoint3 * sample;  ///< Reference to this sample.
    Real separation;         ///< Separation from the source sample.

}; // class NeighboringSample

} // namespace SamplePoint3Internal
} // namespace Algorithms
} // namespace Thea

namespace std {

template <> struct is_trivially_copyable<Thea::Algorithms::SamplePoint3Internal::NeighboringSample> : public std::true_type {};

} // namespace std

namespace Thea {
namespace Algorithms {

/**
 * A class used to store 3D sample points with adjacency information. Encapsulates a sample position, normal and list of
 * neighboring samples.
 */
class SamplePoint3
{
  public:
    typedef SamplePoint3Internal::NeighboringSample Neighbor;  ///< A neighboring sample.
    typedef BoundedSortedArray<Neighbor> NeighborSet;  ///< Set of neighbors, sorted by increasing distance from this sample.

    /** Default constructor. Does not initialize anything. */
    SamplePoint3() {}

    /** Construct with an index and maximum number of neighbors. */
    SamplePoint3(intx index_, intx max_nbrs) : index(index_), nbrs((size_t)max_nbrs) {}

    /** Get the sample index. */
    intx getIndex() const { return index; }

    /** Set the sample index. */
    void setIndex(intx index_) { index = index_; }

    /** Get the position of the sample. */
    Vector3 const & getPosition() const { return p; }

    /** Set the position of the sample. */
    void setPosition(Vector3 const & p_) { p = p_; }

    /** Get the surface normal at the sample position. */
    Vector3 const & getNormal() const { return n; }

    /** Set the surface normal at the sample position. */
    void setNormal(Vector3 const & n_) { n = n_; }

    /** Get the number of neighboring samples. */
    intx numNeighbors() const { return (intx)nbrs.size(); }

    /** Get the set of neighboring samples. */
    NeighborSet const & getNeighbors() const { return nbrs; }

    /** Get a non-const reference to the set of neighboring samples. */
    NeighborSet & getNeighbors() { return nbrs; }

    /** Get the neighboring sample with a given index. */
    Neighbor const & getNeighbor(intx i) const { return nbrs[(size_t)i]; }

    /**
     * Add a new neighboring sample.
     *
     * @return True if the neighbor was successfully added, false otherwise (e.g. if the sample already has enough neighbors at
     *   nearer distances).
     */
    bool addNeighbor(Neighbor const & nbr) { return nbrs.insert(nbr) < nbrs.getCapacity(); }

    /** Clear the set of neighbors of this sample. */
    void clearNeighbors() { nbrs.clear(); }

  private:
    intx index;        ///< Index of the sample.
    Vector3 p;         ///< Sample position.
    Vector3 n;         ///< Sample normal.
    NeighborSet nbrs;  ///< Set of neighboring samples.

}; // class SamplePoint3

template <>
class IsPointN<SamplePoint3, 3>
{
  public:
    static bool const value = true;
};

template <>
class PointTraitsN<SamplePoint3, 3>
{
  public:
    static Vector3 const & getPosition(SamplePoint3 const & sample) { return sample.getPosition(); }
};

template <>
class HasNormalN<SamplePoint3, 3>
{
  public:
    static bool const value = true;
};

template <>
class NormalTraitsN<SamplePoint3, 3>
{
  public:
    static Vector3 const & getNormal(SamplePoint3 const & sample) { return sample.getNormal(); }
};

} // namespace Algorithms
} // namespace Thea

#endif
