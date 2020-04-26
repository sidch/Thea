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

#ifndef __Thea_Algorithms_Filter_hpp__
#define __Thea_Algorithms_Filter_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Algorithms {

/**
 * Interface for all object filters. A filter has an allows() function that takes an object argument and returns a boolean value
 * indicating whether the filter allows the object to pass through or not.
 */
template <typename T>
class Filter
{
  public:
    THEA_DECL_SMART_POINTERS(Filter)

    /** Destructor. */
    virtual ~Filter() {}

    /**
     * Check if the filter allows an object to pass through or not.
     *
     * @return True if the object \a t is allowed through, false if it is blocked.
     */
    virtual bool allows(T const & t) const = 0;

}; // class Filter

/** A filter that allows everything to pass through. */
template <typename T>
class AlwaysPassFilter : public Filter<T>
{
  public:
    bool allows(T const & t) const { return true; }

}; // class AlwaysPassFilter

/** A filter that allows nothing to pass through. */
template <typename T>
class AlwaysBlockFilter : public Filter<T>
{
  public:
    bool allows(T const & t) const { return false; }

}; // class AlwaysPassFilter

} // namespace Algorithms
} // namespace Thea

#endif
