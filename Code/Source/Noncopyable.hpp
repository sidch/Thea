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

#ifndef __Thea_Noncopyable_hpp__
#define __Thea_Noncopyable_hpp__

#include "Platform.hpp"

namespace Thea {

/**
 * A base class for objects that should never be copied. This is achieved by declaring the copy constructor and assignment
 * operator as private members. <b>Never ever</b> try to refer to an object of a derived class using a Noncopyable pointer or
 * reference (in any case this seems semantically weird) -- to ensure this class has zero runtime overhead, the destructor is
 * <b>not virtual</b>.
 */
class THEA_API Noncopyable
{
  protected:
    /** Constructor. */
    Noncopyable() {}

    /** Destructor. */
    ~Noncopyable() {}

  private:
    /**
     * Hidden copy constructor. No body provided since this should never be accessible -- if a linker error occurs then
     * something is seriously wrong.
     */
    THEA_DLL_LOCAL Noncopyable(const Noncopyable &);

    /**
     * Hidden assignment operator. No body provided since this should never be accessible -- if a linker error occurs then
     * something is seriously wrong.
     */
    THEA_DLL_LOCAL Noncopyable const & operator=(Noncopyable const &);

}; // class Noncopyable

} // namespace Thea

#endif
