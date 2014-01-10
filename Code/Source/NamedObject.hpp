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

#ifndef __Thea_NamedObject_hpp__
#define __Thea_NamedObject_hpp__

#include "Common.hpp"

namespace Thea {

/** Interface for an object that has a name. */
class THEA_API AbstractNamedObject
{
  public:
    /** Destructor. */
    virtual ~AbstractNamedObject() {}

    /** Get the name of the object. */
    virtual char const * getName() const = 0;

}; // class AbstractNamedObject

/**
 * An object wrapping a name string.
 *
 * <b>IMPORTANT:</b> Always inherit virtually from this class, i.e. use:
 *
 * \code
 * class MyClass : public virtual NamedObject
 * {
 *   ...
 * };
 * \endcode
 *
 * This also means that a derived class in an inheritance hierarchy involving NamedObject must explicitly call the NamedObject
 * constructor if it wants to initialize the name.
 */
class THEA_API NamedObject : public AbstractNamedObject
{
  public:
    THEA_DEF_POINTER_TYPES(NamedObject, shared_ptr, weak_ptr)

    /** Destructor. */
    virtual ~NamedObject() = 0;

    /** Default constructor. */
    NamedObject() {}

    /** Initializing constructor. */
    NamedObject(std::string const & name_) : name(name_) {}

    /** Copy constructor. */
    NamedObject(NamedObject const & src): name(src.name) {}

    /** Assignment operator. */
    NamedObject & operator=(NamedObject const & src) { name = src.name; return *this; }

    /** Get the name of the object. */
    char const * getName() const { return name.c_str(); }

    /** Set the name of the object. */
    void setName(std::string const & name_) { name = name_; }

  protected:
    /** Access the name string directly, for efficiency. */
    std::string const & getNameStr() const { return name; }

  private:
    std::string name;

}; // class NamedObject

// Inline functions
inline
NamedObject::~NamedObject()
{
  // Pure virtual destructor should have a body
  //
  // http://www.linuxtopia.org/online_books/programming_books/thinking_in_c++/Chapter15_024.html
}

/** Write a short description of a named object to an output stream. */
inline std::ostream &
operator<<(std::ostream & os, NamedObject const & obj)
{
  std::string cls = getClass(obj);

  // Remove template parameters
  size_t last = cls.find('<');

  // Remove namespace qualifiers
  size_t first = cls.rfind(':', last);

  if (last  == std::string::npos) last  = cls.length();
  first = (first == std::string::npos ? 0 : first + 1);

  return os << cls.substr(first, last - first) << " '" << obj.getName() << '\'';
}

} // namespace Thea

THEA_DECL_EXTERN_SMART_POINTERS(Thea::NamedObject)

#endif
