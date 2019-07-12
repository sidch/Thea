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

#ifndef __Thea_Options_hpp__
#define __Thea_Options_hpp__

#include "Common.hpp"
#include "Map.hpp"
#include <boost/any.hpp>

namespace Thea {

/** Abstract interface for key-value options, suitable for passing across DLL boundaries. */
class THEA_API AbstractOptions
{
  public:
    /** Destructor. */
    virtual ~AbstractOptions() {}

    /** Check if an option has been set. */
    virtual int8 hasOption(char const * option_name) const = 0;

    /** Set the value of an integer option. */
    virtual void setInteger(char const * option_name, int64 value) = 0;

    /** Set the value of a floating-point option. */
    virtual void setFloat(char const * option_name, float64 value) = 0;

    /** Set the value of a string option. */
    virtual void setString(char const * option_name, char const * value) = 0;

    /**
     * Get the value of an integer option. If the option has not been set, the default value specified by the last parameter is
     * returned.
     */
    virtual int64 getInteger(char const * option_name, int64 default_value) const = 0;

    /**
     * Get the value of a floating-point option. If the option has not been set, the default value specified by the last
     * parameter is returned.
     */
    virtual float64 getFloat(char const * option_name, float64 default_value) const = 0;

    /**
     * Get the value of a string option. If the option has not been set, the default value specified by the last parameter is
     * returned.
     */
    virtual char const * getString(char const * option_name, char const * default_value) const = 0;

}; // class AbstractOptions

/**
 * A set of options, specified as key-value pairs. Supports a much more general set of value types than the simpler abstract
 * interface it implements.
 */
class THEA_API Options : public AbstractOptions, private Map<std::string, boost::any>
{
  public:
    /** Destructor. */
    ~Options() {}

    // Functions from abstract interface
    int8 hasOption(char const * option_name) const { return find(option_name) != end(); }
    void setInteger(char const * option_name, int64 value) { set<int64>(option_name, value); }
    void setFloat(char const * option_name, float64 value) { set<float64>(option_name, value); }
    void setString(char const * option_name, char const * value) { set<std::string>(option_name, std::string(value)); }
    int64 getInteger(char const * option_name, int64 default_value) const { return get<int64>(option_name, default_value); }
    float64 getFloat(char const * option_name, float64 default_value) const { return get<float64>(option_name, default_value); }
    char const * getString(char const * option_name, char const * default_value) const;

    /** Set the value of an option. */
    template <typename T> void set(char const * option_name, T const & value) { (*this)[option_name] = value; }

    /** Set the value of an option from a C-style string. It is stored and can be read back as a <code>std::string</code>. */
    void set(char const * option_name, char const * value) { (*this)[option_name] = std::string(value); }

    /**
     * Get the value of an option. If the option has not been set, the default value specified by the last parameter is
     * returned.
     */
    template <typename T> T get(char const * option_name, T const & default_value) const
    {
      const_iterator existing = find(option_name);
      if (existing != end())
        return boost::any_cast<T>(existing->second);
      else
        return default_value;
    }

}; // class Options

/**
 * Specialization of get() for C-style strings. In this particular case, a pointer to the C array underlying the held
 * <code>std::string</code> is returned (which may be invalidated by further operations on this Options object). If the
 * option has not been set, the default value specified by the last parameter is returned.
 */
template <>
inline char const *
Options::get<char const *>(char const * option_name, char const * const & default_value) const
{
  const_iterator existing = find(option_name);
  if (existing != end())
  {
    // any_cast with reference type returns handle to held value
    return boost::any_cast<std::string const &>(existing->second).c_str();
  }
  else
    return default_value;
}

// Has to be placed after the specialization of get<char const *>() above.
inline char const *
Options::getString(char const * option_name, char const * default_value) const
{
  return get<char const *>(option_name, default_value);
}

} // namespace Thea

#endif
