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

#ifndef __Thea_Options_hpp__
#define __Thea_Options_hpp__

#include "Common.hpp"
#include "Map.hpp"
#include <boost/any.hpp>

namespace Thea {

/** Interface for key-value options, suitable for passing across DLL boundaries. */
class THEA_API IOptions
{
  public:
    /** Destructor. */
    virtual ~IOptions() {}

    /** Check if an option has been set. */
    virtual int8 THEA_ICALL hasOption(char const * option_name) const = 0;

    /** Set the value of an integer option. */
    virtual void THEA_ICALL setInteger(char const * option_name, int64 value) = 0;

    /** Set the value of a floating-point option. */
    virtual void THEA_ICALL setFloat(char const * option_name, float64 value) = 0;

    /** Set the value of a string option. */
    virtual void THEA_ICALL setString(char const * option_name, char const * value) = 0;

    /**
     * Get the value of an integer option. If the option has not been set, the default value specified by the last parameter is
     * returned.
     */
    virtual int64 THEA_ICALL getInteger(char const * option_name, int64 default_value) const = 0;

    /**
     * Get the value of a floating-point option. If the option has not been set, the default value specified by the last
     * parameter is returned.
     */
    virtual float64 THEA_ICALL getFloat(char const * option_name, float64 default_value) const = 0;

    /**
     * Get the value of a string option. If the option has not been set, the default value specified by the last parameter is
     * returned.
     */
    virtual char const * THEA_ICALL getString(char const * option_name, char const * default_value) const = 0;

}; // class IOptions

/**
 * A set of options, specified as key-value pairs. Supports a much more general set of value types than the simpler abstract
 * interface it implements.
 */
class THEA_API Options : public virtual IOptions, private Map<std::string, boost::any>
{
  public:
    /** Destructor. */
    ~Options() {}

    // Functions from abstract interface
    int8 THEA_ICALL hasOption(char const * option_name) const { return find(option_name) != end(); }
    void THEA_ICALL setInteger(char const * option_name, int64 value) { set<int64>(option_name, value); }
    void THEA_ICALL setFloat(char const * option_name, float64 value) { set<float64>(option_name, value); }
    void THEA_ICALL setString(char const * option_name, char const * value) { set<std::string>(option_name, std::string(value)); }
    int64 THEA_ICALL getInteger(char const * option_name, int64 default_value) const { return get<int64>(option_name, default_value); }
    float64 THEA_ICALL getFloat(char const * option_name, float64 default_value) const { return get<float64>(option_name, default_value); }
    char const * THEA_ICALL getString(char const * option_name, char const * default_value) const;

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
