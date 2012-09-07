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

/** A set of options, specified as key-value pairs. */
class THEA_API Options : private TheaMap<std::string, boost::any>
{
  public:
    /** Check if an option has been set. */
    bool hasOption(std::string const & option_name) const { return find(option_name) != end(); }

    /** Set the value of an option. */
    template <typename T> void set(std::string const & option_name, T const & value) { (*this)[option_name] = value; }

    /** Set the value of an option from a C-style string. It should be read back as a <code>std::string</code>. */
    void set(std::string const & option_name, char const * value) { (*this)[option_name] = std::string(value); }

    /**
     * Get the value of an option. If the option has not been set, the default value specified by the last parameter is
     * returned.
     */
    template <typename T> T get(std::string const & option_name, T const & default_value) const
    {
      const_iterator existing = find(option_name);
      if (existing != end())
        return boost::any_cast<T>(existing->second);
      else
        return default_value;
    }

}; // class Options

} // namespace Thea

#endif
