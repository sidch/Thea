//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_FilePath_hpp__
#define __Thea_FilePath_hpp__

#include "Common.hpp"

namespace Thea {

/**
 * Operations on file paths. These do string manipulation and do <b>not</b> actually access the filesystem.
 *
 * @note Returned paths are in native OS format (e.g. backslashes on Windows, forward slashes on Unix).
 * @note Trailing slashes at the end of directory names, and repeated slashes, are ignored, so "/foo" and "/foo/" are treated
 *   identically, as are "foo/bar" and "foo///bar". In this respect this implementation differs from, say, Boost.Filesystem.
 */
class THEA_API FilePath
{
  public:
    /**
     * Get the basename of the object, consisting of the object name without the path, upto but not including the <i>first</i>
     * period character.
     *
     * @see completeBaseName()
     */
    static std::string baseName(std::string const & path);

    /**
     * Get the complete basename of the object, consisting of the object name without the path, upto but not including the
     * <i>last</i> period character.
     *
     * @see baseName()
     */
    static std::string completeBaseName(std::string const & path);

    /**
     * Get the extension of the object, consisting of all characters after the last period character (or null if no period
     * exists).
     */
    static std::string suffix(std::string const & path);

    /**
     * Get the complete extension of the object, consisting of all characters after the first period character (or null if no
     * period exists).
     */
    static std::string completeSuffix(std::string const & path);

    /** Get the name of the object without the path. Ignores trailing slashes and isolated period characters. */
    static std::string objectName(std::string const & path);

    /** Get the path to the immediate parent of an object. Ignores trailing slashes and isolated period characters. */
    static std::string parent(std::string const & path);

    /**
     * Get the full path to a child object, given the path to the containing directory and the name of the child itself. If the
     * directory name is empty, just the child name is returned.
     */
    static std::string concat(std::string const & parent, std::string const & child);

}; // class FilePath

} // namespace Thea

#endif
