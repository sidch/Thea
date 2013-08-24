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

#ifndef __Thea_FileSystem_hpp__
#define __Thea_FileSystem_hpp__

#include "Common.hpp"

namespace Thea {

/**
 * Filesystem operations. Unlike FilePath, these functions do actually access the filesystem.
 *
 * @note Returned paths are in native OS format (e.g. backslashes on Windows, forward slashes on Unix).
 */
class THEA_API FileSystem
{
  public:
    /** Check if a file or directory exists. */
    static bool exists(std::string const & path);

    /** Check if a file exists, and is indeed a regular file (and not for instance a directory). */
    static bool fileExists(std::string const & path);

    /** Check if a directory exists, and is indeed a directory (and not for instance a file). */
    static bool directoryExists(std::string const & path);

    /** Get the length of a file in bytes. Returns a negative number on failure. */
    static int64 fileSize(std::string const & path);

    /** Resolve a relative path. */
    static std::string resolve(std::string const & path);

    /**
     * Create a directory, including all necessary parents (equivalent to "mkdir -p").
     *
     * @return True on success, false on error.
     */
    static bool createDirectory(std::string const & path);

    /** Get the entire contents of a file as a string. */
    static std::string readWholeFile(std::string const & path)
    {
      std::string s;
      if (!readWholeFile(path, s))
        throw Error("FileSystem: Could not read '" + path + '\'');

      return s;
    }

    /**
     * Get the entire contents of a file as a string.
     *
     * @return True on success, false on error.
     */
    static bool readWholeFile(std::string const & path, std::string & ret);

}; // class FileSystem

} // namespace Thea

#endif
