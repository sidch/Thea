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

/*
-----------------------------------------------------------------------------
This source file is part of OGRE
    (Object-oriented Graphics Rendering Engine)
For the latest info, see http://www.ogre3d.org/

Copyright (c) 2000-2009 Torus Knot Software Ltd

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
-----------------------------------------------------------------------------
*/

#ifndef __Thea_DynLib_hpp__
#define __Thea_DynLib_hpp__

#include "Common.hpp"
#include "NamedObject.hpp"
#include "Map.hpp"

#ifdef THEA_WINDOWS
   struct HINSTANCE__;
   typedef struct HINSTANCE__ * hInstance;
#  define THEA_DYNLIB_HANDLE hInstance
#else
#  define THEA_DYNLIB_HANDLE void *
#endif

namespace Thea {

/** A dynamically loaded library. */
class THEA_API DynLib : public virtual NamedObject
{
  public:
    /**
     * Returns the address of the given symbol from the loaded library.
     *
     * @param sym_name The name of the symbol to search for
     *
     * @return A handle to the symbol on success, NULL on failure.
     */
    void * getSymbol(std::string const & sym_name) const;

  private:
    friend class DynLibManager;

    /** Constructor. */
    DynLib(std::string const & name);

    /** Destructor. */
    ~DynLib();

    /** Load the library. */
    void load();

    /** Unload the library. */
    void unload();

    /** Record a reference to the library. */
    void addRef();

    /** Release a reference to the library. */
    void releaseRef();

    /** Get the number of references to the library. */
    long getRefCount() const;

    /** Get any library error. */
    std::string dynlibError() const;

    THEA_DYNLIB_HANDLE h_inst;  ///< Handle to the loaded library.
    long ref_count;             ///< Number of clients using the library.

}; // class DynLib

/**
 * Manager for dynamically loaded libraries. Keeps track of all open libraries, opens libraries as needed and returns
 * references to already-open libraries.
 */
class THEA_API DynLibManager
{
  public:
    /** Destructor. */
    ~DynLibManager();

    /**
     * Loads a library, if it has not been loaded yet. If it has already been loaded, this function increments its reference
     * count.
     *
     * @param path The path to the library. The extension can be omitted.
     *
     * @return A handle to the opened library.
     *
     * @see unload()
     */
    DynLib * load(std::string const & path);

    /**
     * Unloads a library. This function decrements the reference count of the library. The library will be actually unlinked
     * when there are no more references to it.
     *
     * @param lib The library to be unloaded.
     *
     * @see load()
     */
    void unload(DynLib * lib);

    /** Unload all currently loaded libraries. */
    void unloadAll();

  private:
    typedef TheaMap<std::string, DynLib *> DynLibMap;  ///< Maps paths to the corresponding library objects.

    DynLibMap libs;  ///< Set of loaded libraries.

}; // class DynLibManager

} // namespace Thea

#endif
