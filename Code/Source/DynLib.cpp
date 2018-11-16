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

#include "DynLib.hpp"

#ifdef THEA_WINDOWS
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif

#  if !defined(NOMINMAX) && defined(_MSC_VER)
#    define NOMINMAX  // required to stop windows.h messing up std::min
#  endif

#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

#ifdef THEA_WINDOWS
#  define THEA_DYNLIB_LOAD( a )       LoadLibraryEx( a, NULL, LOAD_WITH_ALTERED_SEARCH_PATH )
#  define THEA_DYNLIB_GETSYM( a, b )  GetProcAddress( a, b )
#  define THEA_DYNLIB_UNLOAD( a )     !FreeLibrary( a )
#else
#  define THEA_DYNLIB_LOAD( a )       dlopen( a, RTLD_LAZY | RTLD_GLOBAL)
#  define THEA_DYNLIB_GETSYM( a, b )  dlsym( a, b )
#  define THEA_DYNLIB_UNLOAD( a )     dlclose( a )
#endif

namespace Thea {

static std::string
DynLib_addExtension(std::string const & name)
{
#if defined(THEA_WINDOWS)
  // Although LoadLibraryEx will add .dll itself when you only specify the library name, if you include a relative path then it
  // does not. So, add it to be sure.
  if (toLower(name.substr(name.length() - 4, 4)) != ".dll")
    return name + ".dll";
#elif defined(THEA_MAC)
  // dlopen() does not add .dylib to the filename, like windows does for .dll
  if (toLower(name.substr(name.length() - 6, 6)) != ".dylib")
    return name + ".dylib";
#else
  // dlopen() does not add .so to the filename, like windows does for .dll
  if (name.substr(name.length() - 3, 3) != ".so")  // assume case-sensitive
     return name + ".so";
#endif

  return name;
}

DynLib::DynLib(std::string const & name)
: NamedObject(name), h_inst(NULL), ref_count(0)
{}

DynLib::~DynLib()
{
  unload();
}

void
DynLib::load()
{
  if (h_inst)
  {
    addRef();
    return;
  }

  std::string name = DynLib_addExtension(getName());
  THEA_LOG << "Loading library '" << getName() << '\'';

  h_inst = (THEA_DYNLIB_HANDLE)THEA_DYNLIB_LOAD(name.c_str());

  if (!h_inst)
    throw Error("Could not load dynamic library '" + name + "' (" + dynlibError() + ')');

  ref_count = 1;
}

void
DynLib::unload()
{
  if (!h_inst) return;

  THEA_LOG << "Unloading library '" << getName() << '\'';

  if (THEA_DYNLIB_UNLOAD(h_inst))
    throw Error("Could not unload dynamic library '" + getNameStr() + "' (" + dynlibError() + ')');

  h_inst = NULL;
  ref_count = 0;
}

void
DynLib::addRef()
{
  ref_count++;
}

void
DynLib::releaseRef()
{
  alwaysAssertM(ref_count > 0, "DynLib: Trying to release non-existent reference");
  ref_count--;

  if (ref_count <= 0)
    unload();
}

long
DynLib::getRefCount() const
{
  return ref_count;
}

void *
DynLib::getSymbol(std::string const & sym_name) const
{
  if (h_inst)
    return (void *)THEA_DYNLIB_GETSYM(h_inst, sym_name.c_str());
  else
    return NULL;
}

std::string
DynLib::dynlibError() const
{
#if defined(THEA_WINDOWS)
  LPVOID lpMsgBuf;
  FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                NULL,
                GetLastError(),
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                (LPTSTR)&lpMsgBuf,
                0,
                NULL);

  std::string ret = (char *)lpMsgBuf;
  LocalFree(lpMsgBuf);  // free the buffer.
  return ret;
#else
  return dlerror();
#endif
}

DynLibManager::~DynLibManager()
{
  unloadAll();
}

DynLib *
DynLibManager::load(std::string const & path)
{
  DynLibMap::iterator loaded = libs.find(path);
  if (loaded != libs.end())
  {
    loaded->second->addRef();
    return loaded->second;
  }
  else
  {
    DynLib * lib = new DynLib(path);
    lib->load();
    libs[path] = lib;
    return lib;
  }
}

void
DynLibManager::unload(DynLib * lib)
{
  if (!lib || lib->getRefCount() <= 0) return;

  lib->releaseRef();

  if (lib->getRefCount() <= 0)  // time to completely unload the lib
  {
    DynLibMap::iterator loaded = libs.find(lib->getName());
    if (loaded != libs.end())
    {
      alwaysAssertM(loaded->second == lib,
                    "DynLibManager: A different library was loaded with the same name (" + std::string(lib->getName()) + ')');
      libs.erase(loaded);
    }

    delete lib;
  }
}

void
DynLibManager::unloadAll()
{
  // Get rid of all libraries
  for (DynLibMap::iterator li = libs.begin(); li != libs.end(); ++li)
    delete li->second;

  libs.clear();
}

} // namespace Thea
