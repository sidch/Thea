//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (c) 2009, Stanford University
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

#ifndef __Thea_Graphics_MeshType_hpp__
#define __Thea_Graphics_MeshType_hpp__

#include "../Common.hpp"

namespace Thea {
namespace Graphics {

//
// The SFINAE tests in this header follow the examples in Vandevoorde and Josuttis, "C++ Templates: The Complete Guide".
//

/**
 * Has boolean member <code>value = true</code> if <code>T</code> is a general mesh, else <code>value = false</code>.
 *
 * @see GeneralMesh
 */
template <typename T>
class IsGeneralMesh
{
  private:
    typedef char One;
    typedef struct { char a[2]; } Two;
    template <typename U> static One test(typename U::GENERAL_MESH_TAG *);
    template <typename U> static Two test(...);

  public:
    static bool const value = (sizeof(test<T>(0)) == 1);
};

/**
 * Has boolean member <code>value = true</code> if <code>T</code> is a DCEL mesh, else <code>value = false</code>.
 *
 * @see DCELMesh
 */
template <typename T>
class IsDCELMesh
{
  private:
    typedef char One;
    typedef struct { char a[2]; } Two;
    template <typename U> static One test(typename U::DCEL_MESH_TAG *);
    template <typename U> static Two test(...);

  public:
    static bool const value = (sizeof(test<T>(0)) == 1);
};

/**
 * Has boolean member <code>value = true</code> if <code>T</code> is a display mesh, else <code>value = false</code>.
 *
 * @see DisplayMesh
 */
template <typename T>
class IsDisplayMesh
{
  private:
    typedef char One;
    typedef struct { char a[2]; } Two;
    template <typename U> static One test(typename U::DISPLAY_MESH_TAG *);
    template <typename U> static Two test(...);

  public:
    static bool const value = (sizeof(test<T>(0)) == 1);
};

/**
 * Has boolean member <code>value = true</code> if <code>T</code> is a CGAL mesh (CGALMesh or CGAL::Polyhedron_3), else
 * <code>value = false</code>.
 *
 * @see CGALMesh, CGAL::Polyhedron_3
 */
template <typename T>
class IsCGALMesh
{
  private:
    typedef char One;
    typedef struct { char a[2]; } Two;
    template <typename U> static One test(typename U::Face_handle *);  // this should be a reasonable differentiator
    template <typename U> static Two test(...);

  public:
    static bool const value = (sizeof(test<T>(0)) == 1);
};

} // namespace Graphics
} // namespace Thea

#endif
