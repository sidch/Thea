//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Thea_MatrixInvert_hpp__
#define __Thea_MatrixInvert_hpp__

// Direct element access via operator()
#define THEA_INVERT_MATRIX_FN              invertMatrix
#define THEA_INVERT_MATRIX_TEMPLATE_TYPE   MatrixT
#define THEA_INVERT_MATRIX_TYPE            MatrixT
#define THEA_MATRIX_GET(m, r, c)           (m)((r), (c))
#define THEA_MATRIX_GET_MUTABLE(m, r, c)   (m)((r), (c))
#define THEA_MATRIX_SET(m, r, c, v)        (m)((r), (c)) = (v)

#include "MatrixInvertTmpl.hpp"

#undef THEA_MATRIX_SET
#undef THEA_MATRIX_GET_MUTABLE
#undef THEA_MATRIX_GET
#undef THEA_INVERT_MATRIX_TYPE
#undef THEA_INVERT_MATRIX_TEMPLATE_TYPE
#undef THEA_INVERT_MATRIX_FN

// Virtual get/set accessors
#define THEA_INVERT_MATRIX_FN              invertAddressableMatrix
#define THEA_INVERT_MATRIX_TEMPLATE_TYPE   T
#define THEA_INVERT_MATRIX_TYPE            AddressableMatrix<T>
#define THEA_MATRIX_GET(m, r, c)           (m).get((r), (c))
#define THEA_MATRIX_GET_MUTABLE(m, r, c)   (m).getMutable((r), (c))
#define THEA_MATRIX_SET(m, r, c, v)        (m).set((r), (c), (v))

#include "MatrixInvertTmpl.hpp"

#undef THEA_MATRIX_SET
#undef THEA_MATRIX_GET_MUTABLE
#undef THEA_MATRIX_GET
#undef THEA_INVERT_MATRIX_TYPE
#undef THEA_INVERT_MATRIX_TEMPLATE_TYPE
#undef THEA_INVERT_MATRIX_FN

#endif
