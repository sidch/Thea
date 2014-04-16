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

#ifndef __Thea_Algorithms_BestFitBox3_hpp__
#define __Thea_Algorithms_BestFitBox3_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Box3.hpp"

namespace Thea {
namespace Algorithms {

/** Approximate best-fit oriented bounding box. */
class THEA_API BestFitBox3
{
  public:
    THEA_DEF_POINTER_TYPES(BestFitBox3, shared_ptr, weak_ptr)

    /** Constructor. */
    BestFitBox3();

    /** Add a point to the set. */
    void addPoint(Vector3 const & point);

    /** Remove all data and (lazily) set the box to null. */
    void clear();

    /** Set the up vector. The computed box will only consider orientations with the given up vector. */
    void setUpVector(Vector3 const & up_) { up = up_.unit(); has_up = true; }

    /** Check if the up vector has been set. */
    bool hasUpVector() const { return has_up; }

    /**
     * Get the up vector, if it has been set.
     *
     * @see hasUpVector();
     */
    Vector3 const & getUpVector() const { return up; }

    /** Clear the up vector. Subsequent alignments will be unconstrained. */
    void clearUpVector() { has_up = false; }

    /** Remove all cached data to free memory, but do <b>not</b> mark the box for recomputation. */
    void releaseMemoryWithoutUpdate();

    /** Get the box object. Will only force a recomputation if new data have been added since the last call to this function. */
    Box3 const & getBox() const;

  private:
    /** Recompute the best-fit oriented bounding box. */
    void update() const;

    TheaArray<Vector3> points;
    bool has_up;
    Vector3 up;

    mutable Box3 box;
    mutable bool updated;

}; // class BestFitBox3

} // namespace Thea
} // namespace Algorithms

#endif
