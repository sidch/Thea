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

#ifndef __Thea_Algorithms_ConvexHull3_hpp__
#define __Thea_Algorithms_ConvexHull3_hpp__

#include "../Common.hpp"
#include "../Graphics/IncrementalMeshBuilder.hpp"
#include "../Array.hpp"
#include "../Vector3.hpp"

namespace Thea {
namespace Algorithms {

/** Convex hull of a set of points in 3-space. */
class THEA_API ConvexHull3
{
  private:

  public:
    THEA_DEF_POINTER_TYPES(ConvexHull3, shared_ptr, weak_ptr)

    /** %Options for computing convex hulls. */
    struct THEA_API Options
    {
      /** %Options for computing approximate convex hulls. */
      struct THEA_API Approx
      {
        long max_vertices_hint;  ///< Hint for maximum number of vertices on the hull.
        Real skin_width;  ///< Amount to expand the hull by, for a little leeway and robustness.

        /** Constructor. */
        Approx(long max_vertices_hint_ = -1, Real skin_width_ = 0)
        : max_vertices_hint(max_vertices_hint_), skin_width(skin_width_) {}

        /** The set of default options. */
        static Approx const & defaults() { static Approx const def; return def; }
      };

      /** %Options for computing exact convex hulls. */
      struct THEA_API Exact
      {
        /** Constructor. */
        Exact() {}

        /** The set of default options. */
        static Exact const & defaults() { static Exact const def; return def; }
      };

      Approx approx;  ///< %Options for computing approximate convex hulls.
      Exact exact;  ///< %Options for computing exact convex hulls.

      /** Constructor. */
      Options(Approx const & approx_ = Approx::defaults(), Exact const & exact_ = Exact::defaults())
      : approx(approx_), exact(exact_) {}

      /** The set of default options. */
      static Options const & defaults() { static Options const def; return def; }
    };

    /** Constructor. */
    ConvexHull3(Options const & options = Options::defaults());

    /** Add a point to the cloud. */
    void addPoint(Vector3 const & point);

    /** Remove all points and (lazily) set the convex hull to null. */
    void clear();

    /** Remove all cached points to free memory, but do <b>not</b> mark the convex hull for recomputation. */
    void releaseMemoryWithoutUpdate();

    /**
     * Compute the approximate convex hull of the data. The polyhedron representing the convex hull is <b><i>appended</i></b> to
     * the supplied mesh.
     */
    template <typename Mesh> void computeApprox(Mesh & mesh) const
    {
      updateApprox();

      typedef Graphics::IncrementalMeshBuilder<Mesh> Builder;
      Builder builder(&mesh);
      builder.begin();

      TheaArray<typename Builder::VertexHandle> approx_vrefs(approx_vertices.size());
      for (size_t i = 0; i < approx_vertices.size(); ++i)
        approx_vrefs[i] = builder.addVertex(approx_vertices[i]);

      typename Builder::VertexHandle face[3];
      for (size_t i = 0; i < approx_indices.size(); i += 3)
      {
        face[0] = approx_vrefs[approx_indices[i    ]];
        face[1] = approx_vrefs[approx_indices[i + 1]];
        face[2] = approx_vrefs[approx_indices[i + 2]];

        builder.addFace(face, face + 3);
      }

      builder.end();
    }

    /** Compute the exact convex hull of the data. */
    template <typename Mesh> void computeExact(Mesh & mesh) const
    { throw Error("ConvexHull3: Exact convex hulls not implemented"); }

  private:
    /** Recompute the approximate convex hull. */
    void updateApprox() const;

    Options options;
    TheaArray<Vector3> points;

    mutable TheaArray<Vector3> approx_vertices;
    mutable TheaArray<size_t> approx_indices;

    mutable bool approx_updated;
    mutable bool exact_updated;

}; // class ConvexHull3

} // namespace Algorithms
} // namespace Thea

#endif
