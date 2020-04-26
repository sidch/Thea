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

#ifndef __Thea_Algorithms_ConvexHull3_hpp__
#define __Thea_Algorithms_ConvexHull3_hpp__

#include "../Common.hpp"
#include "../Graphics/IncrementalMeshBuilder.hpp"
#include "../Array.hpp"
#include "../MatVec.hpp"

namespace Thea {
namespace Algorithms {

/** Convex hull of a set of points in 3-space. */
class THEA_API ConvexHull3
{
  private:

  public:
    THEA_DECL_SMART_POINTERS(ConvexHull3)

    /** %Options for computing convex hulls. */
    struct THEA_API Options
    {
      /** %Options for computing approximate convex hulls. */
      struct THEA_API Approx
      {
        intx max_vertices_hint;  ///< Hint for maximum number of vertices on the hull.
        Real skin_width;  ///< Amount to expand the hull by, for a little leeway and robustness.

        /** Constructor. */
        Approx(intx max_vertices_hint_ = -1, Real skin_width_ = 0)
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

      Array<typename Builder::VertexHandle> approx_vrefs(approx_vertices.size());
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
    Array<Vector3> points;

    mutable Array<Vector3> approx_vertices;
    mutable Array<size_t> approx_indices;

    mutable bool approx_updated;
    mutable bool exact_updated;

}; // class ConvexHull3

} // namespace Algorithms
} // namespace Thea

#endif
