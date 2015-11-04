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

#ifndef __Thea_Algorithms_ImplicitSurfaceMesher_hpp__
#define __Thea_Algorithms_ImplicitSurfaceMesher_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Ball3.hpp"
#include "../Math.hpp"
#include "../UnorderedMap.hpp"
#include "../Graphics/IncrementalMeshBuilder.hpp"
#include "../Graphics/MeshType.hpp"
#include "BloomenthalPolygonizer/polygonizer.h"

#ifdef THEA_ENABLE_CGAL
#  include <CGAL/Surface_mesh_default_triangulation_3.h>
#  include <CGAL/Complex_2_in_triangulation_3.h>
#  include <CGAL/Implicit_surface_3.h>
#  include <CGAL/make_surface_mesh.h>
#endif

namespace Thea {
namespace Algorithms {

/**
 * Mesh an implicit surface defined as the zero level set of a 3D function. The function type <code>F</code> must be
 *
 * \code
 * <real-number-type> (*F)(Vector3 const & p)
 * \endcode
 *
 * or (if it is a functor class) must have the public member function
 *
 * \code
 * <real-number-type> operator()(Vector3 const & p) const
 * \endcode
 */
class THEA_API ImplicitSurfaceMesher
{
  private:
    /** Evaluates function as required by Bloomenthal polygonizer. */
    template <typename FunctorT> struct BloomenthalEval : public BloomenthalPolygonizer::ImplicitFunction
    {
      FunctorT * func;

      /** Constructor. */
      BloomenthalEval(FunctorT * func_) : func(func_)
      {
        alwaysAssertM(func, "ImplicitSurfaceMesher: Surface evaluator function cannot be null");
      }

      /** Evaluate the function. */
      float eval(float x, float y, float z) { return static_cast<float>((*func)(Vector3(x, y, z))); }

    }; // BloomenthalEval

#ifdef THEA_ENABLE_CGAL

    /** Evaluates function as required by CGAL implicit surface wrapper. */
    template <typename FunctorT> struct CGALEval
    {
      typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
      typedef Tr::Geom_traits GT;
      typedef GT::Point_3 Point_3;
      typedef GT::FT FT;

      FunctorT * func;

      /** Constructor. */
      CGALEval(FunctorT * func_) : func(func_) {}

      /** Evaluate the function. */
      FT operator()(Point_3 const & p) const { return static_cast<FT>((*func)(Vector3(p.x(), p.y(), p.z()))); }

    }; // CGALEval

    typedef CGAL::Surface_mesh_default_triangulation_3  Tr;
    typedef CGAL::Complex_2_in_triangulation_3<Tr>      C2t3;

#endif  // THEA_ENABLE_CGAL

  public:
    /** %Options controlling mesh generation by polygonizing the implicit surface. */
    struct THEA_API Options
    {
      /** %Options for meshing the implicit surface via the method of Bloomenthal [1994]. */
      struct THEA_API Bloomenthal
      {
        double cell_size;           ///< Size of the polygonizing cell (negative to select default).
        int max_search_steps;       /**< Limit to how far away we will look for components of the implicit surface (negative to
                                         select default). */
        bool tetrahedralize_cubes;  /**< If true, cubes are divided into tetrahedra and polygonized. Else, cubes are polygonized
                                         directly. (Default: false) */

        /** Constructor. */
        Bloomenthal(double cell_size_ = -1, int max_search_steps_ = -1, bool tetrahedralize_cubes_ = false)
        : cell_size(cell_size_), max_search_steps(max_search_steps_), tetrahedralize_cubes(tetrahedralize_cubes_) {}

        /** Default options. */
        static Bloomenthal const & defaults() { static Bloomenthal const options; return options; }
      };

#ifdef THEA_ENABLE_CGAL

      /**
       * %Options for meshing the implicit surface via the method of Boissonnat and Oudot [2005].
       *
       * @see http://www.cgal.org/Manual/last/doc_html/cgal_manual/Surface_mesher_ref/Class_Surface_mesh_default_criteria_3.html
       */
      struct THEA_API BoissonnatOudot
      {
        double min_facet_angle;        ///< Minimum facet angle, in radians (negative to select default).
        double min_delaunay_radius;    ///< Minimum radius of surface Delaunay balls (negative to select default).
        double min_center_separation;  ///< Minimum center-center distance (negative to select default).

        /** Constructor. */
        BoissonnatOudot(double min_facet_angle_ = -1, double min_delaunay_radius_ = -1, double min_center_separation_ = -1)
        : min_facet_angle(min_facet_angle_),
          min_delaunay_radius(min_delaunay_radius_),
          min_center_separation(min_center_separation_)
        {}

        /** Default options. */
        static BoissonnatOudot const & defaults() { static BoissonnatOudot const options; return options; }
      };

#endif // THEA_ENABLE_CGAL

      /** Constructor. */
      Options(Bloomenthal const & bloomenthal_ = Bloomenthal::defaults()

#ifdef THEA_ENABLE_CGAL
              , BoissonnatOudot const & boissonnat_oudot_ = BoissonnatOudot::defaults()
#endif // THEA_ENABLE_CGAL

      )
      : bloomenthal(bloomenthal_)

#ifdef THEA_ENABLE_CGAL
      , boissonnat_oudot(boissonnat_oudot_)
#endif // THEA_ENABLE_CGAL

      {}

      /** Default options. */
      static Options const & defaults() { static Options const options; return options; }

      Bloomenthal bloomenthal;           ///< %Options for meshing via Bloomenthal [1994].

#ifdef THEA_ENABLE_CGAL
      BoissonnatOudot boissonnat_oudot;  ///< %Options for meshing via Boissonnat and Oudot [2005].
#endif // THEA_ENABLE_CGAL

    }; // struct Options

    /**
     * Polygonize the zero level set of a 3D function to a mesh, using the method of Bloomenthal [1994].
     *
     *   Jules Bloomenthal, "An implicit surface polygonizer", Graphics Gems IV (P. Heckbert, ed.), Academic Press, New York,
     *   1994.
     *
     * @param surface_functor The function whose zero level set defines the surface.
     * @param bounding_ball Bounds the surface.
     * @param pt_near_surface A point on or near the zero level set.
     * @param options Options controlling mesh generation.
     * @param result The output mesh will be stored here (any prior data will <b>not</b> be removed from the mesh).
     */
    template <typename FunctorT, typename MeshT>
    static void meshBloomenthal(FunctorT * surface_functor, Ball3 const & bounding_ball, Vector3 const & pt_near_surface,
                                Options::Bloomenthal const & options, MeshT & result)
    {
      BloomenthalEval<FunctorT> func(surface_functor);
      float size = (float)(options.cell_size < 0 ? bounding_ball.getRadius() / 10.0 : options.cell_size);
      int bounds = options.max_search_steps < 0 ? 10 : options.max_search_steps;
      BloomenthalPolygonizer::Polygonizer polygonizer(&func, size, bounds);
      polygonizer.march(options.tetrahedralize_cubes,
                        (float)pt_near_surface.x(), (float)pt_near_surface.y(), (float)pt_near_surface.z());

      THEA_LOG << "ImplicitSurfaceMesher: " << polygonizer.no_triangles() << " triangles generated via Bloomenthal";

      exportMesh(polygonizer, result);
    }

#ifdef THEA_ENABLE_CGAL

    /**
     * Polygonize the zero level set of a 3D function to a mesh, using the method of Boissonnat and Oudot [2005] as implemented
     * by the CGAL library (http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Surface_mesher/Chapter_main.html)
     *
     *   Jean-Daniel Boissonnat and Steve Oudot, "Provably good sampling and meshing of surfaces", Graphical Models, 67:405-451,
     *   2005.
     *
     * @param surface_functor The function whose zero level set defines the surface.
     * @param bounding_ball Bounds the surface. The functor <em>must</em> evaluate to a negative value at its center.
     * @param options Options controlling mesh generation.
     * @param result The output mesh will be stored here (any prior data will <b>not</b> be removed from the mesh).
     */
    template <typename FunctorT, typename MeshT>
    static void meshBoissonnatOudot(FunctorT * surface_functor, Ball3 const & bounding_ball,
                                    Options::BoissonnatOudot const & options, MeshT & result)
    {
      typedef Tr::Geom_traits                                     GT;
      typedef GT::Sphere_3                                        Sphere_3;
      typedef GT::Point_3                                         Point_3;
      typedef GT::FT                                              FT;
      typedef CGAL::Implicit_surface_3< GT, CGALEval<FunctorT> >  Surface_3;

      Tr tr;                                     // 3D Delaunay triangulation
      C2t3 c2t3(tr);                             // 2D complex in 3D Delaunay triangulation
      CGALEval<FunctorT> func(surface_functor);  // CGAL-style wrapper to evaluate the function

      // Define the surface
      Vector3 const & center = bounding_ball.getCenter();
      Real radius = bounding_ball.getRadius();
      Sphere_3 sph(Point_3(static_cast<FT>(center.x()), static_cast<FT>(center.y()), static_cast<FT>(center.z())),
                   static_cast<FT>(radius * radius));
      Surface_3 surface(func, sph);

      // Define meshing criteria
      CGAL::Surface_mesh_default_criteria_3<Tr> criteria(
          static_cast<FT>(options.min_facet_angle       < 0 ? 30  : Math::radiansToDegrees(options.min_facet_angle)),
          static_cast<FT>(options.min_delaunay_radius   < 0 ? 0.1 : options.min_delaunay_radius),
          static_cast<FT>(options.min_center_separation < 0 ? 0.1 : options.min_center_separation));

      // Generate mesh
      CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

      THEA_LOG << "ImplicitSurfaceMesher: " << c2t3.number_of_facets() << " faces generated via Boissonnat-Oudot";

      exportMesh(c2t3, result);
    }

#endif // THEA_ENABLE_CGAL

  private:
    /** Export from Bloomenthal polygonizer to another mesh type. */
    template <typename MeshT> static void exportMesh(BloomenthalPolygonizer::Polygonizer & polygonizer, MeshT & dst)
    {
      typedef Graphics::IncrementalMeshBuilder<MeshT> Builder;
      Builder builder(&dst);

      TheaArray<typename Builder::VertexHandle> vertices;

      builder.begin();
        for (size_t i = 0; i < polygonizer.no_vertices(); ++i)
        {
          BloomenthalPolygonizer::VERTEX const & bp_p = polygonizer.get_vertex(i);
          BloomenthalPolygonizer::NORMAL const & bp_n = polygonizer.get_normal(i);
          Vector3 p(bp_p.x, bp_p.y, bp_p.z), n(bp_n.x, bp_n.y, bp_n.z);

          typename Builder::VertexHandle vx = builder.addVertex(p, &n);
          if (!vx)
            throw Error("ImplicitSurfaceMesher: Could not add vertex from Bloomenthal polygonizer to mesh");

          vertices.push_back(vx);
        }

        for (size_t i = 0; i < polygonizer.no_triangles(); ++i)
        {
          BloomenthalPolygonizer::TRIANGLE const & tri = polygonizer.get_triangle(i);
          alwaysAssertM(tri.v0 >= 0 && (size_t)tri.v0 < vertices.size()
                     && tri.v1 >= 0 && (size_t)tri.v1 < vertices.size()
                     && tri.v2 >= 0 && (size_t)tri.v2 < vertices.size(),
                        "ImplicitSurfaceMesher: Vertex index in triangle from Bloomenthal polygonizer is out of bounds");

          typename Builder::VertexHandle vx[3] = { vertices[(size_t)tri.v0],
                                                   vertices[(size_t)tri.v1],
                                                   vertices[(size_t)tri.v2] };
          if (!builder.addFace(vx, vx + 3))
            throw Error("ImplicitSurfaceMesher: Could not add triangle from Bloomenthal polygonizer to mesh");
        }
      builder.end();
    }

#ifdef THEA_ENABLE_CGAL

    /** Export C2t3 to another mesh type. */
    template <typename MeshT> static void exportMesh(C2t3 & src, MeshT & dst)
    {
      typedef Tr::Geom_traits::Point_3 Point_3;

      typedef Graphics::IncrementalMeshBuilder<MeshT> Builder;
      Builder builder(&dst);

      typedef TheaUnorderedMap<Tr::Vertex const *, typename Builder::VertexHandle> VertexMap;
      VertexMap vmap;

      builder.begin();
        for (C2t3::Vertex_iterator vi = src.vertices_begin(); vi != src.vertices_end(); ++vi)
        {
          Tr::Vertex const * vin = &(*vi);
          Point_3 const & p = vin->point();

          typename Builder::VertexHandle vout = builder.addVertex(Vector3((Real)p.x(), (Real)p.y(), (Real)p.z()));
          if (!vout)
            throw Error("ImplicitSurfaceMesher: Could not add vertex from Boissonnat-Oudot polygonizer to mesh");

          vmap[vin] = vout;
        }

        typename Builder::VertexHandle face_vertices[3];
        int indices[3];
        for (C2t3::Facet_iterator fi = src.facets_begin(); fi != src.facets_end(); ++fi)
        {
          switch (fi->second)
          {
            case 0: indices[0] = 1; indices[1] = 2; indices[2] = 3; break;
            case 1: indices[0] = 0; indices[1] = 3; indices[2] = 2; break;
            case 2: indices[0] = 3; indices[1] = 0; indices[2] = 1; break;
            case 3: indices[0] = 2; indices[1] = 1; indices[2] = 0; break;
            default:
              throw Error("ImplicitSurfaceMesher: Mesh created by Boissonnat-Oudot polygonizer has invalid face index in cell");
          }

          for (int i = 0; i < 3; ++i)
          {
            Tr::Vertex const * vx = &(*fi->first->vertex(indices[i]));
            typename VertexMap::const_iterator existing = vmap.find(vx);
            if (existing == vmap.end())
              throw Error("ImplicitSurfaceMesher: Mesh created by Boissonnat-Oudot polygonizer refers to unmapped vertex");

            face_vertices[i] = existing->second;
          }

          if (!builder.addFace(face_vertices, face_vertices + 3))
            throw Error("ImplicitSurfaceMesher: Could not add triangle from Boissonnat-Oudot polygonizer to mesh");
        }
      builder.end();
    }

#endif // THEA_ENABLE_CGAL

}; // class ImplicitSurfaceMesher

} // namespace Algorithms
} // namespace Thea

#endif
