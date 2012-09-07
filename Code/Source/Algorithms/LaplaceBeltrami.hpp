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

#ifndef __Thea_Algorithms_LaplaceBeltrami_hpp__
#define __Thea_Algorithms_LaplaceBeltrami_hpp__

#include "../Common.hpp"
#include "../AddressableMatrix.hpp"
#include "../Math.hpp"
#include "../ResizableMatrix.hpp"
#include "../UnorderedMap.hpp"
#include "../Graphics/MeshType.hpp"
#include <boost/type_traits/is_base_of.hpp>
#include <boost/utility/enable_if.hpp>

namespace Thea {
namespace Algorithms {

/**
 * Compute the discrete Laplace-Beltrami operator of a mesh. This only works when the mesh vertices are standard 3-vectors (with
 * dot and cross products defined). The vertices, corresponding to rows (and columns) of the matrix are numbered in sequential
 * order (according to the vertex iterator of the mesh).
 *
 * Two algorithms are provided for computing the L-B operator:
 *
 * - Guoliang Xu, "Discrete Laplace-Beltrami Operator on Sphere and Optimal Spherical Triangulations". Int. J. Comput. Geometry
 *   Appl. 16(1): 75-93 (2006).
 *
 * - M. Belkin, J. Sun, Y. Wang, "Discrete Laplace Operator for Meshed Surfaces", Proc. SoCG (2008).
 */
class THEA_API LaplaceBeltrami
{
  public:
    /** The method used to compute the discrete operator (enum class). */
    struct THEA_API Method
    {
      /** Supported values. */
      enum Value
      {
        /**
         * Guoliang Xu, "Discrete Laplace-Beltrami Operator on Sphere and Optimal Spherical Triangulations". Int. J. Comput.
         * Geometry Appl. 16(1): 75-93 (2006).
         */
        XU_2006,

        /**
         * M. Belkin, J. Sun, Y. Wang, "Discrete Laplace Operator for Meshed Surfaces", Proc. SoCG (2008).
         */
        BSW_2008
      };

      THEA_ENUM_CLASS_BODY(Method)
    };

    /** Compute the discrete Laplace-Beltrami operator for a DCELMesh and store it in the result. */
    template <typename MeshT, typename MatrixT> static void compute(MeshT const & mesh, Method method, MatrixT & result)
    {
      switch (method)
      {
        case Method::XU_2006: computeXu(mesh, result); return;
        case Method::BSW_2008: computeBSW(mesh, result); return;

        default: throw Error("LaplaceBeltrami: Unknown method");
      }
    }

  private:
    /** Shorthand to check that the mesh is of the desired and the matrix is addressable and resizable. */
    template <typename MeshTypeCheckT, typename MatrixT>
    struct CheckTypes
    {
      static bool const value = MeshTypeCheckT::value
                             && boost::is_base_of< AddressableMatrix<typename MatrixT::Value>, MatrixT >::value
                             && boost::is_base_of< ResizableMatrix<typename MatrixT::Value>, MatrixT >::value;
    };

    //==========================================================================================================================
    // Xu 2006
    //==========================================================================================================================

    /**
     * Compute the discrete Laplace-Beltrami operator for a general mesh using the method of [Xu 2006] and store it in the
     * result.
     */
    template <typename MeshT, typename MatrixT>
    static void computeXu(MeshT const & mesh, MatrixT & result,
                          typename boost::enable_if< CheckTypes<Graphics::IsGeneralMesh<MeshT>, MatrixT> >::type * dummy = 0)
    {
      // First sequentially index all vertices of the mesh
      TheaUnorderedMap<typename MeshT::Vertex const *, long> indices;
      long num_vertices = 0;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
        indices[&(*vi)] = num_vertices++;

      // Now compute the discrete Laplace-Beltrami operator
      result.resize(num_vertices, num_vertices);
      result.makeZero();

      typename MeshT::Edge const * first_edge, * ej_next, * ej_prev, * ej;
      typename MeshT::Vertex const * vx, * vj, * vj_prev, * vj_next;
      typename MatrixT::Value denom, cot_a_ij, cot_b_ij, x;
      long i, j, num_visited;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      {
        vx = &(*vi);
        i = indices[vx];

        if (vx->numEdges() > 0)  // not an isolated vertex
        {
          first_edge = *vx->edgesBegin();
          ej_prev = first_edge;
          ej      = ej_prev->nextAroundEndpoint(vx);
          ej_next = ej->nextAroundEndpoint(vx);
          denom = 0;
          num_visited = 0;
          do
          {
            vj      = ej->getOtherEndpoint(vx);
            vj_prev = ej_prev->getOtherEndpoint(vx);
            vj_next = ej_next->getOtherEndpoint(vx);

            cot_a_ij = (typename MatrixT::Value)cot(vx->getPosition(), vj_prev->getPosition(), vj->getPosition());
            cot_b_ij = (typename MatrixT::Value)cot(vx->getPosition(), vj_next->getPosition(), vj->getPosition());
            denom += (cot_a_ij + cot_b_ij) * (vx->getPosition() - vj->getPosition()).length();

            j = indices[vj];
            x = 4 * (cot_a_ij + cot_b_ij);
            result.getMutable(i, j) += x;
            result.getMutable(i, i) -= x;

            ej_prev = ej;
            ej      = ej_next;
            ej_next = ej_next->nextAroundEndpoint(vx);

            if (++num_visited == vx->numEdges() && ej_prev != first_edge)
              throw Error("LaplaceBeltrami: Mesh is not manifold");

          } while (ej_prev != first_edge);

          if (std::abs(denom) > 0)
          {
            for (typename MeshT::Vertex::EdgeConstIterator ei = vi->edgesBegin(); ei != vi->edgesEnd(); ++ei)
            {
              j = indices[(*ei)->getOtherEndpoint(&(*vi))];
              result.getMutable(i, j) /= denom;
            }

            result.getMutable(i, i) /= denom;
          }
        }
      }
    }

    /**
     * Compute the discrete Laplace-Beltrami operator for a DCEL mesh using the method of [Xu 2006] and store it in the result.
     */
    template <typename MeshT, typename MatrixT>
    static void computeXu(MeshT const & mesh, MatrixT & result,
                          typename boost::enable_if< CheckTypes<Graphics::IsDCELMesh<MeshT>, MatrixT> >::type * dummy = 0)
    {
      // First sequentially index all vertices of the mesh
      TheaUnorderedMap<typename MeshT::Vertex const *, long> indices;
      long num_vertices = 0;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
        indices[&(*vi)] = num_vertices++;

      // Now compute the discrete Laplace-Beltrami operator
      result.resize(num_vertices, num_vertices);
      result.makeZero();

      typename MeshT::Halfedge const * he_j_next, * he_j_prev, * he_j;
      typename MeshT::Vertex const * vj, * vj_prev, * vj_next;
      typename MatrixT::Value denom, cot_a_ij, cot_b_ij, x;
      long i, j;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      {
        i = indices[&(*vi)];

        he_j_prev = vi->getHalfedge();
        if (he_j_prev)  // not an isolated vertex
        {
          he_j      = he_j_prev->nextAroundOrigin();
          he_j_next = he_j->nextAroundOrigin();
          denom = 0;
          do
          {
            vj      = he_j->getEnd();
            vj_prev = he_j_prev->getEnd();
            vj_next = he_j_next->getEnd();

            cot_a_ij = (typename MatrixT::Value)cot(vi->getPosition(), vj_prev->getPosition(), vj->getPosition());
            cot_b_ij = (typename MatrixT::Value)cot(vi->getPosition(), vj_next->getPosition(), vj->getPosition());
            denom += (cot_a_ij + cot_b_ij) * (vi->getPosition() - vj->getPosition()).length();

            j = indices[vj];
            x = 4 * (cot_a_ij + cot_b_ij);
            result.getMutable(i, j) += x;
            result.getMutable(i, i) -= x;

            he_j_prev = he_j;
            he_j      = he_j_next;
            he_j_next = he_j_next->nextAroundOrigin();

          } while (he_j_prev != vi->getHalfedge());

          he_j = vi->getHalfedge();
          do
          {
            j = indices[he_j->getEnd()];
            result.getMutable(i, j) /= denom;
            he_j = he_j->nextAroundOrigin();

          } while (he_j != vi->getHalfedge());

          result.getMutable(i, i) /= denom;
        }
      }
    }

#ifdef THEA_ENABLE_CGAL

    /**
     * Compute the discrete Laplace-Beltrami operator for a CGAL mesh using the method of [Xu 2006] and store it in the result.
     */
    template <typename MeshT, typename MatrixT>
    static void computeXu(MeshT const & mesh, MatrixT & result,
                          typename boost::enable_if< CheckTypes<Graphics::IsCGALMesh<MeshT>, MatrixT> >::type * dummy = 0)
    {
      // First sequentially index all vertices of the mesh
      TheaUnorderedMap<typename MeshT::Vertex const *, long> indices;
      long num_vertices = 0;
      for (typename MeshT::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
        indices[&(*vi)] = num_vertices++;

      // Now compute the discrete Laplace-Beltrami operator
      result.resize(num_vertices, num_vertices);
      result.makeZero();

      typename MeshT::Halfedge_around_vertex_const_circulator he_j_next, he_j_prev, he_j;
      typename MeshT::Vertex_const_handle vj, vj_prev, vj_next;
      typename MatrixT::Value denom, cot_a_ij, cot_b_ij, x;
      long i, j;
      for (typename MeshT::Vertex_const_iterator vi = mesh.vertices_begin(); vi != mesh.vertices_end(); ++vi)
      {
        i = indices[&(*vi)];

        denom = 0;
        he_j_next = vi->vertex_begin();
        he_j_prev = he_j_next++;
        he_j      = he_j_next++;
        do
        {
          vj      = he_j->opposite()->vertex();
          vj_prev = he_j_prev->opposite()->vertex();
          vj_next = he_j_next->opposite()->vertex();

          cot_a_ij = (typename MatrixT::Value)cot(vi->point(), vj_prev->point(), vj->point());
          cot_b_ij = (typename MatrixT::Value)cot(vi->point(), vj_next->point(), vj->point());
          denom += (cot_a_ij + cot_b_ij) * (vi->point() - vj->point()).length();

          j = indices[&(*vj)];
          x = 4 * (cot_a_ij + cot_b_ij);
          result.getMutable(i, j) += x;
          result.getMutable(i, i) -= x;

          ++he_j;
          ++he_j_prev;
          ++he_j_next;
        } while (he_j_prev != vi->vertex_begin());

        he_j = vi->vertex_begin();
        do
        {
          j = indices[&(*he_j->opposite()->vertex())];
          result.getMutable(i, j) /= denom;

        } while (++he_j != vi->vertex_begin());

        result.getMutable(i, i) /= denom;
      }
    }

#endif

    /** Cotangent of angle ABC (vertex at B) for vectors in 3-space. The angle less than 180 degrees is chosen. */
    template <typename Vector3T> static double cot(Vector3T const & a, Vector3T const & b, Vector3T const & c)
    {
      Vector3T u = (a - b).unit();
      Vector3T v = (c - b).unit();
      double cos_b = (double)v.dot(u);
      double sin_b = (double)v.cross(u).length();
      return (sin_b < 1e-10) ? 1e+10 : (cos_b / sin_b);
    }

    //==========================================================================================================================
    // Belkin, Sun and Wang 2008
    //==========================================================================================================================

    /**
     * Compute the discrete Laplace-Beltrami operator for a mesh using the method of [Belkin, Sun and Wang 2006] and store it
     * in the result.
     */
    template <typename MeshT, typename MatrixT>
    static void computeBSW(MeshT const & mesh, MatrixT & result)
    {
    }

}; // class LaplaceBeltrami

} // namespace Algorithms
} // namespace Thea

#endif
