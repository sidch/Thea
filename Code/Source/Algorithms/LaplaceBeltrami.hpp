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

#ifndef __Thea_Algorithms_LaplaceBeltrami_hpp__
#define __Thea_Algorithms_LaplaceBeltrami_hpp__

#include "../Common.hpp"
#include "../IAddressableMatrix.hpp"
#include "../Math.hpp"
#include "../UnorderedMap.hpp"
#include "../Graphics/MeshType.hpp"
#include <type_traits>

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

    /** Compute the discrete Laplace-Beltrami operator for a DcelMesh and store it in the result. */
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
                             && std::is_base_of< IAddressableMatrix<typename MatrixT::Value>, MatrixT >::value;
    };

    //==========================================================================================================================
    // Xu 2006
    //==========================================================================================================================

    /**
     * Compute the discrete Laplace-Beltrami operator for a general mesh using the method of [Xu 2006] and store it in the
     * result.
     */
    template < typename MeshT, typename MatrixT,
               typename std::enable_if< CheckTypes<Graphics::IsGeneralMesh<MeshT>, MatrixT>::value, int >::type = 0 >
    static void computeXu(MeshT const & mesh, MatrixT & result)
    {
      // First sequentially index all vertices of the mesh
      UnorderedMap<typename MeshT::Vertex const *, intx> indices;
      intx num_vertices = 0;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
        indices[&(*vi)] = num_vertices++;

      // Now compute the discrete Laplace-Beltrami operator
      if (!result.resize(num_vertices, num_vertices))
        throw Error("LaplaceBeltrami: Could not resize matrix");

      result.setZero();

      typename MeshT::Edge const * first_edge, * ej_next, * ej_prev, * ej;
      typename MeshT::Vertex const * vx, * vj, * vj_prev, * vj_next;
      typename MatrixT::Value denom, cot_a_ij, cot_b_ij, x;
      intx i, j, num_visited;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
      {
        vx = &(*vi);
        i = indices[vx];

        if (vx->numEdges() > 0)  // not an isolated vertex
        {
          // FIXME: Safer and faster to use the fact that faces must be triangles (as cotangent weights are not defined
          // otherwise), and find the third vertex as in SurfaceParametrization. This also avoids the need for the surface to
          // be properly oriented.

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
            denom += (cot_a_ij + cot_b_ij) * (vx->getPosition() - vj->getPosition()).norm();

            j = indices[vj];
            x = 4 * (cot_a_ij + cot_b_ij);
            result.mutableAt(i, j) += x;
            result.mutableAt(i, i) -= x;

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
              result.mutableAt(i, j) /= denom;
            }

            result.mutableAt(i, i) /= denom;
          }
        }
      }
    }

    /**
     * Compute the discrete Laplace-Beltrami operator for a DCEL mesh using the method of [Xu 2006] and store it in the result.
     */
    template < typename MeshT, typename MatrixT,
               typename std::enable_if< CheckTypes<Graphics::IsDcelMesh<MeshT>, MatrixT>::value, int >::type = 0 >
    static void computeXu(MeshT const & mesh, MatrixT & result)
    {
      // First sequentially index all vertices of the mesh
      UnorderedMap<typename MeshT::Vertex const *, intx> indices;
      intx num_vertices = 0;
      for (typename MeshT::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
        indices[&(*vi)] = num_vertices++;

      // Now compute the discrete Laplace-Beltrami operator
      if (!result.resize(num_vertices, num_vertices))
        throw Error("LaplaceBeltrami: Could not resize matrix");

      result.setZero();

      typename MeshT::Halfedge const * he_j_next, * he_j_prev, * he_j;
      typename MeshT::Vertex const * vj, * vj_prev, * vj_next;
      typename MatrixT::Value denom, cot_a_ij, cot_b_ij, x;
      intx i, j;
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
            denom += (cot_a_ij + cot_b_ij) * (vi->getPosition() - vj->getPosition()).norm();

            j = indices[vj];
            x = 4 * (cot_a_ij + cot_b_ij);
            result.mutableAt(i, j) += x;
            result.mutableAt(i, i) -= x;

            he_j_prev = he_j;
            he_j      = he_j_next;
            he_j_next = he_j_next->nextAroundOrigin();

          } while (he_j_prev != vi->getHalfedge());

          he_j = vi->getHalfedge();
          do
          {
            j = indices[he_j->getEnd()];
            result.mutableAt(i, j) /= denom;
            he_j = he_j->nextAroundOrigin();

          } while (he_j != vi->getHalfedge());

          result.mutableAt(i, i) /= denom;
        }
      }
    }

    /** Cotangent of angle ABC (vertex at B) for vectors in 3-space. The angle less than 180 degrees is chosen. */
    template <typename Vector3T> static double cot(Vector3T const & a, Vector3T const & b, Vector3T const & c)
    {
      Vector3T u = (a - b).normalized();
      Vector3T v = (c - b).normalized();
      double cos_b = (double)v.dot(u);
      double sin_b = (double)v.cross(u).norm();
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
      // TODO
      throw FatalError("BSW method not implemented");
    }

}; // class LaplaceBeltrami

} // namespace Algorithms
} // namespace Thea

#endif
