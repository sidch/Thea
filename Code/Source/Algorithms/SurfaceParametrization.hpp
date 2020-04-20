//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2019, Siddhartha Chaudhuri
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

#ifndef __Thea_Algorithms_SurfaceParametrization_hpp__
#define __Thea_Algorithms_SurfaceParametrization_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "../SparseMatVec.hpp"
#include "../UnorderedMap.hpp"
#include "../UnorderedSet.hpp"
#include "IteratorModifiers.hpp"
#include "StdLinearSolver.hpp"
#include <iterator>
#include <type_traits>

namespace Thea {
namespace Algorithms {

/** Project 2-manifolds to a 2D (UV) parameter space. */
class THEA_API SurfaceParametrization
{
  private:
    /** A triplet class satisfying Eigen's requirements, plus allowing triplets to be modified after construction. */
    struct Triplet
    {
      /** Constructor. */
      Triplet(intx r_, intx c_, double val_) : r(r_), c(c_), val(val_) {}

      /** The row of the matrix entry. */
      intx row() const { return r; }

      /** The column of the matrix entry. */
      intx col() const { return c; }

      /** The value of the matrix entry. */
      double value() const { return val; }

      intx r;      ///< The row of the matrix entry.
      intx c;      ///< The column of the matrix entry.
      double val;  ///< The value of the matrix entry.
    };

  public:
    /** Different ways to compute edge weights (enum class). */
    struct WeightType
    {
      /** Supported values. */
      enum Value
      {
        UNIFORM,     ///< Uniform weights.
        COTANGENT,   ///< Cotangent weights.
        MEAN_VALUE,  ///< Mean value coordinates (Floater 2003).
      };

      THEA_ENUM_CLASS_BODY(WeightType)

      THEA_ENUM_CLASS_STRINGS_BEGIN(WeightType)
        THEA_ENUM_CLASS_STRING(UNIFORM,     "uniform")
        THEA_ENUM_CLASS_STRING(COTANGENT,   "cotangent")
        THEA_ENUM_CLASS_STRING(MEAN_VALUE,  "mean_value")
      THEA_ENUM_CLASS_STRINGS_END(WeightType)
    };

    /**
     * Parametrize a collection of mesh faces forming a surface patch.
     *
     * @param faces_begin Iterator to the first face in the sequence.
     * @param faces_end Iterator to (one position beyond) the last face in the sequence.
     * @param weight_type The type of edge weights to use.
     * @param fixed_vertex_params The parameters of the fixed vertices. Can be the same container as \a free_vertex_params.
     * @param free_vertex_params Used to return the parameters of the free vertices, as computed by this function. Can be the
     *   same container as \a fixed_vertex_params.
     *
     * @return True on success, false on error.
     */
    template <typename FaceIteratorT, typename FixedVertexParameterMapT, typename FreeVertexParameterMapT>
    static bool parametrize(FaceIteratorT faces_begin, FaceIteratorT faces_end, WeightType weight_type,
                            FixedVertexParameterMapT const & fixed_vertex_params, FreeVertexParameterMapT & free_vertex_params)
    {
      // Quick check to catch the case when the user has forgotten to provide any boundary parameters. Partial boundary
      // initialization is detected by more expensive means below.
      if (fixed_vertex_params.empty() && faces_begin != faces_end)
      {
        THEA_ERROR << "SurfaceParametrization: Free-boundary parametrization not currently implemented";
        return false;
      }

      typedef typename std::remove_cv<
                  typename std::remove_pointer<
                      typename std::remove_cv< typename FreeVertexParameterMapT::key_type >::type >::type >::type VertexT;
      typedef typename VertexT::Edge                         EdgeT;
      typedef typename VertexT::Face                         FaceT;
      typedef typename FreeVertexParameterMapT::mapped_type  FreeParamT;  // some 2D vector

      // Cache the set of faces so we can quickly check if a face is in the patch or not (and also avoid
      // dereferencing/incrementing the input iterators multiple times, in case this is expensive)
      typedef UnorderedSet<FaceT const *> FaceSet;
      FaceSet faces(makePtrIterator(faces_begin), makePtrIterator(faces_end));

      typedef UnorderedMap<VertexT const *, intx>  VertexIndexMap;
      typedef UnorderedSet<EdgeT const *>          EdgeSet;
      VertexIndexMap interior_verts;
      EdgeSet interior_edges;

      // First build the set of all interior edges
      for (auto fi = faces.begin(); fi != faces.end(); ++fi)
      {
        for (auto fei = (*fi)->edgesBegin(); fei != (*fi)->edgesEnd(); ++fei)
        {
          auto e = *fei;
          if (e->numFaces() > 2)
          {
            THEA_ERROR << "SurfaceParametrization: Mesh is non-manifold (edge is shared by > 2 faces)";
            return false;
          }

          // If this is a boundary edge, make sure all its vertices have fixed parameters
          bool is_boundary = (e->numFaces() < 2);
          if (!is_boundary)
          {
            for (auto efi = e->facesBegin(); efi != e->facesEnd(); ++efi)
              if (faces.find(*efi) == faces.end())  // this is a boundary edge
              {
                is_boundary = true;
                break;
              }
          }

          if (is_boundary)
          {
            if (fixed_vertex_params.find(e->getEndpoint(0)) == fixed_vertex_params.end()
             || fixed_vertex_params.find(e->getEndpoint(1)) == fixed_vertex_params.end())
            {
              THEA_ERROR << "SurfaceParametrization: Boundary vertex lacks fixed parameters";
              return false;
            }
          }
          else
            interior_edges.insert(e);
        }
      }

      // Now build and index the set of all interior vertices
      for (auto e : interior_edges)
        for (intx i = 0; i < 2; ++i)
        {
          auto v = e->getEndpoint(i);
          if (fixed_vertex_params.find(v) == fixed_vertex_params.end()  // already checked every boundary vertex is parametrized
           && interior_verts.find(v) == interior_verts.end())
          {
            intx index = (intx)interior_verts.size();
            interior_verts[v] = index;
          }
        }

      // For every edge (i, j) leaving every interior vertex i, there is a weight w_ij
      intx num_vars = (intx)interior_verts.size();
      Array<Triplet> triplets;
      VectorXd u_consts(num_vars); u_consts.setZero();
      VectorXd v_consts(num_vars); v_consts.setZero();
      bool weight_error = false;
      for (auto v : interior_verts)
      {
        // The diagonal is all 1's
        triplets.push_back(Triplet(v.second, v.second, 1));

        // Add a weight for each neighbor
        double sum_weights = 0;
        size_t tbase = triplets.size();
        for (auto vei = v.first->edgesBegin(); vei != v.first->edgesEnd(); ++vei)
        {
          auto e = *vei;
          auto v1 = e->getOtherEndpoint(v.first);
          double w_ij = getWeight(weight_type, v.first, v1, e, weight_error);
          sum_weights += w_ij;

          auto other = interior_verts.find(v1);
          if (other != interior_verts.end())
            triplets.push_back(Triplet(v.second, other->second, -w_ij));
          else
          {
            auto p = fixed_vertex_params.find(v1);
            debugAssertM(p != fixed_vertex_params.end(), "SurfaceParametrization: Boundary vertex lacks fixed parameters");
            u_consts[v.second] += w_ij * p->second[0];
            v_consts[v.second] += w_ij * p->second[1];
          }
        }

        // Normalization
        for (size_t i = tbase; i < triplets.size(); ++i)
          triplets[i].val /= sum_weights;

        u_consts[v.second] /= sum_weights;
        v_consts[v.second] /= sum_weights;
      }

      // Build the sparse coefficient matrix
      SparseMatrix<double> coeffs(num_vars, num_vars);
      coeffs.setFromTriplets(triplets.begin(), triplets.end());

      // Solve the system for U and V
      StdLinearSolver u_solver(StdLinearSolver::Method::SPARSE_LU);
      StdLinearSolver v_solver(StdLinearSolver::Method::SPARSE_LU);
      if (!u_solver.solve(coeffs, u_consts.data()) || !v_solver.solve(coeffs, v_consts.data()))
      {
        THEA_ERROR << "SurfaceParametrization: Couldn't solve sparse linear system";
        return false;
      }

      float64 const * u_sol = u_solver.getSolution();
      float64 const * v_sol = v_solver.getSolution();
      for (auto v : interior_verts)
      {
        // This is the only place the free vertices are updated, and fixed_vertex_params is not accessed here or later. Hence,
        // the two references can point to the same underlying container.
        free_vertex_params[v.first][0] = static_cast<typename FreeParamT::value_type>(u_sol[v.second]);
        free_vertex_params[v.first][1] = static_cast<typename FreeParamT::value_type>(v_sol[v.second]);
      }

      return true;
    }

  private:
    /**
     * Get the weight of the edge \a e connecting source vertex \a vi to neighboring vertex \a vj. \a error is used to signal if
     * an error occurred.
     */
    template <typename VertexT, typename EdgeT>
    static double getWeight(WeightType weight_type, VertexT const * vi, VertexT const * vj, EdgeT const * e, bool & error)
    {
      error = false;
      switch (weight_type)
      {
        case WeightType::UNIFORM:     return 1;
        case WeightType::COTANGENT:   return getCotangentWeight(vi, vj, e, error);
        case WeightType::MEAN_VALUE:  return getMeanValueWeight(vi, vj, e, error);
        default:
          THEA_ERROR << "SurfaceParametrization: Unsupported weight type";
          error = true; return -1;
      }
    }

    /**
     * Get the cotangent weight of the edge \a e connecting source vertex \a vi to neighboring vertex \a vj. \a error is used to
     * signal if an error occurred.
     */
    template <typename VertexT, typename EdgeT>
    static double getCotangentWeight(VertexT const * vi, VertexT const * vj, EdgeT const * e, bool & error)
    {
      typename EdgeT::Face const * f[2];
      if (getFaces(e, f) != 2) { error = true; return -1; }

      double w = 0;
      for (intx n = 0; n < 2; ++n)
      {
        auto vk = getThirdVertex(f[n], vi, vj);  // assuming triangles makes search faster
        if (!vk) { error = true; return -1; }

        w += cot(vk, vi, vj);
      }

      return w;
    }

    /**
     * Get the mean value weight of the edge \a e connecting source vertex \a vi to neighboring vertex \a vj. \a error is used
     * to signal if an error occurred.
     */
    template <typename VertexT, typename EdgeT>
    static double getMeanValueWeight(VertexT const * vi, VertexT const * vj, EdgeT const * e, bool & error)
    {
      typename EdgeT::Face const * f[2];
      if (getFaces(e, f) != 2) { error = true; return -1; }

      double len_ij = (vj->getPosition() - vi->getPosition()).norm();
      if (Math::fuzzyEq(len_ij, 0.0))
        return 0;

      double w = 0;
      for (intx n = 0; n < 2; ++n)
      {
        auto vk = getThirdVertex(f[n], vi, vj);  // assuming triangles makes search faster
        if (!vk) { error = true; return -1; }

        w += tanHalf(vi, vj, vk, len_ij);
      }

      return w / len_ij;
    }

    /**
     * Get the faces adjacent to an edge \a e. Currently throws an error if called for anything other than an interior, manifold
     * edge with exactly two incident faces.
     */
    template <typename EdgeT, typename FaceT>
    static intx getFaces(EdgeT const * e, FaceT const ** f)
    {
      // Currently, this will only be called for interior edges, which must have exactly two incident faces
      if (e->numFaces() != 2)
      {
        THEA_ERROR << "SurfaceParametrization: Interior edge doesn't have exactly two incident faces";
        return -1;
      }

      auto efi = e->facesBegin();
      f[0] = *efi;
      f[1] = *(++efi);

      return 2;
    };

    /** Get the vertex of a triangular face \a f that is neither of two selected vertices \a vi and \a vj. */
    template <typename FaceT, typename VertexT>
    static VertexT const * getThirdVertex(FaceT const * f, VertexT const * vi, VertexT const * vj)
    {
      if (f->numVertices() != 3)
      {
        THEA_ERROR << "SurfaceParametrization: Selected weight type requires triangular faces";
        return nullptr;
      }

      for (auto fvi = f->verticesBegin(); fvi != f->verticesEnd(); ++fvi)
        if (*fvi != vi && *fvi != vj)
          return *fvi;

      return nullptr;
    }

    /** Get the cotangent of the angle between (\a v0, \a v1) and (\a v0, \a v2). */
    template <typename VertexT>
    static double cot(VertexT const * v0, VertexT const * v1, VertexT const * v2)
    {
      Vector3 const & e01 = v1->getPosition() - v0->getPosition();
      Vector3 const & e02 = v2->getPosition() - v0->getPosition();

      double c = e01.dot(e02);           // NOTE: not normalized
      double s = e01.cross(e02).norm();  // NOTE: not normalized

      return c / s;
    }

    /**
     * Get the tangent of half the angle between (\a v0, \a v1) and (\a v0, \a v2). The length of the first vector is assumed to
     * be precomputed and passed to the function to save a normalization, since mean value coordinates compute this length
     * anyway.
     */
    template <typename VertexT>
    static double tanHalf(VertexT const * v0, VertexT const * v1, VertexT const * v2, double e01_len)
    {
      Vector3 const & e01 =  v1->getPosition() - v0->getPosition();
      Vector3 const & e02 = (v2->getPosition() - v0->getPosition()).normalized();

      // Make use of the formula tan(x/2) = sqrt((1 - cos(x)) / (1 + cos(x))), for x < pi
      double c = e01.dot(e02);  // is larger than cos(x) by a factor of e01_len
      return std::sqrt((e01_len - c) / (e01_len + c));
    }

}; // class SurfaceParametrization

} // namespace Algorithms
} // namespace Thea

#endif
