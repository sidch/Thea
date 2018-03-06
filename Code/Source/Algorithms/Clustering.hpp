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

#ifndef __Thea_Algorithms_Clustering_hpp__
#define __Thea_Algorithms_Clustering_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../CompressedSparseMatrix.hpp"
#include "../Math.hpp"
#include "../Matrix.hpp"
#include "Metrics.hpp"

#ifdef THEA_ENABLE_CLUTO
#  include "cluto.h"
#endif

namespace Thea {

/** Miscellaneous algorithms. */
namespace Algorithms {

/** %Clustering methods. */
class THEA_API Clustering
{
  private:
#ifdef THEA_ENABLE_CLUTO
    typedef CompressedRowMatrix<float, int, int> ClutoMatrix;
#endif

  public:
    THEA_DEF_POINTER_TYPES(Clustering, shared_ptr, weak_ptr)

    /** Flat clustering methods (enum class). */
    struct THEA_API FlatMethod
    {
      /** Supported values. */
      enum Value
      {
        AUTO,                   ///< Automatically choose an appropriate method.

#ifdef THEA_ENABLE_CLUTO
        CLUTO_CLUSTER_DIRECT,   ///< CLUTO direct clustering.
        CLUTO_CLUSTER_RB,       ///< CLUTO clustering via repeated bisection.
        CLUTO_GRAPH_CLUSTER_RB  ///< CLUTO graph-partitioning-based clustering via repeated min-cut bisection.
#endif
      };

      THEA_ENUM_CLASS_BODY(FlatMethod)
    };

    /** Hierarchical clustering methods (enum class). */
    struct THEA_API HierarchicalMethod
    {
      /** Supported values. */
      enum Value
      {
        AUTO,                 ///< Automatically choose an appropriate method.

#ifdef THEA_ENABLE_CLUTO
        CLUTO_CLUSTER,        ///< CLUTO agglomerative clustering.
        CLUTO_CLUSTER_BIASED  ///< CLUTO biased agglomerative clustering.
#endif
      };

      THEA_ENUM_CLASS_BODY(HierarchicalMethod)
    };

    /** Methods to build a hierarchy on top of an existing clustering (enum class). */
    struct THEA_API HierarchyOverlayMethod
    {
      /** Supported values. */
      enum Value
      {
        AUTO,  ///< Automatically choose an appropriate method.

#ifdef THEA_ENABLE_CLUTO
        CLUTO  ///< Use CLUTO.
#endif
      };

      THEA_ENUM_CLASS_BODY(HierarchyOverlayMethod)
    };

    /** Metric for measuring similarity between two points. */
    struct THEA_API SimilarityMeasure
    {
      /** Supported values. */
      enum Value
      {
        AUTO,               ///< Automatically choose an appropriate similarity measure.
        COSINE,             ///< Cosine function.
        CORRELATION_COEFF,  ///< Pearson's correlation coefficient.
        L2_DISTANCE,        ///< Negative euclidean distance.
        EXTENDED_JACCARD    ///< Extended Jaccard coefficient.
      };

      THEA_ENUM_CLASS_BODY(SimilarityMeasure)
    };

    /** How to construct the graph for graph-based clustering. */
    struct THEA_API GraphModel
    {
      /** Supported values. */
      enum Value
      {
        AUTO,  ///< Automatically choose an appropriate graph model.

        /**
         * Add edge between two points if each is in the nearest-neighbor list of the other. The weight of the edge is the
         * similarity of the points.
         */
        SYMMETRIC_DIRECT,

        /**
         * Add edge between two points if each is in the nearest-neighbor list of the other. The weight of the edge is the
         * similarity of the points.
         */
        ASYMMETRIC_DIRECT,

        /**
         * Add edge between two points if each is in the nearest-neighbor list of the other. The weight of the edge is the
         * number of common elements in the nearest neighbor lists of the two points.
         */
        SYMMETRIC_LINK,

        /**
         * Add edge between two points if each is in the nearest-neighbor list of the other. The weight of the edge is the
         * number of common elements in the nearest neighbor lists of the two points.
         */
        ASYMMETRIC_LINK
      };

      THEA_ENUM_CLASS_BODY(GraphModel)
    };

    /** The function to optimize over all clusters. */
    struct THEA_API ClusteringCriterion
    {
      /** Supported values. */
      enum Value
      {
        AUTO,           ///< Automatically choose an appropriate clustering criterion.

#ifdef THEA_ENABLE_CLUTO
        CLUTO_I1,       ///< The I1 from the CLUTO paper.
        CLUTO_I2,       ///< The I2 from the CLUTO paper.
        CLUTO_E1,       ///< The E1 from the CLUTO paper.
        CLUTO_G1,       ///< The G1 from the CLUTO paper.
        CLUTO_G1P,      ///< The G1' from the CLUTO paper.
        CLUTO_H1,       ///< The H1 from the CLUTO paper.
        CLUTO_H2,       ///< The H2 from the CLUTO paper.
        CLUTO_SLINK,    ///< The traditional single-link/MST method.
        CLUTO_SLINK_W,  ///< The traditional single-link/MST method, cluster size weighted.
        CLUTO_CLINK,    ///< The traditional complete-link method.
        CLUTO_CLINK_W,  ///< The traditional complete-link method, cluster size weighted.
        CLUTO_UPGMA,    ///< The traditional UPGMA method.
        CLUTO_CUT,      ///< <em>[Graph clustering only]</em> Edge-cut based.
        CLUTO_RCUT,     ///< <em>[Graph clustering only]</em> Ratio cut.
        CLUTO_NCUT,     ///< <em>[Graph clustering only]</em> Normalized cut.
        CLUTO_MMCUT     ///< <em>[Graph clustering only]</em> Min-Max cut.
#endif
      };

      THEA_ENUM_CLASS_BODY(ClusteringCriterion)
    };

    /** Priority order for choosing the next cluster to split (only for top-down clustering). */
    struct THEA_API SplitPriority
    {
      /** Supported values. */
      enum Value
      {
        AUTO,                    ///< Automatically choose an appropriate splitting priority.
        LARGEST_FIRST,           ///< Split the largest cluster.
        BEST_FIRST,              ///< Split the cluster which will most improve the function being optimized.
        LARGEST_SUBSPACE_FIRST,  ///< Split the cluster whose subspace dimensions will be most reduced by the split.
      };

      THEA_ENUM_CLASS_BODY(SplitPriority)
    };

    /** %Options for flat clustering. */
    struct THEA_API FlatOptions
    {
      FlatMethod           method;
      SimilarityMeasure    similarity_measure;
      GraphModel           graph_model;
      long                 num_nearest_neighbors;
      ClusteringCriterion  clustering_criterion;
      SplitPriority        split_priority;

      /** Constructor. */
      FlatOptions(FlatMethod method_ = FlatMethod::AUTO,
                  SimilarityMeasure similarity_measure_ = SimilarityMeasure::AUTO,
                  GraphModel graph_model_ = GraphModel::AUTO,
                  long num_nearest_neighbors_ = -1,
                  ClusteringCriterion clustering_criterion_ = ClusteringCriterion::AUTO,
                  SplitPriority split_priority_ = SplitPriority::AUTO)
      : method(method_),
        similarity_measure(similarity_measure_),
        graph_model(graph_model_),
        num_nearest_neighbors(num_nearest_neighbors_),
        clustering_criterion(clustering_criterion_),
        split_priority(split_priority_)
      {}

      /** Default options. */
      static FlatOptions const & defaults() { static FlatOptions const options; return options; }
    };

    /**
     * Compute a flat clustering for a set of points. The class <code>T</code> should be a basic numeric type
     * (int/float/double...), Vector2, Vector3, Vector4 or VectorN.
     *
     * @param points Input points.
     * @param labels %Clustering results, each entry is the label of the cluster containing the point.
     * @param options Options for the clustering algorithm.
     * @param num_clusters_hint A hint (which may be ignored) about the likely number of clusters. For k-means this will be used
     *   as the value of k.
     *
     * @return The number of clusters found.
     */
    template <typename T>
    static int computeFlat(TheaArray<T> const & points, TheaArray<int> & labels, FlatOptions options = FlatOptions::defaults(),
                           int num_clusters_hint = 2)
    {
      if (points.size() <= 0)
        return 0;

      int num_clusters_found = 0;

#ifdef THEA_ENABLE_CLUTO
      if (options.method == FlatMethod::AUTO)
        options.method = FlatMethod::CLUTO_GRAPH_CLUSTER_RB;
#endif

      switch (options.method)
      {
#ifdef THEA_ENABLE_CLUTO
        case FlatMethod::CLUTO_CLUSTER_DIRECT:
        case FlatMethod::CLUTO_CLUSTER_RB:
        case FlatMethod::CLUTO_GRAPH_CLUSTER_RB:
        {
          int simfun = 0, grmodel = 0, nnbrs = 0, crfun = 0, cstype = 0;
          toClutoOptions(options, simfun, grmodel, nnbrs, crfun, cstype);

          ClutoMatrix cm;
          toClutoMatrix(points, cm);

          labels.resize(points.size());

          switch (options.method)
          {
            // Refer to the CLUTO manual and to Zhao and Karypis, "Criterion Functions for Document Clustering: Experiments and
            // Analysis", UMN CS 01-040, 2001, for choice of criterion function in each case
            case FlatMethod::CLUTO_CLUSTER_DIRECT:
            {
              CLUTO_VP_ClusterDirect((int)cm.numRows(),                    // nrows
                                     (int)cm.numColumns(),                 // ncols
                                     NULL,                                 // rowptr
                                     NULL,                                 // rowind
                                     &cm.getValues()[0],                   // rowval
                                     simfun,                               // simfun
                                     crfun,                                // crfun
                                     CLUTO_ROWMODEL_NONE,                  // rowmodel
                                     CLUTO_COLMODEL_NONE,                  // colmodel
                                     1.0f,                                 // colprune
                                     10,                                   // ntrials
                                     10,                                   // niter
                                     Random::common().integer(0, 0xFFFF),  // seed
                                     0,                                    // dbglvl
                                     num_clusters_hint,                    // nclusters
                                     &labels[0]);                          // part

              num_clusters_found = num_clusters_hint;
              break;
            }

            case FlatMethod::CLUTO_CLUSTER_RB:
            {
              CLUTO_VP_ClusterRB((int)cm.numRows(),                    // nrows
                                 (int)cm.numColumns(),                 // ncols
                                 NULL,                                 // rowptr
                                 NULL,                                 // rowind
                                 &cm.getValues()[0],                   // rowval
                                 simfun,                               // simfun
                                 crfun,                                // crfun
                                 CLUTO_ROWMODEL_NONE,                  // rowmodel
                                 CLUTO_COLMODEL_NONE,                  // colmodel
                                 1.0f,                                 // colprune
                                 10,                                   // ntrials
                                 10,                                   // niter
                                 Random::common().integer(0, 0xFFFF),  // seed
                                 cstype,                               // cstype
                                 1,                                    // kwayrefine
                                 0,                                    // dbglvl
                                 num_clusters_hint,                    // nclusters
                                 &labels[0]);                          // part

              num_clusters_found = num_clusters_hint;
              break;
            }

            default:  // FlatMethod::CLUTO_GRAPH_CLUSTER_RB
            {
              float cut_value;
              num_clusters_found = CLUTO_VP_GraphClusterRB((int)cm.numRows(),                    // nrows
                                                           (int)cm.numColumns(),                 // ncols
                                                           NULL,                                 // rowptr
                                                           NULL,                                 // rowind
                                                           &cm.getValues()[0],                   // rowval
                                                           simfun,                               // simfun
                                                           CLUTO_ROWMODEL_NONE,                  // rowmodel
                                                           CLUTO_COLMODEL_NONE,                  // colmodel
                                                           1.0f,                                 // colprune
                                                           grmodel,                              // grmodel
                                                           nnbrs,                                // nnbrs
                                                           -1,                                   // edgeprune
                                                           -1,                                   // vtxprune
                                                           5,                                    // mincmp
                                                           10,                                   // ntrials
                                                           Random::common().integer(0, 0xFFFF),  // seed
                                                           cstype,                               // cstype
                                                           0,                                    // dbglvl
                                                           num_clusters_hint,                    // nclusters
                                                           &labels[0],                           // part
                                                           &cut_value);                          // crvalue
            }
          }
          break;
        }
#endif

        default: throw Error("Clustering: Clustering method is not supported");
      }

      return num_clusters_found;
    }

    /**
     * Compute a vertex clustering for a user-specified graph. The class <code>MatrixT</code> should be convertible to
     * CompressedRowMatrix<float, int, int>.
     *
     * @param weights Input adjacency matrix.
     * @param labels %Clustering results, each entry is the label of the cluster containing the point.
     * @param options Options for the clustering algorithm.
     * @param num_clusters_hint A hint (which may be ignored) about the likely number of clusters. For k-means this will be used
     *   as the value of k.
     *
     * @return The number of clusters found.
     */
    template <typename MatrixT>
    static int computeFlat(MatrixT const & weights, TheaArray<int> & labels, FlatOptions options = FlatOptions::defaults(),
                           int num_clusters_hint = 2)
    {
      int num_clusters_found = 0;

#ifdef THEA_ENABLE_CLUTO
      if (options.method == FlatMethod::AUTO)
        options.method = FlatMethod::CLUTO_GRAPH_CLUSTER_RB;
#endif

      switch (options.method)
      {
#ifdef THEA_ENABLE_CLUTO
        case FlatMethod::CLUTO_CLUSTER_DIRECT:
        case FlatMethod::CLUTO_CLUSTER_RB:
        case FlatMethod::CLUTO_GRAPH_CLUSTER_RB:
        {
          alwaysAssertM(weights.numRows() == weights.numColumns(), "Clustering: Adjacency matrix must be square");

          int simfun = 0, grmodel = 0, nnbrs = 0, crfun = 0, cstype = 0;
          toClutoOptions(options, simfun, grmodel, nnbrs, crfun, cstype);

          ClutoMatrix cm;
          toClutoMatrix(weights, cm);
          if (cm.getValues().size() <= 0)  // more zeroes could have been detected during conversion
          {
            THEA_WARNING << "Clustering: Weight matrix is zero, no clusters exist";
            return 0;
          }

          int   * rowptr = (cm.getRowIndices().size()    > 0 ? &cm.getRowIndices()[0]    : NULL);
          int   * rowind = (cm.getColumnIndices().size() > 0 ? &cm.getColumnIndices()[0] : NULL);
          float * rowval = (cm.getValues().size()        > 0 ? &cm.getValues()[0]        : NULL);

          labels.resize((size_t)weights.numRows());

          switch (options.method)
          {
            // Refer to the CLUTO manual and to Zhao and Karypis, "Criterion Functions for Document Clustering: Experiments and
            // Analysis", UMN CS 01-040, 2001, for choice of criterion function in each case
            case FlatMethod::CLUTO_CLUSTER_DIRECT:
            {
              CLUTO_SP_ClusterDirect((int)weights.numRows(),               // nrows
                                     rowptr,                               // rowptr
                                     rowind,                               // rowind
                                     rowval,                               // rowval
                                     crfun,                                // crfun
                                     10,                                   // ntrials
                                     10,                                   // niter
                                     Random::common().integer(0, 0xFFFF),  // seed
                                     0,                                    // dbglvl
                                     num_clusters_hint,                    // nclusters
                                     &labels[0]);                          // part

              num_clusters_found = num_clusters_hint;
              break;
            }

            case FlatMethod::CLUTO_CLUSTER_RB:
            {
              CLUTO_SP_ClusterRB((int)weights.numRows(),               // nrows
                                 rowptr,                               // rowptr
                                 rowind,                               // rowind
                                 rowval,                               // rowval
                                 crfun,                                // crfun
                                 10,                                   // ntrials
                                 10,                                   // niter
                                 Random::common().integer(0, 0xFFFF),  // seed
                                 cstype,                               // cstype
                                 1,                                    // kwayrefine
                                 0,                                    // dbglvl
                                 num_clusters_hint,                    // nclusters
                                 &labels[0]);                          // part

              num_clusters_found = num_clusters_hint;
              break;
            }

            default:  // FlatMethod::CLUTO_GRAPH_CLUSTER_RB
            {
              float cut_value;
              num_clusters_found = CLUTO_SP_GraphClusterRB((int)weights.numRows(),               // nrows
                                                           rowptr,                               // rowptr
                                                           rowind,                               // rowind
                                                           rowval,                               // rowval
                                                           nnbrs,                                // nnbrs
                                                           -1,                                   // edgeprune
                                                           -1,                                   // vtxprune
                                                           5,                                    // mincmp
                                                           10,                                   // ntrials
                                                           Random::common().integer(0, 0xFFFF),  // seed
                                                           cstype,                               // cstype
                                                           0,                                    // dbglvl
                                                           num_clusters_hint,                    // nclusters
                                                           &labels[0],                           // part
                                                           &cut_value);                          // crvalue
            }
          }
          break;
        }
#endif

        default: throw Error("Clustering method is not supported");
      }

      return num_clusters_found;
    }

    /**
     * Compute a hierarchical clustering that respects a given flat clustering. This is currently restricted to producing
     * binary trees.
     */
    template <typename T>
    static void overlayHierarchy(TheaArray<T> const & points, TheaArray<int> const & flat_labels, int num_clusters,
                                 TheaArray<int> & tree_labels, HierarchyOverlayMethod method = HierarchyOverlayMethod::AUTO)
    {
      switch (method)
      {
        default: throw Error("Clustering method is not supported");
      }
    }

    /**
     * Compute a hierarchical clustering that respects a given flat clustering. This is currently restricted to producing
     * binary trees. The class <code>Matrix</code> should be a standard matrix whose entries represent edge weights (0 if the
     * edge is not present).
     */
    template <typename MatrixT>
    static void overlayHierarchy(MatrixT const & weights, TheaArray<int> const & flat_labels, int num_clusters,
                                 TheaArray<int> & tree_labels, HierarchyOverlayMethod method = HierarchyOverlayMethod::AUTO)
    {
      switch (method)
      {
        default: throw Error("Clustering method is not supported");
      }
    }

  private:

#ifdef THEA_ENABLE_CLUTO
    /** Get the number of dimensions of a vector. The dummy argument is needed for template resolution. */
    template <typename T>            static size_t numElems (T const & t);
    template <long N, typename T>  static size_t numElems (VectorN<N, T> const & v);
    template <typename T>            static size_t numElems (TheaArray<T> const & v);

    /** Copy a vector to a float array and return a pointer to the next array position. */
    template <typename T>            static float * copyToFloatArray(T const & t, float * fp);
    template <long N, typename T>  static float * copyToFloatArray(VectorN<N, T> const & v, float * fp);
    template <typename T>            static float * copyToFloatArray(TheaArray<T> const & v, float * fp);

    /** Convert an array of points to a CLUTO matrix. */
    template <typename T> static void toClutoMatrix(TheaArray<T> const & points, ClutoMatrix & cm);

    /**
     * Convert an array of points to a CLUTO matrix, where each point is itself an array (all points must have same dimensions).
     */
    template <typename T> static void toClutoMatrix(TheaArray< TheaArray<T> > const & points, ClutoMatrix & cm);

    /** Convert a dense matrix row-major to a CLUTO dense matrix. */
    template <typename T> static void toClutoMatrix(Matrix<T, MatrixLayout::ROW_MAJOR> const & m, ClutoMatrix & cm);

    /** Convert a dense column-major matrix to a CLUTO dense matrix. */
    template <typename T> static void toClutoMatrix(Matrix<T, MatrixLayout::COLUMN_MAJOR> const & m, ClutoMatrix & cm);

    /** Convert a sparse matrix to a CLUTO matrix. */
    template <typename MatrixT> static void toClutoMatrix(MatrixT const & m, ClutoMatrix & cm);

    /** Convert standard flat clustering options to CLUTO options. */
    static void toClutoOptions(FlatOptions const & options, int & simfun, int & grmodel, int & nnbrs, int & crfun,
                               int & cstype);

    /** Convert a standard similarity measure to its CLUTO equivalent. */
    static int toClutoEnum(SimilarityMeasure sim);

    /** Convert a standard graph model to its CLUTO equivalent. */
    static int toClutoEnum(GraphModel grmodel);

    /** Convert a standard clustering criterion to its CLUTO equivalent. */
    static int toClutoEnum(ClusteringCriterion sim);

    /** Convert a standard splitting priority to its CLUTO equivalent. */
    static int toClutoEnum(SplitPriority sim);
#endif

}; // class Clustering

#ifdef THEA_ENABLE_CLUTO

//=============================================================================================================================
// Get the number of dimensions of a vector. The dummy argument is needed for template resolution.
//=============================================================================================================================

template <typename T>            size_t Clustering::numElems         (T const & t)             { return 1; }
template <> inline               size_t Clustering::numElems<Vector2>(Vector2 const & t)       { return 2; }
template <> inline               size_t Clustering::numElems<Vector3>(Vector3 const & t)       { return 3; }
template <> inline               size_t Clustering::numElems<Vector4>(Vector4 const & t)       { return 4; }
template <long N, typename T>    size_t Clustering::numElems         (VectorN<N, T> const & v) { return (size_t)N; }
template <typename T>            size_t Clustering::numElems         (TheaArray<T> const & v)  { return v.size(); }

//=============================================================================================================================
// Copy a vector to a float array and return a pointer to the next array position.
//=============================================================================================================================

template <typename T> float * Clustering::copyToFloatArray(T const & t, float * fp)
{ *fp = static_cast<float>(t); return fp + 1; }  // plain numbers
template <> inline float * Clustering::copyToFloatArray<Vector2>(Vector2 const & t, float * fp)
{ fp[0] = t[0]; fp[1] = t[1]; return fp + 2; }
template <> inline float * Clustering::copyToFloatArray<Vector3>(Vector3 const & t, float * fp)
{ fp[0] = t[0]; fp[1] = t[1]; fp[2] = t[2]; return fp + 3; }
template <> inline float * Clustering::copyToFloatArray<Vector4>(Vector4 const & t, float * fp)
{ fp[0] = t[0]; fp[1] = t[1]; fp[2] = t[2]; fp[3] = t[3]; return fp + 4; }

template <long N, typename T>
float *
Clustering::copyToFloatArray(VectorN<N, T> const & v, float * fp)
{
  for (size_t i = 0; i < N; ++i)
    fp[i] = static_cast<float>(v[i]);

  return fp + N;
}

template <typename T>
float *
Clustering::copyToFloatArray(TheaArray<T> const & v, float * fp)
{
  for (typename TheaArray<T>::const_iterator vi = v.begin(); vi != v.end(); ++vi, ++fp)
    *fp = static_cast<float>(*vi);

  return fp;
}

//=============================================================================================================================
// Store an array of points in a compressed row matrix.
//=============================================================================================================================

template <typename T>
void
Clustering::toClutoMatrix(TheaArray<T> const & points, ClutoMatrix & cm)
{
  size_t num_points = points.size();
  size_t num_elems  = numElems(points[0]);

  cm.resize((int)num_points, (int)num_elems);
  cm.getRowIndices().clear();
  cm.getValues().resize(num_elems * num_points);

  float * fp = &cm.getValues()[0];
  for (size_t i = 0; i < num_points; ++i)
    fp = copyToFloatArray(points[i], fp);
}

template <typename T>
void
Clustering::toClutoMatrix(TheaArray< TheaArray<T> > const & points, ClutoMatrix & cm)
{
  size_t num_points = points.size();
  size_t num_elems  = points[0].size();

  cm.resize((int)num_points, (int)num_elems);
  cm.getRowIndices().clear();
  cm.getValues().resize(num_elems * num_points);

  float * fp = &cm.getValues()[0];
  for (size_t i = 0; i < num_points; ++i)
  {
    if (points[i].size() != num_elems)
      throw Error("Clustering: Points do not have same dimensions");

    fp = copyToFloatArray(points[i], fp);
  }
}

//=============================================================================================================================
// Convert a dense matrix to a CLUTO dense matrix
//=============================================================================================================================

template <typename T>
void
Clustering::toClutoMatrix(Matrix<T, MatrixLayout::ROW_MAJOR> const & m, ClutoMatrix & cm)
{
  cm.getRowIndices().clear();
  cm.getColumnIndices().clear();
  cm.getValues() = TheaArray<float>(m.begin(), m.end());
}

template <typename T>
void
Clustering::toClutoMatrix(Matrix<T, MatrixLayout::COLUMN_MAJOR> const & m, ClutoMatrix & cm)
{
  cm.getRowIndices().clear();
  cm.getColumnIndices().clear();

  if (!m.isEmpty())
  {
    cm.getValues().resize(m.numRows() * m.numColumns());

    TheaArray<float>::iterator cp = cm.getValues().begin();
    for (int i = 0; i < m.numRows(); ++i)
      for (int j = 0; j < m.numColumns(); ++j, ++cp)
        *cp = m(j, i);
  }
  else
    cm.getValues().clear();
}

//=============================================================================================================================
// Convert a sparse matrix to a CLUTO sparse matrix
//=============================================================================================================================

template <typename MatrixT>
void
Clustering::toClutoMatrix(MatrixT const & m, ClutoMatrix & cm)
{
  cm = ClutoMatrix(m);
}

//=============================================================================================================================
// Convert standard options to CLUTO options
//=============================================================================================================================

inline int
Clustering::toClutoEnum(SimilarityMeasure sim)
{
  switch (sim)
  {
    case SimilarityMeasure::COSINE             :  return CLUTO_SIM_COSINE;
    case SimilarityMeasure::CORRELATION_COEFF  :  return CLUTO_SIM_CORRCOEF;
    case SimilarityMeasure::L2_DISTANCE        :  return CLUTO_SIM_EDISTANCE;
    case SimilarityMeasure::EXTENDED_JACCARD   :  return CLUTO_SIM_EJACCARD;
    default                                    :  throw Error("Clustering: Invalid similarity measure");
  }
}

inline int
Clustering::toClutoEnum(GraphModel grmodel)
{
  switch (grmodel)
  {
    case GraphModel::SYMMETRIC_DIRECT   :  return CLUTO_GRMODEL_SYMETRIC_DIRECT;
    case GraphModel::ASYMMETRIC_DIRECT  :  return CLUTO_GRMODEL_ASYMETRIC_DIRECT;
    case GraphModel::SYMMETRIC_LINK     :  return CLUTO_GRMODEL_SYMETRIC_LINKS;
    case GraphModel::ASYMMETRIC_LINK    :  return CLUTO_GRMODEL_ASYMETRIC_LINKS;
    default                             :  throw Error("Clustering: Invalid graph model");
  }
}

inline int
Clustering::toClutoEnum(ClusteringCriterion cc)
{
  switch (cc)
  {
    case ClusteringCriterion::CLUTO_I1       :  return CLUTO_CLFUN_I1;
    case ClusteringCriterion::CLUTO_I2       :  return CLUTO_CLFUN_I2;
    case ClusteringCriterion::CLUTO_E1       :  return CLUTO_CLFUN_E1;
    case ClusteringCriterion::CLUTO_G1       :  return CLUTO_CLFUN_G1;
    case ClusteringCriterion::CLUTO_G1P      :  return CLUTO_CLFUN_G1P;
    case ClusteringCriterion::CLUTO_H1       :  return CLUTO_CLFUN_H1;
    case ClusteringCriterion::CLUTO_H2       :  return CLUTO_CLFUN_H2;
    case ClusteringCriterion::CLUTO_SLINK    :  return CLUTO_CLFUN_SLINK;
    case ClusteringCriterion::CLUTO_SLINK_W  :  return CLUTO_CLFUN_SLINK_W;
    case ClusteringCriterion::CLUTO_CLINK    :  return CLUTO_CLFUN_CLINK;
    case ClusteringCriterion::CLUTO_CLINK_W  :  return CLUTO_CLFUN_CLINK_W;
    case ClusteringCriterion::CLUTO_UPGMA    :  return CLUTO_CLFUN_UPGMA;
    case ClusteringCriterion::CLUTO_CUT      :  return CLUTO_CLFUN_CUT;
    case ClusteringCriterion::CLUTO_RCUT     :  return CLUTO_CLFUN_RCUT;
    case ClusteringCriterion::CLUTO_NCUT     :  return CLUTO_CLFUN_NCUT;
    case ClusteringCriterion::CLUTO_MMCUT    :  return CLUTO_CLFUN_MMCUT;
    default                                  :  throw Error("Clustering: Invalid clustering criterion");
  }
}

inline int
Clustering::toClutoEnum(SplitPriority sp)
{
  switch (sp)
  {
    case SplitPriority::LARGEST_FIRST           :  return CLUTO_CSTYPE_LARGEFIRST;
    case SplitPriority::BEST_FIRST              :  return CLUTO_CSTYPE_BESTFIRST;
    case SplitPriority::LARGEST_SUBSPACE_FIRST  :  return CLUTO_CSTYPE_LARGESUBSPACEFIRST;
    default                                     :  throw Error("Clustering: Invalid split priority");
  }
}

inline void
Clustering::toClutoOptions(FlatOptions const & options, int & simfun, int & grmodel, int & nnbrs, int & crfun, int & cstype)
{
  debugAssertM(options.method == FlatMethod::CLUTO_CLUSTER_DIRECT
            || options.method == FlatMethod::CLUTO_CLUSTER_RB
            || options.method == FlatMethod::CLUTO_GRAPH_CLUSTER_RB, "Clustering: Options do not specify CLUTO clustering");

  if (options.method != FlatMethod::CLUTO_GRAPH_CLUSTER_RB
   && (options.clustering_criterion == ClusteringCriterion::CLUTO_CUT
    || options.clustering_criterion == ClusteringCriterion::CLUTO_RCUT
    || options.clustering_criterion == ClusteringCriterion::CLUTO_NCUT
    || options.clustering_criterion == ClusteringCriterion::CLUTO_MMCUT))
     throw Error("Clustering: Clustering criterion is only valid for CLUTO graph clustering");

  simfun = (options.similarity_measure == SimilarityMeasure::AUTO ? toClutoEnum(SimilarityMeasure::L2_DISTANCE)
                                                                  : toClutoEnum(options.similarity_measure));

  grmodel = (options.graph_model == GraphModel::AUTO ? toClutoEnum(GraphModel::SYMMETRIC_DIRECT)
                                                     : toClutoEnum(options.graph_model));

  nnbrs = (options.num_nearest_neighbors <= 0 ? 40 : (int)options.num_nearest_neighbors);

  crfun = (options.clustering_criterion == ClusteringCriterion::AUTO ? toClutoEnum(ClusteringCriterion::CLUTO_H2)
                                                                     : toClutoEnum(options.clustering_criterion));

  cstype = (options.split_priority == SplitPriority::AUTO ? toClutoEnum(SplitPriority::BEST_FIRST)
                                                          : toClutoEnum(options.split_priority));
}

#endif

} // namespace Algorithms
} // namespace Thea

#endif
