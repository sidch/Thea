//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_KMeans_hpp__
#define __Thea_Algorithms_KMeans_hpp__

#include "../Common.hpp"
#include "../AtomicInt32.hpp"
#include "../IOStream.hpp"
#include "../Math.hpp"
#include "../Matrix.hpp"
#include "../Random.hpp"
#include "../Serializable.hpp"
#include "../Stopwatch.hpp"
#include "../System.hpp"
#include <boost/thread.hpp>
#include <algorithm>

namespace Thea {
namespace Algorithms {

/** A k-means model for clustering points described by feature vectors. */
class THEA_API KMeans : public Serializable
{
  public:
    THEA_DEF_POINTER_TYPES(KMeans, shared_ptr, weak_ptr)

    /** How to seed the initial centers (enum class). */
    class THEA_API Seeding
    {
      public:
        /** Supported values. */
        enum Value
        {
          RANDOM,             ///< Choose random input points as initial seeds (string: "random").
          K_MEANS_PLUS_PLUS   ///< K-means++ strategy of Arthur/Vassilvitskii '07 (string: "k-means++").
        };

        THEA_ENUM_CLASS_BODY(Seeding)

        THEA_ENUM_CLASS_STRINGS_BEGIN(Seeding)
          THEA_ENUM_CLASS_STRING(RANDOM, "random")
          THEA_ENUM_CLASS_STRING(K_MEANS_PLUS_PLUS, "k-means++")
        THEA_ENUM_CLASS_STRINGS_END(Seeding)

    }; // struct Seeding

    /** Clustering options. */
    /**
     * %Options for the classifier. In most cases, passing a negative value for a normally non-negative parameter auto-selects a
     * suitable value for that parameter.
     */
    class THEA_API Options : public Serializable
    {
      public:
        /** Constructor. Sets default options. */
        Options();

        /* Maximum iterations to seek convergence (ignored if negative, default -1). */
        Options & setMaxIterations(long iters) { max_iterations = iters; return *this; }

        /* Maximum time in seconds to seek convergence (ignored if negative, default -1). */
        Options & setMaxTime(double seconds) { max_time = seconds; return *this; }

        /** How to seed the initial centers (default Seeding::K_MEANS_PLUS_PLUS). */
        Options & setSeeding(Seeding seeding_) { seeding = seeding_; return *this; }

        /** Accelerate computations by parallelization or not (default true). */
        Options & setParallelize(bool value) { parallelize = value; return *this; }

        /** Set whether progress information will be printed to the console or not (default true). */
        Options & setVerbose(bool value) { verbose = value; return *this; }

        /** Load options from a disk file. */
        bool load(std::string const & path);

        /** Save options to a disk file. */
        bool save(std::string const & path) const;

        void deserialize(BinaryInputStream & in, Codec const & codec = Codec_AUTO());
        void serialize(BinaryOutputStream & out, Codec const & codec = Codec_AUTO()) const;
        void deserialize(TextInputStream & in, Codec const & codec = Codec_AUTO());
        void serialize(TextOutputStream & out, Codec const & codec = Codec_AUTO()) const;

        /** Get the set of default options. */
        static Options const & defaults() { static Options const def; return def; }

      private:
        long max_iterations;  ///< Maximum iterations to seek convergence (ignored if negative).
        double max_time;      ///< Maximum time in seconds to seek convergence (ignored if negative).
        Seeding seeding;      ///< How to seed the initial centers.
        bool parallelize;     ///< Accelerate computations by parallelization.
        bool verbose;         ///< Print progress information to the console.

        friend class KMeans;

    }; // class Options

    /** Constructor. */
    KMeans(Options const & options_ = Options::defaults()) : options(options_) {}

    /** Get the current set of options. */
    Options const & getOptions() const { return options; }

    /** Set the current set of options. */
    void setOptions(Options const & options_) { options = options_; }

    /** Get the number of clusters (the value of k). */
    long numClusters() const { return centers.numRows(); }

    /** Get the size of the feature vector of a point or cluster center. */
    long numPointFeatures() const { return centers.numColumns(); }

    /**
     * Get the features of a cluster center. The vector has numPointFeatures() elements.
     *
     * @param i The index of the cluster, in the range [0, numClusters() - 1].
     */
    double const * getClusterCenter(long i) const { return &centers(i, 0); }

    /**
     * Cluster a set of points into \a num_clusters groups, one per row of the input matrix.
     *
     * @param num_clusters Number of clusters to cluster points into (the value of k).
     * @param points The set of points, one vector per row. AddressableMatrixT should have the interface of an AddressableMatrix
     *   of some real-valued scalar type.
     * @param labeling [Optional] Used to return the index of the cluster containing each point. Assumed to be preallocated to
     *   the number of points. Ignored if null.
     * @param squared_distances [Optional] Used to return the square of the distance to the center of the cluster containing
     *   each point. Assumed to be preallocated to the number of points. Ignored if null.
     *
     * @return True if the clustering converged, else false.
     */
    template <typename AddressableMatrixT>
    bool cluster(long num_clusters, AddressableMatrixT const & points, long * labeling = NULL,
                 double * squared_distances = NULL)
    {
      alwaysAssertM(num_clusters > 0, "KMeans: Number of clusters must be at least 1");

      long num_points = points.numRows();
      long num_features = points.numColumns();

      alwaysAssertM(num_points > 0, "KMeans: Number of points must be at least 1");
      alwaysAssertM(num_features > 0, "KMeans: Number of features must be at least 1");

      num_clusters = std::min(num_clusters, num_points);

      if (options.verbose)
        THEA_CONSOLE << "KMeans: Clustering " << num_points << " point(s) into " << num_clusters << " group(s)";

      //=======================================================================================================================
      // Select initial centers
      //=======================================================================================================================

      double start_time = System::time();
      Stopwatch timer;
      if (options.verbose)
        timer.tick();

      if (!selectInitialCenters(num_clusters, points))
        return false;

      if (options.verbose)
      {
        timer.tock();
        THEA_CONSOLE << "KMeans: Selected initial center(s) in " << timer.elapsedTime() << 's';
        timer.tick();
      }

      //=======================================================================================================================
      // Find initial assignment of points to clusters
      //=======================================================================================================================

      TheaArray<long> labeling_local;
      if (!labeling)
      {
        labeling_local.resize((array_size_t)num_points, -1);
        labeling = &labeling_local[0];
      }

      TheaArray<double> sqdist_local;
      if (!squared_distances)
      {
        sqdist_local.resize((array_size_t)num_points);
        squared_distances = &sqdist_local[0];
      }

      mapToClusters(num_clusters, points, labeling, squared_distances);

      //=======================================================================================================================
      // Iterate
      //=======================================================================================================================

      bool converged = false;
      long num_iterations = 0;
      double last_time = System::time();
      while (!converged)
      {
        if (options.max_iterations >= 0 && num_iterations >= options.max_iterations)
          break;

        double curr_time = System::time();
        if (options.max_time >= 0 && curr_time - start_time > options.max_time)
          break;

        // Print an update every 3 seconds or every iteration, whichever is longer
        if (options.verbose && curr_time - last_time > 3)
        {
          THEA_CONSOLE << "KMeans: -- " << num_iterations << " iteration(s)";
          last_time = curr_time;
        }

        updateClusters(points, labeling, squared_distances);  // FIXME: Check for reassignments and || with converged? Unclear.
        converged = !mapToClusters(num_clusters, points, labeling, squared_distances);
        num_iterations++;
      }

      if (options.verbose)
      {
        timer.tock();

        if (converged)
          THEA_CONSOLE << "KMeans: Clusters found after " << num_iterations << " iteration(s) in " << timer.elapsedTime()
                       << 's';
        else
          THEA_CONSOLE << "KMeans: Did not converge after " << num_iterations << " iteration(s) in " << timer.elapsedTime()
                       << 's';
      }

      return converged;
    }

    /**
     * Map a point to its containing cluster. To map many points, use mapToClusters() instead which can parallelize the queries
     * for speed.
     *
     * @param point The feature vector of the point, with numPointFeatures() entries.
     * @param squared_distance [Optional] Used to return the square of the distance to the closest cluster center.
     *
     * @return The index of the cluster containing the point, in the range [0, numClusters() - 1], or a negative value if no
     *   such cluster was found.
     */
    long mapToCluster(double const * point, double * squared_distance = NULL) const
    {
      long label;
      double sqdist;
      mapToCluster(numClusters(), point, label, sqdist);

      if (squared_distance) *squared_distance = sqdist;
      return label;
    }

    /**
     * Map a set of points to their containing clusters. Parallelizes queries if this was enabled in the options.
     *
     * @param points The set of points, one vector per row. AddressableMatrixT should have the interface of an AddressableMatrix
     *   of some real-valued scalar type.
     * @param labeling [Optional] Used to return the index of the cluster containing each point. Assumed to be preallocated to
     *   the number of points. Ignored if null.
     * @param squared_distances [Optional] Used to return the square of the distance to the center of the cluster containing
     *   each point. Assumed to be preallocated to the number of points. Ignored if null.
     */
    template <typename AddressableMatrixT>
    void mapToClusters(AddressableMatrixT const & points, long * labeling, double * squared_distances = NULL) const
    {
      mapToClusters(numClusters(), points, labeling, squared_distances);
    }

    /** Load the k-means model from a disk file. */
    bool load(std::string const & path);

    /** Save the k-means model to a disk file. */
    bool save(std::string const & path) const;

    void deserialize(BinaryInputStream & in, Codec const & codec = Codec_AUTO());
    void serialize(BinaryOutputStream & out, Codec const & codec = Codec_AUTO()) const;
    void deserialize(TextInputStream & in, Codec const & codec = Codec_AUTO());
    void serialize(TextOutputStream & out, Codec const & codec = Codec_AUTO()) const;

  private:
    /** Select initial centers by k-means++. */
    template <typename AddressableMatrixT>
    bool selectInitialCenters(long num_clusters, AddressableMatrixT const & points)
    {
      long num_points = points.numRows();
      long num_features = points.numColumns();

      alwaysAssertM(num_clusters > 0, "KMeans: Must select at least one center");
      alwaysAssertM(num_points >= num_clusters, "KMeans: Cannot select more centers than points");
      alwaysAssertM(num_features > 0, "KMeans: Cannot select centers without any point features");

      if (options.seeding == Seeding::K_MEANS_PLUS_PLUS)
      {
        // Choose initial cluster centers by k-means++ [Arthur/Vassilvitskii '07]:
        // 1. Choose one center uniformly at random from among the data points.
        // 2. For each data point x, compute D(x), the distance between x and the nearest center that has already been chosen.
        // 3. Choose one new data point at random as a new center, using a weighted probability distribution where a point x is
        //    chosen with probability proportional to D(x)^2.
        // 4. Repeat steps 2 and 3 until k centers have been chosen.

        centers.resize(num_clusters, num_features);
        centers.fill(0);

        // First center is randomly chosen
        long index = Random::common().integer(0, num_points - 1);
        addPointToCenter(points, index, 0);

        // Subsequent centers by k-means++
        TheaArray<double> sqdist((array_size_t)num_points);
        double start_time = System::time();
        for (long i = 1; i < num_clusters; ++i)
        {
          mapToClusters(i, points, NULL, &sqdist[0]);

          // Sample next center from points with probability proportional to sqdist
          double sum_sqdist = 0;
          for (array_size_t j = 0; j < sqdist.size(); ++j)
            sum_sqdist += sqdist[j];

          double r = Random::common().uniform(0, sum_sqdist);
          sum_sqdist = 0;
          index = num_points - 1;  // to compensate for numerical error when r is approximately = sum_sqdist
          for (array_size_t j = 0; j < sqdist.size(); ++j)
          {
            sum_sqdist += sqdist[j];
            if (sum_sqdist >= r)
            {
              index = (long)j;
              break;
            }
          }

          addPointToCenter(points, index, i);

          double curr_time = System::time();
          if (curr_time - start_time > 3)  // print an update every 3 seconds or every iteration, whichever is longer
          {
            THEA_CONSOLE << "KMeans: -- selected " << i + 1 << " center(s)";
            start_time = curr_time;
          }
        }
      }
      else if (options.seeding == Seeding::RANDOM)
      {
        centers.resize(num_clusters, num_features);
        centers.fill(0);

        TheaArray<long> indices((array_size_t)num_points);
        for (long i = 0; i < num_points; ++i)
          indices[(array_size_t)i] = i;

        // Select num_clusters random points as initial centers
        Random::common().randomShuffle((int32)num_points, (int32)num_clusters, &indices[0]);

        for (long i = 0; i < num_clusters; ++i)
          addPointToCenter(points, indices[(array_size_t)i], i);
      }
      else
      {
        THEA_ERROR << "KMeans: Unsupported seeding strategy " << options.seeding.toString();
        return false;
      }

      return true;
    }

    /** Worker class for parallelizing the mapping of points to their closest centers. */
    template <typename AddressableMatrixT> class ClusterMapper
    {
      public:
        /** Constructor. */
        ClusterMapper(KMeans const * parent_, long num_clusters_, AddressableMatrixT const * points_, long points_begin_,
                     long points_end_, long * cluster_indices_, double * cluster_sqdists_)
        : parent(parent_), num_clusters(num_clusters_), points(points_), points_begin(points_begin_), points_end(points_end_),
          cluster_indices(cluster_indices_), cluster_sqdists(cluster_sqdists_)
        {}

        /** Main function, called once per thread. */
        void operator()()
        {
          long index = -1;
          double sqdist = -1;
          bool changed = false;
          for (long i = points_begin; i < points_end; ++i)
          {
            parent->mapToCluster(num_clusters, *points, i, index, sqdist);

            if (cluster_indices)
            {
              changed = changed || (cluster_indices[i] != index);
              cluster_indices[i] = index;
            }

            if (cluster_sqdists)
              cluster_sqdists[i] = sqdist;
          }

          if (changed)
            parent->flag.increment();
        }

      private:
        KMeans const * parent;
        long num_clusters;
        AddressableMatrixT const * points;
        long points_begin;
        long points_end;
        long * cluster_indices;
        double * cluster_sqdists;

    }; // class ClusterMapper

    template <typename AddressableMatrixT> friend class ClusterMapper;

    /**
     * Map each point to its closest center.
     *
     * @return True if \a cluster_indices is non-null and the assignment to centers changed as a result of a call to this
     *   function, else false.
     */
    template <typename AddressableMatrixT>
    bool mapToClusters(long num_clusters, AddressableMatrixT const & points, long * cluster_indices,
                      double * cluster_sqdists = NULL) const
    {
      long num_points = points.numRows();
      unsigned int concurrency = boost::thread::hardware_concurrency();
      flag = 0;

      if (options.parallelize && concurrency > 1 && num_points > (long)(2 * concurrency))
      {
        boost::thread_group pool;
        double points_per_thread = num_points / (double)concurrency;

        long points_begin = 0;
        for (unsigned int i = 0; i < concurrency; ++i)
        {
          long points_end = std::min((long)Math::round(points_begin + points_per_thread), num_points);

          pool.add_thread(new boost::thread(ClusterMapper<AddressableMatrixT>(this,
                                                                             num_clusters,
                                                                             &points,
                                                                             points_begin,
                                                                             points_end,
                                                                             cluster_indices,
                                                                             cluster_sqdists)));
          points_begin = points_end;
        }

        pool.join_all();
      }
      else
      {
        ClusterMapper<AddressableMatrixT> mapper(this, num_clusters, &points, 0, num_points, cluster_indices, cluster_sqdists);
        mapper();
      }

      return (flag.value() > 0);
    }

    /** Map a point to its nearest center. */
    template <typename AddressableMatrixT>
    void mapToCluster(long num_clusters, AddressableMatrixT const & points, long point_index, long & cluster_index,
                     double & cluster_sqdist) const
    {
      long num_features = centers.numColumns();
      fvec.resize((array_size_t)num_features);
      points.getRow(point_index, &fvec[0]);

      mapToCluster(num_clusters, &fvec[0], cluster_index, cluster_sqdist);
    }

    /** Map a point to its nearest center. */
    void mapToCluster(long num_clusters, double const * point, long & cluster_index, double & cluster_sqdist) const
    {
      long num_features = centers.numColumns();
      cluster_index = -1;
      cluster_sqdist = -1;
      for (long i = 0; i < num_clusters; ++i)
      {
        // Compute squared distance to this center
        double sqdist = 0;
        for (long j = 0; j < num_features; ++j)
        {
          double diff = point[j] - centers(i, j);
          sqdist += (diff * diff);
        }

        if (i == 0 || sqdist < cluster_sqdist)
        {
          cluster_index = i;
          cluster_sqdist = sqdist;
        }
      }
    }

    /**
     * Update each center to be the centroid of its cluster.
     *
     * @return True if points were reassigned to fix empty clusters, else false.
     */
    template <typename AddressableMatrixT>
    bool updateClusters(AddressableMatrixT const & points, long * cluster_indices, double * cluster_sqdists)
    {
      long num_clusters = centers.numRows();
      long num_points = points.numRows();
      long num_features = points.numColumns();

      centers.fill(0);

      TheaArray<long> num_assigned((array_size_t)num_clusters, 0);
      for (long i = 0; i < num_points; ++i)
      {
        long cc_index = cluster_indices[i];

        addPointToCenter(points, i, cc_index);
        num_assigned[(array_size_t)cc_index]++;
      }

      // If a cluster is empty then:
      //   1. Find the biggest cluster.
      //   2. Find the farthest from the center point in the biggest cluster.
      //   3. Exclude the farthest point from the biggest cluster and form a new 1-point cluster.
      //
      // (follows OpenCV)

      bool reassigned = false;

      for (long i = 0; i < num_clusters; ++i)
      {
        if (num_assigned[(array_size_t)i] <= 0)
        {
          long max_cluster = (long)(std::max_element(num_assigned.begin(), num_assigned.end()) - num_assigned.begin());
          alwaysAssertM(num_assigned[(array_size_t)max_cluster] > 1, "KMeans: Maximum cluster has < 2 points");

          long farthest = -1;
          for (long j = 0; j < num_points; ++j)
            if (cluster_indices[j] == max_cluster && (farthest < 0 || cluster_sqdists[j] > cluster_sqdists[farthest]))
              farthest = j;

          alwaysAssertM(farthest >= 0, "KMeans: Farthest point in maximum cluster not found");

          // Reassign the point
          cluster_indices[farthest] = i;
          cluster_sqdists[farthest] = 0.0;

          addPointToCenter(points, farthest, i);
          subtractPointFromCenter(points, farthest, max_cluster);

          num_assigned[(array_size_t)i]++;
          num_assigned[(array_size_t)max_cluster]--;

          reassigned = true;
        }
      }

      for (long i = 0; i < num_clusters; ++i)
      {
        long n = num_assigned[(array_size_t)i];
        if (n > 0)
        {
          for (long j = 0; j < num_features; ++j)
            centers(i, j) /= n;
        }
      }

      return reassigned;
    }

    /** Add the feature vector of a point to the coordinates of a center. */
    template <typename AddressableMatrixT>
    void addPointToCenter(AddressableMatrixT const & points, long point_index, long cluster_index)
    {
      fvec.resize((array_size_t)centers.numColumns());
      points.getRow(point_index, &fvec[0]);

      for (array_size_t i = 0; i < fvec.size(); ++i)
        centers(cluster_index, (long)i) += fvec[i];
    }

    /** Subtract the feature vector of a point from the coordinates of a center. */
    template <typename AddressableMatrixT>
    void subtractPointFromCenter(AddressableMatrixT const & points, long point_index, long cluster_index)
    {
      fvec.resize((array_size_t)centers.numColumns());
      points.getRow(point_index, &fvec[0]);

      for (array_size_t i = 0; i < fvec.size(); ++i)
        centers(cluster_index, (long)i) -= fvec[i];
    }

    Matrix<double> centers;
    Options options;

    mutable TheaArray<double> fvec;
    mutable AtomicInt32 flag;

}; // class KMeans

} // namespace Algorithms
} // namespace Thea

#endif
