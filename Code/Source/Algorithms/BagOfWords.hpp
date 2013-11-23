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

#ifndef __Thea_Algorithms_BagOfWords_hpp__
#define __Thea_Algorithms_BagOfWords_hpp__

#include "../Common.hpp"
#include "../AtomicInt32.hpp"
#include "../Math.hpp"
#include "../Matrix.hpp"
#include "../Random.hpp"
#include "../Stopwatch.hpp"
#include "../System.hpp"
#include <boost/thread.hpp>

namespace Thea {
namespace Algorithms {

/** A bag-of-words model for classifying objects described as a set of points, each point having a feature vector. */
class THEA_API BagOfWords
{
  public:
    THEA_DEF_POINTER_TYPES(BagOfWords, shared_ptr, weak_ptr)

    /** Constructor. */
    BagOfWords(bool use_threads_ = true) : use_threads(use_threads_) {}

    /** Get the number of words in the vocabulary. */
    long numWords() const { return centers.numRows(); }

    /** Get the size of the feature vector of a point. */
    long numPointFeatures() const { return centers.numColumns(); }

    /**
     * Train the Bag-of-Words model from a set of training points, one per row of the input matrix. The points are clustered
     * into a vocabulary of \a num_words words during training.
     *
     * @param training_points The set of training points, one vector per row. AddressableMatrixT should have the interface of an
     *   AddressableMatrix of some real-valued scalar type.
     * @param num_words Number of words to cluster points into.
     */
    template <typename AddressableMatrixT> void train(long num_words, AddressableMatrixT const & training_points)
    {
      alwaysAssertM(num_words > 0, "BagOfWords: Number of words must be at least 1");

      long num_points = training_points.numRows();
      long num_features = training_points.numColumns();

      alwaysAssertM(num_points > 0, "BagOfWords: Number of points must be at least 1");
      alwaysAssertM(num_features > 0, "BagOfWords: Number of features must be at least 1");

      THEA_CONSOLE << "BagOfWords: Training model with " << num_words << " word(s) from " << num_points << " point(s)";

      Stopwatch timer;
      timer.tick();

      selectInitialCenters(num_words, training_points);

      timer.tock();
      THEA_CONSOLE << "BagOfWords: Selected initial center(s) in " << timer.elapsedTime() << 's';
      timer.tick();

      TheaArray<long> labeling;
      mapToCenters(num_words, training_points, &labeling);

      bool changed = true;
      long num_iterations = 0;
      double start_time = System::time();
      do
      {
        updateCenters(training_points, labeling);
        changed = mapToCenters(num_words, training_points, &labeling);
        num_iterations++;

        double curr_time = System::time();
        if (curr_time - start_time > 3)  // print an update every 3 seconds or every iteration, whichever is longer
        {
          THEA_CONSOLE << "BagOfWords: -- " << num_iterations << " iteration(s)";
          start_time = curr_time;
        }

      } while (changed);

      timer.tock();
      THEA_CONSOLE << "BagOfWords: Words found after " << num_iterations << " iteration(s) in " << timer.elapsedTime() << 's';
    }

    /**
     * Compute the histogram of word frequencies for a set of points. Each point is assigned to its closest word, and the
     * histogram bin corresponding to the word incremented by 1. Alternatively, a set of point weights may be specified, in
     * which case the histogram bin is incremented by the weight of the point.
     *
     * @param points The points for which to build the word histogram. AddressableMatrixT should have the interface of an
     *   AddressableMatrix of some real-valued scalar type.
     * @param histogram The vector of word frequencies, assumed to be preallocated to numWords() entries. U should be a
     *   real-valued scalar type (e.g. <code>long</code> or <code>double</code>).
     * @param point_weights [Optional] The contribution of each point to the histogram (1 if null).
     */
    template <typename AddressableMatrixT, typename U>
    void computeWordFrequencies(AddressableMatrixT const & points, U * histogram, double const * point_weights = NULL) const
    {
      alwaysAssertM(centers.numRows() > 0, "BagOfWords: No words found -- model has not been trained");

      long num_features = centers.numColumns();
      alwaysAssertM(num_features == points.numColumns(),
                    "BagOfWords: Number of point features of test object doesn't match training data");

      long num_words = centers.numRows();
      for (long i = 0; i < num_words; ++i)
        histogram[i] = 0;

      long num_points = points.numRows();
      long center_index = -1;
      double center_sqdist = -1;
      for (long i = 0; i < num_points; ++i)
      {
        mapToCenter(num_words, points, i, center_index, center_sqdist);

        debugAssertM(center_index >= 0 && center_index < num_words, "BagOfWords: Invalid closest center index");

        histogram[center_index] += static_cast<U>(point_weights ? point_weights[i] : 1);
      }
    }

  private:
    /** Select initial centers by k-means++. */
    template <typename AddressableMatrixT>
    void selectInitialCenters(long num_words, AddressableMatrixT const & training_points)
    {
      long num_points = training_points.numRows();
      long num_features = training_points.numColumns();

      // Choose initial cluster centers by k-means++ [Arthur/Vassilvitskii '07]:
      // 1. Choose one center uniformly at random from among the data points.
      // 2. For each data point x, compute D(x), the distance between x and the nearest center that has already been chosen.
      // 3. Choose one new data point at random as a new center, using a weighted probability distribution where a point x is
      //    chosen with probability proportional to D(x)^2.
      // 4. Repeat steps 2 and 3 until k centers have been chosen.

      centers.resize(num_words, num_features);
      centers.fill(0);

      // First center is randomly chosen
      long index = Random::common().integer(0, num_points - 1);
      addPointToCenter(training_points, index, 0);

      // Subsequent centers by k-means++
      TheaArray<double> sqdist;
      double start_time = System::time();
      for (long i = 1; i < num_words; ++i)
      {
        mapToCenters(i, training_points, NULL, &sqdist);

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

        addPointToCenter(training_points, index, i);

        double curr_time = System::time();
        if (curr_time - start_time > 3)  // print an update every 3 seconds or every iteration, whichever is longer
        {
          THEA_CONSOLE << "BagOfWords: -- selected " << i + 1 << " center(s)";
          start_time = curr_time;
        }
      }
    }

    /** Worker class for parallelizing the mapping of points to their closest centers. */
    template <typename AddressableMatrixT> class CenterMapper
    {
      public:
        /** Constructor. */
        CenterMapper(BagOfWords const * parent_, long num_centers_, AddressableMatrixT const * points_, long points_begin_,
                     long points_end_, TheaArray<long> * center_indices_, TheaArray<double> * center_sqdists_)
        : parent(parent_), num_centers(num_centers_), points(points_), points_begin(points_begin_), points_end(points_end_),
          center_indices(center_indices_), center_sqdists(center_sqdists_)
        {}

        /** Main function, called once per thread. */
        void operator()()
        {
          long index = -1;
          double sqdist = -1;
          bool changed = false;
          for (long i = points_begin; i < points_end; ++i)
          {
            parent->mapToCenter(num_centers, *points, i, index, sqdist);

            if (center_indices)
            {
              changed = changed || ((*center_indices)[(array_size_t)i] != index);
              (*center_indices)[(array_size_t)i] = index;
            }

            if (center_sqdists)
              (*center_sqdists)[(array_size_t)i] = sqdist;
          }

          if (changed)
            parent->flag.increment();
        }

      private:
        BagOfWords const * parent;
        long num_centers;
        AddressableMatrixT const * points;
        long points_begin;
        long points_end;
        TheaArray<long> * center_indices;
        TheaArray<double> * center_sqdists;

    }; // class CenterMapper

    template <typename AddressableMatrixT> friend class CenterMapper;

    /**
     * Map each point to its closest center.
     *
     * @return True if \a center_indices is non-null and the assignment to centers changed as a result of a call to this
     *   function, else false.
     */
    template <typename AddressableMatrixT>
    bool mapToCenters(long num_centers, AddressableMatrixT const & points, TheaArray<long> * center_indices,
                      TheaArray<double> * center_sqdists = NULL) const
    {
      long num_points = points.numRows();

      if (center_indices) center_indices->resize((array_size_t)num_points, -1);
      if (center_sqdists) center_sqdists->resize((array_size_t)num_points);

      flag = 0;

      unsigned int concurrency = boost::thread::hardware_concurrency();
      if (use_threads && concurrency > 1 && num_points > (long)(2 * concurrency))
      {
        boost::thread_group pool;
        double points_per_thread = num_points / (double)concurrency;

        long points_begin = 0;
        for (unsigned int i = 0; i < concurrency; ++i)
        {
          long points_end = (long)Math::round(points_begin + points_per_thread);

          pool.add_thread(new boost::thread(CenterMapper<AddressableMatrixT>(this,
                                                                             num_centers,
                                                                             &points,
                                                                             points_begin,
                                                                             points_end,
                                                                             center_indices,
                                                                             center_sqdists)));
          points_begin = points_end;
        }

        pool.join_all();
      }
      else
      {
        CenterMapper<AddressableMatrixT> mapper(this, num_centers, &points, 0, num_points, center_indices, center_sqdists);
        mapper();
      }

      return (flag.value() > 0);
    }

    /** Map a point to its nearest center. */
    template <typename AddressableMatrixT>
    void mapToCenter(long num_centers, AddressableMatrixT const & points, long point_index, long & center_index,
                     double & center_sqdist) const
    {
      long num_features = centers.numColumns();

      fvec.resize((array_size_t)num_features);
      points.getRow(point_index, &fvec[0]);

      center_index = -1;
      center_sqdist = -1;
      for (long i = 0; i < num_centers; ++i)
      {
        // Compute squared distance to this center
        double sqdist = 0;
        for (long j = 0; j < num_features; ++j)
        {
          double diff = fvec[(array_size_t)j] - centers(i, j);
          sqdist += (diff * diff);
        }

        if (i == 0 || sqdist < center_sqdist)
        {
          center_index = i;
          center_sqdist = sqdist;
        }
      }
    }

    /** Update each center to be the centroid of its cluster */
    template <typename AddressableMatrixT>
    void updateCenters(AddressableMatrixT const & points, TheaArray<long> const & center_indices)
    {
      long num_words = centers.numRows();
      long num_points = points.numRows();
      long num_features = points.numColumns();

      centers.fill(0);

      TheaArray<long> num_assigned((array_size_t)num_words, 0);
      for (long i = 0; i < num_points; ++i)
      {
        long cc_index = center_indices[(array_size_t)i];

        addPointToCenter(points, i, cc_index);
        num_assigned[(array_size_t)cc_index]++;
      }

      for (long i = 0; i < num_words; ++i)
      {
        long n = num_assigned[(array_size_t)i];
        if (n > 0)
        {
          for (long j = 0; j < num_features; ++j)
            centers(i, j) /= n;
        }
      }
    }

    /** Add the feature vector of a point to the coordinates of a center. */
    template <typename AddressableMatrixT>
    void addPointToCenter(AddressableMatrixT const & points, long point_index, long center_index)
    {
      fvec.resize((array_size_t)centers.numColumns());
      points.getRow(point_index, &fvec[0]);

      for (array_size_t i = 0; i < fvec.size(); ++i)
        centers(center_index, (long)i) += fvec[i];
    }

    Matrix<double> centers;
    bool use_threads;

    mutable TheaArray<double> fvec;
    mutable AtomicInt32 flag;

}; // class BagOfWords

} // namespace Algorithms
} // namespace Thea

#endif
