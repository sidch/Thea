//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2015, Siddhartha Chaudhuri
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

#ifndef __Thea_Algorithms_Histogram_hpp__
#define __Thea_Algorithms_Histogram_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "../Math.hpp"

namespace Thea {
namespace Algorithms {

/** A 1D histogram of real values. */
class THEA_API Histogram
{
  public:
    THEA_DEF_POINTER_TYPES(Histogram, shared_ptr, weak_ptr)

    /**
     * Constructor to wrap a pre-allocated array of bins.
     *
     * @param num_bins_ The number of bins in the histogram.
     * @param bins_ A preallocated array of \a num_bins_ bins.
     * @param min_value_ The lower limit of the histogram range.
     * @param max_value_ The upper limit of the histogram range.
     * @param set_zero If true, the bins are set to zero by the constructor.
     */
    Histogram(long num_bins_, double * bins_, double min_value_ = 0, double max_value_ = 1, bool set_zero = true)
    : owns_memory(false), num_bins(num_bins_), bins(bins_), sum_values(0)
    {
      alwaysAssertM(num_bins_ > 0, "Histogram: Number of bins must be positive");

      setRange(min_value_, max_value_);

      if (set_zero)
        setZero();
      else
        recomputeSum();
    }

    /**
     * Constructor to allocate an internal array for the histogram, with a specified number of bins.
     *
     * @param num_bins_ The number of bins in the histogram.
     * @param min_value_ The lower limit of the histogram range.
     * @param max_value_ The upper limit of the histogram range.
     */
    Histogram(long num_bins_, double min_value_ = 0, double max_value_ = 1)
    : owns_memory(true), num_bins(num_bins_), sum_values(0)
    {
      alwaysAssertM(num_bins_ > 0, "Histogram: Number of bins must be positive");
      bins = new double[num_bins_];

      setRange(min_value_, max_value_);
      setZero();
    }

    /**
     * Copy constructor. Always allocates a copy of the source bin array, regardless of whether it was internally allocated or
     * not.
     */
    Histogram(Histogram const & src)
    : owns_memory(false), num_bins(src.num_bins), min_value(src.min_value), max_value(src.max_value),
      bin_size(src.bin_size), sum_values(src.sum_values)
    {
      bins = new double[num_bins];
      for (long i = 0; i < num_bins; ++i)
        bins[i] = src.bins[i];
    }

    /** Destructor. */
    ~Histogram() { if (owns_memory) delete [] bins; }

    /** Get the number of bins in the histogram. */
    long numBins() const { return num_bins; }

    /** Get the array of bins in the histogram. */
    double * getBins() { return bins; }

    /** Get the array of bins in the histogram. */
    double const * getBins() const { return bins; }

    /** Get the value of a specific bin. */
    double getBin(long index) const
    {
      debugAssertM(index >= 0 && index < num_bins, format("Histogram: Index %ld out of range [%ld, %ld)", index, 0L, num_bins));
      return bins[index];
    }

    /** Get a mutable reference to the value of a specific bin. */
    double & getBin(long index)
    {
      debugAssertM(index >= 0 && index < num_bins, format("Histogram: Index %ld out of range [%ld, %ld)", index, 0L, num_bins));
      return bins[index];
    }

    /** Get the bin corresponding to a particular value, without inserting the latter. */
    long getBinForValue(double value) const
    {
      return Math::clamp((long)std::floor((value - min_value) / bin_size), 0, num_bins - 1);
    }

    /** Get the sum of values in the histogram. */
    double sumValues() const { return sum_values; }

    /** Get the lower limit of the range of the histogram. */
    double minValue() const { return min_value; }

    /** Get the upper limit of the range of the histogram. */
    double maxValue() const { return max_value; }

    /** Set the range of the histogram. */
    void setRange(double min_value_, double max_value_)
    {
      alwaysAssertM(min_value_ < max_value_, "Histogram: Range must be non-empty");

      min_value = min_value_;
      max_value = max_value_;
      bin_size = (max_value_ - min_value_) / num_bins;
    }

    /** Set all bins to zero. */
    void setZero()
    {
      for (long i = 0; i < num_bins; ++i)
        bins[i] = 0.0;

      sum_values = 0.0;
    }

    /** Add a value to the histogram and return the index of the bin in which it was inserted. */
    long insert(double value)
    {
      long bin = getBinForValue(value);
      bins[bin] += 1.0;
      sum_values += 1.0;

      return bin;
    }

    /** Remove a value from the histogram and return the index of the bin from which it was removed. */
    long remove(double value)
    {
      long bin = getBinForValue(value);
      bins[bin] -= 1.0;
      sum_values -= 1.0;

      return bin;
    }

    /** Add a set of values to the histogram. */
    template <typename IteratorT> void insert(IteratorT begin, IteratorT end)
    {
      for (IteratorT ii = begin; ii != end; ++ii)
        insert(static_cast<double>(*ii));
    }

    /** Remove a set of values from the histogram. */
    template <typename IteratorT> void remove(IteratorT begin, IteratorT end)
    {
      for (IteratorT ii = begin; ii != end; ++ii)
        remove(static_cast<double>(*ii));
    }

    /**
     * Insert the contents of the bins of another histogram into this one. The histograms must have the same number of bins and
     * (if \a match_ranges is true) the same range.
     *
     * @see remove(Histogram const &)
     */
    void insert(Histogram const & h, bool match_ranges = true)
    {
      alwaysAssertM(num_bins == h.num_bins, "Histogram: Base and increment histograms must have same bin counts");

      if (match_ranges)
      {
        alwaysAssertM(Math::fuzzyEq(min_value, h.min_value) && Math::fuzzyEq(max_value, h.max_value),
                      "Histogram: Base and increment histograms must have same ranges");
      }

      for (long i = 0; i < num_bins; ++i)
        bins[i] += h.bins[i];

      sum_values += h.sum_values;
    }

    /**
     * Subtract the contents of the bins of another histogram from this one. The histograms must have the same number of bins
     * and (if \a match_ranges is true) the same range.
     *
     * @see insert(Histogram const &)
     */
    void remove(Histogram const & h, bool match_ranges = true)
    {
      alwaysAssertM(num_bins == h.num_bins, "Histogram: Base and decrement histograms must have same bin counts");

      if (match_ranges)
      {
        alwaysAssertM(Math::fuzzyEq(min_value, h.min_value) && Math::fuzzyEq(max_value, h.max_value),
                      "Histogram: Base and decrement histograms must have same ranges");
      }

      for (long i = 0; i < num_bins; ++i)
        bins[i] -= h.bins[i];

      sum_values -= h.sum_values;
    }

    /**
     * Normalize the histogram by scaling the bin values so that they sum to 1. If they initially sum to zero, they are left
     * unchanged.
     */
    void normalize()
    {
      if (sum_values > 0)
      {
        for (long i = 0; i < num_bins; ++i)
          bins[i] /= sum_values;

        sum_values = 1.0;
      }
    }

  private:
    /** Recompute the sum of bin values. */
    void recomputeSum()
    {
      sum_values = 0.0;
      for (long i = 0; i < num_bins; ++i)
        sum_values += bins[i];
    }

    bool owns_memory;    ///< True if this object internally allocated the array of histogram bins.
    long num_bins;       ///< The number of histogram bins.
    double * bins;       ///< The array of histogram bins.
    double min_value;    ///< The upper limit of binned values.
    double max_value;    ///< The lower limit of binned values.
    double bin_size;     ///< The size of a histogram bin.
    double sum_values;   ///< The sum of values in the histogram.

}; // class Histogram

} // namespace Algorithms
} // namespace Thea

#endif
