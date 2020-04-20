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
#include "KMeans.hpp"
#include "../AbstractAddressableMatrix.hpp"
#include "../Serializable.hpp"
#include <type_traits>

namespace Thea {
namespace Algorithms {

/** A bag-of-words model for classifying objects described as a set of points, each point having a feature vector. */
class THEA_API BagOfWords : public Serializable
{
  public:
    THEA_DECL_SMART_POINTERS(BagOfWords)

    /** Options for building the model. */
    typedef KMeans::Options Options;

    /** Constructor. */
    BagOfWords(Options const & options_ = Options::defaults()) : vocabulary(options_) {}

    /** Get the current set of options. */
    Options const & getOptions() const { return vocabulary.getOptions(); }

    /** Set the current set of options. */
    void setOptions(Options const & options_) { vocabulary.setOptions(options_); }

    /** Get the number of words in the vocabulary. */
    intx numWords() const { return vocabulary.numClusters(); }

    /** Get the size of the feature vector of a point. */
    intx numPointFeatures() const { return vocabulary.numPointFeatures(); }

    /**
     * Train the Bag-of-Words model from a set of training points, one per row of the input matrix. The points are clustered
     * into a vocabulary of \a num_words words during training.
     *
     * @param num_words Number of words to cluster points into.
     * @param training_points The set of training points, one vector per row. AddressableMatrixT should have the interface of an
     *   AddressableMatrix of some real-valued scalar type.
     *
     * @return True if the training converged, else false.
     */
    template < typename AddressableMatrixT,
               typename std::enable_if< std::is_base_of< AbstractAddressableMatrix<typename AddressableMatrixT::Value>,
                                                         AddressableMatrixT >::value, int >::type = 0 >
    bool train(intx num_words, AddressableMatrixT const & training_points)
    {
      return vocabulary.cluster(num_words, training_points);
    }

    /**
     * Compute the histogram of word frequencies for a set of points. Each point is assigned to its closest word, and the
     * histogram bin corresponding to the word incremented by 1. Alternatively, a set of point weights may be specified, in
     * which case the histogram bin is incremented by the weight of the point.
     *
     * @param points The points for which to build the word histogram. AddressableMatrixT should have the interface of an
     *   AddressableMatrix of some real-valued scalar type.
     * @param histogram The vector of word frequencies, assumed to be preallocated to numWords() entries. U should be a
     *   real-valued scalar type (e.g. <code>intx</code> or <code>double</code>).
     * @param point_weights [Optional] The contribution of each point to the histogram (1 if null).
     */
    template < typename AddressableMatrixT, typename U,
               typename std::enable_if< std::is_base_of< AbstractAddressableMatrix<typename AddressableMatrixT::Value>,
                                                         AddressableMatrixT >::value, int >::type = 0 >
    void computeWordFrequencies(AddressableMatrixT const & points, U * histogram, double const * point_weights = nullptr) const
    {
      intx num_points = points.rows();
      Array<intx> labeling((size_t)num_points);
      vocabulary.mapToClusters(points, &labeling[0]);

      intx num_words = numWords();
      for (intx i = 0; i < num_words; ++i)
        histogram[i] = 0;

      for (size_t i = 0; i < labeling.size(); ++i)
      {
        if (labeling[i] >= 0)
          histogram[labeling[i]] += static_cast<U>(point_weights ? point_weights[i] : 1);
      }
    }

    /** Load the bag-of-words model from a disk file. */
    bool load(std::string const & path) { return vocabulary.load(path); }

    /** Save the bag-of-words model to a disk file. */
    bool save(std::string const & path) const { return vocabulary.save(path); }

    void read(BinaryInputStream & in, Codec const & codec = Codec_AUTO(), bool read_block_header = false)
    { vocabulary.read(in, codec, read_block_header); }

    void write(BinaryOutputStream & out, Codec const & codec = Codec_AUTO(), bool write_block_header = false) const
    { vocabulary.write(out, codec, write_block_header); }

    void read(TextInputStream & in, Codec const & codec = Codec_AUTO()) { vocabulary.read(in, codec); }
    void write(TextOutputStream & out, Codec const & codec = Codec_AUTO()) const { vocabulary.write(out, codec); }

  private:
    KMeans vocabulary;  ///< The set of learned words.

}; // class BagOfWords

} // namespace Algorithms
} // namespace Thea

#endif
