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
// First version: 2013
//
//============================================================================

#ifndef __Thea_Algorithms_BagOfWords_hpp__
#define __Thea_Algorithms_BagOfWords_hpp__

#include "../Common.hpp"
#include "KMeans.hpp"
#include "../IAddressableMatrix.hpp"
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
               typename std::enable_if< std::is_base_of< IAddressableMatrix<typename AddressableMatrixT::Value>,
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
               typename std::enable_if< std::is_base_of< IAddressableMatrix<typename AddressableMatrixT::Value>,
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

    void read(BinaryInputStream & in, Codec const & codec = CodecAuto(), bool read_block_header = false)
    { vocabulary.read(in, codec, read_block_header); }

    void write(BinaryOutputStream & out, Codec const & codec = CodecAuto(), bool write_block_header = false) const
    { vocabulary.write(out, codec, write_block_header); }

    void read(TextInputStream & in, Codec const & codec = CodecAuto()) { vocabulary.read(in, codec); }
    void write(TextOutputStream & out, Codec const & codec = CodecAuto()) const { vocabulary.write(out, codec); }

  private:
    KMeans vocabulary;  ///< The set of learned words.

}; // class BagOfWords

} // namespace Algorithms
} // namespace Thea

#endif
