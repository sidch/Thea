//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2020, Siddhartha Chaudhuri
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

#ifndef __Thea_MatrixIO_hpp__
#define __Thea_MatrixIO_hpp__

#include "IOStream.hpp"
#include "MatVec.hpp"
#include "SparseMatVec.hpp"
#include <algorithm>
#include <type_traits>

namespace Thea {

/** Base class for codecs that read/write matrices to/from streams. */
template <typename MatrixT>
class MatrixCodec : public Codec
{
  public:
    /** Read the matrix from the input stream. */
    virtual void readMatrix(MatrixT & m, BinaryInputStream & input, BlockHeader const * block_header) const = 0;

}; // class MatrixCodec

/**
 * Codec for reading/writing matrices to/from comma-separated values (CSV) text files. Instead of the comma, a different
 * separating character may be specified (typically tab or space).
 *
 * For dense matrices, each line is a row of the matrix, with elements separated by the separating character.
 *
 * For sparse matrices, the format (assuming the field separator is ',') is a sequence of comma-separated (row, column, value)
 * triplets, preceded by a line containing the dimensions of the matrix:
 *
 * \code
 * #rows,#cols
 * R0,C1,V0
 * R1,C1,V1
 * ...
 * RN, CN, VN
 * \endcode
 *
 * (This is not a standardized format, just a simple custom adaptation of CSV for sparse matrices.)
 *
 * If a block header is specified for readMatrix(), it is used to set the separating character (deserialized from the custom
 * field of the header) for this read, and the character specified when constructing the codec is ignored.
 */
template <typename MatrixT, typename Enable = void>
class CodecCSV : public MatrixCodec<MatrixT> {};

namespace CodecCSVInternal {

// If a block header is present, use it to create a new input stream wrapping the matrix data block, and also get the field
// separation character.
inline BinaryInputStream::Ptr
parseHeader(BinaryInputStream & input, Codec::BlockHeader const * block_header, char & sep)
{
  if (block_header)
  {
    sep = (char)(block_header->custom & 0xFF);
    return std::make_shared<BinaryInputStream>(input, (int64)block_header->data_size);
  }

  return BinaryInputStream::Ptr();
}

} // namespace CodecCSVInternal

// Specialization of CodecCSV for dense matrices.
template <typename MatrixT>
class CodecCSV< MatrixT, typename std::enable_if< std::is_base_of< Eigen::DenseBase<MatrixT>, MatrixT >::value >::type >
: public MatrixCodec<MatrixT>
{
  public:
    /** Construct the codec with a given field separation character (typically comma, tab or space). */
    CodecCSV(char separator_ = ',') : separator(separator_) {}

    char const * getName() const { static char const * my_name = "CSV"; return my_name; }
    Codec::MagicString const & getMagic() const { static Codec::MagicString const magic = Codec::toMagic("CSV"); return magic; }

    void readMatrix(MatrixT & m, BinaryInputStream & input, Codec::BlockHeader const * block_header) const
    {
      BinaryInputStream * in = &input;
      char sep = separator;
      auto tmp_in = CodecCSVInternal::parseHeader(input, block_header, sep);
      if (tmp_in) { in = tmp_in.get(); }

      Array<std::string> lines;
      while (in->hasMore())
        lines.push_back(in->readLine());

      if (lines.empty())
      {
        m.resize(0, 0);
        return;
      }

      intx num_rows = (intx)lines.size();
      intx num_cols = (intx)std::count(lines[0].begin(), lines[0].end(), sep) + 1;
      m.resize(num_rows, num_cols);

      std::string field;
      for (intx i = 0; i < num_rows; ++i)
      {
        std::istringstream line_in(lines[(size_t)i]);
        for (intx j = 0; j < num_cols; ++j)
        {
          if (!std::getline(line_in, field, sep))
            throw Error(format("%s: Could not read number from row %ld, column %ld of stream '%s'",
                               getName(), (long)i, (long)j, input.getName()));

          m(i, j) = std::stod(field);
        }
      }
    }

  private:
    char separator;  ///< Field separator character.

}; // class CodecCSV< Eigen::DenseBase<MatrixT> >

// Specialization of CodecCSV for sparse matrices.
template <typename MatrixT>
class CodecCSV< MatrixT, typename std::enable_if< std::is_base_of< Eigen::SparseMatrixBase<MatrixT>, MatrixT >::value >::type >
: public MatrixCodec<MatrixT>
{
  public:
    /** Construct the codec with a given field separation character (typically comma, tab or space). */
    CodecCSV(char separator_ = ',') : separator(separator_) {}

    char const * getName() const { static char const * my_name = "CSV"; return my_name; }
    Codec::MagicString const & getMagic() const { static Codec::MagicString const magic = Codec::toMagic("CSV"); return magic; }

    void readMatrix(MatrixT & m, BinaryInputStream & input, Codec::BlockHeader const * block_header) const
    {
      BinaryInputStream * in = &input;
      char sep = separator;
      auto tmp_in = CodecCSVInternal::parseHeader(input, block_header, sep);
      if (tmp_in) { in = tmp_in.get(); }

      if (!in->hasMore())
      {
        m.resize(0, 0);
        m.data().squeeze();
        return;
      }

      // Read dimensions of matrix
      std::string rstr, cstr;
      {
        std::istringstream line_in(in->readLine());
        if (!std::getline(line_in, rstr, sep) || !std::getline(line_in, cstr, sep))
          throw Error(format("%s: Could not read matrix dimensions from stream '%s'", getName(), input.getName()));

        intx nrows = (intx)std::stol(rstr), ncols = (intx)std::stol(cstr);
        if (nrows < 0 || ncols < 0)
          throw Error(format("%s: Negative matrix dimension in stream '%s'", getName(), input.getName()));

        m.resize(nrows, ncols);
        m.data().squeeze();
      }

      // Read triplets and initialize matrix
      {
        typedef typename MatrixT::value_type Value;
        Array< Eigen::Triplet<Value> > triplets;
        std::string vstr;
        while (in->hasMore())
        {
          std::istringstream line_in(in->readLine());
          if (!std::getline(line_in, rstr, sep) || !std::getline(line_in, cstr, sep) || !std::getline(line_in, vstr, sep))
            throw Error(format("%s: Could not read triplet on line %ld of stream '%s'",
                               getName(), (long)triplets.size() + 2, input.getName()));

          intx row = (intx)std::stol(rstr), col = (intx)std::stol(cstr);
          long double value = std::stold(vstr);
          if (row < 0 || row >= m.rows() || col < 0 || col >= m.cols())
            throw Error(format("%s: Out-of-bounds matrix position (%ld, %ld) in stream '%s'",
                               getName(), (long)row, (long)col, input.getName()));

          triplets.push_back(Eigen::Triplet<Value>(row, col, static_cast<Value>(value)));
        }

        m.setFromTriplets(triplets.begin(), triplets.end());
      }
    }

  private:
    char separator;  ///< Field separator character.

}; // class CodecCSV

template <typename MatrixT>
void
BinaryInputStream::readMatrix(MatrixT & m, Codec const & codec, bool read_block_header)
{
  Codec::BlockHeader bh;
  if (read_block_header)
    bh.read(*this);

  if (codec == Codec_AUTO())
  {
    // Default codecs
    static CodecCSV<MatrixT> const codec_csv;
    static MatrixCodec<MatrixT> const * default_codecs[] = { &codec_csv };

    if (!read_block_header)
      throw Error(getNameStr() + ": Cannot autodetect codec without a block header");

    for (size_t i = 0; i < sizeof(default_codecs) / sizeof(MatrixCodec<MatrixT> const *); ++i)
      if (bh.magic == default_codecs[i]->getMagic())
      {
        default_codecs[i]->readMatrix(m, *this, &bh);
        return;
      }

    throw Error(getNameStr() + ": Could not auto-detect a suitable codec to read the matrix");
  }
  else
  {
    MatrixCodec<MatrixT> const * mat_codec = dynamic_cast< MatrixCodec<MatrixT> const * >(&codec);
    if (!mat_codec)
      throw Error(getNameStr() + ": Codec is not a matrix codec");

    if (read_block_header && bh.magic != mat_codec->getMagic())
      throw Error(getNameStr() + ": Magic string mismatch");

    mat_codec->readMatrix(m, *this, (read_block_header ? &bh : NULL));
  }
}

} // namespace Thea

#endif
