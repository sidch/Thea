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
// First version: 2020
//
//============================================================================

#ifndef __Thea_MatrixIO_hpp__
#define __Thea_MatrixIO_hpp__

#include "BinaryInputStream.hpp"
#include "BinaryOutputStream.hpp"
#include "FilePath.hpp"
#include "MatVec.hpp"
#include "SparseMatVec.hpp"
#include <algorithm>
#include <cctype>
#include <type_traits>

namespace Thea {

/** Abstract base class for codecs that read/write matrices to/from streams. */
template <typename MatrixT>
class MatrixCodec : public Codec
{
  public:
    /**
     * Given an input stream and an optional (already extracted) block header, check if the codec is suitable for decoding the
     * stream. Among other things, this function can check the magic string in the block header (if present), whether the stream
     * name has an appropriate extension (if no block header is present), and whether the serialized and target matrix types are
     * similarly sparse or dense.
     */
    virtual bool matches(BinaryInputStream & input, Codec::BlockHeader const * block_header) const = 0;

    /**
     * Read the matrix from an input stream. Any block header, if present, is assumed to have already been extracted from the
     * stream and stored in the object pointed to by \a block_header. Else, \a block_header should be a null pointer.
     */
    virtual void readMatrix(MatrixT & m, BinaryInputStream & input, BlockHeader const * block_header) const = 0;

    /** Write the matrix to an output stream, optionally preceded by a block header. */
    virtual void writeMatrix(MatrixT const & m, BinaryOutputStream & output, bool write_block_header) const = 0;

}; // class MatrixCodec

namespace MatrixIOInternal {

// From the block header, return whether the matrix stored in a stream is sparse or dense.
inline bool
isSparse(Codec::BlockHeader const & block_header)
{
  return (bool)(block_header.custom & 0x00FF);  // all standard codecs use least-significant byte to store sparsity flag
}

} // namespace MatrixIOInternal

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

// CSV parameters that are stored in the custom field of the block header.
struct Params
{
  bool is_sparse;
  char separator;
};

// If a block header is present, use it to create a new input stream wrapping the matrix data block, and also get the parameters
// stored in the custom field. Else, leave the parameters unchanged, and return a null pointer.
inline BinaryInputStream::Ptr
parseHeader(BinaryInputStream & input, Codec::BlockHeader const * block_header, Params & params)
{
  if (block_header)
  {
    params.is_sparse = MatrixIOInternal::isSparse(*block_header);
    params.separator = (char)(block_header->custom & 0xFF00);
    return std::make_shared<BinaryInputStream>(input, (int64)block_header->data_size);
  }

  return BinaryInputStream::Ptr();
}

inline char
detectSeparatorFromLine(std::string const & line, char default_sep)
{
  // Look for separators in the order comma, tab, space
  if (line.find_first_of(',')  != std::string::npos) { return ',';  }
  if (line.find_first_of('\t') != std::string::npos) { return '\t'; }
  if (line.find_first_of(' ')  != std::string::npos) { return ' ';  }
  return default_sep;
}

// Check from the content of a CSV stream (not the block header) whether it is likely to represent a sparse matrix or a dense
// matrix
inline bool
isSparse(BinaryInputStream & input, char sep)
{
  std::string line1, line2;
  auto curr_pos = input.getPosition();
    if (input.hasMore()) line1 = input.readLine();
    if (input.hasMore()) line2 = input.readLine();
  input.setPosition(curr_pos);

  if (line1.empty() || line2.empty())
    return false;  // can't really guess, so assume it is dense

  if (sep < 0) sep = detectSeparatorFromLine(line1, ',');
  if (sep < 0) return false;  // couldn't detect separator, again default to dense

  auto n1 = std::count(line1.begin(), line1.end(), sep);
  auto n2 = std::count(line2.begin(), line2.end(), sep);

  return (n1 == 1 && n2 == 2);  // a non-empty sparse matrix has "#rows,#cols" on the first line, and a triplet on the second
}

} // namespace CodecCSVInternal

// Specialization of CodecCSV for dense matrices.
template <typename MatrixT>
class CodecCSV< MatrixT, typename std::enable_if< std::is_base_of< Eigen::DenseBase<MatrixT>, MatrixT >::value >::type >
: public MatrixCodec<MatrixT>
{
  public:
    /**
     * Constructor.
     *
     * @param separator_ Field separation character, typically comma, tab or space. If negative, it will be auto-detected
     *   (independently) from each input stream.
     */
    CodecCSV(char separator_ = -1) : separator(separator_) {}

    char const * getName() const { static char const * my_name = "CSV"; return my_name; }
    Codec::MagicString const & getMagic() const { static Codec::MagicString const magic = Codec::toMagic("CSV"); return magic; }

    bool matches(BinaryInputStream & input, Codec::BlockHeader const * block_header) const
    {
      return block_header ? (block_header->magic == getMagic() && !MatrixIOInternal::isSparse(*block_header))
                          : (toLower(FilePath::extension(input.getName())) == "csv"
                          && !CodecCSVInternal::isSparse(input, separator));
    }

    void readMatrix(MatrixT & m, BinaryInputStream & input, Codec::BlockHeader const * block_header) const
    {
      BinaryInputStream * in = &input;
      CodecCSVInternal::Params params{ /* is_sparse = */ false, /* separator = */ separator };
      auto tmp_in = CodecCSVInternal::parseHeader(input, block_header, params);
      if (tmp_in) { in = tmp_in.get(); }

      if (params.is_sparse)
        throw Error(format("%s: Matrix in stream '%s' is sparse, but a dense target matrix was given",
                           getName(), input.getName()));

      Array<std::string> lines;
      while (in->hasMore())
      {
        lines.push_back(in->readLine());

        // If the separator is still not set, try to autodetect it from the first line
        if (lines.size() == 1 && params.separator < 0)
        {
          params.separator = CodecCSVInternal::detectSeparatorFromLine(lines[0], '\n');  // assume entire line is one field
          if (params.separator < 0)
            throw Error(format("%s: Could not autodetect the field separator from the stream '%s'",
                               getName(), input.getName()));
        }
      }

      if (lines.empty())
      {
        m.resize(0, 0);
        return;
      }

      intx num_rows = (intx)lines.size();
      intx num_cols = (intx)std::count(lines[0].begin(), lines[0].end(), params.separator) + 1;
      m.resize(num_rows, num_cols);

      std::string field;
      for (intx i = 0; i < num_rows; ++i)
      {
        std::istringstream line_in(lines[(size_t)i]);
        for (intx j = 0; j < num_cols; ++j)
        {
          if (!std::getline(line_in, field, params.separator))
            throw Error(format("%s: Could not read matrix element (%ld, %ld) from stream '%s'",
                               getName(), (long)i, (long)j, input.getName()));

          m(i, j) = fromString<typename MatrixT::value_type>(field);
        }
      }
    }

    void writeMatrix(MatrixT const & m, BinaryOutputStream & output, bool write_block_header) const
    {
      Codec::BlockHeader bh(this->getMagic());
      if (write_block_header)
        bh.markAndSkip(output);

      // No need to set endianness, CSV is a purely text-based format

      char sep = (separator < 0 ? ',' : separator);
      for (intx i = 0; i < m.rows(); ++i)
      {
        std::ostringstream row;
        for (intx j = 0; j < m.cols(); ++j)
        {
          if (j > 0) row << sep;
          row << m(i, j);
        }
        row << '\n';
        output.writeBytes((int64)row.str().length(), row.str().data());
      }

      if (write_block_header)
        bh.calcAndWrite(output);
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
    /**
     * Constructor.
     *
     * @param separator_ Field separation character, typically comma, tab or space. If negative, it will be auto-detected
     *   (independently) from each input stream.
     */
    CodecCSV(char separator_ = -1) : separator(separator_) {}

    char const * getName() const { static char const * my_name = "CSV"; return my_name; }
    Codec::MagicString const & getMagic() const { static Codec::MagicString const magic = Codec::toMagic("CSV"); return magic; }

    bool matches(BinaryInputStream & input, Codec::BlockHeader const * block_header) const
    {
      return block_header ? (block_header->magic == getMagic() && MatrixIOInternal::isSparse(*block_header))
                          : (toLower(FilePath::extension(input.getName())) == "csv"
                          && CodecCSVInternal::isSparse(input, separator));
    }

    void readMatrix(MatrixT & m, BinaryInputStream & input, Codec::BlockHeader const * block_header) const
    {
      BinaryInputStream * in = &input;
      CodecCSVInternal::Params params{ /* is_sparse = */ true, /* separator = */ separator };
      auto tmp_in = CodecCSVInternal::parseHeader(input, block_header, params);
      if (tmp_in) { in = tmp_in.get(); }

      if (!params.is_sparse)
        throw Error(format("%s: Matrix in stream '%s' is dense, but a sparse target matrix was given",
                           getName(), input.getName()));

      if (!in->hasMore())
      {
        m.resize(0, 0);
        m.data().squeeze();
        return;
      }

      // Read dimensions of matrix
      std::string rstr, cstr;
      {
        auto line = in->readLine();

        // If the separator is still not set, try to autodetect it from the first line
        if (params.separator < 0)
        {
          params.separator = CodecCSVInternal::detectSeparatorFromLine(line, -1);
          if (params.separator < 0)
            throw Error(format("%s: Could not autodetect the field separator from the stream '%s'",
                               getName(), input.getName()));
        }

        std::istringstream line_in(line);
        if (!std::getline(line_in, rstr, params.separator) || !std::getline(line_in, cstr, params.separator))
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
          if (!std::getline(line_in, rstr, params.separator)
           || !std::getline(line_in, cstr, params.separator)
           || !std::getline(line_in, vstr, params.separator))
            throw Error(format("%s: Could not read triplet on line %ld of stream '%s'",
                               getName(), (long)triplets.size() + 2, input.getName()));

          intx row = (intx)std::stol(rstr), col = (intx)std::stol(cstr);
          auto value = fromString<typename MatrixT::value_type>(vstr);
          if (row < 0 || row >= m.rows() || col < 0 || col >= m.cols())
            throw Error(format("%s: Out-of-bounds matrix position (%ld, %ld) in stream '%s'",
                               getName(), (long)row, (long)col, input.getName()));

          triplets.push_back(Eigen::Triplet<Value>(row, col, value));
        }

        m.setFromTriplets(triplets.begin(), triplets.end());
      }
    }

    void writeMatrix(MatrixT const & m, BinaryOutputStream & output, bool write_block_header) const
    {
      Codec::BlockHeader bh(this->getMagic());
      if (write_block_header)
        bh.markAndSkip(output);

      // No need to set endianness, CSV is a purely text-based format

      char sep = (separator < 0 ? ',' : separator);
      output.printf("%ld%c%ld\n", (long)m.rows(), sep, (long)m.cols());

      for (intx i = 0; i < m.outerSize(); ++i)
      {
        std::ostringstream oss;
        for (typename MatrixT::InnerIterator it(m, i); it; ++it)
        {
          oss << it.row() << sep << it.col() << sep << it.value() << '\n';
          output.writeBytes((int64)oss.str().length(), oss.str().data());
        }
      }

      if (write_block_header)
        bh.calcAndWrite(output);
    }

  private:
    char separator;  ///< Field separator character.

}; // class CodecCSV

template <typename MatrixT>
intx
BinaryInputStream::readMatrix(bool read_block_header, MatrixT & m, Codec const & codec)
{
  Codec::BlockHeader bh_obj, * bh = nullptr;
  if (read_block_header)
  {
    bh_obj.read(*this);
    bh = &bh_obj;
  }

  if (codec == CodecAuto())
  {
    // Default codecs
    static CodecCSV<MatrixT> const codec_csv;
    static MatrixCodec<MatrixT> const * default_codecs[] = { &codec_csv };

    for (size_t i = 0; i < sizeof(default_codecs) / sizeof(MatrixCodec<MatrixT> const *); ++i)
      if (default_codecs[i]->matches(*this, bh))
      {
        default_codecs[i]->readMatrix(m, *this, bh);
        return 0;  // always
      }

    throw Error(getNameStr() + ": Could not autodetect a suitable codec to read the matrix");
  }
  else
  {
    MatrixCodec<MatrixT> const * mat_codec = dynamic_cast< MatrixCodec<MatrixT> const * >(&codec);
    if (!mat_codec)
      throw Error(getNameStr() + ": Codec is not a matrix codec for this matrix type");

    if (bh && bh->magic != mat_codec->getMagic())
      throw Error(getNameStr() + ": Magic string mismatch");

    mat_codec->readMatrix(m, *this, bh);
  }

  return 0;  // always
}

inline intx
BinaryInputStream::readMatrixHelper(intx index, bool read_block_header)
{
  throw Error(getNameStr() + ": Could not autodetect a suitable codec and matrix type to read the input");
}

template <typename CandidateMatrixT, typename... MoreMatrixTypes>
intx
BinaryInputStream::readMatrixHelper(intx index, bool read_block_header, CandidateMatrixT & m, MoreMatrixTypes & ... rest)
{
  // Try reading the first candidate matrix type
  auto init_pos = getPosition();
  try
  {
    readMatrix(read_block_header, m);
    return index;
  }
  catch (...)
  {
    setPosition(init_pos);  // read failed, rewind silently
  }

  return readMatrixHelper(index + 1, read_block_header, rest...);  // recurse for the other types
}

template <typename MatrixT>
void
BinaryOutputStream::writeMatrix(MatrixT const & m, Codec const & codec, bool write_block_header)
{
  if (codec == CodecAuto())
    throw Error(getNameStr() + ": You must explicitly choose a codec for writing matrices");

  MatrixCodec<MatrixT> const * mat_codec = dynamic_cast< MatrixCodec<MatrixT> const * >(&codec);
  if (!mat_codec)
    throw Error(getNameStr() + ": Codec is not a matrix codec for this matrix type");

  mat_codec->writeMatrix(m, *this, write_block_header);
}

} // namespace Thea

#endif
