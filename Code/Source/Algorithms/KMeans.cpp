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

#include "KMeans.hpp"

namespace Thea {
namespace Algorithms {

KMeans::Options::Options()
: max_iterations(-1),
  max_time(-1),
  seeding(Seeding::K_MEANS_PLUS_PLUS),
  parallelize(true),
  verbose(true)
{
}

bool
KMeans::Options::load(std::string const & path)
{
  try
  {
    TextInputStream in(path, Serializable::configReadSettings());
    read(in);
  }
  THEA_CATCH(return false;, ERROR, "KMeans: Could not load options from input file '%s'", path.c_str())

  return true;
}

bool
KMeans::Options::save(std::string const & path) const
{
  try
  {
    TextOutputStream out(path, Serializable::configWriteSettings());
    write(out);
  }
  THEA_CATCH(return false;, ERROR, "KMeans: Could not save options to output file '%s'", path.c_str())

  return true;
}

void
KMeans::Options::read(BinaryInputStream & input, Codec const & codec, bool read_block_header)
{
  (void)read_block_header;  // ignored

  { BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);

    max_iterations = input.readInt64();
    max_time = input.readFloat64();
    seeding = Seeding(input.readAlignedString(4));
    parallelize = (input.readInt32() != 0);
    verbose = (input.readInt32() != 0);
  }
}

void
KMeans::Options::read(TextInputStream & input, Codec const & codec)
{
  *this = defaults();

  while (input.hasMore())
  {
    std::string field = input.readSymbol();
    if (field == "END_OPTIONS")
      break;

    input.readSymbol("=");

    if (field == "max_iterations")
      max_iterations = (intx)input.readNumber();
    else if (field == "max_time")
      max_time = input.readNumber();
    else if (field == "seeding")
      seeding = Seeding(input.readString());
    else if (field == "parallelize")
      parallelize = input.readBoolean();
    else if (field == "verbose")
      verbose = input.readBoolean();
  }
}

void
KMeans::Options::write(BinaryOutputStream & output, Codec const & codec, bool write_block_header) const
{
  (void)write_block_header;  // ignored

  { BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);

    output.writeInt64(max_iterations);
    output.writeFloat64(max_time);
    output.writeAlignedString(seeding.toString(), 4);
    output.writeInt32(parallelize ? 1 : 0);
    output.writeInt32(verbose ? 1 : 0);
  }
}

void
KMeans::Options::write(TextOutputStream & output, Codec const & codec) const
{
  output.printf("max_iterations = %ld\n", max_iterations);
  output.printf("max_time = %lg\n", max_time);
  output.printf("seeding = \"%s\"\n", seeding.toString().c_str());
  output.printf("parallelize = %s\n", (parallelize ? "true" : "false"));
  output.printf("verbose = %s\n", (verbose ? "true" : "false"));
  output.printf("END_OPTIONS\n");
}

bool
KMeans::load(std::string const & path)
{
  try
  {
    TextInputStream in(path.c_str(), Serializable::configReadSettings());
    read(in);
  }
  THEA_CATCH(return false;, ERROR, "KMeans: Could not load model from input file '%s'", path.c_str())

  return true;
}

bool
KMeans::save(std::string const & path) const
{
  try
  {
    TextOutputStream out(path.c_str(), Serializable::configWriteSettings());
    write(out);
  }
  THEA_CATCH(return false;, ERROR, "KMeans: Could not load model from input file '%s'", path.c_str())

  return true;
}

void
KMeans::read(BinaryInputStream & input, Codec const & codec, bool read_block_header)
{
  if (read_block_header)
    input.skip(Codec::BlockHeader::SERIALIZED_LENGTH);  // not used

  options.read(input);

  { BinaryInputStream::EndiannessScope scope(input, Endianness::LITTLE);

    intx num_clusters = input.readInt64();
    if (num_clusters < 0) throw Error("KMeans: Negative number of clusters read");

    intx num_point_features = input.readInt64();
    if (num_point_features < 0) throw Error("KMeans: Negative number of point features read");

    centers.resize(num_clusters, num_point_features);
    for (intx i = 0; i < num_clusters; ++i)
      for (intx j = 0; j < num_point_features; ++j)
        centers(i, j) = input.readFloat64();
  }
}

void
KMeans::read(TextInputStream & input, Codec const & codec)
{
  options.read(input);

  intx num_clusters = (intx)input.readNumber();
  if (num_clusters < 0) throw Error("KMeans: Negative number of clusters read");

  intx num_point_features = (intx)input.readNumber();
  if (num_point_features < 0) throw Error("KMeans: Negative number of point features read");

  centers.resize(num_clusters, num_point_features);
  for (intx i = 0; i < num_clusters; ++i)
    for (intx j = 0; j < num_point_features; ++j)
      centers(i, j) = input.readNumber();
}

void
KMeans::write(BinaryOutputStream & output, Codec const & codec, bool write_block_header) const
{
  Codec::BlockHeader header("KMeans");
  if (write_block_header)
    header.markAndSkip(output);

  options.write(output);

  { BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);

    intx num_clusters = numClusters();
    intx num_point_features = numPointFeatures();

    output.writeInt64(num_clusters);
    output.writeInt64(num_point_features);

    for (intx i = 0; i < num_clusters; ++i)
      for (intx j = 0; j < num_point_features; ++j)
        output.writeFloat64(centers(i, j));
  }

  if (write_block_header)
    header.calcAndWrite(output);
}

void
KMeans::write(TextOutputStream & output, Codec const & codec) const
{
  options.write(output);
  output.writeNewline();

  intx num_clusters = numClusters();
  intx num_point_features = numPointFeatures();

  output.printf("%ld %ld\n", num_clusters, num_point_features);

  for (intx i = 0; i < num_clusters; ++i)
  {
    for (intx j = 0; j < num_point_features; ++j)
    {
      if (j == num_point_features - 1)
        output.printf("%lg\n", centers(i, j));
      else
        output.printf("%lg ", centers(i, j));
    }
  }
}

} // namespace Algorithms
} // namespace Thea
