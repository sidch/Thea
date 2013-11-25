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
    deserialize(in);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "KMeans: Could not load options from input file '%s'", path.c_str())

  return true;
}

bool
KMeans::Options::save(std::string const & path) const
{
  try
  {
    TextOutputStream out(path, Serializable::configWriteSettings());
    serialize(out);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "KMeans: Could not save options to output file '%s'", path.c_str())

  return true;
}

void
KMeans::Options::deserialize(BinaryInputStream & input, Codec const & codec)
{
  max_iterations = input.readInt64();
  max_time = input.readFloat64();
  seeding = Seeding(input.readAlignedString(4));
  parallelize = (input.readInt32() != 0);
  verbose = (input.readInt32() != 0);
}

void
KMeans::Options::deserialize(TextInputStream & input, Codec const & codec)
{
  *this = defaults();

  while (input.hasMore())
  {
    std::string field = input.readSymbol();
    if (field == "END_OPTIONS")
      break;

    input.readSymbol("=");

    if (field == "max_iterations")
      max_iterations = (long)input.readNumber();
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
KMeans::Options::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  output.writeInt64(max_iterations);
  output.writeFloat64(max_time);
  output.writeAlignedString(seeding.toString(), 4);
  output.writeInt32(parallelize ? 1 : 0);
  output.writeInt32(verbose ? 1 : 0);
}

void
KMeans::Options::serialize(TextOutputStream & output, Codec const & codec) const
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
    options.deserialize(in);
    deserialize(in);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "KMeans: Could not load model from input file '%s'", path.c_str())

  return true;
}

bool
KMeans::save(std::string const & path) const
{
  try
  {
    TextOutputStream out(path.c_str(), Serializable::configWriteSettings());
    options.serialize(out);
    out.writeNewline();
    serialize(out);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "KMeans: Could not load model from input file '%s'", path.c_str())

  return true;
}

void
KMeans::deserialize(BinaryInputStream & input, Codec const & codec)
{
  long num_clusters = input.readInt64();
  if (num_clusters < 0) throw Error("KMeans: Negative number of clusters read");

  long num_point_features = input.readInt64();
  if (num_point_features < 0) throw Error("KMeans: Negative number of point features read");

  centers.resize(num_clusters, num_point_features);
  for (long i = 0; i < num_clusters; ++i)
    for (long j = 0; j < num_point_features; ++j)
      centers(i, j) = input.readFloat64();
}

void
KMeans::deserialize(TextInputStream & input, Codec const & codec)
{
  long num_clusters = (long)input.readNumber();
  if (num_clusters < 0) throw Error("KMeans: Negative number of clusters read");

  long num_point_features = (long)input.readNumber();
  if (num_point_features < 0) throw Error("KMeans: Negative number of point features read");

  centers.resize(num_clusters, num_point_features);
  for (long i = 0; i < num_clusters; ++i)
    for (long j = 0; j < num_point_features; ++j)
      centers(i, j) = input.readNumber();
}

void
KMeans::serialize(BinaryOutputStream & output, Codec const & codec) const
{
  long num_clusters = numClusters();
  long num_point_features = numPointFeatures();

  output.writeInt64(num_clusters);
  output.writeInt64(num_point_features);

  for (long i = 0; i < num_clusters; ++i)
    for (long j = 0; j < num_point_features; ++j)
      output.writeFloat64(centers(i, j));
}

void
KMeans::serialize(TextOutputStream & output, Codec const & codec) const
{
  long num_clusters = numClusters();
  long num_point_features = numPointFeatures();

  output.printf("%ld %ld\n", num_clusters, num_point_features);

  for (long i = 0; i < num_clusters; ++i)
  {
    for (long j = 0; j < num_point_features; ++j)
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
