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

#include "Codec.hpp"
#include "Iostream.hpp"

namespace Thea {

namespace CodecInternal {

intx const RESERVED_LENGTH = 8;  // additional space allocated but currently unused in a serialized BlockHeader

} // namespace CodecInternal

Codec::BlockHeader::BlockHeader(BinaryInputStream & in)
: magic(zeroMagic()), data_size(0), header_pos(0)
{
  read(in);
}

void
Codec::BlockHeader::read(BinaryInputStream & in)
{
  BinaryInputStream::EndiannessScope scope(in, Endianness::LITTLE);
  in.readBytes(MAGIC_LENGTH, magic.data());
  data_size = in.readUInt64();
  in.skip(CodecInternal::RESERVED_LENGTH);
}

void
Codec::BlockHeader::write(BinaryOutputStream & out) const
{
  BinaryOutputStream::EndiannessScope scope(out, Endianness::LITTLE);  // headers are always little endian
  out.writeBytes(MAGIC_LENGTH, magic.data());
  out.writeUInt64(data_size);

  static uint8 const RESERVED_ZEROES[CodecInternal::RESERVED_LENGTH] = {};  // zero initialization
  out.writeBytes(CodecInternal::RESERVED_LENGTH, RESERVED_ZEROES);
}

int64
Codec::BlockHeader::markAndSkip(BinaryOutputStream & out)
{
  header_pos = out.getPosition();
  out.skip(SERIALIZED_LENGTH);
  return header_pos;
}

void
Codec::BlockHeader::calcAndWrite(BinaryOutputStream & out)
{
  int64 final_pos = out.getPosition();
  data_size = (uint64)(final_pos - header_pos - SERIALIZED_LENGTH);
  out.setPosition(header_pos);
  write(out);
  out.setPosition(final_pos);
}

} // namespace Thea
