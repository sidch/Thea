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

#include "Codec.hpp"
#include "IOStream.hpp"

namespace Thea {

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
}

void
Codec::BlockHeader::write(BinaryOutputStream & out) const
{
  BinaryOutputStream::EndiannessScope scope(out, Endianness::LITTLE);  // headers are always little endian
  out.writeBytes(MAGIC_LENGTH, magic.data());
  out.writeUInt64(data_size);
}

int64
Codec::BlockHeader::markAndSkip(BinaryOutputStream & out)
{
  header_pos = out.getPosition();
  out.skip(BLOCK_HEADER_LENGTH);
  return header_pos;
}

void
Codec::BlockHeader::calcAndWrite(BinaryOutputStream & out)
{
  int64 final_pos = out.getPosition();
  data_size = (uint64)(final_pos - header_pos - BLOCK_HEADER_LENGTH);
  out.setPosition(header_pos);
  write(out);
  out.setPosition(final_pos);
}

} // namespace Thea
