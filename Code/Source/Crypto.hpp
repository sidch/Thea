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

#ifndef __Thea_Crypto_hpp__
#define __Thea_Crypto_hpp__

#include "Common.hpp"

namespace Thea {

/** Cryptographic functions. */
class THEA_API Crypto
{
  public:
    /** Get the CRC32 hash of a sequence of bytes. */
    static uint32 crc32(void const * byte, size_t num_bytes);

}; // class Crypto

} // namespace Thea

#endif
