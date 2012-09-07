//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#include "Zernike2.hpp"

#ifndef THEA_NO_ZERNIKE

#include "../Math.hpp"

namespace Thea {
namespace Algorithms {

Zernike2::Zernike2(Options const & opts_)
: opts(opts_), lut_generated(false)
{
  alwaysAssertM(opts.angular_steps > 0 && opts.radial_steps > 0 && opts.lut_radius > 0,
                "Zernike2: Lookup dimensions must be greater than zero");
}

void
Zernike2::generateBasisLUT() const
{
  if (lut_generated) return;

  int lut_size = 2 * opts.lut_radius + 1;
  lut.resize(boost::extents[opts.angular_steps][opts.radial_steps][lut_size + 1][lut_size + 1]);  // one row/column of padding

  int max_radius = opts.lut_radius;
  double angle, temp, radius;

  for (int y = 0; y < lut_size; ++y)
    for (int x = 0; x < lut_size; ++x)
    {
      radius = std::sqrt((double)(Math::square(x - max_radius) + Math::square(y - max_radius)));

      if(radius < max_radius)
      {
        angle = std::atan2((double)(y - max_radius), (double)(x - max_radius));

        for (int p = 0; p < opts.angular_steps; ++p)
          for (int r = 0; r < opts.radial_steps; ++r)
          {
            temp = std::cos(radius * Math::pi() * r / max_radius);
            lut[p][r][x][y] = temp * std::complex<double>(std::cos(angle * p), std::sin(angle * p));
          }
      }
      else
      {
        for (int p = 0; p < opts.angular_steps; ++p)
          for (int r = 0; r < opts.radial_steps; ++r)
            lut[p][r][x][y] = std::complex<double>(0.0, 0.0);
      }
    }

  // Add some padding in case we index beyond the border
  for (int p = 0; p < opts.angular_steps; ++p)
    for (int r = 0; r < opts.radial_steps; ++r)
      for (int i = 0; i <= lut_size; ++i)
        lut[p][r][i][lut_size] = lut[p][r][lut_size][i] = std::complex<double>(0.0, 0.0);

  lut_generated = true;
}

} // namespace Algorithms
} // namespace Thea

#endif // !defined(THEA_NO_ZERNIKE)
