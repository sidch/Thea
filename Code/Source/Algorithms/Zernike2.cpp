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
// First version: 2011
//
//============================================================================

#include "Zernike2.hpp"
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
  lut.resize(opts.angular_steps, opts.radial_steps, lut_size + 1, lut_size + 1);  // one row/column of padding

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
            lut(p, r, x, y) = temp * std::complex<double>(std::cos(angle * p), std::sin(angle * p));
          }
      }
      else
      {
        for (int p = 0; p < opts.angular_steps; ++p)
          for (int r = 0; r < opts.radial_steps; ++r)
            lut(p, r, x, y) = std::complex<double>(0.0, 0.0);
      }
    }

  // Add some padding in case we index beyond the border
  for (int p = 0; p < opts.angular_steps; ++p)
    for (int r = 0; r < opts.radial_steps; ++r)
      for (int i = 0; i <= lut_size; ++i)
        lut(p, r, i, lut_size) = lut(p, r, lut_size, i) = std::complex<double>(0.0, 0.0);

  lut_generated = true;
}

} // namespace Algorithms
} // namespace Thea
