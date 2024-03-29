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
// First version: 2009
//
//============================================================================

#include "Math.hpp"
#include "Array.hpp"
#include <algorithm>

namespace Thea {
namespace Math {

intx
treeDepth(intx num_elems, int max_elems_in_leaf, Real max_child_fraction)
{
  alwaysAssertM(num_elems >= 0, "Math: Can't compute tree depth for negative number of elements");
  alwaysAssertM(max_elems_in_leaf > 0, "Math: Can't compute tree depth for non-positive number of elements at leaf");
  alwaysAssertM(max_child_fraction > 0 && max_child_fraction < 1,
                "Math: Max fraction of elements assigned to a child must be in range (0, 1)");

  if (num_elems <= 0) return 0;

  double log_inv_mcf = -std::log(max_child_fraction);
  intx est_depth = (intx)std::ceil(std::log(num_elems / (double)max_elems_in_leaf) / log_inv_mcf);

  return est_depth < 0 ? 0 : est_depth;  // check shouldn't be necessary but do it just in case
}

namespace MathInternal {

inline double
cbrt(double x)
{
  return x < 0 ? -std::pow(-x, 1.0 / 3.0) : std::pow(x, 1.0 / 3.0);
}

double const SOLVER_EPSILON = 1e-30;

} // namespace MathInternal

// Root of linear equation c0 + c1 * x = 0
int
solveLinear(double c0, double c1, double * root)
{
  using namespace MathInternal;

  if (std::fabs(c1) < SOLVER_EPSILON)
    return 0;

  *root = -c0 / c1;
  return 1;
}

// Distinct real roots x1, x2 of quadratic equation c0 + c1 * x + c2 * x^2 = 0
int
solveQuadratic(double c0, double c1, double c2, double * roots)
{
  using namespace MathInternal;

  if (std::fabs(c2) < SOLVER_EPSILON)
    return solveLinear(c1, c0, roots);
  else
  {
    double d = c1 * c1 - 4 * c2 * c0;

    if (std::fabs(d) < SOLVER_EPSILON)
    {
      *roots = -c1 / (2 * c2);
      return 1;
    }
    else if (d > 0)
    {
      double sqrt_d = std::sqrt(d);
      roots[0] = (-c1 + sqrt_d) / (2 * c2);
      roots[1] = (-c1 - sqrt_d) / (2 * c2);
      return 2;
    }
  }

  return 0;
}

// Real roots of cubic equation c0 + c1 * x + c2 * x^2 + c3 * x^3 = 0
int
solveCubic(double c0, double c1, double c2, double c3, double * roots)
{
  using namespace MathInternal;

  static double const RAD_120 = degreesToRadians(120.0);
  static double const RAD_240 = degreesToRadians(240.0);

  int num_roots = 0;

  if (std::fabs(c3) < SOLVER_EPSILON)
    num_roots = solveQuadratic(c0, c1, c2, roots);
  else if (std::fabs(c2) < SOLVER_EPSILON && std::fabs(c1) < SOLVER_EPSILON)
  {
    *roots = std::pow(-c0 / c3, 1.0 / 3.0);
    num_roots = 1;
  }
  else if (std::fabs(c0) < SOLVER_EPSILON)
  {
    *roots = 0;
    num_roots = 1;

    double quad_roots[2];
    int num_quad_roots = solveQuadratic(c1, c2, c3, quad_roots);
    for (int i = 0; i < num_quad_roots; ++i)
      if (std::fabs(quad_roots[i]) > SOLVER_EPSILON)
        roots[num_roots++] = quad_roots[i];
  }
  else
  {
    double k = 0, p = 0, q = 0;
    if (std::fabs(c2) > SOLVER_EPSILON)
    {
      k = (-c2) / (3 * c3);
      p = (3 * c3 * c1 - c2 * c2) / (-3 * c3 * c3);
      q = (2 * c2 * c2 * c2 - 9 * c3 * c2 * c1 + 27 * c3 * c3 * c0) / (-27 * c3 * c3 * c3);
    }
    else
    {
      k = 0;
      p = -c1 / c3;
      q = -c0 / c3;
    }

    double p_d3_e3 = p * p * p / 27;
    double w = q * q / 4 - p_d3_e3;

    if (w < 0)
    {
      double cos3a = q / (2 * std::sqrt(p_d3_e3));
      double a = std::acos(cos3a) / 3;
      double t = std::sqrt(p * (4.0 / 3.0));

      roots[0] = t * std::cos(a) + k;
      roots[1] = t * std::cos(a + RAD_120) + k;
      roots[2] = t * std::cos(a + RAD_240) + k;
      num_roots = 3;
    }
    else
    {
      if (w < SOLVER_EPSILON)
        *roots = 2 * MathInternal::cbrt(q / 2) + k;
      else
        *roots = MathInternal::cbrt(q / 2 + std::sqrt(w)) + MathInternal::cbrt(q / 2 - std::sqrt(w)) + k;

      num_roots = 1;
    }
  }

  return num_roots;
}

// Real roots of quartic equation c0 + c1 * x + c2 * x^2 + c3 * x^3 + c4 * x^4 = 0
int
solveQuartic(double c0, double c1, double c2, double c3, double c4, double * roots)
{
  using namespace MathInternal;

  int num_roots = 0;

  if (std::fabs(c4) < SOLVER_EPSILON)
    num_roots = solveCubic(c0, c1, c2, c3, roots);
  else if (std::fabs(c3) < SOLVER_EPSILON && std::fabs(c2) < SOLVER_EPSILON && std::fabs(c1) < SOLVER_EPSILON)
  {
    double x = -c0 / c4;
    if (std::fabs(x) < SOLVER_EPSILON)
    {
      *roots = 0;
      num_roots = 1;
    }
    else if (x > 0)
    {
      roots[0] = std::sqrt(std::sqrt(x));
      roots[1] = -roots[0];
      num_roots = 2;
    }
    else
      num_roots = 0;
  }
  else if (std::fabs(c0) < SOLVER_EPSILON)
  {
    *roots = 0;
    num_roots = 1;

    double cubic_roots[3];
    int num_cubic_roots = solveCubic(c1, c2, c3, c4, cubic_roots);
    for (int i = 0; i < num_cubic_roots; ++i)
      if (std::fabs(cubic_roots[i]) > SOLVER_EPSILON)
        roots[num_roots++] = cubic_roots[i];
  }
  else
  {
    double a3 = c3 / c4;
    double a2 = c2 / c4;
    double a1 = c1 / c4;
    double a0 = c0 / c4;

    double cubic_roots[3];
    int num_cubic_roots = solveCubic(a1 * a1 - a0 * (4 * a2 - a3 * a3), 8 * a0 - 2 * a1 * a3, 4 * a2, -8, cubic_roots);
    if (num_cubic_roots <= 0)
      num_roots = 0;
    else
    {
      double hi_root = cubic_roots[0];
      if (num_cubic_roots > 1) hi_root = std::max(hi_root, cubic_roots[1]);
      if (num_cubic_roots > 2) hi_root = std::max(hi_root, cubic_roots[2]);

      num_roots = 0;

      double h2 = hi_root * hi_root;
      if (h2 > a0)
      {
        double b = std::sqrt(h2 - a0);
        double m = (a3 * hi_root - a1) / (2 * b);

        double quad_roots[2][2] = { {0, 0}, {0, 0} };
        int nqr[2] = {
          solveQuadratic(hi_root + b, (a3 / 2) + m, 1, (double *)quad_roots[0]),
          solveQuadratic(hi_root - b, (a3 / 2) - m, 1, (double *)quad_roots[1])
        };

        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < nqr[i]; ++j)
            roots[num_roots++] = quad_roots[i][j];
      }
    }
  }

  return num_roots;
}

} // namespace Math
} // namespace Thea
