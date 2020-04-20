#include "../BezierN.hpp"

using namespace std;
using namespace Thea;

int
main(int argc, char * argv[])
{
  BezierN<3, Real> b(3);

  b.setControl(0, Vector3(0, 0, 0));
  b.setControl(1, Vector3(1, 0, 0));
  b.setControl(2, Vector3(1, 1, 0));
  b.setControl(3, Vector3(2, 1, 0));

  THEA_CONSOLE << "b(0.0) = " << toString(b.getPoint(0.0));
  THEA_CONSOLE << "b(0.5) = " << toString(b.getPoint(0.5));
  THEA_CONSOLE << "b(1.0) = " << toString(b.getPoint(1.0));

  Array<Vector3> even_p(10);
  Array<Real> even_u(10);
  b.getEvenlySpacedPoints(10, &even_p[0], &even_u[0]);
  THEA_CONSOLE << "even = [";
  for (size_t i = 0; i < even_p.size(); ++i)
    THEA_CONSOLE << "  " << toString(even_p[i]) << ", " << even_u[i] << (i < even_p.size() - 1 ? "," : "");
  THEA_CONSOLE << "];";

  THEA_CONSOLE << "tangents = [";
  for (size_t i = 0; i < even_u.size(); ++i)
  {
    Vector3 t = b.getTangent(even_u[i]);
    THEA_CONSOLE << "  " << toString(t) << (i < even_p.size() - 1 ? "," : "");
  }
  THEA_CONSOLE << "];";

  THEA_CONSOLE << "normals = [";
  for (size_t i = 0; i < even_u.size(); ++i)
  {
    Vector3 n = b.getNormal(even_u[i]);
    THEA_CONSOLE << "  " << toString(n) << (i < even_p.size() - 1 ? "," : "");
  }
  THEA_CONSOLE << "];";

  THEA_CONSOLE << "binormals = [";
  for (size_t i = 0; i < even_u.size(); ++i)
  {
    Vector3 bn = b.getBinormal(even_u[i]);
    THEA_CONSOLE << "  " << toString(bn) << (i < even_p.size() - 1 ? "," : "");
  }
  THEA_CONSOLE << "];";

  static Vector3 const d[] = {
    Vector3(0.0, 0.0, 0.0),
    Vector3(0.0, 0.5, 0.0),
    Vector3(1.1, 1.4, 0.0),
    Vector3(2.1, 1.6, 0.0),
    Vector3(3.2, 1.1, 0.0),
    Vector3(4.0, 0.2, 0.0),
    Vector3(4.0, 0.0, 0.0),
  };
  static size_t const dsize = sizeof(d) / sizeof(Vector3);

  Array<Real> u(dsize);

  for (int reparam_iters = 0; reparam_iters <= 10; ++reparam_iters)
  {
    if (reparam_iters == 10)
      reparam_iters = 100;

    double sqerr = b.fitToPoints(&d[0], &d[0] + dsize, nullptr, &u[0], true, reparam_iters);
    THEA_CONSOLE << "sqerr" << reparam_iters << " = " << sqerr << ';';

    THEA_CONSOLE << "points" << reparam_iters << " = [";
    for (long i = 0; i < b.numControls(); ++i)
    {
      Vector3 c = b.getControl(i);
      THEA_CONSOLE << "  " << toString(c) << (i < b.numControls() - 1 ? "," : "");
    }
    THEA_CONSOLE << "];";
  }

  b.getEvenlySpacedPoints(10, &even_p[0], &even_u[0]);
  THEA_CONSOLE << "even = [";
  for (size_t i = 0; i < even_p.size(); ++i)
    THEA_CONSOLE << "  " << toString(even_p[i]) << ", " << even_u[i] << (i < even_p.size() - 1 ? "," : "");
  THEA_CONSOLE << "];";

  return 0;
}
