#include "../Common.hpp"
#include "../Algorithms/MetricL2.hpp"
#include "../AffineTransform2.hpp"
#include "../AffineTransform3.hpp"
#include "../Math.hpp"
#include "../RigidTransform2.hpp"
#include "../RigidTransform3.hpp"
#include "../Transformable.hpp"
#include <iostream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Algorithms;

void testMetrics(int argc, char * argv[]);

int
main(int argc, char * argv[])
{
  try
  {
    testMetrics(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  // Hooray, all tests passed
  cout << "Metrics: Test completed" << endl;
  return 0;
}

void
testMetrics(int argc, char * argv[])
{
  cout << MetricL2::distance<2, Real>(Vector2::zero(), Vector2(1, 1)) << endl;
  cout << MetricL2::distance<3, Real>(Vector3::zero(), Vector3(1, 1, 1)) << endl;
  cout << MetricL2::distance<4, Real>(Vector4::zero(), Vector4(1, 1, 1, 1)) << endl;

  AxisAlignedBox3 aabb(Vector3::zero(), Vector3(1, 1, 1));
  Vector3 p(2, 2, 2);

  cout << MetricL2::distance<3, Real>(aabb, p) << endl;
  cout << MetricL2::distance<3, Real>(p, aabb) << endl;
  cout << MetricL2::distance<3, Real>(aabb, aabb) << endl;

  Ball3 ball(Vector3(-2, -2, -2), 1);
  cout << MetricL2::distance<3, Real>(ball, p) << endl;
  cout << MetricL2::distance<3, Real>(p, ball) << endl;
  cout << MetricL2::distance<3, Real>(aabb, ball) << endl;
  cout << MetricL2::distance<3, Real>(ball, aabb) << endl;

  LocalTriangle3 tri(Vector3(-1, -1, -1), Vector3(0.5, 0.5, 0), Vector3(0, 0.5, 0.5));
  cout << MetricL2::distance<3, Real>(tri, p) << endl;
  cout << MetricL2::distance<3, Real>(p, tri) << endl;
  cout << MetricL2::distance<3, Real>(tri, ball) << endl;
  cout << MetricL2::distance<3, Real>(ball, tri) << endl;

  Matrix4 m1 = Matrix4(Matrix3::rotationEulerAnglesXYZ(Math::degreesToRadians(30),
                                                       Math::degreesToRadians(30),
                                                       Math::degreesToRadians(30)), Vector3::zero());
  Matrix4 m2 = Matrix4::homTranslation(Vector3(1, 1, 1));
  RigidTransform3 rt = RigidTransform3::translation(Vector3(5, 5, 5))
                     * RigidTransform3::rotationAxisAngle(Vector3(-1, 1, -1), Math::degreesToRadians(45));
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&p, &m1),     makeTransformedObject(&p, &m2)) << endl;
  cout << MetricL2::distance<3, Real>(p,                                  makeTransformedObject(&p, &m2)) << endl;
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&p, &m1),     p) << endl;
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&p, &m1),     makeTransformedObject(&ball, &rt)) << endl;
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&ball, &rt),  makeTransformedObject(&p, &m2)) << endl;
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&p, &m1),     makeTransformedObject(&tri, &m2)) << endl;
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&tri, &m1),   makeTransformedObject(&p, &m2)) << endl;
  cout << MetricL2::distance<3, Real>(makeTransformedObject(&tri, &m1),   makeTransformedObject(&tri, &m2)) << endl;

#if 0
  cout << MetricL2::distance<VectorN<10, ot>, VectorN<10, T> >
  cout << MetricL2::distance< AxisAlignedBoxN<N, T>, VectorN<N, T> >
  cout << MetricL2::distance< VectorN<N, T>, AxisAlignedBoxN<N, T> >
  cout << MetricL2::distance< AxisAlignedBoxN<N, T>, AxisAlignedBoxN<N, T> >
  cout << MetricL2::distance< Triangle3<VertexTripleType>, Vector3 >
  cout << MetricL2::distance< Vector3, Triangle3<VertexTripleType> >
  cout << MetricL2::distance< Triangle3<VertexTripleType1>, Triangle3<VertexTripleType2> >
  cout << MetricL2::distance< Triangle3<VertexTripleType>, Ball3 >
  cout << MetricL2::distance< Ball3, Triangle3<VertexTripleType> >
#endif
}
