#include "../Common.hpp"
#include "../Algorithms/IntersectionTester.hpp"
#include "../Algorithms/BvhN.hpp"
#include "../Algorithms/MetricL2.hpp"
#include "../Algorithms/PointTraitsN.hpp"
#include "../Algorithms/RayIntersectionTester.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Ball3.hpp"
#include "../BoundedSortedArrayN.hpp"
#include "../Triangle3.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace Thea;
using namespace Algorithms;

void testPointBvh();
void testTriangleBvh();

int
main(int argc, char * argv[])
{
  try
  {
    testPointBvh();
    cout << endl;
    testTriangleBvh();
  }
  THEA_CATCH(return -1;, ERROR, "%s", "An error occurred")

  // Hooray, all tests passed
  cout << "BvhN: Test completed" << endl;
  return 0;
}

// A custom point class.
struct MyCustomPoint
{
  MyCustomPoint() {}
  MyCustomPoint(string const & name_, Vector3 const & position_) : name(name_), position(position_) {}
  MyCustomPoint(MyCustomPoint const & src) : name(src.name), position(src.position) {}

  string name;
  Vector3 position;
};

// A custom class that stores the positions of a triangle's three vertices, plus the name of the triangle (you can use this to
// store the ID of the triangle instead, or other info)
struct MyCustomTriangleVertexTriple
{
  // Default constructor
  MyCustomTriangleVertexTriple() {}

  // Initializing constructor
  MyCustomTriangleVertexTriple(string const & name_, Vector3 const & v0_, Vector3 const & v1_, Vector3 const & v2_)
  : name(name_)
  {
    v[0] = v0_;
    v[1] = v1_;
    v[2] = v2_;
  }

  // Copy constructor
  MyCustomTriangleVertexTriple(MyCustomTriangleVertexTriple const & src)
  : name(src.name)
  {
    v[0] = src.v[0];
    v[1] = src.v[1];
    v[2] = src.v[2];
  }

  // Get the i'th vertex. This is the operative function that a vertex triple needs to have.
  Vector3 const & getVertex(int i) const { return v[i]; }

  string name;
  Vector3 v[3];
};

// Our custom triangle class just wraps the vertex triple above. We don't need to specify any additional traits classes because
// the BVH already knows how to handle any specialization of Triangle<3, ...>. The vertex triple tells Triangle how to get the
// three vertices of the triangle, and the rest is automatically set up. To access the vertex triple (and any custom info inside
// it), use triangle.getVertices(), which returns a reference to the wrapped vertex triple.
typedef TriangleN<3, MyCustomTriangleVertexTriple, Real> MyCustomTriangle;

namespace Thea {
namespace Algorithms {

// Tell the BVH that a MyCustomPoint object is a 3D point.
template <>
struct IsPointN<MyCustomPoint, 3>
{
  static bool const value = true;
};

// A specialization of the traits class, to obtain the position of a MyCustomPoint. The BVH requires this traits class to
// get the position of an arbitrary point type.
template <>
struct PointTraitsN<MyCustomPoint, 3, Real>
{
  static Vector3 const & getPosition(MyCustomPoint const & np) { return np.position; }
};

} // namespace Algorithms
} // namespace Thea

// This function is called by the BVH for every point encountered in a range. The args are the index of the point (in the array
// returned by bvh.getElements()) and a reference to the point itself (as cached by the BVH).
bool printPoint(intx index, MyCustomPoint & np)
{
  cout << "  Found point '" << np.name << "' at position " << toString(np.position) << endl;
  return false;  // the range query stops when this function returns true -- here we don't want that to happen
}

void
testPointBvh()
{
  cout << "=========================\n"
       << "Testing BVH on points\n"
       << "=========================" << endl;

  //============================================================================================================================
  // Generate input data
  //============================================================================================================================

  // Create a bunch of named points
  vector<MyCustomPoint> points;
  static int NUM_POINTS = 100;
  for (int i = 0; i < NUM_POINTS; ++i)
  {
    ostringstream oss; oss << "Point " << i;
    Vector3 random_pt(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX);
    MyCustomPoint p(oss.str(), random_pt);
    points.push_back(p);
  }
  cout << "Generated " << points.size() << " random points" << endl;

  //============================================================================================================================
  // Create BVH for the data
  //============================================================================================================================

  // Create a BVH for all the points
  typedef BvhN<MyCustomPoint, 3> Bvh;
  Bvh bvh(points.begin(), points.end());  // To reinitialize the tree later, call bvh.init(begin, end). For fast
                                                // reinitialization, set the deallocate_previous_memory arg of init() to false.
  cout << "Created BVH for points" << endl;

  //============================================================================================================================
  // Three ways of doing a range query
  //
  // IMPORTANT: Range queries can return the same element twice, both in the functor-based approach and when returning vectors!
  //
  //============================================================================================================================

  // A range query that prints all the points in a ball. This functor-based approach is the most general way of doing a range
  // query
  Ball3 ball(Vector3(0.5f, 0.5f, 0.5f), 0.25f);
  cout << "\nLooking for all points in ball " << ball.toString() << ':' << endl;
  bvh.processRange<IntersectionTester>(ball, printPoint);  // processes all points in the ball until the functor returns true

  // Another way of writing a range query, that explicitly returns all the elements in the ball. This might be slower because of
  // memory allocation by push_back in vector.
  vector<MyCustomPoint> elems_in_range;
  bvh.rangeQuery<IntersectionTester>(ball, elems_in_range);
  cout << "\nRange query returned " << elems_in_range.size() << " elements (with possible duplications)" << endl;

  // Yet another way of writing a range query, that returns the indices of all the elements in the ball. This might also be
  // slower because of memory allocation by push_back in vector.
  vector<intx> indices_of_elems_in_range;
  bvh.rangeQueryIndices<IntersectionTester>(ball, indices_of_elems_in_range);
  cout << "\nRange query returned " << elems_in_range.size() << " indices (with possible duplications)" << endl;

  //============================================================================================================================
  // Using an axis-aligned box as the query range instead
  //============================================================================================================================

  // A range query that prints all the points in a box
  AxisAlignedBox3 box(Vector3(0.25f, 0.25f, 0.25f), Vector3(0.75f, 0.75f, 0.75f));
  cout << "\nLooking for all points in box " << box.toString() << ':' << endl;
  bvh.processRange<IntersectionTester>(box, printPoint);

  //============================================================================================================================
  // Finding the nearest neighbor of a query point
  //============================================================================================================================

  // Find the point nearest to a query point, using the L2 norm
  Vector3 query(0.5f, 0.5f, 0.5f);
  double dist_bound = 0.25;  // we'll limit the search to all points within a distance of 0.25; passing -1 turns this off
  double dist = 0;  // this will contain the distance to the returned point
  intx nn_index = bvh.closestElement<MetricL2>(query, dist_bound, UniversalCompatibility(), &dist);
  if (nn_index >= 0)
    cout << "\nThe point nearest the query " << toString(query) << " is " << bvh.getElements()[nn_index].name
         << " at distance " << dist << endl;
  else
    cout << "\nNo nearest neighbor found" << endl;

  //============================================================================================================================
  // Finding the k-nearest neighbors of a query point
  //============================================================================================================================

  // Find the 3 points nearest to the query point defined above, using the L2 norm and the same upper bound on the distance
  BoundedSortedArrayN<3, Bvh::NeighborPair> nbrs;
  intx num_nbrs = bvh.kClosestPairs<MetricL2>(query, nbrs, dist_bound, UniversalCompatibility());
  if (num_nbrs > 0)
  {
    cout << '\n' << num_nbrs << " neighbors (max 3) found for query " << toString(query) << ':' << endl;
    for (size_t i = 0; i < nbrs.size(); ++i)
    {
      cout << "  " << bvh.getElements()[nbrs[i].getTargetIndex()].name << " at distance "
           << nbrs[i].getDistance<MetricL2>() << endl;
    }
  }
  else
    cout << "\nNo nearest neighbor found (max 3)" << endl;

  //============================================================================================================================
  // Finding the nearest neighbors between two sets of points
  //============================================================================================================================

  // Generate a second set of points, and find the nearest pair of points between the two sets
  vector<MyCustomPoint> new_points;
  for (int i = 0; i < NUM_POINTS; ++i)
  {
    ostringstream oss; oss << "NewPoint " << i;
    Vector3 random_pt(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX);
    MyCustomPoint p(oss.str(), random_pt);
    new_points.push_back(p);
  }
  Bvh new_bvh(new_points.begin(), new_points.end());

  typedef Bvh::NeighborPair NeighborPair;

  // -1 means there's no limit on the maximum allowed separation. Setting a (small) positive value can make this query run
  // *much* faster.
  NeighborPair nn_pair = bvh.closestPair<MetricL2>(new_bvh, -1, UniversalCompatibility(), true);
  if (nn_pair.isValid())
  {
    cout << "\nThe nearest neighbors are "
         << new_bvh.getElements()[nn_pair.getQueryIndex()].name << ' ' << toString(nn_pair.getQueryPoint()) << " and "
         << bvh.getElements()[nn_pair.getTargetIndex()].name << ' ' << toString(nn_pair.getTargetPoint())
         << " at separation " << nn_pair.getDistance<MetricL2>() << endl;

    cout << "    Query point is at " << toString(new_bvh.getElements()[nn_pair.getQueryIndex()].position) << endl;
    cout << "    Target point is at " << toString(bvh.getElements()[nn_pair.getTargetIndex()].position) << endl;
  }
  else
    cout << "\nNo nearest neighbors found between the two BVHs" << endl;

  //============================================================================================================================
  // Finding the k-nearest neighbors of a set of query points
  //============================================================================================================================

  // Find the 3 points nearest to the query point defined above, using the L2 norm and the same upper bound on the distance
  nbrs.clear();  // not necessary but let's do this anyway
  num_nbrs = bvh.kClosestPairs<MetricL2>(new_bvh, nbrs, -1, UniversalCompatibility(), true);
  if (num_nbrs > 0)
  {
    cout << '\n' << num_nbrs << " pairs of nearest neighbors (max 3) found for query point set:" << endl;
    for (intx i = 0; i < nbrs.size(); ++i)
    {
      cout << "  (" << new_bvh.getElements()[nbrs[i].getQueryIndex()].name << ' ' << toString(nbrs[i].getQueryPoint()) << ", "
                    << bvh.getElements()[nbrs[i].getTargetIndex()].name << ' ' << toString(nbrs[i].getTargetPoint())
                    << ") at distance " << nbrs[i].getDistance<MetricL2>() << endl;

      cout << "      Query point is at " << toString(new_bvh.getElements()[nbrs[i].getQueryIndex()].position) << endl;
      cout << "      Target point is at " << toString(bvh.getElements()[nbrs[i].getTargetIndex()].position) << endl;
    }
  }
  else
    cout << "\nNo nearest neighbor found for the query point set (max 3)" << endl;
}

void
testTriangleBvh()
{
  cout << "============================\n"
       << "Testing BVH on triangles\n"
       << "============================" << endl;

  //============================================================================================================================
  // Generate input data
  //============================================================================================================================

  // Create a bunch of triangles
  vector<MyCustomTriangle> triangles;
  static int NUM_TRIANGLES = 100;
  for (int i = 0; i < NUM_TRIANGLES; ++i)
  {
    ostringstream oss; oss << "Triangle " << i;
    Vector3 v[3] = { Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX),
                     Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX),
                     Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX) };

    MyCustomTriangle tri(MyCustomTriangleVertexTriple(oss.str(), v[0], v[1], v[2]));
    triangles.push_back(tri);
  }
  cout << "Generated " << triangles.size() << " random triangles" << endl;

  //============================================================================================================================
  // Create BVH for the data
  //============================================================================================================================

  // Create a BVH for all the triangles
  typedef BvhN<MyCustomTriangle, 3> Bvh;
  Bvh bvh(triangles.begin(), triangles.end());  // To reinitialize the tree later, call bvh.init(begin, end). For fast
                                                      // reinitialization, set the deallocate_previous_memory arg of init() to
                                                      // false.
  bvh.enableNearestNeighborAcceleration();
  cout << "Created BVH for triangles" << endl;

  //============================================================================================================================
  // Finding the nearest neighbor of a query point
  //============================================================================================================================

  // Find the triangle nearest to a query point, using the L2 norm
  Vector3 query(0.5f, 0.5f, 0.5f);
  double dist_bound = 0.25;  // we'll limit the search to all triangles within a distance of 0.25; passing -1 turns this off
  double dist = 0;  // this will contain the distance to the nearest triangle
  Vector3 closest_point;  // this will contain the closest point on the nearest triangle
  intx nn_index = bvh.closestElement<MetricL2>(query, dist_bound, UniversalCompatibility(), &dist, &closest_point);
  if (nn_index >= 0)
    cout << "\nThe triangle nearest the query " << query.transpose() << " is "
         << bvh.getElements()[nn_index].getVertices().name
         << " at distance " << dist << ", with closest point" << toString(closest_point) << endl;
  else
    cout << "\nNo nearest neighbor of the query point found" << endl;

  //============================================================================================================================
  // Finding the nearest neighbors between two sets of triangles. Slow, because of triangle-triangle distance tests. Is
  // accelerated by setting an upper bound on the maximum separation between the nearest neighbors (this replaces the -1 in the
  // call to closestPair()).
  //============================================================================================================================

  // Generate a second set of triangles, and find the nearest pair of triangles between the two sets
  vector<MyCustomTriangle> new_triangles;
  for (int i = 0; i < NUM_TRIANGLES; ++i)
  {
    ostringstream oss; oss << "NewTriangle " << i;
    Vector3 v[3] = { Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX),
                     Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX),
                     Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX) };

    MyCustomTriangle tri(MyCustomTriangleVertexTriple(oss.str(), v[0], v[1], v[2]));
    new_triangles.push_back(tri);
  }
  Bvh new_bvh(new_triangles.begin(), new_triangles.end());
  new_bvh.enableNearestNeighborAcceleration();

  typedef Bvh::NeighborPair NeighborPair;

  NeighborPair nn_pair = bvh.closestPair<MetricL2>(new_bvh, -1, UniversalCompatibility(), true);
  if (nn_pair.isValid())
    cout << "\nThe nearest neighbors are "
         << new_bvh.getElements()[nn_pair.getQueryIndex()].getVertices().name << " and "
         << bvh.getElements()[nn_pair.getTargetIndex()].getVertices().name
         << " at separation " << nn_pair.getDistance<MetricL2>()
         << ", with closest points " << toString(nn_pair.getQueryPoint()) << " and " << toString(nn_pair.getTargetPoint())
         << " respectively" << endl;
  else
    cout << "\nNo nearest neighbors found between the two BVHs" << endl;

  //============================================================================================================================
  // Ray-triangle intersection. There are three types of intersection queries, each giving strictly more information than the
  // last. They should all take roughly the same time for BVHs of triangles.
  //============================================================================================================================

  // Generate a random ray from the origin into the positive quadrant (the ray's direction vector need not be a unit vector,
  // but here we'll use one)
  Ray3 ray(Vector3::Zero(), Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX).normalized());

  // For now we won't limit the maximum hit time of the ray
  Real max_time = -1;

  cout << "\nTesting intersection with ray " << ray.toString() << endl;

  // First type of query: does the ray intersect any triangle in the tree?
  if (bvh.rayIntersects<RayIntersectionTester>(ray, max_time))
    cout << "Ray intersects a triangle in the BVH" << endl;
  else
    cout << "Ray does not intersect any triangle in the BVH" << endl;

  // Second type of query: what is the hit time of the ray?
  Real hit_time = bvh.rayIntersectionTime<RayIntersectionTester>(ray, max_time);
  if (hit_time >= 0)
    cout << "Ray intersects a triangle in the BVH after time " << hit_time
         << " (at point " << toString(ray.getPoint(hit_time)) << ')' << endl;
  else
    cout << "Ray does not intersect any triangle in the BVH" << endl;

  // Third type of query: full info about the intersection point
  RayStructureIntersection3 isec = bvh.rayStructureIntersection<RayIntersectionTester>(ray, max_time);
  if (isec.isValid())
  {
    cout << "Ray intersects a triangle in the BVH:\n"
         << "    hit time = " << isec.getTime() << '\n'
         << "    intersected triangle = " << bvh.getElements()[isec.getElementIndex()].getVertices().name << '\n'
         << "    intersection point = " << toString(ray.getPoint(isec.getTime())) << '\n'
         << "    normal at intersection point = " << toString(isec.getNormal()) << endl;
  }
  else
    cout << "Ray does not intersect any triangle in the BVH" << endl;
}
