#include "../Common.hpp"
#include "../Algorithms/KDTreeN.hpp"
#include "../Algorithms/MetricL2.hpp"
#include "../Algorithms/IntersectionTester.hpp"
#include "../Algorithms/RayIntersectionTester.hpp"
#include "../AxisAlignedBox3.hpp"
#include "../Ball3.hpp"
#include "../BoundedSortedArrayN.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace Thea;
using namespace Algorithms;

void testPointKDTree();
void testTriangleKDTree();

int
main(int argc, char * argv[])
{
  try
  {
    testPointKDTree();
    cout << endl;
    testTriangleKDTree();
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  // Hooray, all tests passed
  cout << "KDTreeN: Test completed" << endl;
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
// the kd-tree already knows how to handle any specialization of Triangle3<...>. The vertex triple tells Triangle3 how to get
// the three vertices of the triangle, and the rest is automatically set up. To access the vertex triple (and any custom info
// inside it), use triangle.getVertices(), which returns a reference to the wrapped vertex triple.
typedef Triangle3<MyCustomTriangleVertexTriple> MyCustomTriangle;

namespace Thea {
namespace Algorithms {

// Tell the kd-tree that a MyCustomPoint object is a 3D point.
template <>
struct IsPointN<MyCustomPoint, 3>
{
  static bool const value = true;
};

// A specialization of the traits class, to obtain the position of a MyCustomPoint. The kd-tree requires this traits class to
// get the position of an arbitrary point type.
template <>
struct PointTraitsN<MyCustomPoint, 3>
{
  static Vector3 const & getPosition(MyCustomPoint const & np) { return np.position; }
};

} // namespace Algorithms
} // namespace Thea

// This function is called by the kd-tree for every point encountered in a range. The args are the index of the point (in the
// array returned by kdtree.getElements()) and a reference to the point itself (as cached by the kd-tree).
bool printPoint(long index, MyCustomPoint & np)
{
  cout << "  Found point '" << np.name << "' at position " << np.position.toString() << endl;
  return false;  // the range query stops when this function returns true -- here we don't want that to happen
}

void
testPointKDTree()
{
  cout << "=========================\n"
       << "Testing kd-tree on points\n"
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
  // Create kd-tree for the data
  //============================================================================================================================

  // Create a kd-tree for all the points
  typedef KDTreeN<MyCustomPoint, 3> KDTree;
  KDTree kdtree(points.begin(), points.end());  // To reinitialize the tree later, call kdtree.init(begin, end). For fast
                                                // reinitialization, set the deallocate_previous_memory arg of init() to false.
  cout << "Created kd-tree for points" << endl;

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
  kdtree.processRangeUntil<IntersectionTester>(ball, printPoint);  // processes all points in the ball until the functor returns
                                                                   // true

  // Another way of writing a range query, that explicitly returns all the elements in the ball. This might be slower because of
  // memory allocation by push_back in vector.
  vector<MyCustomPoint> elems_in_range;
  kdtree.rangeQuery<IntersectionTester>(ball, elems_in_range);
  cout << "\nRange query returned " << elems_in_range.size() << " elements (with possible duplications)" << endl;

  // Yet another way of writing a range query, that returns the indices of all the elements in the ball. This might also be
  // slower because of memory allocation by push_back in vector.
  vector<long> indices_of_elems_in_range;
  kdtree.rangeQueryIndices<IntersectionTester>(ball, indices_of_elems_in_range);
  cout << "\nRange query returned " << elems_in_range.size() << " indices (with possible duplications)" << endl;

  //============================================================================================================================
  // Using an axis-aligned box as the query range instead
  //============================================================================================================================

  // A range query that prints all the points in a box
  AxisAlignedBox3 box(Vector3(0.25f, 0.25f, 0.25f), Vector3(0.75f, 0.75f, 0.75f));
  cout << "\nLooking for all points in box " << box.toString() << ':' << endl;
  kdtree.processRangeUntil<IntersectionTester>(box, printPoint);

  //============================================================================================================================
  // Finding the nearest neighbor of a query point
  //============================================================================================================================

  // Find the point nearest to a query point, using the L2 norm
  Vector3 query(0.5f, 0.5f, 0.5f);
  double dist_bound = 0.25;  // we'll limit the search to all points within a distance of 0.25; passing -1 turns this off
  double dist = 0;  // this will contain the distance to the returned point
  long nn_index = kdtree.closestElement<MetricL2>(query, dist_bound, &dist);
  if (nn_index >= 0)
    cout << "\nThe point nearest the query " << query.toString() << " is " << kdtree.getElements()[nn_index].name
         << " at distance " << dist << endl;
  else
    cout << "\nNo nearest neighbor found" << endl;

  //============================================================================================================================
  // Finding the k-nearest neighbors of a query point
  //============================================================================================================================

  // Find the 3 points nearest to the query point defined above, using the L2 norm and the same upper bound on the distance
  BoundedSortedArrayN<3, KDTree::NeighborPair> nbrs;
  long num_nbrs = kdtree.kClosestPairs<MetricL2>(query, nbrs, dist_bound);
  if (num_nbrs > 0)
  {
    cout << '\n' << num_nbrs << " neighbors (max 3) found for query " << query.toString() << ':' << endl;
    for (long i = 0; i < nbrs.size(); ++i)
    {
      cout << "  " << kdtree.getElements()[nbrs[i].getTargetIndex()].name << " at distance "
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
  KDTree new_kdtree(new_points.begin(), new_points.end());

  typedef KDTree::NeighborPair NeighborPair;
  NeighborPair nn_pair = kdtree.closestPair<MetricL2>(new_kdtree, -1, true);  // -1 means there's no limit on the maximum
                                                                              // allowed separation. Setting a (small) positive
                                                                              // value can make this query run *much* faster.
  if (nn_pair.isValid())
  {
    cout << "\nThe nearest neighbors are "
         << new_kdtree.getElements()[nn_pair.getQueryIndex()].name << ' ' << nn_pair.getQueryPoint() << " and "
         << kdtree.getElements()[nn_pair.getTargetIndex()].name << ' ' << nn_pair.getTargetPoint()
         << " at separation " << nn_pair.getDistance<MetricL2>() << endl;

    cout << "    Query point is at " << new_kdtree.getElements()[nn_pair.getQueryIndex()].position << endl;
    cout << "    Target point is at " << kdtree.getElements()[nn_pair.getTargetIndex()].position << endl;
  }
  else
    cout << "\nNo nearest neighbors found between the two kd-trees" << endl;

  //============================================================================================================================
  // Finding the k-nearest neighbors of a set of query points
  //============================================================================================================================

  // Find the 3 points nearest to the query point defined above, using the L2 norm and the same upper bound on the distance
  nbrs.clear();  // not necessary but let's do this anyway
  num_nbrs = kdtree.kClosestPairs<MetricL2>(new_kdtree, nbrs, -1, true);
  if (num_nbrs > 0)
  {
    cout << '\n' << num_nbrs << " pairs of nearest neighbors (max 3) found for query point set:" << endl;
    for (long i = 0; i < nbrs.size(); ++i)
    {
      cout << "  (" << new_kdtree.getElements()[nbrs[i].getQueryIndex()].name << ' ' << nbrs[i].getQueryPoint() << ", "
                    << kdtree.getElements()[nbrs[i].getTargetIndex()].name << ' ' << nbrs[i].getTargetPoint()
                    << ") at distance " << nbrs[i].getDistance<MetricL2>() << endl;

      cout << "      Query point is at " << new_kdtree.getElements()[nbrs[i].getQueryIndex()].position << endl;
      cout << "      Target point is at " << kdtree.getElements()[nbrs[i].getTargetIndex()].position << endl;
    }
  }
  else
    cout << "\nNo nearest neighbor found for the query point set (max 3)" << endl;
}

void
testTriangleKDTree()
{
  cout << "============================\n"
       << "Testing kd-tree on triangles\n"
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
  // Create kd-tree for the data
  //============================================================================================================================

  // Create a kd-tree for all the triangles
  typedef KDTreeN<MyCustomTriangle, 3> KDTree;
  KDTree kdtree(triangles.begin(), triangles.end());  // To reinitialize the tree later, call kdtree.init(begin, end). For fast
                                                      // reinitialization, set the deallocate_previous_memory arg of init() to
                                                      // false.
  kdtree.enableNearestNeighborAcceleration();
  cout << "Created kd-tree for triangles" << endl;

  //============================================================================================================================
  // Finding the nearest neighbor of a query point
  //============================================================================================================================

  // Find the triangle nearest to a query point, using the L2 norm
  Vector3 query(0.5f, 0.5f, 0.5f);
  double dist_bound = 0.25;  // we'll limit the search to all triangles within a distance of 0.25; passing -1 turns this off
  double dist = 0;  // this will contain the distance to the nearest triangle
  Vector3 closest_point;  // this will contain the closest point on the nearest triangle
  long nn_index = kdtree.closestElement<MetricL2>(query, dist_bound, &dist, &closest_point);
  if (nn_index >= 0)
    cout << "\nThe triangle nearest the query " << query.toString() << " is "
         << kdtree.getElements()[nn_index].getVertices().name
         << " at distance " << dist << ", with closest point" << closest_point.toString() << endl;
  else
    cout << "\nNo nearest neighbor of the query point found" << endl;

  //============================================================================================================================
  // Finding the nearest neighbors between two sets of triangles. Slow, because of triangle-triangle distance tests. Can be
  // accelerated by setting an upper bound on the maximum separation between the nearest neighbors (this replaces the -1 in the
  // call to closestPair()). This upper bound can be computed by finding the distance between two auxiliary kd-trees, each on
  // point samples from the corresponding set of triangles, e.g. the triangle vertices. This is basically the approach taken by
  // CGAL's AABB_Tree::accelerate_distance_queries(). However, I haven't added this to my KD-tree implementation yet, so you'll
  // need to do this explicitly if it's really required.
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
  KDTree new_kdtree(new_triangles.begin(), new_triangles.end());
  new_kdtree.enableNearestNeighborAcceleration();

  typedef KDTree::NeighborPair NeighborPair;
  NeighborPair nn_pair = kdtree.closestPair<MetricL2>(new_kdtree, -1, true);  // -1 means there's no limit on the maximum
                                                                              // allowed separation. Setting a (small) positive
                                                                              // value can make this query run *much* faster.
  if (nn_pair.isValid())
    cout << "\nThe nearest neighbors are "
         << new_kdtree.getElements()[nn_pair.getQueryIndex()].getVertices().name << " and "
         << kdtree.getElements()[nn_pair.getTargetIndex()].getVertices().name
         << " at separation " << nn_pair.getDistance<MetricL2>()
         << ", with closest points " << nn_pair.getQueryPoint().toString() << " and " << nn_pair.getTargetPoint().toString()
         << " respectively" << endl;
  else
    cout << "\nNo nearest neighbors found between the two kd-trees" << endl;

  //============================================================================================================================
  // Ray-triangle intersection. There are three types of intersection queries, each giving strictly more information than the
  // last. They all take roughly the same time for kd-trees of triangles, so there's no special reason to favour
  //============================================================================================================================

  // Generate a random ray from the origin into the positive quadrant (the ray's direction vector need not be a unit vector,
  // but here we'll use one)
  Ray3 ray(Vector3::zero(), Vector3(rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX, rand() / (Real)RAND_MAX).unit());

  // For now we won't limit the maximum hit time of the ray
  Real max_time = -1;

  cout << "\nTesting intersection with ray " << ray.toString() << endl;

  // First type of query: does the ray intersect any triangle in the tree?
  if (kdtree.rayIntersects<RayIntersectionTester>(ray, max_time))
    cout << "Ray intersects a triangle in the kd-tree" << endl;
  else
    cout << "Ray does not intersect any triangle in the kd-tree" << endl;

  // Second type of query: what is the hit time of the ray?
  Real hit_time = kdtree.rayIntersectionTime<RayIntersectionTester>(ray, max_time);
  if (hit_time >= 0)
    cout << "Ray intersects a triangle in the kd-tree after time " << hit_time
         << " (at point " << ray.getPoint(hit_time).toString() << ')' << endl;
  else
    cout << "Ray does not intersect any triangle in the kd-tree" << endl;

  // Third type of query: full info about the intersection point
  RayStructureIntersection3 isec = kdtree.rayStructureIntersection<RayIntersectionTester>(ray, max_time);
  if (isec.isValid())
  {
    cout << "Ray intersects a triangle in the kd-tree:\n"
         << "    hit time = " << isec.getTime() << '\n'
         << "    intersected triangle = " << kdtree.getElements()[isec.getElementIndex()].getVertices().name << '\n'
         << "    intersection point = " << ray.getPoint(isec.getTime()).toString() << '\n'
         << "    normal at intersection point = " << isec.getNormal().toString() << endl;
  }
  else
    cout << "Ray does not intersect any triangle in the kd-tree" << endl;
}
