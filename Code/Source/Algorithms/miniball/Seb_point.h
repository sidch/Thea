// Synopsis: Simple point class
//
// Authors: Martin Kutz <kutz@math.fu-berlin.de>,
//          Kaspar Fischer <kf@iaeth.ch>

#ifndef SEB_POINT_H
#define SEB_POINT_H

#include "Seb_configure.h"

namespace SEB_NAMESPACE {

  template<unsigned int D, typename Float>
  class Point
  // A simple class representing a d-dimensional point.
  {
  public: // types:
    typedef Float const * Const_iterator;
    typedef Float * Iterator;

  public: // construction and destruction:

    Point()
    // Constructs a d-dimensional point with undefined coordinates.
    {
    }

    template<typename InputIterator>
    Point(InputIterator first)
    // Constructs a d-dimensional point with Cartesian center
    // coordinates [first,first+d).
    {
      for (unsigned int i = 0; i < D; ++i, ++first)
        c[i] = *first;
    }

  public: // access:

    const Float& operator[](unsigned int i) const
    // Returns a const-reference to the i-th coordinate.
    {
      SEB_ASSERT(0 <= i && i < D);
      return c[i];
    }

    Float& operator[](unsigned int i)
    // Returns a reference to the i-th coordinate.
    {
      SEB_ASSERT(0 <= i && i < D);
      return c[i];
    }

    Const_iterator begin() const
    // Returns a const-iterator to the first of the d Cartesian coordinates.
    {
      return c.begin();
    }

    Const_iterator end() const
    // Returns the past-the-end iterator corresponding to begin().
    {
      return c.end();
    }

    unsigned int dims() const
    // Returns the number of dimensions of the point
    {
      return D;
    }

  private: // member fields:
    Float c[D];       // Cartesian center coordinates
  };

} // namespace SEB_NAMESPACE

#endif // SEB_POINT_H
