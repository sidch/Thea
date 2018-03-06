#include "../Common.hpp"
#include "../Array.hpp"
#include "../Image.hpp"
#include "../ImageMatrix.hpp"
#include "../Math.hpp"
#include "../Vector3.hpp"
#include "../Vector4.hpp"
#include "../VectorN.hpp"
#include "../Algorithms/Zernike2.hpp"
#include <complex>
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

double
mag(std::complex<double> const & c)
{
  return abs(c);
}

template <long N, typename T>
VectorN<N, T> mag(VectorN< N, std::complex<T> > const & c)
{
  VectorN<N, T> result;
  for (long i = 0; i < N; ++i)
    result[i] = abs(c[i]);

  return result;
}

template <typename U, typename T>
int
zernike(ImageMatrix<T> const & m)
{
  Zernike2 zernike;

  TheaArray<U> moments;
  zernike.compute(m, 0.5 * (m.numColumns() - 1), 0.5 * (m.numRows() - 1),
                  0.5 * sqrt(Math::square(m.numColumns()) + Math::square(m.numRows())), moments);

  for (size_t i = 0; i < moments.size(); ++i)
  {
    cout << "Moment magnitude [" << i << "] = " << mag(moments[i]) << endl;
  }

  return 0;
}

int
main(int argc, char * argv[])
{
  try
  {
    if (argc < 2)
    {
      cerr << "Usage: " << argv[0] << " <image-file>" << endl;
      return 0;
    }

    Image image(argv[1]);
    switch (image.getType())
    {
      case Image::Type::LUMINANCE_8U:
        return zernike< std::complex<double> >(ImageMatrix<uint8>(&image));

      case Image::Type::RGB_8U:
        return zernike< VectorN< 3, std::complex<double> > >(ImageMatrix< VectorN<3, uint8> >(&image));

      case Image::Type::RGBA_8U:
        return zernike< VectorN< 4, std::complex<double> > >(ImageMatrix< VectorN<4, uint8> >(&image));

      default:
        cerr << "Unsupported image format" << endl;
        return -1;
    }
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  cout << "Zernike: Test completed" << endl;
  return 0;
}
