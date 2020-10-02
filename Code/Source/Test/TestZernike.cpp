#include "../Common.hpp"
#include "../Array.hpp"
#include "../Image.hpp"
#include "../ImageMatrix.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "../Algorithms/Zernike2.hpp"
#include <complex>
#include <iostream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

template <typename T>
int
zernike(ImageMatrix<T> const & m)
{
  Zernike2 zernike;
  Zernike2::MomentMatrix<3, Real> moments;
  zernike.compute(m, 0.5 * (m.cols() - 1), 0.5 * (m.rows() - 1),
                     0.5 * sqrt(Math::square(m.cols()) + Math::square(m.rows())), moments);

  for (intx i = 0; i < moments.cols(); ++i)
    THEA_CONSOLE<< "Moment magnitude [" << i << "] = " << moments.col(i).norm();

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
        return zernike(ImageMatrix<ColorL8>(&image));

      case Image::Type::RGB_8U:
        return zernike(ImageMatrix<ColorRgb8>(&image));

      case Image::Type::RGBA_8U:
        return zernike(ImageMatrix<ColorRgba8>(&image));

      default:
        cerr << "Unsupported image format" << endl;
        return -1;
    }
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "An error occurred")

  cout << "Zernike: Test completed" << endl;
  return 0;
}
