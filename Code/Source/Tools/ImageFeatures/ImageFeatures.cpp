#include "../../Common.hpp"
#include "../../Algorithms/ImageFeatures/ShapeContext.hpp"
#include "../../Array.hpp"
#include "../../Image.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool computeShapeContext(Image const & image, long num_radial_bins, long num_polar_bins, TheaArray<Real> & values);

int
main(int argc, char * argv[])
{
  string image_path;
  string out_path;

  bool invert = false;

  int curr_opt = 0;
  for (int i = 1; i < argc; ++i)
  {
    if (!beginsWith(argv[i], "--"))
    {
      switch (curr_opt)
      {
        case 0: image_path = argv[i]; break;
        case 1: out_path = argv[i]; break;
      }

      curr_opt++;
    }
    else if (string(argv[i]) == "--invert")
      invert = true;
  }

  if (curr_opt < 2)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " <image> <outfile> [<feature0> <feature1> ...]";
    THEA_CONSOLE << "    <featureN> must be one of:";
    THEA_CONSOLE << "        --sc[:num-radial-bins,num-polar-bins]";
    THEA_CONSOLE << "        --invert (not a feature, inverts pixel values before computing features)";

    return 0;
  }

  Image image;
  try
  {
    image.load(image_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load image %s", image_path.c_str())

  if (invert && !image.invert())
  {
    THEA_ERROR << "Could not invert image";
    return -1;
  }

  int width  = image.getWidth();
  int height = image.getHeight();

  // Compute features
  TheaArray< TheaArray<Real> > features((array_size_t)(width * height));
  TheaArray<string> feat_names;

  for (int i = 1; i < argc; ++i)
  {
    if (!beginsWith(argv[i], "--"))
      continue;

    string feat = string(argv[i]).substr(2);
    if (feat == "sc" || beginsWith(feat, "sc:"))
    {
      long num_radial_bins = 5;
      long num_polar_bins = 12;

      if (feat != "sc")
      {
        if (sscanf(feat.c_str(), "sc:%ld,%ld", &num_radial_bins, &num_polar_bins) != 2)
        {
          THEA_ERROR << "Could not parse number of shape context bins";
          return -1;
        }
      }

      TheaArray<Real> values;
      if (!computeShapeContext(image, num_radial_bins, num_polar_bins, values))
        return -1;

      alwaysAssertM(values.size() == num_radial_bins * num_polar_bins * features.size(),
                    "Number of shape contexts don't match number of pixels");

      Real const * entry = &values[0];
      long entry_size = num_radial_bins * num_polar_bins;
      for (array_size_t j = 0; j < features.size(); ++j, entry += entry_size)
        features[j].insert(features[j].end(), entry, entry + entry_size);
    }
    else if (feat == "invert")
      continue;
    else
    {
      THEA_WARNING << "Ignoring unsupported feature type: " << feat;
      continue;
    }

    feat_names.push_back(feat);
  }

  ostringstream feat_str;
  for (array_size_t i = 0; i < feat_names.size(); ++i)
  {
    if (i > 0) feat_str << ", ";
    feat_str << feat_names[i];
  }

  THEA_CONSOLE << "Computed " << feat_names.size() << " feature(s): " << feat_str.str();

  // Write features to file
  ofstream out(out_path.c_str());
  if (!out)
  {
    THEA_ERROR << "Could not open output file " << out_path << " for writing features";
    return -1;
  }

  long num_written = 0;
  for (int i = 0; i < height; ++i)
  {
    for (int j = 0; j < width; ++j)
    {
      TheaArray<Real> const & entry = features[(array_size_t)(i * width + j)];

      bool all_zero = true;
      for (array_size_t k = 0; k < entry.size(); ++k)
        if (entry[k] > 0)
        {
          all_zero = false;
          break;
        }

      if (all_zero)
        continue;

      out << i << ' ' << j;

      for (array_size_t k = 0; k < entry.size(); ++k)
        out << ' ' << entry[k];

      out << '\n';

      num_written++;
    }
  }

  THEA_CONSOLE << "Wrote " << num_written << '/' << features.size() << " feature vectors to " << out_path;

  return 0;
}

bool
computeShapeContext(Image const & image, long num_radial_bins, long num_polar_bins, TheaArray<Real> & values)
{
  THEA_CONSOLE << "Computing shape contexts";

  values.clear();

  try
  {
    ImageFeatures::ShapeContext sc(image);
    sc.compute(num_radial_bins, num_polar_bins, values);
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "Could not compute shape context")

  THEA_CONSOLE << "  -- done";

  return true;
}
