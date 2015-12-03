#include "../../Common.hpp"
#include "../../Algorithms/MeshFeatures/Global/DistanceHistogram.hpp"
#include "../../Algorithms/MeshFeatures/Local/Curvature.hpp"
#include "../../Algorithms/MeshFeatures/Local/ShapeDiameter.hpp"
#include "../../Algorithms/BestFitSphere3.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Graphics/GeneralMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Array.hpp"
#include "../../IOStream.hpp"
#include "../../Vector3.hpp"
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace Thea;
using namespace Graphics;
using namespace Algorithms;

typedef GeneralMesh<> Mesh;
typedef MeshGroup<Mesh> MG;

struct MeshScaleType
{
  enum Value
  {
    BSPHERE,
    BBOX,
    AVG_DIST,
  };

  THEA_ENUM_CLASS_BODY(MeshScaleType)

  string toString() const
  {
    switch (*this)
    {
      case BSPHERE: return "bounding sphere";
      case BBOX: return "bounding box";
      case AVG_DIST: return "average distance from centroid";
      default: return "unknown";
    }
  }
};

MeshScaleType mesh_scale_type = MeshScaleType::BSPHERE;
double mesh_scale = 1;
bool is_oriented = false;  // all normals point outwards
bool abs_values = false;

int usage(int argc, char * argv[]);
double meshScale(MG & mg, MeshScaleType mesh_scale_type);
bool computeDistanceHistogram(MG const & mg, long num_bins, long num_samples, DistanceType dist_type, double max_distance,
                              double reduction_ratio, TheaArray<double> & values);
bool computeCurvatureHistogram(MG const & mg, long num_bins, long num_samples, double reduction_ratio,
                               TheaArray<double> & values);
bool computeSDFHistogram(MG const & mg, long num_bins, long num_samples, TheaArray<double> & values);

int
main(int argc, char * argv[])
{
  string mesh_path;
  string out_path;

  bool shift_to_01 = false;
  bool feat_scale = false;
  double feat_scale_factor = 1;
  bool binary = false;

  int curr_opt = 0;
  for (int i = 1; i < argc; ++i)
  {
    string arg = argv[i];
    if (!beginsWith(arg, "--"))
    {
      switch (curr_opt)
      {
        case 0: mesh_path = arg; break;
        case 1: out_path = arg; break;
        default:
        {
          THEA_ERROR << "Too many positional arguments";
          return -1;
        }
      }

      curr_opt++;
    }
    else if (arg == "--abs")
    {
      abs_values = true;
    }
    else if (arg == "--binary")
    {
      binary = true;
    }
    else if (beginsWith(arg, "--featscale="))
    {
      if (sscanf(arg.c_str(), "--featscale=%lf", &feat_scale_factor) == 1)
        feat_scale = true;
      else
      {
        THEA_ERROR << "Couldn't parse feat_scale factor";
        return -1;
      }
    }
    else if (arg == "--is-oriented")
    {
      is_oriented = true;
    }
    else if (beginsWith(arg, "--meshscale="))
    {
      string type = arg.substr(strlen("--meshscale="));
      if (type == "bsphere")
        mesh_scale_type = MeshScaleType::BSPHERE;
      else if (type == "bbox")
        mesh_scale_type = MeshScaleType::BBOX;
      else if (type == "avgdist")
        mesh_scale_type = MeshScaleType::AVG_DIST;
      else
      {
        THEA_ERROR << "Unsupported mesh scale type";
        return -1;
      }
    }
    else if (arg == "--shift01")
    {
      shift_to_01 = true;
    }
    else
      continue;

    argv[i][0] = 0;  // zero out the argument so it won't be considered as a feature later
  }

  if (curr_opt < 2)
    return usage(argc, argv);

  MG mg;
  try
  {
    mg.load(mesh_path);
    mesh_scale = meshScale(mg, mesh_scale_type);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load mesh %s", mesh_path.c_str())

  THEA_CONSOLE << "Loaded mesh from " << mesh_path << " with scale " << mesh_scale << " (based on "
               << mesh_scale_type.toString() << ')';

  // Compute features
  TheaArray<double> features;
  TheaArray<string> feat_names;

  for (int i = 1; i < argc; ++i)
  {
    string feat = string(argv[i]);
    if (beginsWith(feat, "--ch="))
    {
      long num_bins, num_samples;
      double reduction_ratio;

      long num_params = sscanf(feat.c_str(), "--ch=%ld,%ld,%lf", &num_bins, &num_samples, &reduction_ratio);
      if (num_params < 1)
      {
        THEA_ERROR << "Couldn't parse curvature histogram parameters";
        return -1;
      }
      else if (num_params < 3)
      {
        THEA_WARNING << "Sample reduction ratio not specified for curvature histogram, using default value";
        reduction_ratio = -1;

        if (num_params < 2)
        {
          THEA_WARNING << "Number of samples not specified for curvature histogram, using default value";
          num_samples = -1;
        }
      }

      if (num_bins <= 0)
      {
        THEA_ERROR << "Number of histogram bins must be > 0";
        return -1;
      }

      TheaArray<double> values;
      if (!computeCurvatureHistogram(mg, num_bins, num_samples, reduction_ratio, values))
        return -1;

      features.insert(features.end(), values.begin(), values.end());
    }
    else if (beginsWith(feat, "--dh="))
    {
      char dist_str[260];
      DistanceType dist_type;
      long num_bins, num_samples;
      double max_distance;
      double reduction_ratio;

      long num_params = sscanf(feat.c_str(), "--dh=%256[^,],%ld,%ld,%lf,%lf",
                               dist_str, &num_bins, &num_samples, &max_distance, &reduction_ratio);
      if (num_params < 2)
      {
        THEA_ERROR << "Couldn't parse distance histogram parameters";
        return -1;
      }

      if (!dist_type.fromString(dist_str))
      {
        THEA_ERROR << "Unknown distance metric: " << dist_str;
        return -1;
      }

      if (num_params < 5)
      {
        THEA_WARNING << "Sample reduction ratio not specified for distance histogram, using default value";
        reduction_ratio = -1;

        if (num_params < 4)
        {
          THEA_WARNING << "Distance limit for not specified for distance histogram, using default of mesh scale";
          max_distance = -1;

          if (num_params < 3)
          {
            THEA_WARNING << "Number of samples not specified for distance histogram, using default value";
            num_samples = -1;
          }
        }
      }

      if (num_bins <= 0)
      {
        THEA_ERROR << "Number of histogram bins must be > 0";
        return -1;
      }

      TheaArray<double> values;
      if (!computeDistanceHistogram(mg, num_bins, num_samples, dist_type, max_distance, reduction_ratio, values))
        return -1;

      features.insert(features.end(), values.begin(), values.end());
    }
    else if (beginsWith(feat, "--sdf="))
    {
      long num_bins, num_samples;

      long num_params = sscanf(feat.c_str(), "--sdf=%ld,%ld", &num_bins, &num_samples);
      if (num_params < 1)
      {
        THEA_ERROR << "Couldn't parse SDF histogram parameters";
        return -1;
      }
      else if (num_params < 2)
      {
        THEA_WARNING << "Number of samples not specified for SDF histogram, using default value";
        num_samples = -1;
      }

      if (num_bins <= 0)
      {
        THEA_ERROR << "Number of histogram bins must be > 0";
        return -1;
      }

      TheaArray<double> values;
      if (!computeSDFHistogram(mg, num_bins, num_samples, values))
        return -1;

      features.insert(features.end(), values.begin(), values.end());
    }
    else
    {
      if (!feat.empty())
        THEA_WARNING << "Ignoring unsupported option: " << feat;

      continue;
    }

    feat_names.push_back(feat.substr(2));
  }

  ostringstream feat_str;
  for (array_size_t i = 0; i < feat_names.size(); ++i)
  {
    if (i > 0) feat_str << ", ";
    feat_str << feat_names[i];
  }

  THEA_CONSOLE << "Computed " << feat_names.size() << " feature(s): " << feat_str.str();

  // Write features to file
  if (binary)
  {
    BinaryOutputStream out(out_path, Endianness::LITTLE);
    if (!out.ok())
    {
      THEA_ERROR << "Could not open output file " << out_path << " for writing features";
      return -1;
    }

    out.writeInt64((int64)features.size());
    for (array_size_t i = 0; i < features.size(); ++i)
    {
      double f = features[i];

      if (feat_scale) f *= feat_scale_factor;
      if (shift_to_01) f = 0.5 * (1.0 + f);
      if (abs_values) f = fabs(f);

      out.writeFloat32((float32)f);
    }
  }
  else
  {
    ofstream out(out_path.c_str());
    if (!out)
    {
      THEA_ERROR << "Could not open output file " << out_path << " for writing features";
      return -1;
    }

    for (array_size_t i = 0; i < features.size(); ++i)
    {
      double f = features[i];

      if (feat_scale) f *= feat_scale_factor;
      if (shift_to_01) f = 0.5 * (1.0 + f);
      if (abs_values) f = fabs(f);

      if (i > 0) out << ' ';
      out << f;
    }

    out << '\n';
  }

  THEA_CONSOLE << "Wrote " << features.size() << " feature(s) to " << out_path;

  return 0;
}

int
usage(int argc, char * argv[])
{
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Usage: " << argv[0] << " <mesh> <outfile> [<feature0> <feature1> ...]";
  THEA_CONSOLE << "    <featureN> must be one of:";
  THEA_CONSOLE << "        --ch=<num-bins>[,<num-samples>[,<reduction-ratio>]]";
  THEA_CONSOLE << "        --dh=<metric>,<num-bins>[,<num-samples>[,<max_distance>[,<reduction-ratio>]]]";
  THEA_CONSOLE << "        --sdf=<num-bins>[,<num-samples>]";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "    The following options may also be specified:";
  THEA_CONSOLE << "        --abs (uses the absolute value of every feature)";
  THEA_CONSOLE << "        --binary (outputs features in binary format)";
  THEA_CONSOLE << "        --featscale=<factor> (scales feature values by the factor)";
  THEA_CONSOLE << "        --is-oriented (assumes mesh normals consistently point outward)";
  THEA_CONSOLE << "        --meshscale={bsphere|bbox|avgdist} (used to set neighborhood scales)";
  THEA_CONSOLE << "        --normalize (rescale mesh so --meshscale == 1)";
  THEA_CONSOLE << "        --shift01 (maps features in [-1, 1] to [0, 1])";
  THEA_CONSOLE << "";

  return -1;
}

double
meshScale(MG & mg, MeshScaleType mesh_scale_type)
{
  switch (mesh_scale_type)
  {
    case MeshScaleType::BBOX:
    {
      mg.updateBounds();
      return mg.getBounds().getExtent().length();
    }

    case MeshScaleType::AVG_DIST:
    {
      MeshSampler<Mesh> sampler(mg);
      TheaArray<Vector3> samples;
      sampler.sampleEvenlyByArea(50000, samples);

      if (samples.size() <= 1)
        return 0.0;

      Vector3 centroid = CentroidN<Vector3, 3>::compute(samples.begin(), samples.end());

      double sum_dists = 0;
      for (array_size_t i = 0; i < samples.size(); ++i)
        sum_dists += (samples[i] - centroid).length();

      return (3.0 * sum_dists) / samples.size();  // 3 instead of 2 since the avg is in general smaller than the max
    }

    default:
    {
      BestFitSphere3 bsphere;
      PointCollectorN<BestFitSphere3, 3>(&bsphere).addMeshVertices(mg);
      return bsphere.getDiameter();
    }
  }
}

bool
computeDistanceHistogram(MG const & mg, long num_bins, long num_samples, DistanceType dist_type, double max_distance,
                         double reduction_ratio, TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing " << dist_type.toString() << " distance histogram";

  MeshFeatures::Global::DistanceHistogram<> dh(mg, num_samples, (Real)mesh_scale);

  values.resize((array_size_t)num_bins);
  Histogram histogram(num_bins, &values[0]);
  dh.compute(histogram, dist_type, (Real)max_distance, (Real)reduction_ratio);

  histogram.normalize();

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeCurvatureHistogram(MG const & mg, long num_bins, long num_samples, double reduction_ratio, TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing curvature histogram";

  values.resize((array_size_t)num_bins);
  Histogram histogram(num_bins, &values[0], (abs_values ? 0.0 : -1.0), 1.0);

  MeshFeatures::Local::Curvature<> projcurv(mg, num_samples, (Real)mesh_scale);

  if (num_samples < 0)
    num_samples = projcurv.numSamples();

  if (reduction_ratio < 0)
    reduction_ratio = 5000.0 / num_samples;

  long num_queries = Math::clamp((long)std::ceil(reduction_ratio * num_samples), 0L, num_samples - 1);
  TheaArray<int32> query_indices((array_size_t)num_queries);
  Random::common().sortedIntegers(0, (int32)num_samples - 1, (int32)num_queries, &query_indices[0]);

  for (array_size_t i = 0; i < query_indices.size(); ++i)
  {
    long index = query_indices[i];
    Vector3 p = projcurv.getSamplePosition(index);
    Vector3 n = projcurv.getSampleNormal(index);

    double curv = projcurv.computeProjectedCurvature(p, n);
    if (abs_values)
      curv = fabs(curv);

    histogram.insert(curv);
  }

  histogram.normalize();

  THEA_CONSOLE << "  -- done";

  return true;
}

bool
computeSDFHistogram(MG const & mg, long num_bins, long num_samples, TheaArray<double> & values)
{
  THEA_CONSOLE << "Computing SDF histogram";

  values.resize((array_size_t)num_bins);
  Histogram histogram(num_bins, &values[0], 0.0, 1.0);

  MeshFeatures::Local::ShapeDiameter<Mesh> sdf(mg, (Real)mesh_scale);

  if (num_samples < 0)
    num_samples = 5000;

  MeshSampler<Mesh> sampler(mg);
  TheaArray<Vector3> positions, normals;
  sampler.sampleEvenlyByArea(num_samples, positions, &normals);

  for (long i = 0; i < num_samples; ++i)
  {
    Vector3 p = positions[(array_size_t)i];
    Vector3 n = normals[(array_size_t)i];

    double v0 = sdf.compute(p, n, true);
    if (v0 < 0)
      v0 = sdf.compute(p, n, false);

    if (!is_oriented)
    {
      double v1 = sdf.compute(p, -n, true);
      if (v1 < 0)
        v1 = sdf.compute(p, -n, false);

      v0 = ((v1 < 0 || (v0 >= 0 && v0 < v1)) ? v0 : v1);
    }

    histogram.insert(v0);
  }

  histogram.normalize();

  THEA_CONSOLE << "  -- done";

  return true;
}
