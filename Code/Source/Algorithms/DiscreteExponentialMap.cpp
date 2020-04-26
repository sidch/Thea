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
// First version: 2014
//
//============================================================================

#include "DiscreteExponentialMap.hpp"
#include "SampleGraph.hpp"
#include "ShortestPaths.hpp"
#include "../AffineTransform3.hpp"
#include "../Map.hpp"
#include "../Math.hpp"
#include "../Plane3.hpp"
#include <functional>

namespace Thea {
namespace Algorithms {

namespace DiscreteExponentialMapInternal {

class Impl
{
  public:
    typedef DiscreteExponentialMap::Options Options;
    typedef DiscreteExponentialMap::ParameterMap ParameterMap;

  private:
    struct ParamData
    {
      ParamData const * pred_data;
      Vector2 uv;
      AffineTransform3 uv_transform;
      Vector3 proj_p;
    };

    typedef SampleGraph::SurfaceSample SurfaceSample;
    typedef ShortestPaths<SampleGraph> Geodesics;
    typedef Map<Geodesics::VertexHandle, ParamData> ParamDataMap;  // FIXME: Can we use UnorderedMap?

  public:
    Impl(Options const & options_) : options(options_), radius(0) {}

    void parametrize(SampleGraph const & sample_graph, intx origin_index_, Vector3 const & u_axis_, Vector3 const & v_axis_,
                     Real radius_)
    {
      clear();

      SurfaceSample * origin_sample = const_cast<SurfaceSample *>(&sample_graph.getSamples()[(size_t)origin_index_]);
      origin = origin_sample->getPosition();
      u_axis = u_axis_.normalized();  // renormalize to be safe
      v_axis = v_axis_.normalized();
      tangent_plane = Plane3::fromPointAndNormal(origin, u_axis.cross(v_axis));
      radius = radius_;
      blend_bandwidth_squared = (options.blendUpwind() ? Math::square(3 * sample_graph.getAverageSeparation()) : 0);

      geodesics.dijkstraWithCallback(const_cast<SampleGraph &>(sample_graph), origin_sample, std::ref(*this), radius);
    }

    Vector2 getParameters(intx sample_index, bool & has_parameters) const
    {
      ParameterMap::const_iterator existing = params.find(sample_index);
      if (existing != params.end())
      {
        has_parameters = true;
        return existing->second;
      }
      else
      {
        has_parameters = false;
        return Vector2::Zero();
      }
    }

    ParameterMap const & getParameterMap() const
    {
      return params;
    }

    CoordinateFrame3 getTangentFrame() const
    {
      return CoordinateFrame3::_fromAffine(AffineTransform3(u_axis, v_axis, u_axis.cross(v_axis), origin));
    }

    Real getRadius() const
    {
      return radius;
    }

    void clear()
    {
      params.clear();
      param_data.clear();
    }

    // This class also acts as the callback during Dijkstra search. VertexHandle is a pointer to a sample.
    bool operator()(Geodesics::VertexHandle vertex, double distance, bool has_pred, Geodesics::VertexHandle pred)
    {
      ParamData curr_data;

      if (has_pred)
      {
        ParamDataMap::const_iterator existing_pred = param_data.find(pred);
        alwaysAssertM(existing_pred != param_data.end(), "SampleGraph: No parametrization data associated with predecessor");
        ParamData const & pred_data = existing_pred->second;

        Vector3 p = vertex->getPosition();
        Vector3 sum_p = pred_data.uv_transform * p;
        Real sum_weights = 1;

        if (options.blendUpwind())
        {
          // Blend in the positions of nearby upwind points (the predecessor's visited neighbors) to make the estimate a bit
          // more robust
          SurfaceSample::NeighborSet const & pred_nbrs = pred->getNeighbors();
          for (int i = 0; i < pred_nbrs.size(); ++i)
          {
            SurfaceSample::Neighbor const & pred_nbr = pred_nbrs[i];
            ParamDataMap::const_iterator existing_nbr = param_data.find(pred_nbr.getSample());
            if (existing_nbr != param_data.end())  // already assigned parameters
            {
              ParamData const & pred_nbr_data = existing_nbr->second;

              Vector3 offset = pred_nbr.getSample()->getPosition() - p;
              Real weight = kernelFastGaussianSqDistUnscaled(offset.squaredNorm(), blend_bandwidth_squared);

              sum_p += weight * (pred_nbr_data.uv_transform * p);
              sum_weights += weight;
            }
          }
        }

        // Now do the DEM unwinding based on the averaged position
        unwind(tangent_plane, sum_p / sum_weights, &pred_data, curr_data);
        Vector3 offset = curr_data.proj_p - origin;
        curr_data.uv = Vector2(offset.dot(u_axis), offset.dot(v_axis));
        curr_data.pred_data = &pred_data;  // hopefully the map class guarantees this will be preserved

        params[vertex->getIndex()] = (options.normalize() ? curr_data.uv / radius : curr_data.uv);

        // THEA_CONSOLE << "Assigned params " << curr_data.uv;
      }
      else
      {
        unwind(tangent_plane, vertex->getPosition(), nullptr, curr_data);
        Vector3 offset = curr_data.proj_p - origin;
        curr_data.uv = Vector2(offset.dot(u_axis), offset.dot(v_axis));
        curr_data.pred_data = nullptr;

        params[vertex->getIndex()] = (options.normalize() ? curr_data.uv / radius : curr_data.uv);

        // THEA_CONSOLE << "Assigned params " << curr_data.uv << " (no pred)";
      }

      param_data[vertex] = curr_data;

      return false;
    }

  private:
    Real kernelFastGaussianSqDistUnscaled(Real squared_dist, Real squared_bandwidth)
    {
      return Math::fastMinusExp((float)(squared_dist / squared_bandwidth));
    }

    // Parametrize the point as an increment from its predecessor. Sets curr.proj_p and curr.uv_transform.
    void
    unwind(Plane3 const & tangent_plane, Vector3 pos_in_pred_frame, ParamData const * pred, ParamData & curr)
    {
      if (pred)
      {
        ParamData const * pred_pred = pred->pred_data;
        if (pred_pred)
        {
          Vector3 edge = pos_in_pred_frame - pred->proj_p;
          Vector3 prev_edge = pred->proj_p - pred_pred->proj_p;

          Vector3 right = prev_edge.cross(tangent_plane.getNormal()).normalized();
          Vector3 v = edge - (edge.dot(right) * right);
          if (v.squaredNorm() > 1.0e-10f)  // TODO: Is this a generally applicable epsilon?
          {
            AffineTransform3 pred_inc = AffineTransform3::translation(pred->proj_p)
                                      * AffineTransform3::rotationArc(v.normalized(), prev_edge.normalized(), false)
                                      * AffineTransform3::translation(-pred->proj_p);
            curr.proj_p = pred_inc * pos_in_pred_frame;
            curr.uv_transform = pred_inc * pred->uv_transform;
          }
          else
          {
            curr.proj_p = tangent_plane.closestPoint(pos_in_pred_frame);
            curr.uv_transform = AffineTransform3::translation(curr.proj_p - pos_in_pred_frame) * pred->uv_transform;
          }
        }
        else  // the predecessor is the source
        {
          Vector3 edge = pos_in_pred_frame - pred->proj_p;
          Vector3 tp = tangent_plane.closestPoint(pos_in_pred_frame);
          Vector3 v = tp - pred->proj_p;
          if (v.squaredNorm() > 1.0e-10f)  // TODO: Is this a generally applicable epsilon?
          {
            AffineTransform3 pred_inc = AffineTransform3::translation(pred->proj_p)
                                      * AffineTransform3::rotationArc(edge.normalized(), v.normalized(), false)
                                      * AffineTransform3::translation(-pred->proj_p);
            curr.proj_p = pred_inc * pos_in_pred_frame;
            curr.uv_transform = pred_inc * pred->uv_transform;
          }
          else  // hope this never happens, because we're going to squash the edge to zero length
          {
            curr.proj_p = tp;
            curr.uv_transform = AffineTransform3::translation(curr.proj_p - pos_in_pred_frame) * pred->uv_transform;
          }
        }
      }
      else  // this is the source
      {
        // Map to origin
        Vector3 tp = tangent_plane.closestPoint(pos_in_pred_frame);  // just in case, but it should be on the plane anyway
        curr.proj_p = tp;
        curr.uv_transform = AffineTransform3::translation(tp - pos_in_pred_frame);
      }
    }

    Options options;
    ParameterMap params;
    ParamDataMap param_data;
    Geodesics geodesics;
    Plane3 tangent_plane;
    Vector3 origin;
    Vector3 u_axis, v_axis;
    Real radius;
    Real blend_bandwidth_squared;

}; // class Impl

} // namespace DiscreteExponentialMapInternal

DiscreteExponentialMap::DiscreteExponentialMap(Options const & options_)
: impl(new DiscreteExponentialMapInternal::Impl(options_))
{}

DiscreteExponentialMap::~DiscreteExponentialMap()
{
  delete impl;
}

void
DiscreteExponentialMap::parametrize(SampleGraph const & sample_graph, intx origin_index, Vector3 const & u_axis,
                                    Vector3 const & v_axis, Real radius)
{
  impl->parametrize(sample_graph, origin_index, u_axis, v_axis, radius);
}

Vector2
DiscreteExponentialMap::getParameters(intx sample_index, bool & has_parameters) const
{
  return impl->getParameters(sample_index, has_parameters);
}

DiscreteExponentialMap::ParameterMap const &
DiscreteExponentialMap::getParameterMap() const
{
  return impl->getParameterMap();
}

CoordinateFrame3
DiscreteExponentialMap::getTangentFrame() const
{
  return impl->getTangentFrame();
}

Real
DiscreteExponentialMap::getRadius() const
{
  return impl->getRadius();
}

void
DiscreteExponentialMap::clear()
{
  impl->clear();
}

} // namespace Algorithms
} // namespace Thea
