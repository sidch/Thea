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
// First version: 2015
//
//============================================================================

#include "Mesh.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Math.hpp"

namespace Browse3D {

Mesh::IndexVertexMap Mesh::index_to_vertex;
Mesh::IndexFaceMap Mesh::index_to_face;

MeshGroup *
Mesh::getAncestor(intx generations) const
{
  if (generations < 1)
  {
    THEA_ERROR << getName() << ": Ancestor generation gap must be >= 1";
    return nullptr;
  }

  MeshGroup * anc = parent;
  if (!anc)
      return anc;

  while (--generations > 0 && anc->getParent())
    anc = anc->getParent();

  return anc;
}

bool
Mesh::hasAncestor(MeshGroup const * anc) const
{
  MeshGroup const * a = parent;
  while (a)
  {
    if (anc == a)
      return true;

    a = a->getParent();
  }

  return false;
}

void
Mesh::updateFeatures() const
{
  if (valid_features)
    return;

  static size_t const NUM_BINS = 64;

  features.resize(NUM_BINS + 1);
  std::fill(features.begin(), features.end(), 0.0);
  valid_features = true;

  if (numVertices() <= 0)
  {
    THEA_CONSOLE << getName() << ": No vertices, zero feature vector";
    return;
  }

  Vector3 centroid = Algorithms::CentroidN<Vertex, 3>::compute(verticesBegin(), verticesEnd());

  Array<Real> all_dists;
  for (VertexConstIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
    all_dists.push_back((vi->getPosition() - centroid).norm());

  Real dmax = *std::max_element(all_dists.begin(), all_dists.end());
  if (dmax <= 1e-20)
  {
    THEA_CONSOLE << getName() << ": Coincident vertices, zero feature vector";
    return;
  }

  double bin_size = dmax / NUM_BINS;
  for (size_t i = 0; i < all_dists.size(); ++i)
  {
    double d = all_dists[i] / bin_size;
    int bin = Math::clamp((int)std::floor(d), 0, (int)NUM_BINS - 1);

    // The middle of the bin is at bin + 0.5. We split the contribution of the value to both this bin and its neighbor.
    double d_curr = 1.0, d_prev = 0.0, d_next = 0.0;
    if (d < bin + 0.5 && bin > 0)
    {
      double w = (bin + 0.5) - d;
      d_curr = 1 - w;
      d_prev = w;
    }
    else if (d > bin + 0.5 && bin < (int)NUM_BINS - 1)
    {
      double w = d - (bin + 0.5);
      d_curr = 1 - w;
      d_next = w;
    }

    features[(size_t)bin] += d_curr;
    if (bin > 0)                 { features[(size_t)bin - 1] += d_prev; }
    if (bin < (int)NUM_BINS - 1) { features[(size_t)bin + 1] += d_next; }
  }

  features.back() = dmax;  // keep track of the overall scale of the mesh

  THEA_DEBUG << getName() << ": Features: [" << stringJoin(features, ", ") << ']';
}

namespace MeshInternal {

bool
areSimilarFeatureVectors(Array<double> const & f0, intx nv0, Array<double> const & f1, intx nv1)
{
  if (nv0 != nv1)
    return false;

  alwaysAssertM(f0.size() == f1.size(), "Feature vectors have different sizes");

  // Compare histograms
  double diff = 0;
  for (size_t i = 0; i + 1 < f0.size(); ++i)
    diff += std::fabs(f0[i] - f1[i]);

  static double const HIST_THRESHOLD = 1e-2;
  if (diff > HIST_THRESHOLD * nv0)
    return false;

  // Compare scales
  static double const SCALE_THRESHOLD = 1e-2;
  if (!Math::fuzzyEq(f0.back(), f1.back(), SCALE_THRESHOLD * (std::fabs(f0.back()) + 1)))
    return false;

  return true;
}

void
countVertices(MeshGroup const & mg, intx & num_vertices)
{
  for (MeshGroup::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
    num_vertices += (*mi)->numVertices();

  for (MeshGroup::GroupConstIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
    countVertices(**ci, num_vertices);
}

// A rather hacky way of forming a joint descriptor, suitable only for exact matches like we want here
void
accumGroupFeatures(MeshGroup const & mg, Array<double> & features)
{
  for (MeshGroup::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
  {
    Array<double> const & mf = (*mi)->getFeatures();
    if (features.empty())
      features.resize(mf.size());

    alwaysAssertM(mf.size() == features.size(), "Feature vectors have different sizes");

    for (size_t j = 0; j < features.size(); ++j)
      features[j] += mf[j];
  }

  for (MeshGroup::GroupConstIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
    accumGroupFeatures(**ci, features);
}

} // namespace MeshInternal

bool
isSimilarTo(Mesh const & lhs, Mesh const & rhs)
{
  if (&lhs == &rhs)
    return true;

  return MeshInternal::areSimilarFeatureVectors(lhs.getFeatures(), lhs.numVertices(), rhs.getFeatures(), rhs.numVertices());
}

bool
isSimilarTo(MeshGroup const & lhs, MeshGroup const & rhs)
{
  using namespace MeshInternal;

  if (&lhs == &rhs)
    return true;

  intx nv0 = 0, nv1 = 0;
  countVertices(lhs, nv0);
  countVertices(rhs, nv1);

  if (nv0 != nv1)
    return false;

  Array<double> f0, f1;
  accumGroupFeatures(lhs, f0);
  accumGroupFeatures(rhs, f1);

  return areSimilarFeatureVectors(f0, nv0, f1, nv1);
}

bool
isSimilarTo(Mesh const & lhs, MeshGroup const & rhs)
{
  using namespace MeshInternal;

  intx nv0 = lhs.numVertices();
  intx nv1 = 0;
  countVertices(rhs, nv1);

  if (nv0 != nv1)
    return false;

  Array<double> const & f0 = lhs.getFeatures();
  Array<double> f1;
  accumGroupFeatures(rhs, f1);

  return areSimilarFeatureVectors(f0, nv0, f1, nv1);
}

bool
isSimilarTo(MeshGroup const & lhs, Mesh const & rhs)
{
  return isSimilarTo(rhs, lhs);
}

} // namespace Browse3D
