//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2015, Siddhartha Chaudhuri
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holders nor the names of contributors
// to this software may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//============================================================================

#include "Mesh.hpp"
#include "../../Algorithms/CentroidN.hpp"
#include "../../Math.hpp"

namespace Browse3D {

Mesh::IndexVertexMap Mesh::index_to_vertex;
Mesh::IndexFaceMap Mesh::index_to_face;

MeshGroup *
Mesh::getAncestor(long generations) const
{
  if (generations < 1)
  {
    THEA_ERROR << getName() << ": Ancestor generation gap must be >= 1";
    return NULL;
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

  TheaArray<Real> all_dists;
  for (VertexConstIterator vi = verticesBegin(); vi != verticesEnd(); ++vi)
    all_dists.push_back((vi->getPosition() - centroid).length());

  Real dmax = *std::max_element(all_dists.begin(), all_dists.end());
  if (dmax <= 1e-20)
  {
    THEA_CONSOLE << getName() << ": Coincident vertices, zero feature vector";
    return;
  }

  double bin_size = dmax / NUM_BINS;
  for (size_t i = 0; i < all_dists.size(); ++i)
  {
    int bin = Math::clamp((int)std::floor(all_dists[i] / bin_size), 0, (int)NUM_BINS - 1);
    double bin_max = (bin + 1) * bin_size;
    features[(size_t)bin] += (all_dists[i] / bin_max);  // make it more discriminative by not binning 1
  }

  features.back() = dmax;  // keep track of the overall scale of the mesh

  // THEA_CONSOLE << getName() << ": Features: [" << seqStr(features.begin(), features.end()) << ']';
}

namespace MeshInternal {

bool
areSimilarFeatureVectors(TheaArray<double> const & f0, long nv0, TheaArray<double> const & f1, long nv1)
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
countVertices(MeshGroup const & mg, long & num_vertices)
{
  for (MeshGroup::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
    num_vertices += (*mi)->numVertices();

  for (MeshGroup::GroupConstIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
    countVertices(**ci, num_vertices);
}

// A rather hacky way of forming a joint descriptor, suitable only for exact matches like we want here
void
accumGroupFeatures(MeshGroup const & mg, TheaArray<double> & features)
{
  for (MeshGroup::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
  {
    TheaArray<double> const & mf = (*mi)->getFeatures();
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

  long nv0 = 0, nv1 = 0;
  countVertices(lhs, nv0);
  countVertices(rhs, nv1);

  if (nv0 != nv1)
    return false;

  TheaArray<double> f0, f1;
  accumGroupFeatures(lhs, f0);
  accumGroupFeatures(rhs, f1);

  return areSimilarFeatureVectors(f0, nv0, f1, nv1);
}

bool
isSimilarTo(Mesh const & lhs, MeshGroup const & rhs)
{
  using namespace MeshInternal;

  long nv0 = lhs.numVertices();
  long nv1 = 0;
  countVertices(rhs, nv1);

  if (nv0 != nv1)
    return false;

  TheaArray<double> const & f0 = lhs.getFeatures();
  TheaArray<double> f1;
  accumGroupFeatures(rhs, f1);

  return areSimilarFeatureVectors(f0, nv0, f1, nv1);
}

bool
isSimilarTo(MeshGroup const & lhs, Mesh const & rhs)
{
  return isSimilarTo(rhs, lhs);
}

} // namespace Browse3D
