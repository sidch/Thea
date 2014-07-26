//============================================================================
//
// This file is part of the Thea project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2013, Siddhartha Chaudhuri/Princeton University
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

#ifndef __Thea_Algorithms_MeshSampler_hpp__
#define __Thea_Algorithms_MeshSampler_hpp__

#include "../Common.hpp"
#include "../Graphics/MeshGroup.hpp"
#include "MeshTriangles.hpp"
#include "../Random.hpp"
#include "../Vector3.hpp"

namespace Thea {
namespace Algorithms {

/** Sample points from a mesh. */
template <typename MeshT>
class MeshSampler
{
  public:
    typedef MeshT                         Mesh;       ///< The mesh class.
    typedef MeshTriangles<Mesh>           Triangles;  ///< A set of mesh triangles.
    typedef typename Triangles::Triangle  Triangle;   ///< A triangle belonging to a face of the mesh.

    /**
     * Constructs the object to compute sample points on a given mesh. The mesh must persist as long as this object does.
     * Initializes internal data structures that do not need to be recomputed for successive calls to functions that generate
     * samples on this mesh.
     */
    MeshSampler(Mesh const & mesh) : num_external_tris(0), external_tris(NULL)
    {
      tris.add(const_cast<Mesh &>(mesh));
    }

    /**
     * Constructs the object to compute sample points on a given mesh group. The mesh group must persist as long as this object
     * does. Initializes internal data structures that do not need to be recomputed for successive calls to functions that
     * generate samples on this mesh.
     */
    MeshSampler(Graphics::MeshGroup<Mesh> const & mesh_group) : num_external_tris(0), external_tris(NULL)
    {
      tris.add(const_cast<Graphics::MeshGroup<Mesh> &>(mesh_group));
    }

    /**
     * Constructs the object to compute sample points on a mesh with a precomputed triangulation. The triangles must persist as
     * long as this object does.
     */
    MeshSampler(long num_tris, Triangle const * tris) : num_external_tris(num_tris), external_tris(tris)
    {
      alwaysAssertM(num_tris <= 0 || tris, "MeshSampler: Triangle list cannot be null");
    }

    /** Destructor. */
    ~MeshSampler() {}

    /**
     * Samples points from the mesh evenly by area.
     *
     * @param approx_num_samples The approximate number of samples to compute. This need <em>not</em> be exactly the number
     *   actually computed (and returned by the function).
     * @param positions Used to return the sample positions.
     * @param face_normals If not null, used to return the normals of the parent faces of the samples.
     * @param triangles If not null, used to return the mesh triangles from which the samples were selected.
     *
     * @return The number of samples computed.
     */
    long sampleEvenlyByArea(long approx_num_samples,
                            TheaArray<Vector3> & positions,
                            TheaArray<Vector3> * face_normals = NULL,
                            TheaArray<Triangle const *> * triangles = NULL) const
    {
      positions.clear();
      if (face_normals) face_normals->clear();
      if (triangles) triangles->clear();

      Triangle const * tarray;
      long tcount;
      if (num_external_tris > 0)
      {
        tarray = external_tris;
        tcount = num_external_tris;
      }
      else
      {
        if (tris.isEmpty())
          return 0;
        else
        {
          tarray = &tris.getTriangles()[0];
          tcount = tris.numTriangles();
        }
      }

      double total_area = 0;
      for (long i = 0; i < tcount; ++i)
        total_area += tarray[i].getArea();

      if (total_area <= 0)
        return 0;

      for (long i = 0; i < tcount; ++i)
      {
        double frac_samples = (approx_num_samples * tarray[i].getArea()) / total_area;
        while (frac_samples > 0)
        {
          // Last (possible) sample?
          if (frac_samples < 1 && Random::common().uniform01() > frac_samples)
            break;

          positions.push_back(tarray[i].randomPoint());
          if (face_normals) face_normals->push_back(tarray[i].getNormal());
          if (triangles) triangles->push_back(&tarray[i]);

          frac_samples -= 1.0;
        }
      }

      return (long)positions.size();
    }

  private:
    Triangles tris;  ///< An internally computed triangulation of the mesh.
    long num_external_tris;  ///< Number of externally supplied, precomputed mesh triangles.
    Triangle const * external_tris;  ///< Array of externally supplied, precomputed mesh triangles.

}; // class MeshSampler

} // namespace Algorithms
} // namespace Thea

#endif
