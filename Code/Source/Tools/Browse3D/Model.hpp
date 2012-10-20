//============================================================================
//
// This file is part of the Browse3D project.
//
// This software is covered by the following BSD license, except for portions
// derived from other works which are covered by their respective licenses.
// For full licensing information including reproduction of these external
// licenses, see the file LICENSE.txt provided in the documentation.
//
// Copyright (C) 2011, Siddhartha Chaudhuri/Stanford University
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

#ifndef __Browse3D_Model_hpp__
#define __Browse3D_Model_hpp__

#include "Common.hpp"
#include "GraphicsWidget.hpp"
#include "MeshFwd.hpp"
#include "../../Algorithms/KDTree3.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Algorithms/RayQueryStructure3.hpp"
#include <QObject>

class QMouseEvent;

namespace Browse3D {

class PointCloud;
typedef shared_ptr<PointCloud> PointCloudPtr;

/**
 * An indexed vertex of a mesh suitable for kd-tree storage, saves memory over DisplayMeshVertex and allows forward-declaration
 * of the mesh.
 */
struct IndexedMeshVertex
{
  IndexedMeshVertex() {}
  IndexedMeshVertex(Mesh * mesh_, long index_, Vector3 const & position_) : mesh(mesh_), index(index_), position(position_) {}

  Mesh * mesh;
  long index;
  Vector3 position;

}; // struct IndexedMeshVertex

} // namespace Browse3D

namespace Thea {
namespace Algorithms {

// Specify that a mesh vertex is a logical 3D point. */
template <>
class IsPointN<Browse3D::IndexedMeshVertex, 3>
{
  public:
    static bool const value = true;
};

// Map a mesh vertex to its 3D position. */
template <>
struct PointTraitsN<Browse3D::IndexedMeshVertex, 3>
{
  static Vector3 const & getPosition(Browse3D::IndexedMeshVertex const & t) { return t.position; }
};

} // namespace Algorithms
} // namespace Thea

namespace Browse3D {

/** The model manipulated by the user. */
class Model : public QObject, public GraphicsWidget
{
    Q_OBJECT

  public:
    typedef Thea::Algorithms::MeshKDTree<Mesh> KDTree;  ///< A kd-tree on mesh triangles.
    typedef Algorithms::RayStructureIntersection3 RayStructureIntersection3;  /**< Intersection of a ray with an acceleration
                                                                                   structure. */
    typedef Thea::Algorithms::KDTree3<IndexedMeshVertex> VertexKDTree;  ///< A kd-tree on mesh vertices.

    /** A sample point on the surface. */
    struct Sample
    {
      Sample() {}

      Sample(Mesh * mesh_, long face_index_, Vector3 const & position_, QString const & type_ = "", QString const & label_ = "")
      : mesh(mesh_), face_index(face_index_), position(position_), type(type_), label(label_)
      {}

      Mesh * mesh;
      long face_index;
      Vector3 position;
      QString type;
      QString label;

    }; // struct Sample

    //========================================================================================================================
    // Basic functions
    //========================================================================================================================

    /** Constructor. */
    Model(QString const & initial_mesh = "");

    /** Destructor. */
    ~Model();

    /** Get the name of the model. */
    QString getName() const;

    /** Get the filename of the currently loaded model. */
    QString const & getFilename() const { return filename; }

    /** Is the model empty? */
    bool isEmpty() const;

    /** Clear the model and invalidate all associated structures. */
    void clear();

    /** Invalidate all associated structures. */
    void invalidateAll();

    //========================================================================================================================
    // KD-trees on mesh triangles and vertices
    //========================================================================================================================

    /** Invalidate the kd-tree of the model. The kd-tree will be lazily recomputed. */
    void invalidateKDTree();

    /** Invalidate the kd-tree on the vertices of the model. The kd-tree will be lazily recomputed. */
    void invalidateVertexKDTree();

    /** Compute the kd-tree of the model, if the current kd-tree is invalid. */
    void updateKDTree() const;

    /** Compute the kd-tree on the vertices of the model, if the current kd-tree is invalid. */
    void updateVertexKDTree() const;

    /** Get the kd-tree for the model. By default, it will be recomputed if it is currently invalid. */
    KDTree const & getKDTree(bool recompute_if_invalid = true) const;

    /** Get the kd-tree on the vertices of the model. By default, it will be recomputed if it is currently invalid. */
    VertexKDTree const & getVertexKDTree(bool recompute_if_invalid = true) const;

    //========================================================================================================================
    // Ray-shooting and nearest neighbor queries
    //========================================================================================================================

    /** Check if a ray intersects the model, optionally restricting the search to a maximum hit time. */
    bool rayIntersects(Ray3 const & ray, Real max_time = -1) const;

    /**
     * Compute the time after which a ray will intersect the model, optionally restricting the search to a maximum hit time. If
     * the ray does not intersect the model, a negative value is returned.
     */
    Real rayIntersectionTime(Ray3 const & ray, Real max_time = -1) const;

    /**
     * Compute the intersection of a ray with the model, optionally restricting the search to a maximum hit time. If the ray does
     * not intersect the model, an invalid intersection is returned.
     */
    RayStructureIntersection3 rayIntersection(Ray3 const & ray, Real max_time = -1) const;

    /**
     * Get the point on the model closest to a query point.
     *
     * @param query The query point.
     * @param distance_bound Maximum allowed distance to closest point. Pass a negative value to disable this limit.
     * @param min_dist Used to return distance to closest point.
     * @param closest_pt Used to return the closest point itself.
     * @param closest_pt_normal Used to return the surface normal at the closest point.
     * @param accelerate_with_vertices If true, the function first computes an approximate nearest neighbor using the set of
     *   mesh vertices, to make the exact computation faster. This parameter should normally always be true.
     *
     * @return The index of the kd-tree triangle containing the closest point, if one is found within the distance bound, else
     *   a negative number.
     */
    long closestPoint(Vector3 const & query, Real distance_bound = -1, Real * min_dist = NULL, Vector3 * closest_pt = NULL,
                      Vector3 * closest_pt_normal = NULL, bool accelerate_with_vertices = true) const;

    /**
     * Select the nearest point on the model along a ray.
     *
     * @return The distance, in multiples of ray length, to the picked point.
     */
    Real pick(Ray3 const & ray);

    /** Invalidate the last picked point. */
    void invalidatePick();

    //========================================================================================================================
    // Interaction
    //========================================================================================================================

    /** [Qt] Called when a mouse button is pressed. */
    void mousePressEvent(QMouseEvent * event);

    /** [Qt] Called when the mouse is moved. */
    void mouseMoveEvent(QMouseEvent * event);

    /** [Qt] Called when a mouse button is released. */
    void mouseReleaseEvent(QMouseEvent * event);

    //========================================================================================================================
    // Samples
    //========================================================================================================================

    /** Get the number of sample points. */
    long numSamples() const { return (long)samples.size(); }

    /** Get the set of sample points. */
    TheaArray<Sample> const & getSamples() const { return samples; }

    /** Add a sample point. */
    void addSample(Sample const & sample);

    /** Add the currently picked point as a sample. */
    bool addPickedSample(QString const & label, bool snap_to_vertex);

    /** Remove a sample. */
    void removeSample(long index);

    /** Select a particular sample. */
    void selectSample(long index);

    /** Load samples from a file. */
    bool loadSamples(QString const & filename_);

    /** Save samples to a file. */
    bool saveSamples(QString const & filename_) const;

    /** Get the path to the file in which samples are stored. */
    QString getSamplesFilename() const;

    //========================================================================================================================
    // Bounding boxes
    //========================================================================================================================

    AxisAlignedBox3 const & getBounds() const;

    void updateBounds();

    //========================================================================================================================
    // Display
    //========================================================================================================================

    void uploadToGraphicsSystem(Graphics::RenderSystem & render_system);

    void draw(Graphics::RenderSystem & render_system,
              Graphics::RenderOptions const & options = Graphics::RenderOptions::defaults()) const;

  public slots:
    /**
     * Load the model from a disk file.
     *
     * @return True if the model was successfully loaded, else false.
     */
    bool load(QString const & filename_);

    /**
     * Select a file via a file dialog and load the model from it.
     *
     * @return True if the model was successfully loaded, else false.
     */
    bool selectAndLoad();

  signals:
    /**
     * Emitted when the model changes its geometry.
     *
     * @param model Always points to the emitting model.
     */
    void geometryChanged(Model const * model);

    /**
     * Emitted when the filename of the model changes.
     *
     * @param current_filename The current filename of the model.
     */
    void filenameChanged(QString const & current_filename);

    /**
     * Emitted when the model needs to be redrawn.
     *
     * @param model Always points to the emitting model.
     */
    void needsRedraw(Model const * model);

    /**
     * Emitted when a displayed list of samples needs to be synced with the model.
     *
     * @param model Always points to the emitting model.
     */
    void needsSyncSamples(Model const * model);

  private:
    /** Clear the model mesh. */
    void clearMesh();

    /** Clear the point cloud. */
    void clearPoints();

    MeshGroupPtr mesh_group;
    PointCloudPtr point_cloud;
    QString filename;

    AxisAlignedBox3 bounds;

    TheaArray<Sample> samples;
    bool valid_pick;
    Sample picked_sample;
    long selected_sample;

    mutable bool valid_kdtree;
    mutable KDTree * kdtree;

    mutable bool valid_vertex_kdtree;
    mutable VertexKDTree * vertex_kdtree;

}; // class Model

} // namespace Browse3D

#endif
