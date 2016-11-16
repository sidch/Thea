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
#include "Segment.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MeshKDTree.hpp"
#include "../../Algorithms/PointTraitsN.hpp"
#include "../../Algorithms/RayQueryStructureN.hpp"
#include "../../AffineTransform3.hpp"
#include "../../Transformable.hpp"
#include <wx/event.h>

class wxMouseEvent;

namespace Browse3D {

class ModelDisplay;
class PointCloud;
typedef shared_ptr<PointCloud> PointCloudPtr;

} // namespace Browse3D

wxDECLARE_EVENT(EVT_MODEL_PATH_CHANGED,         wxCommandEvent);
wxDECLARE_EVENT(EVT_MODEL_GEOMETRY_CHANGED,     wxCommandEvent);
wxDECLARE_EVENT(EVT_MODEL_NEEDS_REDRAW,         wxCommandEvent);
wxDECLARE_EVENT(EVT_MODEL_NEEDS_SYNC_SAMPLES,   wxCommandEvent);
wxDECLARE_EVENT(EVT_MODEL_NEEDS_SYNC_SEGMENTS,  wxCommandEvent);

namespace Browse3D {

/** The model manipulated by the user. */
class Model : public GraphicsWidget, public Transformable<AffineTransform3>, public wxEvtHandler
{
  private:
    typedef Transformable<AffineTransform3> TransformableBaseT;

  public:
    typedef Thea::Algorithms::MeshKDTree<Mesh> KDTree;  ///< A kd-tree on mesh triangles.
    typedef Algorithms::RayStructureIntersection3 RayStructureIntersection3;  /**< Intersection of a ray with an acceleration
                                                                                   structure. */
    typedef Thea::Algorithms::KDTreeN<MeshVertex const *, 3> VertexKDTree;  ///< A kd-tree on mesh vertices.

    /** A sample point on the surface. */
    struct Sample
    {
      Sample() {}

      Sample(Mesh * mesh_, long face_index_, Vector3 const & position_, std::string const & type_ = "",
             std::string const & label_ = "")
      : mesh(mesh_), face_index(face_index_), position(position_), type(type_), label(label_)
      {}

      Mesh * mesh;
      long face_index;
      Vector3 position;
      std::string type;
      std::string label;

    }; // struct Sample

    //========================================================================================================================
    // Basic functions
    //========================================================================================================================

    /** Constructor. */
    Model(std::string const & initial_mesh = "");

    /** Destructor. */
    ~Model();

    /** Get the name of the model. */
    std::string getName() const;

    /** Get the path of the currently loaded model. */
    std::string const & getPath() const { return path; }

    /** Is the model empty? */
    bool isEmpty() const;

    /** Clear the model and invalidate all associated structures. */
    void clear();

    /** Invalidate all associated structures. */
    void invalidateAll();

    void setTransform(AffineTransform3 const & trans_);

    void clearTransform();

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
     * Compute the intersection of a ray with the model, optionally restricting the search to a maximum hit time. If the ray
     * does not intersect the model, an invalid intersection is returned.
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

    /** [wxWidgets] Called when a mouse button is pressed. */
    void mousePressEvent(wxMouseEvent & event);

    /** [wxWidgets] Called when the mouse is moved. */
    void mouseMoveEvent(wxMouseEvent & event);

    /** [wxWidgets] Called when a mouse button is released. */
    void mouseReleaseEvent(wxMouseEvent & event);

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
    bool addPickedSample(std::string const & label, bool snap_to_vertex);

    /** Remove a sample. */
    void removeSample(long index);

    /** Select a particular sample. */
    void selectSample(long index);

    /** Load samples from a file. */
    bool loadSamples(std::string const & path_);

    /** Save samples to a file. */
    bool saveSamples(std::string const & path_) const;

    /** Get the path to the file in which samples are stored. */
    std::string getSamplesPath() const;

    //========================================================================================================================
    // Segments
    //========================================================================================================================

    /**
     * Flip the selection status of the submesh along a ray.
     *
     * @param ray The ray to use for picking.
     * @param extend_to_similar If true, all meshes similar to the picked mesh have their selection states set to match the
     *   latter.
     *
     * @return The distance, in multiples of ray length, to the picked point.
     */
    Real togglePickMesh(Ray3 const & ray, bool extend_to_similar = false);

    /** Expand/contract the selection by \a offset steps in the mesh group hierarchy. */
    void promotePickedSegment(long offset = 1);

    /** Check if a mesh is currently selected by picking. */
    bool isPicked(Mesh const * mesh) const { return picked_segment.hasMesh(mesh); }

    /** Invalidate the currently picked segment selection. */
    void invalidatePickedSegment() { picked_segment.clear(); }

    /** Get the number of labeled segments. */
    long numSegments() const { return (long)segments.size(); }

    /** Get the set of labeled segments. */
    TheaArray<Segment> const & getSegments() const { return segments; }

    /** Add a segment. */
    void addSegment(Segment const & segment);

    /** Add the currently picked segment to the set of labeled segments. */
    bool addPickedSegment(std::string const & label);

    /** Remove a segment from the list of labeled segments. */
    void removeSegment(long index);

    /** Get the segment containing a given mesh, or null if there is no such segment. */
    Segment const * getSegment(Mesh const * mesh, int * index = NULL) const
    { return const_cast<Model *>(this)->getSegment(mesh, index); }

    /** Get the segment containing a given mesh, or null if there is no such segment. */
    Segment * getSegment(Mesh const * mesh, int * index = NULL);

    /** Select a particular segment. */
    void selectSegment(long index);

    /** Load labeled segments from a file. */
    bool loadSegments(std::string const & path_);

    /** Save labeled segments to a file. */
    bool saveSegments(std::string const & path_) const;

    /** Get the path to the file in which labeled segments are stored. */
    std::string getSegmentsPath() const;

    //========================================================================================================================
    // Features
    //========================================================================================================================

    /** Load features from a file. */
    bool loadFeatures(std::string const & path_);

    /** Get the path of the currently loaded features. */
    std::string const & getFeaturesPath() const { return features_path; }

    /** Check if the model has currently loaded features. */
    bool hasFeatures() const { return has_features; }

    //========================================================================================================================
    // Face labels
    //========================================================================================================================

    /** Load face labels from a file. */
    bool loadFaceLabels(std::string const & path_);

    /** Get the path of the currently loaded face labels. */
    std::string const & getFaceLabelsPath() const { return face_labels_path; }

    /** Check if the model has currently loaded face labels. */
    bool hasFaceLabels() const { return has_face_labels; }

    //========================================================================================================================
    // Bounding boxes
    //========================================================================================================================

    /** Get the bounding box in object coordinates. */
    AxisAlignedBox3 const & getBounds() const;

    /** Get the bounding box in world coordinates. */
    AxisAlignedBox3 getTransformedBounds() const;

    void updateBounds();

    //========================================================================================================================
    // Display
    //========================================================================================================================

    /** Register a display window to receive events from this model. */
    void registerDisplay(ModelDisplay * display);

    /** Deregister a display window from receiving events from this model. */
    void deregisterDisplay(ModelDisplay * display);

    /** Get the default color of the model. */
    ColorRGBA const & getColor() const { return color; }

    /** Set the default color of the model. */
    void setColor(ColorRGBA const & color_) { color = color_; }

    void uploadToGraphicsSystem(Graphics::RenderSystem & render_system);

    void draw(Graphics::RenderSystem & render_system,
              Graphics::RenderOptions const & options = Graphics::RenderOptions::defaults()) const;

    //========================================================================================================================
    // GUI callbacks
    //========================================================================================================================

    /**
     * Load the model from a disk file.
     *
     * @return True if the model was successfully loaded, else false.
     */
    bool load(std::string path_);

    /**
     * Select a file via a file dialog and load the model from it.
     *
     * @return True if the model was successfully loaded, else false.
     */
    bool selectAndLoad();

  private:
    /** Clear the model mesh. */
    void clearMesh();

    /** Clear the point cloud. */
    void clearPoints();

    /** Get the default path to the file in which features are stored. */
    std::string getDefaultFeaturesPath() const;

    /** Get the default path to the file in which the face labels are stored. */
    std::string getDefaultFaceLabelsPath() const;

    /** Draw the mesh group colored by segment. */
    void drawSegmentedMeshGroup(MeshGroupPtr mesh_group, int depth, int & node_index, Graphics::RenderSystem & render_system,
                                Graphics::RenderOptions const & options) const;

    MeshGroupPtr mesh_group;
    PointCloudPtr point_cloud;
    std::string path;

    std::string features_path;
    bool has_features;

    std::string face_labels_path;
    bool has_face_labels;

    ColorRGBA color;
    AxisAlignedBox3 bounds;

    TheaArray<Sample> samples;
    bool valid_pick;
    Sample picked_sample;
    long selected_sample;

    TheaArray<Segment> segments;
    Segment picked_segment;
    long segment_depth_promotion;
    long selected_segment;

    mutable bool valid_kdtree;
    mutable KDTree * kdtree;

    mutable bool valid_vertex_kdtree;
    mutable VertexKDTree * vertex_kdtree;

}; // class Model

} // namespace Browse3D

#endif
