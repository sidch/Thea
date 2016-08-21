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

#ifndef __Browse3D_MainWindow_hpp__
#define __Browse3D_MainWindow_hpp__

#include "Common.hpp"
#include <wx/frame.h>

class wxCheckBox
class wxListBox;
class wxNotebook;
class wxTextCtrl;

namespace Browse3D {

class Model;
class ModelDisplay;

/** The main application window. */
class MainWindow : public wxFrame
{
    typedef wxFrame BaseType;

  public:
    /** Constructor. */
    explicit MainWindow(wxWindow * parent = NULL);

    /** Destructor. */
    ~MainWindow();

    /** Initialize the main window. Should be called at application startup. */
    void init();

    /** Get the display widget showing the model. */
    ModelDisplay * getRenderDisplay();

    /** Get the active model. */
    Model const * getModel() const { return model; }

    /** Get the active model. */
    Model * getModel() { return model; }

    /** Get the number of overlay models. */
    long numOverlays() const { return (long)overlays.size(); }

    /** Get the set of overlay models. */
    Model const * const * getOverlays() const { return overlays.empty() ? NULL : &overlays[0]; }

    /** Check if point-picking is on. */
    bool pickPoints() const;

    /** Check if segment-selection is on. */
    bool pickSegments() const;

    /** Select and load a model. */
    void selectAndLoadModel();

    /** Load the previous model in the directory. */
    void loadPreviousModel();

    /** Load the next model in the directory. */
    void loadNextModel();

    /** Load the previous set of features in the features directory. */
    void loadPreviousFeatures();

    /** Load the next set of features in the features directory. */
    void loadNextFeatures();

    /** Add the currently picked point to the set of samples. */
    void addPickedSample();

    /** Select the sample indicated by the table selection. */
    void selectSample();

    /** Remove the selected sample. */
    void removeSelectedSample();

    /** Sync the displayed list of samples with the model. */
    void syncSamples();

    /** Add the currently picked segment to the set of segments. */
    void addPickedSegment();

    /** Expand the picked segment by one level in the hierarchy. */
    void expandPickedSegment();

    /** Contract the picked segment by one level in the hierarchy. */
    void contractPickedSegment();

    /** Select the segment indicated by the table selection. */
    void selectSegment();

    /** Remove the selected segment. */
    void removeSelectedSegment();

    /** Sync the displayed list of segments with the model. */
    void syncSegments();

    /** Show/hide the toolbox. */
    void setShowToolbox(bool value);

    /** Turn point-picking on/off. */
    void setPickPoints(bool value);

    /** Turn segment-picking on/off. */
    void setPickSegments(bool value);

    //=========================================================================================================================
    // GUI callbacks etc
    //=========================================================================================================================

    /** Set the window title. */
    void SetTitle(wxString const & title);

    /** Called when the window is closed. */
    void OnExit(wxCommandEvent & event);

  private:
    /** Get rid of all overlay models. */
    void clearOverlays();

    // Models
    Model * model;
    TheaArray<Model *> overlays;

    // Widgets
    ModelDisplay * model_display;
    wxNotebook * toolbox;

    wxListBox * points_table;
    wxTextCtrl * point_label;
    wxCheckBox * pick_points_snap_to_vertex;

    wxListBox * segments_table;
    wxTextCtrl * segment_label;

}; // class MainWindow

} // namespace Browse3D

#endif
