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
// First version: 2011
//
//============================================================================

#ifndef __Browse3D_MainWindow_hpp__
#define __Browse3D_MainWindow_hpp__

#include "Common.hpp"
#include <wx/frame.h>

class wxCheckBox;
class wxListBox;
class wxNotebook;
class wxSplitterWindow;
class wxTextCtrl;
class wxUpdateUIEvent;

namespace Browse3D {

class Model;
class ModelDisplay;

/** Holds MainWindow UI elements. */
struct MainWindowUI
{
  wxSplitterWindow * main_splitter;

  ModelDisplay * model_display;
  wxNotebook * toolbox;

  wxListBox * points_table;
  wxTextCtrl * point_label;
  wxCheckBox * point_snap_to_vertex;

  wxListBox * segments_table;
  wxTextCtrl * segment_label;

  /** Default constructor. Sets all pointers to null. */
  MainWindowUI();

}; // struct MainWindow

/** The main application window. */
class MainWindow : public wxFrame
{
  private:
    typedef wxFrame BaseType;

  public:
    /** Custom item IDs. */
    enum ItemID
    {
      ID_VIEW_SHADED,
      ID_VIEW_WIREFRAME,
      ID_VIEW_SHADED_WIREFRAME,
      ID_VIEW_TWO_SIDED,
      ID_VIEW_FLAT_SHADED,
      ID_VIEW_FIT,
      ID_VIEW_PRINT_CAMERA,
      ID_GO_PREV,
      ID_GO_NEXT,
      ID_GO_PREV_FEATURES,
      ID_GO_NEXT_FEATURES,
      ID_TOOLS_SCREENSHOT,
      ID_TOOLS_TOOLBOX,

      ID_SEGMENT_LABEL,
      ID_SEGMENT_ADD,
      ID_SEGMENT_REMOVE,
      ID_SEGMENT_EXPAND,
      ID_SEGMENT_CONTRACT,

      ID_POINT_LABEL,
      ID_POINT_ADD,
      ID_POINT_REMOVE,

    }; // enum EventID

    /** Constructor. */
    explicit MainWindow(wxWindow * parent = nullptr);

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
    intx numOverlays() const { return (intx)overlays.size(); }

    /** Get the set of overlay models. */
    Model const * const * getOverlays() const { return overlays.empty() ? nullptr : &overlays[0]; }

    /** Check if point-picking is on. */
    bool pickPoints() const;

    /** Check if segment-selection is on. */
    bool pickSegments() const;

    /** [wxWidgets] Set the window title. */
    void SetTitle(wxString const & title);

    //=========================================================================================================================
    // GUI callbacks
    //=========================================================================================================================

    /** Called when program exits. */
    void OnExit(wxEvent & event = DUMMY_EVENT);

    /** Set the window title. */
    void setTitle(wxEvent & event = DUMMY_EVENT);

    /** Select and load a model. */
    void selectAndLoadModel(wxEvent & event = DUMMY_EVENT);

    /** Load the previous model in the directory. */
    void loadPreviousModel(wxEvent & event = DUMMY_EVENT);

    /** Load the next model in the directory. */
    void loadNextModel(wxEvent & event = DUMMY_EVENT);

    /** Load the previous set of features in the features directory. */
    void loadPreviousFeatures(wxEvent & event = DUMMY_EVENT);

    /** Load the next set of features in the features directory. */
    void loadNextFeatures(wxEvent & event = DUMMY_EVENT);

    /** Add the currently picked point to the set of samples. */
    void addPickedSample(wxEvent & event = DUMMY_EVENT);

    /** Select the sample indicated by the table selection. */
    void selectSample(wxEvent & event = DUMMY_EVENT);

    /** Remove the selected sample. */
    void removeSelectedSample(wxEvent & event = DUMMY_EVENT);

    /** Sync the displayed list of samples with the model. */
    void syncSamples(wxEvent & event = DUMMY_EVENT);

    /** Add the currently picked segment to the set of segments. */
    void addPickedSegment(wxEvent & event = DUMMY_EVENT);

    /** Expand the picked segment by one level in the hierarchy. */
    void expandPickedSegment(wxEvent & event = DUMMY_EVENT);

    /** Contract the picked segment by one level in the hierarchy. */
    void contractPickedSegment(wxEvent & event = DUMMY_EVENT);

    /** Select the segment indicated by the table selection. */
    void selectSegment(wxEvent & event = DUMMY_EVENT);

    /** Remove the selected segment. */
    void removeSelectedSegment(wxEvent & event = DUMMY_EVENT);

    /** Sync the displayed list of segments with the model. */
    void syncSegments(wxEvent & event = DUMMY_EVENT);

    /** Show/hide the toolbox. */
    void setShowToolbox(wxCommandEvent & event);

    /** Called when the toolbox panel changes. */
    void toolboxPanelChanged(wxEvent & event = DUMMY_EVENT);

    /** Synchronize states of menu and toolbar buttons etc. */
    void updateUI(wxUpdateUIEvent & event);

    /** Refresh the model display. */
    void refreshDisplay(wxEvent & event = DUMMY_EVENT);

  private:
    /** Get rid of all overlay models. */
    void clearOverlays();

    /** Show or hide the toolbox. */
    void setShowToolbox(bool value);

    // Models
    Model * model;
    Array<Model *> overlays;

    // Widgets
    MainWindowUI ui;

}; // class MainWindow

} // namespace Browse3D

#endif
