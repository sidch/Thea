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

#include "MainWindow.hpp"
#include "App.hpp"
#include "Model.hpp"
#include "ModelDisplay.hpp"
#include "Util.hpp"
#include "../../Application.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include <wx/accel.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/listbox.h>
#include <wx/menu.h>
#include <wx/notebook.h>
#include <wx/panel.h>
#include <wx/sizer.h>
#include <wx/splitter.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>

namespace Browse3D {

// Definition of object declared in common header
wxCommandEvent DUMMY_EVENT;

static int const SEGMENTS_TAB_INDEX  =     0;
static int const POINTS_TAB_INDEX    =     1;
static Real const MIN_SPLIT_SIZE     =   300;

MainWindowUI::MainWindowUI()
: model_display(NULL),
  toolbox(NULL),
  points_table(NULL),
  point_label(NULL),
  point_snap_to_vertex(NULL),
  segments_table(NULL),
  segment_label(NULL)
{}

MainWindow::MainWindow(wxWindow * parent)
: BaseType(parent, wxID_ANY, "Browse3D", wxDefaultPosition, wxSize(800, 600)),
  model(NULL)
{
  init();
}

void
MainWindow::init()
{
  //==========================================================================================================================
  // Menu
  //==========================================================================================================================

  // Set up the main menu
  wxMenuBar * menubar = new wxMenuBar();

  // File Menu
  wxMenu * file_menu = new wxMenu();
  file_menu->Append(wxID_OPEN,    "&Open");
  file_menu->Append(wxID_SAVEAS,  "&Save");
  file_menu->AppendSeparator();
  file_menu->Append(wxID_EXIT,    "&Quit");
  menubar->Append(file_menu, "&File");

  // View menu
  wxMenu * view_menu = new wxMenu();
  wxMenu * rendering_menu = new wxMenu();
    rendering_menu->AppendRadioItem(ID_VIEW_SHADED,            "&Shaded");
    rendering_menu->AppendRadioItem(ID_VIEW_WIREFRAME,         "&Wireframe");
    rendering_menu->AppendRadioItem(ID_VIEW_SHADED_WIREFRAME,  "S&haded + wireframe");
    rendering_menu->AppendSeparator();
    rendering_menu->AppendCheckItem(ID_VIEW_TWO_SIDED,         "&Two-sided lighting");
    rendering_menu->AppendCheckItem(ID_VIEW_FLAT_SHADING,      "&Flat shading");
  view_menu->AppendSubMenu(rendering_menu,  "&Rendering");
  view_menu->Append(ID_VIEW_FIT,            "&Fit view to model");
  menubar->Append(view_menu, "&View");

  // Go menu
  wxMenu * go_menu = new wxMenu();
  go_menu->Append(ID_GO_PREV,           "&Previous model");
  go_menu->Append(ID_GO_NEXT,           "&Next model");
  go_menu->AppendSeparator();
  go_menu->Append(ID_GO_PREV_FEATURES,  "Previous features");
  go_menu->Append(ID_GO_NEXT_FEATURES,  "Next features");
  menubar->Append(go_menu, "&Go");

  // Tools menu
  wxMenu * tools_menu = new wxMenu();
  tools_menu->Append(ID_TOOLS_SCREENSHOT,        "&Save screenshot");
  tools_menu->AppendCheckItem(ID_TOOLS_TOOLBOX,  "&Toolbox");
  menubar->Append(tools_menu, "&Tools");

  // About menu
  wxMenu * help_menu = new wxMenu();
  help_menu->Append(wxID_ABOUT,  "&About");
  menubar->Append(help_menu, "&Help");

  SetMenuBar(menubar);

  //==========================================================================================================================
  // Toolbar
  //==========================================================================================================================

// #define SHOW_TOOLBAR
#ifdef SHOW_TOOLBAR
  wxToolBar * toolbar = CreateToolBar(wxTB_HORIZONTAL | wxTB_TEXT | wxTB_NOICONS, wxID_ANY, "Main toolbar");
  toolbar->AddTool(wxID_OPEN, "Open", wxNullBitmap, "Open a file");
  toolbar->AddTool(ID_GO_PREV, "Previous", wxNullBitmap, "Go to the previous model");
  toolbar->AddTool(ID_GO_NEXT, "Next", wxNullBitmap, "Go to the next model");
  toolbar->AddSeparator();
  toolbar->AddTool(ID_VIEW_FIT, "Fit", wxNullBitmap, "Fit view to model");
  toolbar->AddSeparator();
  toolbar->AddRadioTool(ID_VIEW_SHADED, "S", wxNullBitmap, wxNullBitmap, "Shaded polygons");
  toolbar->AddRadioTool(ID_VIEW_WIREFRAME, "W", wxNullBitmap, wxNullBitmap, "Wireframe");
  toolbar->AddRadioTool(ID_VIEW_SHADED_WIREFRAME, "SW", wxNullBitmap, wxNullBitmap, "Shading + wireframe");
  toolbar->AddSeparator();
  toolbar->AddTool(ID_TOOLS_TOOLBOX, "Toolbox", wxNullBitmap, "Show/hide toolbox");

  toolbar->Realize();
#endif

  //==========================================================================================================================
  // Main layout
  //==========================================================================================================================

  static Real const SASH_GRAVITY = 0.67;
  ui.main_splitter = new wxSplitterWindow(this, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                                          wxSP_3D | wxSP_BORDER | wxSP_PERMIT_UNSPLIT | wxSP_LIVE_UPDATE);
  ui.main_splitter->SetSashGravity(SASH_GRAVITY);
  ui.main_splitter->SetMinimumPaneSize(MIN_SPLIT_SIZE);
  ui.main_splitter->SetSashInvisible(false);

  // Create the model
  model = new Model;

  // An OpenGL display box for the model
  ui.model_display = new ModelDisplay(ui.main_splitter, model);

  // A tabbed pane for the toolbox
  ui.toolbox = new wxNotebook(ui.main_splitter, wxID_ANY);

  // Segment picking interface
  wxPanel * segments_panel = new wxPanel(ui.toolbox);
  wxBoxSizer * segments_sizer = new wxBoxSizer(wxVERTICAL);
  segments_panel->SetSizer(segments_sizer);
  ui.segments_table = new wxListBox(segments_panel, wxID_ANY);
  segments_sizer->Add(ui.segments_table, 1, wxEXPAND, 0);

  wxStaticText * segment_label_caption = new wxStaticText(segments_panel, wxID_ANY, "Label");
  ui.segment_label = new wxTextCtrl(segments_panel, ID_SEGMENT_LABEL);
  wxBoxSizer * segments_label_sizer = new wxBoxSizer(wxHORIZONTAL);
  segments_label_sizer->Add(segment_label_caption, 0, wxEXPAND | wxRIGHT, 5);
  segments_label_sizer->Add(ui.segment_label, 1, wxEXPAND | wxLEFT, 5);
  segments_sizer->Add(segments_label_sizer, 0, wxEXPAND | wxTOP, 5);

  wxButton * segment_add_btn       =  new wxButton(segments_panel, ID_SEGMENT_ADD,       "Add segment");
  wxButton * segment_remove_btn    =  new wxButton(segments_panel, ID_SEGMENT_REMOVE,    "Remove segment");
  wxButton * segment_expand_btn    =  new wxButton(segments_panel, ID_SEGMENT_EXPAND,    "Expand selection");
  wxButton * segment_contract_btn  =  new wxButton(segments_panel, ID_SEGMENT_CONTRACT,  "Contract selection");

  wxBoxSizer * segments_add_remove_sizer = new wxBoxSizer(wxHORIZONTAL);
  segments_add_remove_sizer->Add(segment_add_btn, 1, wxEXPAND | wxRIGHT, 5);
  segments_add_remove_sizer->Add(segment_remove_btn, 1, wxEXPAND | wxLEFT, 5);
  segments_sizer->Add(segments_add_remove_sizer, 0, wxEXPAND | wxTOP, 5);

  wxBoxSizer * segments_exp_contr_sizer = new wxBoxSizer(wxHORIZONTAL);
  segments_exp_contr_sizer->Add(segment_expand_btn, 1, wxEXPAND | wxRIGHT, 5);
  segments_exp_contr_sizer->Add(segment_contract_btn, 1, wxEXPAND | wxLEFT, 5);
  segments_sizer->Add(segments_exp_contr_sizer, 0, wxEXPAND | wxTOP, 5);

  ui.toolbox->AddPage(segments_panel, "Segments");

  // Point picking interface
  wxPanel * points_panel = new wxPanel(ui.toolbox);
  wxBoxSizer * points_sizer = new wxBoxSizer(wxVERTICAL);
  points_panel->SetSizer(points_sizer);
  ui.points_table = new wxListBox(points_panel, wxID_ANY);
  points_sizer->Add(ui.points_table, 1, wxEXPAND, 0);

  wxStaticText * point_label_caption = new wxStaticText(points_panel, wxID_ANY, "Label");
  ui.point_label = new wxTextCtrl(points_panel, ID_POINT_LABEL);
  wxBoxSizer * points_label_sizer = new wxBoxSizer(wxHORIZONTAL);
  points_label_sizer->Add(point_label_caption, 0, wxEXPAND | wxRIGHT, 5);
  points_label_sizer->Add(ui.point_label, 1, wxEXPAND | wxLEFT, 5);
  points_sizer->Add(points_label_sizer, 0, wxEXPAND | wxTOP, 5);

  wxButton * point_add_btn     =  new wxButton(points_panel, ID_POINT_ADD,     "Add point");
  wxButton * point_remove_btn  =  new wxButton(points_panel, ID_POINT_REMOVE,  "Remove point");
  wxBoxSizer * points_btn_sizer = new wxBoxSizer(wxHORIZONTAL);
  points_btn_sizer->Add(point_add_btn, 1, wxEXPAND | wxRIGHT, 5);
  points_btn_sizer->Add(point_remove_btn, 1, wxEXPAND | wxLEFT, 5);
  points_sizer->Add(points_btn_sizer, 0, wxEXPAND | wxTOP, 5);

  ui.point_snap_to_vertex = new wxCheckBox(points_panel, wxID_ANY, "Snap point to vertex");
  points_sizer->Add(ui.point_snap_to_vertex, 0, wxEXPAND | wxTOP, 5);
  ui.toolbox->AddPage(points_panel, "Points");

  // The main UI is split into two panes
  ui.main_splitter->SplitVertically(ui.model_display, ui.toolbox, -MIN_SPLIT_SIZE);

  // Do the top-level layout
  wxBoxSizer * main_sizer = new wxBoxSizer(wxVERTICAL);
  SetSizer(main_sizer);
  main_sizer->Add(ui.main_splitter, 1, wxEXPAND);

  //==========================================================================================================================
  // Callbacks
  //==========================================================================================================================

  Bind(wxEVT_MENU, &MainWindow::selectAndLoadModel, this, wxID_OPEN);
  Bind(wxEVT_MENU, &MainWindow::OnExit, this, wxID_EXIT);

  Bind(wxEVT_MENU, &ModelDisplay::renderShaded, ui.model_display, ID_VIEW_SHADED);
  Bind(wxEVT_MENU, &ModelDisplay::renderWireframe, ui.model_display, ID_VIEW_WIREFRAME);
  Bind(wxEVT_MENU, &ModelDisplay::renderShadedWireframe, ui.model_display, ID_VIEW_SHADED_WIREFRAME);
  Bind(wxEVT_MENU, &ModelDisplay::setTwoSided, ui.model_display, ID_VIEW_TWO_SIDED);
  Bind(wxEVT_MENU, &ModelDisplay::setFlatShading, ui.model_display, ID_VIEW_FLAT_SHADING);
  Bind(wxEVT_MENU, &ModelDisplay::fitViewToModel, ui.model_display, ID_VIEW_FIT);

  Bind(wxEVT_MENU, &MainWindow::loadPreviousModel, this, ID_GO_PREV);
  Bind(wxEVT_MENU, &MainWindow::loadNextModel, this, ID_GO_NEXT);
  Bind(wxEVT_MENU, &MainWindow::loadPreviousFeatures, this, ID_GO_PREV_FEATURES);
  Bind(wxEVT_MENU, &MainWindow::loadNextFeatures, this, ID_GO_NEXT_FEATURES);

  Bind(wxEVT_MENU, &ModelDisplay::saveScreenshot, ui.model_display, ID_TOOLS_SCREENSHOT);
  Bind(wxEVT_MENU, &MainWindow::toggleToolboxVisible, this, ID_TOOLS_TOOLBOX);

  Bind(wxEVT_BUTTON, &MainWindow::expandPickedSegment, this, ID_SEGMENT_EXPAND);
  Bind(wxEVT_BUTTON, &MainWindow::contractPickedSegment, this, ID_SEGMENT_CONTRACT);
  Bind(wxEVT_BUTTON, &MainWindow::addPickedSegment, this, ID_SEGMENT_ADD);
  Bind(wxEVT_BUTTON, &MainWindow::removeSelectedSegment, this, ID_SEGMENT_REMOVE);
  ui.segments_table->Bind(wxEVT_LISTBOX, &MainWindow::selectSegment, this);

  // Since deselection events are not generated (known wx bug #15603)
  ui.segments_table->Bind(wxEVT_LEFT_UP, &MainWindow::selectSegment, this);
  ui.segments_table->Bind(wxEVT_RIGHT_UP, &MainWindow::selectSegment, this);
  ui.segments_table->Bind(wxEVT_MIDDLE_UP, &MainWindow::selectSegment, this);

  Bind(wxEVT_BUTTON, &MainWindow::addPickedSample, this, ID_POINT_ADD);
  Bind(wxEVT_BUTTON, &MainWindow::removeSelectedSample, this, ID_POINT_REMOVE);
  ui.points_table->Bind(wxEVT_LISTBOX, &MainWindow::selectSample, this);

  // Since deselection events are not generated (known wx bug #15603)
  ui.points_table->Bind(wxEVT_LEFT_UP, &MainWindow::selectSample, this);
  ui.points_table->Bind(wxEVT_RIGHT_UP, &MainWindow::selectSample, this);
  ui.points_table->Bind(wxEVT_MIDDLE_UP, &MainWindow::selectSample, this);

  ui.toolbox->Bind(wxEVT_NOTEBOOK_PAGE_CHANGED, &MainWindow::refreshDisplay, this);

  model->Bind(EVT_MODEL_PATH_CHANGED, &MainWindow::setTitle, this);
  model->Bind(EVT_MODEL_NEEDS_SYNC_SAMPLES, &MainWindow::syncSamples, this);
  model->Bind(EVT_MODEL_NEEDS_SYNC_SEGMENTS, &MainWindow::syncSegments, this);

  Bind(wxEVT_UPDATE_UI, &MainWindow::updateUI, this);  // synchronize menu and toolbar buttons

  //==========================================================================================================================
  // Keyboard shortcuts for menu items
  //==========================================================================================================================

  wxAcceleratorEntry accel[] = {
    wxAcceleratorEntry(wxACCEL_CTRL, (int)'0', ID_VIEW_FIT),
    wxAcceleratorEntry(wxACCEL_CTRL, (int)',', ID_GO_PREV),
    wxAcceleratorEntry(wxACCEL_CTRL, (int)'.', ID_GO_NEXT),
    wxAcceleratorEntry(wxACCEL_CTRL, (int)'G', ID_TOOLS_SCREENSHOT),
    wxAcceleratorEntry(wxACCEL_CTRL, (int)'T', ID_TOOLS_TOOLBOX),
    wxAcceleratorEntry(wxACCEL_CTRL, (int)'[', ID_GO_PREV_FEATURES),
    wxAcceleratorEntry(wxACCEL_CTRL, (int)']', ID_GO_NEXT_FEATURES),
  };
  int num_accel = (int)(sizeof(accel) / sizeof(wxAcceleratorEntry));
  wxAcceleratorTable accel_table(num_accel, accel);
  SetAcceleratorTable(accel_table);

  // Load the initial model, if any
  bool loaded = model->load(app().options().model);
  if (loaded)
  {
    model->setTransform(app().options().model_transform);

    // Load overlays
    overlays.clear();
    for (size_t i = 0; i < app().options().overlays.size(); ++i)
    {
      Model * overlay = new Model;
      loaded = overlay->load(app().options().overlays[i]);
      if (loaded)
      {
        overlay->setTransform(app().options().overlay_transforms[i]);
        overlay->setColor(getPaletteColor((long)i));
        overlays.push_back(overlay);
      }
      else
      {
        delete overlay;
        clearOverlays();
        break;
      }
    }
  }

  //==========================================================================================================================
  // Initial view
  //==========================================================================================================================

  // We have to both set the menu item and call the function since wxEVT_MENU is not generated without actually clicking
  tools_menu->FindItem(ID_TOOLS_TOOLBOX)->Check(false);   setToolboxVisible(false);
  rendering_menu->FindItem(ID_VIEW_SHADED)->Check(true);  ui.model_display->renderShaded();

  rendering_menu->FindItem(ID_VIEW_TWO_SIDED)->Check(app().options().two_sided);
  ui.model_display->setTwoSided(app().options().two_sided);

  rendering_menu->FindItem(ID_VIEW_FLAT_SHADING)->Check(app().options().flat);
  ui.model_display->setFlatShading(app().options().flat);

/*
  setPickSegments(false);
  setPickPoints(false);
  ui->actionToolsToolbox->setChecked(false);
  ui->pickPointsSnapToVertex->setChecked(false);
*/
}

MainWindow::~MainWindow()
{
  model->Unbind(EVT_MODEL_PATH_CHANGED, &MainWindow::setTitle, this);
  model->Unbind(EVT_MODEL_NEEDS_SYNC_SAMPLES, &MainWindow::syncSamples, this);
  model->Unbind(EVT_MODEL_NEEDS_SYNC_SEGMENTS, &MainWindow::syncSegments, this);

  ui.model_display->setModel(NULL);  // this is necessary else we get a segfault when the base class destructor is called after
                                     // this, and can't find the model to deregister callbacks when destroying the display.
  clearOverlays();
  delete model;
}

ModelDisplay *
MainWindow::getRenderDisplay()
{
  return ui.model_display;
}

void
MainWindow::SetTitle(wxString const & title)
{
  if (title.empty())
    BaseType::SetTitle("Browse3D");
  else
  {
    std::string filename = FilePath::objectName(title.ToStdString());
    BaseType::SetTitle(filename + " - Browse3D (" + title + ")");
  }
}

//=============================================================================================================================
// GUI callbacks
//=============================================================================================================================

void
MainWindow::setTitle(wxEvent & event)
{
  SetTitle(model->getPath());
}

void
MainWindow::selectAndLoadModel(wxEvent & event)
{
  if (model->selectAndLoad())
    clearOverlays();
}

void
getMeshPatterns(TheaArray<std::string> & patterns)
{
  patterns.clear();
  patterns.push_back("*.3ds");
  patterns.push_back("*.obj");
  patterns.push_back("*.off");
  patterns.push_back("*.off.bin");
  patterns.push_back("*.ply");
  patterns.push_back("*.pts");
}

void
getFeaturePatterns(TheaArray<std::string> & patterns)
{
  patterns.clear();
  patterns.push_back("*.arff");
  patterns.push_back("*.arff.*");
  patterns.push_back("*.features");
  patterns.push_back("*.features.*");
}

long
fileIndex(TheaArray<std::string> const & files, std::string const & file, TheaArray<std::string> const * patterns = NULL)
{
  std::string fname = FilePath::objectName(file);
  for (size_t i = 0; i < files.size(); ++i)
    if (fname == FilePath::objectName(files[i]))
      return (long)i;

  return -1;
}

long
fileIndex(std::string const & dir, std::string const & file, TheaArray<std::string> & files,
          TheaArray<std::string> const * patterns = NULL)
{
  files.clear();

  if (FileSystem::fileExists(dir))
  {
    files.push_back(dir);
  }
  else
  {
    std::string pat = (patterns ? stringJoin(*patterns, ' ') : "");
    if (FileSystem::getDirectoryContents(dir, files, FileSystem::ObjectType::FILE, pat,
                                         false /* recursive */,
                                         true /* ignore_case */) <= 0)
      return -1;
  }

  return fileIndex(files, file, patterns);
}

void
MainWindow::loadPreviousModel(wxEvent & event)
{
  TheaArray<std::string> patterns;
  getMeshPatterns(patterns);

  TheaArray<std::string> files;
  long index = fileIndex(FilePath::parent(FileSystem::resolve(model->getPath())), model->getPath(), files, &patterns);
  if (files.empty())
    return;

  if (index < 0 || index >= (long)files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  clearOverlays();

  if (index == 0)
    model->load(files[files.size() - 1]);
  else
    model->load(files[index - 1]);
}

void
MainWindow::loadNextModel(wxEvent & event)
{
  TheaArray<std::string> patterns;
  getMeshPatterns(patterns);

  TheaArray<std::string> files;
  long index = fileIndex(FilePath::parent(FileSystem::resolve(model->getPath())), model->getPath(), files, &patterns);
  if (files.empty())
    return;

  if (index < 0 || index >= (long)files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  clearOverlays();

  if (index == (long)files.size() - 1)
    model->load(files[0]);
  else
    model->load(files[index + 1]);
}

void
MainWindow::loadPreviousFeatures(wxEvent & event)
{
  TheaArray<std::string> patterns;
  getFeaturePatterns(patterns);

  TheaArray<std::string> files;
  long index = fileIndex(app().options().features, model->getFeaturesPath(), files, &patterns);
  if (files.empty())
    return;

  if (index < 0 || index >= (long)files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  if (index == 0)
    model->loadFeatures(files[files.size() - 1]);
  else
    model->loadFeatures(files[index - 1]);

  THEA_CONSOLE << "Loaded features " << model->getFeaturesPath();
}

void
MainWindow::loadNextFeatures(wxEvent & event)
{
  TheaArray<std::string> patterns;
  getFeaturePatterns(patterns);

  TheaArray<std::string> files;
  long index = fileIndex(app().options().features, model->getFeaturesPath(), files, &patterns);
  if (files.empty())
    return;

  if (index < 0 || index >= (long)files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  if (index == (long)files.size() - 1)
    model->loadFeatures(files[0]);
  else
    model->loadFeatures(files[index + 1]);

  THEA_CONSOLE << "Loaded features " << model->getFeaturesPath();
}

void
MainWindow::clearOverlays()
{
  for (size_t i = 0; i < overlays.size(); ++i)
    delete overlays[i];

  overlays.clear();
}

void
MainWindow::addPickedSample(wxEvent & event)
{
  std::string label = trimWhitespace(ui.point_label->GetValue().ToStdString());
  if (label.empty())
    return;

  THEA_CONSOLE << "Adding sample with label " << label;

  model->addPickedSample(label, ui.point_snap_to_vertex->GetValue());
  model->invalidatePick();

  ui.points_table->Append(label);
}

void
MainWindow::syncSamples(wxEvent & event)
{
  TheaArray<Model::Sample> const & samples = model->getSamples();
  wxArrayString labels;
  for (size_t i = 0; i < samples.size(); ++i)
  {
    THEA_CONSOLE << "Adding sample with label " << samples[i].label;
    labels.Add(samples[i].label);
  }

  ui.points_table->Set(labels);
}

void
MainWindow::removeSelectedSample(wxEvent & event)
{
  wxArrayInt sel;
  ui.points_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    return;

  int index = sel.Item(0);
  model->selectSample(index);
  THEA_CONSOLE << "Removing sample " << index << " with label " << ui.points_table->GetString(index);

  model->removeSample(index);
  ui.points_table->Delete(index);
}

void
MainWindow::selectSample(wxEvent & event)
{
  wxArrayInt sel;
  ui.points_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    model->selectSample(-1);
  else
    model->selectSample(sel.Item(0));
}

bool
MainWindow::pickPoints() const
{
  return ui.toolbox->IsShown() && ui.toolbox->GetSelection() == POINTS_TAB_INDEX;
}

void
MainWindow::addPickedSegment(wxEvent & event)
{
  std::string label = trimWhitespace(ui.segment_label->GetValue().ToStdString());
  if (label.empty())
    return;

  THEA_CONSOLE << "Adding segment with label " << label;

  model->addPickedSegment(label);
  model->invalidatePickedSegment();

  ui.segments_table->Append(label);
}

void
MainWindow::expandPickedSegment(wxEvent & event)
{
  model->promotePickedSegment(1);
}

void
MainWindow::contractPickedSegment(wxEvent & event)
{
  model->promotePickedSegment(-1);
}

void
MainWindow::syncSegments(wxEvent & event)
{
  TheaArray<Segment> const & segments = model->getSegments();
  wxArrayString labels;
  for (size_t i = 0; i < segments.size(); ++i)
  {
    THEA_CONSOLE << "Adding segment with label " << segments[i].getLabel();
    labels.Add(segments[i].getLabel());
  }

  ui.segments_table->Set(labels);
}

void
MainWindow::removeSelectedSegment(wxEvent & event)
{
  wxArrayInt sel;
  ui.segments_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    return;

  int index = sel.Item(0);
  model->selectSegment(index);
  THEA_CONSOLE << "Removing segment " << index << " with label " << ui.segments_table->GetString(index);

  model->removeSegment(index);
  ui.segments_table->Delete(index);
}

void
MainWindow::selectSegment(wxEvent & event)
{
  wxArrayInt sel;
  ui.segments_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    model->selectSegment(-1);
  else
    model->selectSegment(sel.Item(0));
}

bool
MainWindow::pickSegments() const
{
  return ui.toolbox->IsShown() && ui.toolbox->GetSelection() == SEGMENTS_TAB_INDEX;
}

void
MainWindow::toggleToolboxVisible(wxEvent & event)
{
  setToolboxVisible(!ui.toolbox->IsShown());
}

void
MainWindow::setToolboxVisible(wxCommandEvent & event)
{
  setToolboxVisible(event.IsChecked());
}

void
MainWindow::setToolboxVisible(bool value)
{
  if (ui.toolbox->IsShown() == value)
    return;

  ui.toolbox->Show(value);
  if (value)
    ui.main_splitter->SplitVertically(ui.model_display, ui.toolbox, -MIN_SPLIT_SIZE);
  else
    ui.main_splitter->Unsplit(ui.toolbox);
}

void
MainWindow::updateUI(wxUpdateUIEvent & event)
{
  // TODO
}

void
MainWindow::refreshDisplay(wxEvent & event)
{
  ui.model_display->Refresh();
}

void
MainWindow::OnExit(wxEvent & event)
{
  Close(true);
}

} // namespace Browse3D
