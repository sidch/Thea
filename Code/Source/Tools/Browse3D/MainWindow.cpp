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

namespace Browse3D {

static int const SEGMENTS_TAB_INDEX  =  0;
static int const POINTS_TAB_INDEX    =  1;

MainWindowUI::MainWindowUI()
: model_display(NULL),
  toolbox(NULL),
  points_table(NULL),
  point_label(NULL),
  pick_points_snap_to_vertex(NULL),
  wxListBox * segments_table(NULL),
  wxTextCtrl * segment_label(NULL)
{}

MainWindow::MainWindow(wxWindow * parent)
: BaseType(parent, wxID_ANY),
  model(NULL),
  ui.model_display(NULL)
{}

void
MainWindow::init()
{
  // Set the window title
  setWindowTitle("Browse3D");

  // Set up the main menu
  menubar = new wxMenuBar();

  // File Menu
  wxMenu * file_menu = new wxMenu();
  file_menu->Append(wxID_OPEN,    "&Open"));
  file_menu->Append(wxID_SAVEAS,  "&Save"));
  file_menu->AppendSeparator();
  file_menu->Append(wxID_EXIT,    "&Quit"));
  menubar->Append(file_menu, "&File"));

  // View menu
  wxMenu * view_menu = new wxMenu();
  wxMenu * rendering_menu = new wxMenu();
    rendering_menu->Append(ID_VIEW_SHADED,            "&Shaded");
    rendering_menu->Append(ID_VIEW_WIREFRAME,         "&Wireframe");
    rendering_menu->Append(ID_VIEW_SHADED_WIREFRAME,  "S&haded + wireframe");
    rendering_menu->AppendSeparator();
    rendering_menu->Remove(rendering_menu->Append(-1, wxEmptyString));  // resets radio grouping
    rendering_menu->Append(ID_VIEW_TWO_SIDED,         "&Two-sided lighting");
    rendering_menu->Append(ID_VIEW_FLAT_SHADING,      "&Flat shading");
  view_menu->AppendSubMenu(rendering_menu,  "&Rendering");
  view_menu->Append(ID_VIEW_FIT,            "&Fit view to model");
  menubar->Append(help_menu, "&View");

  // Go menu
  wxMenu * go_menu = new wxMenu();
  go_menu->Append(ID_GO_PREV,           "&Previous model"));
  go_menu->Append(ID_GO_NEXT,           "&Next model"));
  go_menu->AppendSeparator();
  go_menu->Append(ID_GO_PREV_FEATURES,  "Previous features"));
  go_menu->Append(ID_GO_NEXT_FEATURES,  "Next features"));
  menubar->Append(go_menu, "&Go"));

  // Tools menu
  wxMenu * tools_menu = new wxMenu();
  tools_menu->Append(ID_TOOLS_SCREENSHOT,  "&Save screenshot");
  tools_menu->Append(ID_TOOLS_TOOLBOX,     "&Toolbox");
  menubar->Append(tools_menu, "&Tools");

  // About menu
  wxMenu * help_menu = new wxMenu();
  help_menu->Append(wxID_ABOUT,  "&About");
  menubar->Append(help_menu, "&Help");

  SetMenuBar(menubar);

/*
  // Icons for menu options/buttons
  ui->actionFileOpen->setIcon(QIcon::fromTheme("document-open",
        QIcon(toQString(Application::getFullResourcePath("Icons/Tango/scalable/document-open.svg")))));
  ui->actionFileSaveAs->setIcon(QIcon::fromTheme("document-save-as",
        QIcon(toQString(Application::getFullResourcePath("Icons/Tango/scalable/document-save-as.svg")))));

  ui->actionViewFitViewToModel->setIcon(QIcon::fromTheme("zoom-fit-best",
        QIcon(toQString(Application::getFullResourcePath("Icons/Tango/scalable/zoom-fit-best.svg")))));
  ui->actionViewWireframe->setIcon(QIcon::fromTheme("view-wireframe",
        QIcon(toQString(Application::getFullResourcePath("Icons/Custom/48x48/view-wireframe.png")))));
  ui->actionViewShaded->setIcon(QIcon::fromTheme("view-shaded",
        QIcon(toQString(Application::getFullResourcePath("Icons/Custom/48x48/view-shaded.png")))));
  ui->actionViewShadedWireframe->setIcon(QIcon::fromTheme("view-shaded-wireframe",
        QIcon(toQString(Application::getFullResourcePath("Icons/Custom/48x48/view-shaded-wireframe.png")))));

  ui->actionGoPrevious->setIcon(QIcon::fromTheme("go-previous",
        QIcon(toQString(Application::getFullResourcePath("Icons/Tango/scalable/go-previous.svg")))));
  ui->actionGoNext->setIcon(QIcon::fromTheme("go-next",
        QIcon(toQString(Application::getFullResourcePath("Icons/Tango/scalable/go-next.svg")))));

  view_type_action_group = new QActionGroup(this);
  ui->actionViewShaded->setActionGroup(view_type_action_group);
  ui->actionViewWireframe->setActionGroup(view_type_action_group);
  ui->actionViewShadedWireframe->setActionGroup(view_type_action_group);

  // Shortcuts for menu options
  ui->actionFileOpen->setShortcuts(QKeySequence::Open);
  ui->actionFileSaveAs->setShortcuts(QKeySequence::SaveAs);
  ui->actionFileQuit->setShortcuts(QKeySequence::Quit);
*/

  wxSplitterWindow * main_splitter = new wxSplitterWindow(this, wxID_ANY);
  main_splitter->SetSashGravity(0.67);
  main_splitter->SetMinimumPaneSize(50);

  // Create the model
  model = new Model;

  // An OpenGL display box for the model
  ui.model_display = new ModelDisplay(main_splitter, model);
  wxBoxSizer * model_sizer = new wxBoxSizer(wxVERTICAL);
  model_sizer->Add(ui.model_display, 1, wxEXPAND, 0);

  // A tabbed pane for the toolbox
  ui.toolbox = new wxNotebook(main_splitter);

  // Segment picking interface
  wxPanel * segments_panel = new wxPanel(ui.toolbox);
  wxBoxSizer * segments_sizer = new wxBoxSizer(wxVERTICAL);
  ui.segments_table = new wxListBox(segments_panel);
  ui.toolbox->Add(ui.segments_table);

  // Point picking interface
  wxPanel * points_panel = new wxPanel(ui.toolbox);
  wxBoxSizer * points_sizer = new wxBoxSizer(wxVERTICAL);
  ui.points_table = new wxListBox(points_panel);
  ui.toolbox->Add(ui.points_table);

  // Setup signal/slot connections
/*
  connect(ui->actionFileOpen, SIGNAL(triggered(bool)), this, SLOT(selectAndLoadModel()));
  connect(ui->actionFileQuit, SIGNAL(triggered(bool)), this, SLOT(close()));

  connect(ui->actionViewFitViewToModel, SIGNAL(triggered(bool)), ui.model_display, SLOT(fitViewToModel()));
  connect(ui->actionViewWireframe, SIGNAL(triggered(bool)), ui.model_display, SLOT(renderWireframe()));
  connect(ui->actionViewShaded, SIGNAL(triggered(bool)), ui.model_display, SLOT(renderShaded()));
  connect(ui->actionViewShadedWireframe, SIGNAL(triggered(bool)), ui.model_display, SLOT(renderShadedWireframe()));
  connect(ui->actionViewTwoSidedLighting, SIGNAL(toggled(bool)), ui.model_display, SLOT(setTwoSided(bool)));
  connect(ui->actionViewFlatShading, SIGNAL(toggled(bool)), ui.model_display, SLOT(setFlatShading(bool)));

  connect(ui->actionGoPrevious, SIGNAL(triggered(bool)), this, SLOT(loadPreviousModel()));
  connect(ui->actionGoNext,     SIGNAL(triggered(bool)), this, SLOT(loadNextModel()));

  connect(ui->actionGoPreviousFeatures, SIGNAL(triggered(bool)), this, SLOT(loadPreviousFeatures()));
  connect(ui->actionGoNextFeatures,     SIGNAL(triggered(bool)), this, SLOT(loadNextFeatures()));

  connect(ui->actionToolsSaveScreenshot, SIGNAL(triggered(bool)), ui.model_display, SLOT(saveScreenshot()));
  connect(ui->actionToolsToolbox, SIGNAL(toggled(bool)), this, SLOT(setShowToolbox(bool)));

  connect(model, SIGNAL(filenameChanged(QString const &)), this, SLOT(setWindowTitle(QString const &)));
  connect(model, SIGNAL(needsSyncSamples(Model const *)), this, SLOT(syncSamples()));
  connect(model, SIGNAL(needsSyncSegments(Model const *)), this, SLOT(syncSegments()));

  connect(ui.toolbox, SIGNAL(currentChanged(int)), this, SLOT(update()));

  connect(ui->buttonExpandSegment, SIGNAL(clicked()), this, SLOT(expandPickedSegment()));
  connect(ui->buttonContractSegment, SIGNAL(clicked()), this, SLOT(contractPickedSegment()));
  connect(ui->buttonAddSegment, SIGNAL(clicked()), this, SLOT(addPickedSegment()));
  connect(ui->buttonRemoveSegment, SIGNAL(clicked()), this, SLOT(removeSelectedSegment()));
  connect(ui.segments_table, SIGNAL(itemSelectionChanged()), this, SLOT(selectSegment()));

  connect(ui->buttonAddPoint, SIGNAL(clicked()), this, SLOT(addPickedSample()));
  connect(ui->buttonRemovePoint, SIGNAL(clicked()), this, SLOT(removeSelectedSample()));
  connect(ui->pointsTable, SIGNAL(itemSelectionChanged()), this, SLOT(selectSample()));

  // Set/sync default toggle values
  ui->actionViewShaded->trigger();

  ui->actionViewTwoSidedLighting->setChecked(app().options().two_sided);
  ui.model_display->setTwoSided(ui->actionViewTwoSidedLighting->isChecked());

  ui->actionViewFlatShading->setChecked(app().options().flat);
  ui.model_display->setFlatShading(ui->actionViewFlatShading->isChecked());

  setPickSegments(false);
  setPickPoints(false);
  ui->actionToolsToolbox->setChecked(false);
  ui->pickPointsSnapToVertex->setChecked(false);
*/

  // Load the initial model, if any
  bool loaded = model->load(app().options().model);
  if (loaded)
  {
    model->setTransform(app().options().model_transform);

    // Load overlays
    overlays.clear();
    for (int i = 0; i < app().options().overlays.size(); ++i)
    {
      Model * overlay = new Model;
      loaded = overlay->load(app().options().overlays[i]);
      if (loaded)
      {
        overlay->setTransform(app().options().overlay_transforms[i]);
        overlay->setColor(getPaletteColor(i));
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
}

MainWindow::~MainWindow()
{
  delete model;
}

ModelDisplay *
MainWindow::getRenderDisplay()
{
  return ui.model_display;
}

void
MainWindow::selectAndLoadModel()
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
fileIndex(std::string const & dir, std::string const & file, TheaArray<std::string> const * patterns = NULL)
{
  TheaArray<std::string> files;
  if (FileSystem::getDirectoryContents(dir, files, FileSystem::ObjectType::FILE, stringJoin(patterns, ' '), false) <= 0)
    return;

  std::string fname = FilePath::objectName(file);
  for (array_size_t i = 0; i < files.size(); ++i)
    if (fname == FilePath::objectName(files[i]))
      return (long)i;

  return -1;
}

void
MainWindow::loadPreviousModel()
{
  TheaArray<std::string> patterns;
  getMeshPatterns(patterns);

  long index = fileIndex(app().options().features, model->getPath(), &patterns);
  if (index < 0)  // maybe the file was deleted recently?
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
MainWindow::loadNextModel()
{
  TheaArray<std::string> patterns;
  getMeshPatterns(patterns);

  long index = fileIndex(app().options().features, model->getPath(), &patterns);
  if (index < 0)  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  clearOverlays();

  if (index == files.size() - 1)
    model->load(files[0]);
  else
    model->load(files[index + 1]);
}

void
MainWindow::loadPreviousFeatures()
{
  TheaArray<std::string> patterns;
  getFeaturePatterns(patterns);

  long index = fileIndex(app().options().features, model->getFeaturesPath(), &patterns);
  if (index < 0)  // maybe the file was deleted recently?
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
MainWindow::loadNextFeatures()
{
  TheaArray<std::string> patterns;
  getFeaturePatterns(patterns);

  long index = fileIndex(app().options().features, model->getFeaturesPath(), &patterns);
  if (index < 0)  // maybe the file was deleted recently?
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
  for (array_size_t i = 0; i < overlays.size(); ++i)
    delete overlays[i];

  overlays.clear();
}

void
MainWindow::addPickedSample()
{
  std::string label = ui.point_label->GetValue().ToStdString();
  THEA_CONSOLE << "Adding sample with label" << label;

  model->addPickedSample(label, pick_points_snap_to_vertex->GetValue());
  model->invalidatePick();

  ui.points_table->Append(label);
}

void
MainWindow::syncSamples()
{
  TheaArray<Model::Sample> const & samples = model->getSamples();
  std::vector<std::string> labels((size_t)samples.size());
  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    THEA_CONSOLE << "Adding sample with label" << samples[i].label;
    labels[(size_t)i] = segments[i].label;
  }

  ui.points_table->Set(labels);
}

void
MainWindow::removeSelectedSample()
{
  wxArrayInt sel;
  ui.points_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    return;

  int index = sel.Item(0);
  model->selectSample(index);
  THEA_CONSOLE << "Removing sample " << item << " with label " << ui.points_table->GetString(index);

  model->removeSample(index);
  ui.points_table->Delete(index);
}

void
MainWindow::selectSample()
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
MainWindow::setPickPoints(bool value)
{
  ui.toolbox->Show(value);
  if (value)
    ui.toolbox->SetSelection(POINTS_TAB_INDEX);
}

void
MainWindow::addPickedSegment()
{
  std::string label = ui.segment_label->GetValue().ToStdString();
  THEA_CONSOLE << "Adding segment with label" << label;

  model->addPickedSegment(label);
  model->invalidatePickedSegment();

  ui.segments_table->Append(label);
}

void
MainWindow::expandPickedSegment()
{
  model->promotePickedSegment(1);
}

void
MainWindow::contractPickedSegment()
{
  model->promotePickedSegment(-1);
}

void
MainWindow::syncSegments()
{
  TheaArray<Segment> const & segments = model->getSegments();
  std::vector<std::string> labels((size_t)segments.size());
  for (array_size_t i = 0; i < segments.size(); ++i)
  {
    THEA_CONSOLE << "Adding segment with label" << segments[i].getLabel();
    labels[(size_t)i] = segments[i].getLabel();
  }

  ui.segments_table->Set(labels);
}

void
MainWindow::removeSelectedSegment()
{
  wxArrayInt sel;
  ui.segments_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    return;

  int index = sel.Item(0);
  model->selectSegment(index);
  THEA_CONSOLE << "Removing segment " << item << " with label " << ui.segments_table->GetString(index);

  model->removeSegment(index);
  ui.segments_table->Delete(index);
}

void
MainWindow::selectSegment()
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
MainWindow::setPickSegments(bool value)
{
  ui.toolbox->Show(value);
  if (value)
    ui.toolbox->SetSelection(SEGMENTS_TAB_INDEX);
}

void
MainWindow::setShowToolbox(bool value)
{
  ui.toolbox->Show(value);
}

//=============================================================================================================================
// GUI callbacks etc
//=============================================================================================================================

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

void
MainWindow::OnExit(wxCommandEvent & event)
{
  Close(true);
}

} // namespace Browse3D
