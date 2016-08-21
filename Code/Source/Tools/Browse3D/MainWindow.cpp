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

MainWindow::MainWindow(wxWindow * parent)
: BaseType(parent, wxID_ANY),
  model(NULL),
  model_display(NULL)
{
}

void
MainWindow::init()
{
  setWindowTitle("Browse3D");

  view_type_action_group = new QActionGroup(this);
  ui->actionViewShaded->setActionGroup(view_type_action_group);
  ui->actionViewWireframe->setActionGroup(view_type_action_group);
  ui->actionViewShadedWireframe->setActionGroup(view_type_action_group);

  // Shortcuts for menu options
  ui->actionFileOpen->setShortcuts(QKeySequence::Open);
  ui->actionFileSaveAs->setShortcuts(QKeySequence::SaveAs);
  ui->actionFileQuit->setShortcuts(QKeySequence::Quit);

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

  // Create the model
  model = new Model;

  // An OpenGL display box for the model
  model_display = new ModelDisplay(this, model);
  ui->modelLayout->addWidget(model_display);
  model_display->show();

  // Segment picking interface
  ui->segments_table->setColumnCount(1);
  ui->segments_table->horizontalHeader()->setStretchLastSection(true);

  // Point picking interface
  ui->pointsTable->setColumnCount(1);
  ui->pointsTable->horizontalHeader()->setStretchLastSection(true);

  // Setup signal/slot connections
  connect(ui->actionFileOpen, SIGNAL(triggered(bool)), this, SLOT(selectAndLoadModel()));
  connect(ui->actionFileQuit, SIGNAL(triggered(bool)), this, SLOT(close()));

  connect(ui->actionViewFitViewToModel, SIGNAL(triggered(bool)), model_display, SLOT(fitViewToModel()));
  connect(ui->actionViewWireframe, SIGNAL(triggered(bool)), model_display, SLOT(renderWireframe()));
  connect(ui->actionViewShaded, SIGNAL(triggered(bool)), model_display, SLOT(renderShaded()));
  connect(ui->actionViewShadedWireframe, SIGNAL(triggered(bool)), model_display, SLOT(renderShadedWireframe()));
  connect(ui->actionViewTwoSidedLighting, SIGNAL(toggled(bool)), model_display, SLOT(setTwoSided(bool)));
  connect(ui->actionViewFlatShading, SIGNAL(toggled(bool)), model_display, SLOT(setFlatShading(bool)));

  connect(ui->actionGoPrevious, SIGNAL(triggered(bool)), this, SLOT(loadPreviousModel()));
  connect(ui->actionGoNext,     SIGNAL(triggered(bool)), this, SLOT(loadNextModel()));

  connect(ui->actionGoPreviousFeatures, SIGNAL(triggered(bool)), this, SLOT(loadPreviousFeatures()));
  connect(ui->actionGoNextFeatures,     SIGNAL(triggered(bool)), this, SLOT(loadNextFeatures()));

  connect(ui->actionToolsSaveScreenshot, SIGNAL(triggered(bool)), model_display, SLOT(saveScreenshot()));
  connect(ui->actionToolsToolbox, SIGNAL(toggled(bool)), this, SLOT(setShowToolbox(bool)));

  connect(model, SIGNAL(filenameChanged(QString const &)), this, SLOT(setWindowTitle(QString const &)));
  connect(model, SIGNAL(needsSyncSamples(Model const *)), this, SLOT(syncSamples()));
  connect(model, SIGNAL(needsSyncSegments(Model const *)), this, SLOT(syncSegments()));

  connect(ui->toolbox, SIGNAL(currentChanged(int)), this, SLOT(update()));

  connect(ui->buttonExpandSegment, SIGNAL(clicked()), this, SLOT(expandPickedSegment()));
  connect(ui->buttonContractSegment, SIGNAL(clicked()), this, SLOT(contractPickedSegment()));
  connect(ui->buttonAddSegment, SIGNAL(clicked()), this, SLOT(addPickedSegment()));
  connect(ui->buttonRemoveSegment, SIGNAL(clicked()), this, SLOT(removeSelectedSegment()));
  connect(ui->segments_table, SIGNAL(itemSelectionChanged()), this, SLOT(selectSegment()));

  connect(ui->buttonAddPoint, SIGNAL(clicked()), this, SLOT(addPickedSample()));
  connect(ui->buttonRemovePoint, SIGNAL(clicked()), this, SLOT(removeSelectedSample()));
  connect(ui->pointsTable, SIGNAL(itemSelectionChanged()), this, SLOT(selectSample()));

  // Set/sync default toggle values
  ui->actionViewShaded->trigger();

  ui->actionViewTwoSidedLighting->setChecked(app().options().two_sided);
  model_display->setTwoSided(ui->actionViewTwoSidedLighting->isChecked());

  ui->actionViewFlatShading->setChecked(app().options().flat);
  model_display->setFlatShading(ui->actionViewFlatShading->isChecked());

  setPickSegments(false);
  setPickPoints(false);
  ui->actionToolsToolbox->setChecked(false);
  ui->pickPointsSnapToVertex->setChecked(false);

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
  return model_display;
}

void
MainWindow::selectAndLoadModel()
{
  if (model->selectAndLoad())
    clearOverlays();
}

void
getMeshPatterns(QStringList & patterns)
{
  patterns << "*.3ds" << "*.obj" << "*.off" << "*.off.bin" << "*.pts";
}

void
getFeaturePatterns(QStringList & patterns)
{
  patterns << "*.arff" << "*.arff.*" << "*.features" << "*.features.*";
}

QStringList
getDirFiles(QString const & filename, QStringList const & patterns)
{
  QFileInfo info(filename);
  QDir dir = info.dir();
  if (dir.exists())
    return dir.entryList(patterns, QDir::Files | QDir::Readable, QDir::Name | QDir::IgnoreCase);
  else
    return QStringList();
}

void
MainWindow::loadPreviousModel()
{
  QStringList patterns; getMeshPatterns(patterns);
  QStringList files = getDirFiles(model->getFilename(), patterns);
  if (files.isEmpty())
    return;

  QFileInfo info(model->getFilename());
  int index = files.indexOf(info.fileName());
  if (index < 0 || index >= files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  clearOverlays();

  if (index == 0)
    model->load(info.dir().filePath(files.last()));
  else
    model->load(info.dir().filePath(files[index - 1]));
}

void
MainWindow::loadNextModel()
{
  QStringList patterns; getMeshPatterns(patterns);
  QStringList files = getDirFiles(model->getFilename(), patterns);
  if (files.isEmpty())
    return;

  QFileInfo info(model->getFilename());
  int index = files.indexOf(info.fileName());
  if (index < 0 || index >= files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  clearOverlays();

  if (index == files.size() - 1)
    model->load(info.dir().filePath(files.first()));
  else
    model->load(info.dir().filePath(files[index + 1]));
}

void
MainWindow::loadPreviousFeatures()
{
  QFileInfo fdir(app().options().features);
  if (!fdir.isDir())
    return;

  QStringList patterns; getFeaturePatterns(patterns);
  QStringList files = getDirFiles(app().options().features, patterns);
  if (files.isEmpty())
    return;

  QFileInfo info(model->getFeaturesFilename());
  int index = files.indexOf(info.fileName());
  if (index < 0 || index >= files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  if (index == 0)
  {
    model->loadFeatures(fdir.dir().filePath(files.last()));
    THEA_CONSOLE << info.dir().filePath(files.last());
  }
  else
  {
    model->loadFeatures(fdir.dir().filePath(files[index - 1]));
    THEA_CONSOLE << info.dir().filePath(files.last());
  }

  THEA_CONSOLE << "Loaded features " << model->getFeaturesFilename();
}

void
MainWindow::loadNextFeatures()
{
  QFileInfo fdir(app().options().features);
  if (!fdir.isDir())
    return;

  QStringList patterns; getFeaturePatterns(patterns);
  QStringList files = getDirFiles(app().options().features, patterns);
  if (files.isEmpty())
    return;

  QFileInfo info(model->getFeaturesFilename());
  int index = files.indexOf(info.fileName());
  if (index < 0 || index >= files.size())  // maybe the file was deleted recently?
    index = 0;
  else if (files.size() == 1)
    return;

  if (index == files.size() - 1)
    model->loadFeatures(fdir.dir().filePath(files.first()));
  else
    model->loadFeatures(fdir.dir().filePath(files[index + 1]));

  THEA_CONSOLE << "Loaded features " << model->getFeaturesFilename();
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
  std::string label = point_label->GetValue().ToStdString();
  THEA_CONSOLE << "Adding sample with label" << label;

  model->addPickedSample(label, pick_points_snap_to_vertex->GetValue());
  model->invalidatePick();

  points_table->Append(label);
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

  points_table->Set(labels);
}

void
MainWindow::removeSelectedSample()
{
  wxArrayInt sel;
  points_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    return;

  int index = sel.Item(0);
  model->selectSample(index);
  THEA_CONSOLE << "Removing sample " << item << " with label " << points_table->GetString(index);

  model->removeSample(index);
  points_table->Delete(index);
}

void
MainWindow::selectSample()
{
  wxArrayInt sel;
  points_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    model->selectSample(-1);
  else
    model->selectSample(sel.Item(0));
}

bool
MainWindow::pickPoints() const
{
  return toolbox->IsShown() && toolbox->GetSelection() == POINTS_TAB_INDEX;
}

void
MainWindow::setPickPoints(bool value)
{
  toolbox->Show(value);
  if (value)
    toolbox->SetSelection(POINTS_TAB_INDEX);
}

void
MainWindow::addPickedSegment()
{
  std::string label = segment_label->GetValue().ToStdString();
  THEA_CONSOLE << "Adding segment with label" << label;

  model->addPickedSegment(label);
  model->invalidatePickedSegment();

  segments_table->Append(label);
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

  ui->segments_table->Set(labels);
}

void
MainWindow::removeSelectedSegment()
{
  wxArrayInt sel;
  segments_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    return;

  int index = sel.Item(0);
  model->selectSegment(index);
  THEA_CONSOLE << "Removing segment " << item << " with label " << segments_table->GetString(index);

  model->removeSegment(index);
  segments_table->Delete(index);
}

void
MainWindow::selectSegment()
{
  wxArrayInt sel;
  segments_table->GetSelections(sel);
  if (sel.GetCount() <= 0)
    model->selectSegment(-1);
  else
    model->selectSegment(sel.Item(0));
}

bool
MainWindow::pickSegments() const
{
  return toolbox->IsShown() && toolbox->GetSelection() == SEGMENTS_TAB_INDEX;
}

void
MainWindow::setPickSegments(bool value)
{
  toolbox->Show(value);
  if (value)
    toolbox->SetSelection(SEGMENTS_TAB_INDEX);
}

void
MainWindow::setShowToolbox(bool value)
{
  toolbox->Show(value);
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
    std::string filename = FilePath::nodeName(title.ToStdString());
    BaseType::SetTitle(filename + " - Browse3D (" + title + ")");
  }
}

void
MainWindow::OnExit(wxCommandEvent & event)
{
  Close(true);
}

} // namespace Browse3D
