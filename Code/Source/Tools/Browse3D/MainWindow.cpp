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
#include "ui_MainWindow.h"
#include "Util.hpp"
#include "../../Application.hpp"
#include <QActionGroup>
#include <QDir>
#include <QFileInfo>
#include <QList>
#include <QStringList>
#include <QTableWidget>
#include <QTableWidgetItem>

namespace Browse3D {

MainWindow::MainWindow(QWidget * parent)
: QMainWindow(parent),
  ui(new Ui::MainWindow),
  model_display(NULL)
{
}

void
MainWindow::init()
{
  ui->setupUi(this);
  setWindowTitle("");

  view_type_action_group = new QActionGroup(this);
  ui->actionViewShaded->setActionGroup(view_type_action_group);
  ui->actionViewWireframe->setActionGroup(view_type_action_group);
  ui->actionViewShadedWireframe->setActionGroup(view_type_action_group);

  // Shortcuts for menu options
  ui->actionFileOpen->setShortcuts(QKeySequence::Open);
  ui->actionFileSaveAs->setShortcuts(QKeySequence::SaveAs);
  ui->actionFileQuit->setShortcuts(QKeySequence::Quit);
  ui->actionGoPrevious->setShortcuts(QKeySequence::MoveToPreviousPage);
  ui->actionGoNext->setShortcuts(QKeySequence::MoveToNextPage);

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

  // Point picked interface
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

  connect(ui->actionGoPrevious, SIGNAL(triggered(bool)), this, SLOT(loadPreviousModel()));
  connect(ui->actionGoNext,     SIGNAL(triggered(bool)), this, SLOT(loadNextModel()));

  connect(ui->actionToolsSaveScreenshot, SIGNAL(triggered(bool)), model_display, SLOT(saveScreenshot()));
  connect(ui->actionToolsPickPoints, SIGNAL(toggled(bool)), this, SLOT(setPickPoints(bool)));

  connect(model, SIGNAL(filenameChanged(QString const &)), this, SLOT(setWindowTitle(QString const &)));
  connect(model, SIGNAL(needsSyncSamples(Model const *)), this, SLOT(syncSamples()));

  connect(ui->buttonAddPoint, SIGNAL(clicked()), this, SLOT(addPickedSample()));
  connect(ui->buttonRemovePoint, SIGNAL(clicked()), this, SLOT(removeSelectedSample()));
  connect(ui->pointsTable, SIGNAL(itemSelectionChanged()), this, SLOT(selectSample()));

  // Set/sync default toggle values
  ui->actionViewShaded->trigger();

  ui->actionViewTwoSidedLighting->setChecked(app().options().two_sided);
  model_display->setTwoSided(ui->actionViewTwoSidedLighting->isChecked());

  setPickPoints(false);
  ui->actionToolsPickPoints->setChecked(false);
  ui->pickPointsSnapToVertex->setChecked(true);

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
  delete ui;
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

QStringList
getDirFiles(QString const & filename)
{
  QFileInfo info(filename);
  QDir dir = info.dir();
  if (dir.exists())
  {
    QStringList name_filter;
    name_filter << "*.3ds" << "*.obj" << "*.off" << "*.off.bin" << "*.pts";

    return dir.entryList(name_filter, QDir::Files | QDir::Readable, QDir::Name | QDir::IgnoreCase);
  }
  else
    return QStringList();
}

void
MainWindow::loadPreviousModel()
{
  QStringList files = getDirFiles(model->getFilename());
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
  QStringList files = getDirFiles(model->getFilename());
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
MainWindow::clearOverlays()
{
  for (array_size_t i = 0; i < overlays.size(); ++i)
    delete overlays[i];

  overlays.clear();
}

void
MainWindow::setWindowTitle(QString const & title)
{
  if (title.isEmpty())
    BaseType::setWindowTitle("Browse3D");
  else
  {
    QFileInfo info(title);
    BaseType::setWindowTitle(info.fileName() + " - Browse3D (" + title + ")");
  }
}

void
MainWindow::addPickedSample()
{
  QString label = ui->pointLabel->text();
  model->addPickedSample(label, ui->pickPointsSnapToVertex->isChecked());
  model->invalidatePick();

  qDebug() << "Adding sample with label" << label;

  QTableWidgetItem * item = new QTableWidgetItem(label);
  item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

  ui->pointsTable->setRowCount((int)model->numSamples());
  ui->pointsTable->setItem((int)model->numSamples() - 1, 0, item);
}

void
MainWindow::syncSamples()
{
  ui->pointsTable->setRowCount((int)model->numSamples());
  TheaArray<Model::Sample> const & samples = model->getSamples();
  for (array_size_t i = 0; i < samples.size(); ++i)
  {
    qDebug() << "Adding sample with label" << samples[i].label;

    QTableWidgetItem * item = new QTableWidgetItem(samples[i].label);
    item->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    ui->pointsTable->setItem((int)i, 0, item);
  }
}

void
MainWindow::removeSelectedSample()
{
  QList<QTableWidgetItem *> sel = ui->pointsTable->selectedItems();
  if (sel.isEmpty())
    return;

  QTableWidgetItem * item = sel.front();
  qDebug() << "Removing sample" << item->row() << "with label" << item->text();

  model->removeSample(item->row());
  ui->pointsTable->removeRow(item->row());
}

void
MainWindow::selectSample()
{
  QList<QTableWidgetItem *> sel = ui->pointsTable->selectedItems();
  if (sel.isEmpty())
    model->selectSample(-1);
  else
  {
    QTableWidgetItem * item = sel.front();
    model->selectSample(item->row());  // ignore column
  }
}

bool
MainWindow::pickPoints() const
{
  return ui->actionToolsPickPoints->isChecked();
}

void
MainWindow::setPickPoints(bool value)
{
  ui->pickPointsPane->setVisible(value);
}

void
MainWindow::closeEvent(QCloseEvent * event)
{
}

} // namespace Browse3D
