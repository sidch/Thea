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
#include <QMainWindow>

namespace Ui {

class MainWindow;

}

class QActionGroup;

namespace Browse3D {

class Model;
class ModelDisplay;

/** The main application window. */
class MainWindow : public QMainWindow
{
    Q_OBJECT

    typedef QMainWindow BaseType;

  public:
    /** Constructor. */
    explicit MainWindow(QWidget * parent = NULL);

    /** Destructor. */
    ~MainWindow();

    /** Initialize the main window. Should be called at application startup. */
    void init();

    /** Get the display widget showing the model. */
    ModelDisplay * getRenderDisplay();

    /** Check if point-picking is on. */
    bool pickPoints() const;

  public slots:
    /** Select and load a model. */
    void selectAndLoadModel();

    /** Load the previous model in the directory. */
    void loadPreviousModel();

    /** Load the next model in the directory. */
    void loadNextModel();

    /** Set the window title. */
    void setWindowTitle(QString const & title);

    /** Add the currently picked point to the set of samples. */
    void addPickedSample();

    /** Select the sample indicated by the table selection. */
    void selectSample();

    /** Remove the selected sample. */
    void removeSelectedSample();

    /** Sync the displayed list of samples with the model. */
    void syncSamples();

    /** Turn point-picking on/off. */
    void setPickPoints(bool value);

  protected:
    /** Called just before the window is closed. */
    void closeEvent(QCloseEvent * event);

  private:
    Ui::MainWindow * ui;
    QActionGroup * view_type_action_group;
    Model * model;
    ModelDisplay * model_display;

}; // class MainWindow

} // namespace Browse3D

#endif
