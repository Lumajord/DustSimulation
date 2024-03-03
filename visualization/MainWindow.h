 /***************************************************************************
 *   Simulation of particle collisions/agglomeration						*
 *   Copyright (C) 2009 Alexander Seizinger									*
 *	 alexander.seizinger[...]gmx.net										*
 *																			*
 *   This program is free software; you can redistribute it and/or modify	*
 *   it under the terms of the GNU General Public License as published by	*
 *   the Free Software Foundation; either version 2 of the License, or		*
 *   (at your option) any later version.									*
 *																			*
 *   This program is distributed in the hope that it will be useful,		*
 *   but WITHOUT ANYs WARRANTY; without even the implied warranty of		*
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the			*
 *   GNU General Public License for more details.							*
 *																			*
 *   You should have received a copy of the GNU General Public License		*
 *   along with this program; if not, write to the							*
 *   Free Software Foundation, Inc.,										*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.				*
 ***************************************************************************/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "OpenGLWidget.h"

#include <QMainWindow>
#include <QMessageBox>
#include <QAction>
#include <QMenuBar>
#include <QMenu>
#include <QInputDialog>

#include "Constants.h"
#include "StartSimWidget.h"

#include "CollideAgglomeratesWidget.h"
#include "InitAggregateWidget.h"
#include "InitBoxWidget.h"
#include "InitCakeWidget.h"
#include "InitChainWidget.h"
#include "InitWallWidget.h"
#include "AddParticlesWidget.h"
#include "RotationWidget.h"
#include "SliceWidget.h"
#include "DuplicationWidget.h"
#include "ModifyWallWidget.h"
#include "DisplayOptionsWidget.h"
#include "VisualizationOptionsWidget.h"
#include "AgglomerateInfoWidget.h"
#include "MaterialInfoWidget.h"
#include "SimulationThread.h"
#include "SimulationLib.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(SimulationThread *sim_thread);
	~MainWindow();

	// open gl window used for visualization
	OpenGLWidget *glView;

	StartSimWidget *startSimWidget;

	InitChainWidget *initChainWidget;

	InitAggregateWidget *initAggregateWidget;

	InitCakeWidget *initCakeWidget;

	CollideAgglomeratesWidget *collideAgglomeratesWidget;

	WallCollisionWidget *wallCollisionWidget;

	InitBoxWidget *initBoxWidget;

	AddParticlesWidget *addParticlesWidget;

	RotationWidget *rotationWidget;

	SliceWidget *sliceWidget;

	DuplicationWidget *duplicationWidget;

	ModifyWallWidget *modifyWallWidget;

	MaterialInfoWidget *materialInfoWidget;

	AgglomerateInfoWidget *agglomerateInfoWidget;

	VisualizationOptionsWidget *visualizationOptionsWidget;

	DisplayOptionsWidget *displayOptionsWidget;

	// thread running the simulation
	SimulationThread *sim_thread;

	// true if simulation is paused
	bool pause;

	// if true, a screenshot will be taken after the specified number of timesteps
	bool video_mode;

	//
	float cms[3];

	char screenshot_path[200];

	int screenshot_counter;

protected:
	void resizeEvent(QResizeEvent *event);
	void keyPressEvent(QKeyEvent *event);
	void closeEvent(QCloseEvent *event);

private:
    QMenu *fileMenu;
	QMenu *initMenu;
	QMenu *modifyMenu;
	QMenu *aggregateMenu;
	QMenu *materialMenu;
	QMenu *displayMenu;
	QMenu *viewMenu;
	QMenu *helpMenu;
	
	QAction *runSimAct;
	QAction *loadSimAct;
	QAction *saveSimAct;
	QAction *importAct;
	QAction *savePositionsAct;
	QAction *screenshotAct;
    QAction *folderToScreenshotsAct;
    QAction *posFolderToScreenshotsAct;
	QAction *enterSeedAct;
	QAction *exitAct;

	QAction *initChainAct;
	QAction *initAggregateAct;
	QAction *initCakeAct;
	QAction *collideAgglomeratesAct;
	QAction *collideWallAct;
	QAction *initBoxAct;
	QAction *initPullStrengthTestAct;
	QAction *initShearStrengthTestAct;

	QAction *addParticlesAct;
	QAction *rotateAct;
	QAction *sliceAct;
	QAction *duplicateAct;
	QAction *filterFragmentsAct;
	QAction *removeWallsAct;
	QAction *modifyWallsAct;
	QAction *resetContactsAct;
	QAction *disturbParticlesAct;

	QAction *displayCrossSectionAct;
	QAction *moveToOriginAct;
	QAction *showAgglomerateInfoAct;
	QAction *printContactHistogramAct;
	QAction *printFragmentVelocitiesAct;
	QAction *printElasticChargeAct;
	QAction *printFillingFactorProfileAct;
    QAction *printFillingFactorProfileSphereAct;

	QAction *printFractalDimensionAct;

	QAction *loadMaterialAct;
	
	QAction *aboutAct;
	QAction *calcAct;
	QAction *centerViewAct;
	QAction *toggleGyrationRadiusAct;
	QAction *showVisualizationOptionsAct;
	QAction *showDisplayOptionsAct;
	QAction *showMaterialInfoAct;
	QAction *enterInfoTextAct;

public slots:
	void showLoadWidget();
    void showFolderToScreenshotsWidget();
    void showPosFolderToScreenshotsWidget();
	void showSaveWidget();
	void showImportWidget();
	void showSavePositionsWidget();
	void showRunSimWidget();
    void enterSeed();

	void showInitAggregateWidget();
	void showInitCakeWidget();
	void showInitChainWidget();
	void showInitWallWidget();
	void showInitBoxWidget();
	void initPullStrengthTest();
	void initShearStrengthTest();
	void showCollideAgglomeratesWidget();
	void showAddParticlesWidget();
	void showRotationWidget();
	void showSliceWidget();
	void showDuplicationWidget();
	void showModifyWallsWidget();
	void disturbParticles();
	void showAgglomerateInfo();
	void printContactHistogram();
	void printFragmentVelocities();
	void printElasticCharge();
	void printFillingFactorProfile();
    void printFillingFactorProfileSphere();
	void printFractalDimension();
	void showVisualizationOptions();
	void showDisplayOptions();
	void showAboutMessage();
	void showErrorMessage(ErrorCode);

	void initChain();
	void loadMaterial();
	void showMaterialInfo();
	void enterInfoText();
	void displayCrossSection();

	void updateVisualization();
	void updateAgglomerateInfo();
	void centerView();
	void followCMS();
	void customScreenshot();
	void autoScreenshot();
	void numberOfParticlesChanged(int);

#ifdef DRAW_CONTACT_POINTERS
	void contactPointerUpdateDone();
#endif

signals:
	void signalInitSimChain(int, int, int, double, double);
	void signalImpactChainOnAgglomerate(int, double, double);
	void signalInitPullStrengthTest(double);
	void signalInitShearStrengthTest(double);
	void signalDisturbParticles(double);
	
	void signalLoadDataFromFile(char*);
	void signalImportDataFromFile(char*);
	void signalSaveDataToFile(char*);
	void signalLoadMaterial(char*);
	void signalSavePositionsToFile(char*);
    void signalFolderToScreenshots(char*);
    void signalPosFolderToScreenshots(char*);

	void signalStartSim(double, double, double, int, int, int, bool);
	void signalPauseSimulation();
	void signalContinueSimulation();

#ifdef DRAW_CONTACT_POINTERS
	void signalGetContactPointers(multimap<int, vec3>*);
#endif
};

#endif
