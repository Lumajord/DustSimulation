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

#ifndef SIMULATIONTHREAD_H
#define SIMULATIONTHREAD_H


#include "SimulationVisualization.h"
#include <QtCore>

class SimulationThread : public QThread
{
	Q_OBJECT

public:
	SimulationThread();
	~SimulationThread();
    
	void run();

	// resizes arrays where data of drawn objects is stored
	void resizeArrays(int new_number_of_particles, int new_number_of_walls);

	void updateGyrationRadius();

	void getCMS(float *cms);

	float getGyrationRadius(float *cms);

	// the actual simulation (+ visualization routines)
	SimulationVis sim;

	// clock to determine if updating positions/colors is necessary 
	QTime last_visualization_update;

	// true if simulation has been paused
	bool pause;

	// false while particle positions & colors and walls are updated (to prevent open gl from drawing)
	bool can_draw;

	int update_counter;

	// gyration radius
	bool draw_gyration_radius;
	std::vector<GyrationRadiusInfo> gyration_radii;

	// position/color arrays
	int number_of_particles;
	double *particle_positions;
	float *particle_colors;

	// walls
	int number_of_walls;
	double *wall_positions;
	double *wall_normals;
	float *alpha_values;

	// true if screenshots are saved after screenshot_interval integration steps
	bool video_mode;
	int screenshot_interval;

public slots:
	void loadFromFile(char*);
    void importFromFile(char*);
	void saveToFile(char*);
	void savePositionsToFile(char*);
	void loadMaterial(char*);

	void initChain(int, int, int, double, double);
	void initChainBox(int, double, double);
	void impactChainOnAgglomerate(int, double, double);
	void initCluster(int, int, double, double, double);
	void initBAMAggregate(unsigned int, double, double, BAMSelectionMethod);
	void initFractalAggregate(unsigned int, double, double, double, unsigned int, unsigned int);
	void initTreeAggregate(unsigned int, double, double, double, double, double, double, double, double);
	void initSimpleCubicLattice(unsigned int, unsigned int, unsigned int, double);
	void initHexagonalLattice(unsigned int, unsigned int, unsigned int, double);
	void initRBDCake(unsigned int, double, double, double, double);
	void initBAMCake(unsigned int, double, double, double, double, double, double, BAMSelectionMethod);
	void initNeedleDeposition(unsigned int);
	void collideAgglomerates(const char*, const char*, bool, bool, double, double, double, bool);
	void hitAndStick(const char*, const char*, bool, bool, double);
	void impactMultiProjectiles(const char*, const char*, unsigned int, double, double, bool, double, bool);
	void impactMultiPCAProjectiles(const char*, unsigned int, unsigned int, double, double, bool, double, bool);
	void collideAgglomerateWithWall(const char*, double, int, double, bool);
	void initCompression(const char*, bool, bool, bool, double, double, double);
	void initCompressionRelaxation(const char*, bool, double, double, double);
	void initShockwave(const char*, bool, double, double, double);
	void initDynamicCompression(const char*, bool, double, double, double);
	void initOpenBox(const char*, bool, bool, double);
	void initPullStrengthTest(double);
	void initShearStrengthTest(double);

	void addBAMParticles(unsigned int, double, double, const char*, BAMSelectionMethod);
	void addFractalChains(const char*, unsigned int, unsigned int, unsigned int);
	void rotateAggregate(double, double, double, double);
	void sliceBox(const char*, double, double, double);
	void sliceSphere(const char*, double, bool);
	void sliceCylinder(const char*, double);
	void sliceTop(const char*, double);
	void sliceBottom(const char*, double);
	void duplicateAggregate(int, int, int, bool, bool);
	void filterFragments();
	void setWallProperties(int, double, double, double, double, double, double);
	void removeWalls();
	void resetContacts();
	void disturbParticles(double);

	void moveToOrigin();
	void getCrossSection(double*, double*);

	void startSimulation(double, double, int, int, int, int, bool, bool);
	void pauseSimulation();
	void continueSimulation();

	void toggleGyrationRadius();

	// updates particle positions & colors and walls and tells the GL window to redraw scene
	void updateDrawnObjects();

	// updates text info such as filling factor, coordination number etc.
	void updateInfo();

	void updateAgglomerateInfo(AgglomerateInfo *agg_info, bool only_inner_particles);
	void printContactHistogram(const char *filename, bool only_inner_particles);
	void printElasticCharge(const char *filename);
	void printFragmentVelocities(const char *filename);
	void printFillingFactorProfile(const char *filename);
    void printFillingFactorProfileSphere(const char *filename);
	void printFractalDimension(const char *filename);

#ifdef DRAW_CONTACT_POINTERS
	void StoreContactPointers(multimap<int,vec3> *);
#endif

signals:
	void signalInitFailed(ErrorCode);
	void signalSaveFailed();
	void signalLoadMaterialFailed();
	void signalMaterialLoaded();

	void signalNumberOfParticlesChanged(int);
	void signalNumberOfWallsChanged(int);

	void signalUpdateGLView();

	void signalTakeScreenshot();

	void signalParticleUpdateDone(int);
	void signalWallsUpdateDone();

#ifdef DRAW_CONTACT_POINTERS
	void ContactPointerUpdateDone();
#endif
};

#endif
