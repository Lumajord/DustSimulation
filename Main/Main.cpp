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

#ifdef WIN32
	#pragma warning(disable: 4005 4018 4244 4267 4996) // signed/unsigned, loss of precision, deprecation...
#endif



// Simulation
#include <stdio.h>
#include <stdlib.h>

// GUI
#include "MainWindow.h"
#include <QApplication>
#include "SimulationThread.h"

#ifndef TRACK_PARTICLE_ORIENTATION
	#error Tracking of particle orientation not enabled!
#endif

#ifndef TRACK_DISSIPATED_ENERGY
	#error Tracking of dissipated energy not enabled!
#endif

#ifndef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	#error Tracking of dissipated energy per particle not enabled!
#endif

#ifndef ENABLE_GRAVITY
	#error Gravity not enabled!
#endif

// for sctreenshots
QWaitCondition wait_for_screenshot;
QMutex mutex;

int main(int argc, char **argv)
{
	SimulationThread sim_thread;

	QApplication app(argc, argv);

	QGLFormat glf = QGLFormat::defaultFormat(); 
	glf.setSampleBuffers(true); 
	glf.setSamples(2); 
	QGLFormat::setDefaultFormat(glf); 

	MainWindow window(&sim_thread);
	window.show();

	// stop sim thread if still running
	sim_thread.pause = true;
	sim_thread.wait();

    return app.exec();
}
