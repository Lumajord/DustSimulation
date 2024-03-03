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

#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H

#include "GL/glew.h"
#if defined (_WIN32)
#include "GL/wglew.h"
#endif


#include "GL/glu.h"
#include <vector>
#include <QtOpenGL/QGLWidget>
#include <QOpenGLFunctions>

#include <QMessageBox>

#include "Constants.h"
#include "Wall.h"
#include "Shaders.h"
#include "OpenGLNavigationWidget.h"

const float FIELD_OF_VIEW = 45.0;

class SimulationThread;

struct char16{ char data[16]; };

class OpenGLWidget : public OpenGLNavigationWidget
{
    Q_OBJECT

public:

#ifdef DRAW_CONTACT_POINTERS
	OpenGLWidget(QWidget *parent, vector<WallInfo> **walls, multimap<int, vec3> **contact_pointers);
#else
    OpenGLWidget(QWidget *parent, SimulationThread *sim_thread);
#endif 

    ~OpenGLWidget();

	void setNumberOfParticles(int number_of_particles);

	GLuint compileProgram(const char *vsource, const char *fsource);
	
    GLuint program;

	int particle_array_size;

    GLuint vbPosition;
    GLuint vbColor;

	SimulationThread *sim_thread;

	// true if initialization went fine
	bool initialized;

	// quadric used for walls
	GLUquadricObj *quadric;
	GLUquadricObj *quadric2;

	float point_size_scale_factor;

	// texture storing noise
	static const int noise_text_size = 128;
	GLuint noise_texture;
	GLubyte *noise_texture_ptr;

	// stores contact pointers
	std::multimap<int, vec3> **contact_pointers;

	// used to determine text width/height
	QFont font;
	QFont font_script;
	QFontMetrics *font_metrics;

	// stores tics for the key
    std::vector< char16 > key_tic_labels;
	unsigned int key_tics;
	bool key_log_scale;
	double key_min_value;
	double key_max_value;

public slots:
	void selectShader();
	void updateVertexBuffers();
	void updateKey();
    void paintGL();

protected:
    void initializeGL();
	void renderKey();
    void resizeGL(int width, int height);
};

// precalculated values for gravity vector
static const int LOD = 17;
static const double cos_table[LOD] = {cos(0.0), cos(0.125*M_PI), cos(0.25*M_PI), cos(0.375*M_PI), cos(0.5*M_PI), cos(0.625*M_PI), cos(0.75*M_PI), cos(0.875*M_PI), cos(M_PI), cos(1.125*M_PI), cos(1.25*M_PI), cos(1.375*M_PI), cos(1.5*M_PI), cos(1.625*M_PI), cos(1.75*M_PI), cos(1.875*M_PI), cos(2.0*M_PI)};
static const double sin_table[LOD] = {sin(0.0), sin(0.125*M_PI), sin(0.25*M_PI), sin(0.375*M_PI), sin(0.5*M_PI), sin(0.625*M_PI), sin(0.75*M_PI), sin(0.875*M_PI), sin(M_PI), sin(1.125*M_PI), sin(1.25*M_PI), sin(1.375*M_PI), sin(1.5*M_PI), sin(1.625*M_PI), sin(1.75*M_PI), sin(1.875*M_PI), sin(2.0*M_PI)};
//static const double cos_table[LOD] = {cos(0.0), cos(0.25*M_PI), cos(0.5*M_PI), cos(0.75*M_PI), cos(M_PI), cos(1.25*M_PI), cos(1.5*M_PI), cos(1.75*M_PI), cos(2.0*M_PI)};
//static const double sin_table[LOD] = {sin(0.0), sin(0.25*M_PI), sin(0.5*M_PI), sin(0.75*M_PI), sin(M_PI), sin(1.25*M_PI), sin(1.5*M_PI), sin(1.75*M_PI), sin(2.0*M_PI)};

#endif
