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

#include "InitBoxWidget.h"

extern double box_wall_speed;
extern double box_dynamic_pressure;
extern double box_stop_filling_factor;
extern double box_stop_dissipation_factor;
extern double box_stop_pressure;
extern double box_side_wall_mod;
extern int box_mode;
extern Qt::CheckState box_random_orientation;
extern Qt::CheckState box_moving_bottom_wall;
extern Qt::CheckState box_side_walls;
extern QString file_path;

InitBoxWidget::InitBoxWidget()
{
	resize(300, 510);
	setMinimumWidth(300);
	setMaximumWidth(300);
	setMinimumHeight(510);
	setMaximumHeight(510);

	setWindowTitle (tr("Init Box"));

	selectAgglomerateButton = new QPushButton(tr("Select Agglomerate"), this);
	selectAgglomerateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectAgglomerateButton->setGeometry(10, 10, 280, 30);
	connect(selectAgglomerateButton, SIGNAL(clicked(bool)), this, SLOT(selectAgglomerate()));

	useCurrentAgglomerateCheckBox = new QCheckBox("Use current agglomerate", this);
	useCurrentAgglomerateCheckBox->setGeometry(10, 50, 280, 30);
	useCurrentAgglomerateCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	useCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	connect(useCurrentAgglomerateCheckBox, SIGNAL(stateChanged(int)), this, SLOT(useCurrentAgglomerateCheckBoxChanged(int)));

	randomOrientationCheckBox = new QCheckBox("Random orientation of agglomerates", this);
	randomOrientationCheckBox->setGeometry(10, 80, 280, 30);
	randomOrientationCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomOrientationCheckBox->setCheckState(box_random_orientation);
	connect(randomOrientationCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRandomOrientation(int)));

	sideWallsCheckBox = new QCheckBox("Side walls", this);
	sideWallsCheckBox->setGeometry(10, 110, 280, 30);
	sideWallsCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	sideWallsCheckBox->setCheckState(box_side_walls);
	connect(sideWallsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setSideWalls(int)));

	wallSpeedEdit = new QDoubleSpinBox(this);
	wallSpeedEdit->setGeometry(210, 150, 80, 30); 
	wallSpeedEdit->setRange(0, 100000);
	wallSpeedEdit->setValue(box_wall_speed);
	wallSpeedEdit->setSingleStep(10);
	connect(wallSpeedEdit, SIGNAL(valueChanged(double)), this, SLOT(setWallSpeed(double)));

	movingBottomWallCheckBox = new QCheckBox("Move bottom wall", this);
	movingBottomWallCheckBox->setGeometry(10, 180, 280, 30);
	movingBottomWallCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	movingBottomWallCheckBox->setCheckState(box_moving_bottom_wall);
	movingBottomWallCheckBox->setToolTip(tr("If checked the bottom wall is moving with the same velocity as the top wall"));
	connect(movingBottomWallCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setMovingBottomWall(int)));

	sideWallModEdit = new QDoubleSpinBox(this);
	sideWallModEdit->setGeometry(210, 220, 80, 30); 
	sideWallModEdit->setRange(0, 100000);
	sideWallModEdit->setDecimals(4);
	sideWallModEdit->setValue(box_side_wall_mod);
	sideWallModEdit->setSingleStep(0.001);
	sideWallModEdit->setToolTip(tr("Modifies the strength of the rolling/sliding interaction of particles with the side walls"));
	connect(sideWallModEdit, SIGNAL(valueChanged(double)), this, SLOT(setSideWallMod(double)));

	dynamicPressureEdit = new QDoubleSpinBox(this);
	dynamicPressureEdit->setGeometry(210, 260, 80, 30); 
	dynamicPressureEdit->setRange(0, 100000000);
	dynamicPressureEdit->setValue(box_dynamic_pressure);
	dynamicPressureEdit->setSingleStep(1);
	dynamicPressureEdit->setToolTip(tr("Mass of the top wall will adjusted to exert the desired pressure\n(only affects dynamic compression mode)"));
	connect(dynamicPressureEdit, SIGNAL(valueChanged(double)), this, SLOT(setDynamicPressure(double)));

	stopFillingFactorEdit = new QDoubleSpinBox(this);
	stopFillingFactorEdit->setGeometry(210, 300, 80, 30); 
	stopFillingFactorEdit->setRange(0.0, 0.75);
	stopFillingFactorEdit->setValue(box_stop_filling_factor);
	stopFillingFactorEdit->setSingleStep(0.01);
	stopFillingFactorEdit->setToolTip(tr("Compression will be stopped when this filling fator is reached"));
	connect(stopFillingFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(setStopFillingFactor(double)));

	stopDissipationFactorEdit = new QDoubleSpinBox(this);
	stopDissipationFactorEdit->setGeometry(210, 340, 80, 30); 
	stopDissipationFactorEdit->setRange(0.0, 1.0);
	stopDissipationFactorEdit->setValue(box_stop_dissipation_factor);
	stopDissipationFactorEdit->setSingleStep(0.01);
	stopDissipationFactorEdit->setEnabled(false);
	stopDissipationFactorEdit->setToolTip(tr("After stopping the wall, the simulation will continue until this ammount of kinetic energy has been dissipated\n(only affects relaxation mode)"));
	connect(stopDissipationFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(setStopDissipationFactor(double)));

	stopPressureEdit = new QDoubleSpinBox(this);
	stopPressureEdit->setGeometry(210, 380, 80, 30); 
	stopPressureEdit->setRange(0, 10000000);
	stopPressureEdit->setValue(box_stop_pressure);
	stopPressureEdit->setSingleStep(10);
	stopPressureEdit->setToolTip(tr("Simulation is stopped when the pressure on the bottom wall exceeds this value\n(only affects shockwave mode)"));
	connect(stopPressureEdit, SIGNAL(valueChanged(double)), this, SLOT(setPenetrationDepth(double)));

	modeSelectionBox = new QComboBox(this);
	modeSelectionBox->setGeometry(160, 420, 130, 30);
	modeSelectionBox->addItem("Static Compression");
	modeSelectionBox->addItem("Relaxation");
	modeSelectionBox->addItem("Shockwave");
	modeSelectionBox->addItem("Dynamic Compression");
	modeSelectionBox->addItem("Open Box");
	modeSelectionBox->setCurrentIndex(box_mode);
	modeSelectionBox->setToolTip(tr("Static Compression: Top wall moves downward at constant speed until stop filling factor is reached.\nRelaxation: Similar to static compression but simulation is continued until the specified ammount of kinetic energy is dissipated\nDynamic Compression: Wall speed influenced by particle<->wall interaction\nOpen Box: No top wall"));
	connect(modeSelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(modeChanged(int)));

	initButton = new QPushButton(tr("Init Compression Box"), this);
	initButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initButton->setGeometry(10, 470, 280, 30);
	connect(initButton, SIGNAL(clicked(bool)), this, SLOT(init()));

	strcpy(filename, "");
	modeChanged(box_mode);
}

InitBoxWidget::~InitBoxWidget(void)
{
}

void InitBoxWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 170, tr("Wall speed (in cm/s):"));
	painter.drawText(10, 240, tr("Side wall modifier:"));
	painter.drawText(10, 280, tr("Dynamic pressure (in kPa):"));
	painter.drawText(10, 320, tr("Stop filling factor):"));
	painter.drawText(10, 360, tr("Stop dissipation factor:"));
	painter.drawText(10, 400, tr("Stop bottom pressure (in kPa):"));
	painter.drawText(10, 440, tr("Select mode:"));
}

void InitBoxWidget::setWallSpeed(double value)
{
	box_wall_speed = value;
}

void InitBoxWidget::setDynamicPressure(double value)
{
	box_dynamic_pressure = value;
}

void InitBoxWidget::setStopFillingFactor(double value)
{
	box_stop_filling_factor = value;
}

void InitBoxWidget::setStopDissipationFactor(double value)
{
	box_stop_dissipation_factor = value;
}

void InitBoxWidget::setPenetrationDepth(double value)
{
	box_stop_pressure = value;
}

void InitBoxWidget::setSideWallMod(double value)
{
	box_side_wall_mod = value;
}

void InitBoxWidget::setRandomOrientation(int value)
{
	box_random_orientation = (Qt::CheckState)value;
}

void InitBoxWidget::setMovingBottomWall(int value)
{
	box_moving_bottom_wall = (Qt::CheckState)value;
}

void InitBoxWidget::setSideWalls(int value)
{
	box_side_walls = (Qt::CheckState)value;

	if(box_side_walls == Qt::Checked)
		sideWallModEdit->setEnabled(true);
	else
		sideWallModEdit->setEnabled(false);
}

void InitBoxWidget::useCurrentAgglomerateCheckBoxChanged(int value)
{
	if(value == Qt::Unchecked)
	{
		if( strcmp(filename, "") == 0)
			useCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	}
}

void InitBoxWidget::selectAgglomerate()
{
	QString temp = QFileDialog::getOpenFileName(this, tr("Select Agglomerate"), file_path, tr("Particle Files (*.dat)"));

	if(!temp.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(temp).path();
        strcpy(filename, temp.toLatin1().data());

		selectAgglomerateButton->setText(QFileInfo(temp).fileName());
		useCurrentAgglomerateCheckBox->setCheckState(Qt::Unchecked);
	}
}

void InitBoxWidget::init()
{
	if(useCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
	{
		if(box_mode == 0)
			emit signalInitCompression(NULL, box_random_orientation, box_side_walls, box_moving_bottom_wall, box_wall_speed, box_stop_filling_factor, box_side_wall_mod);
		else if(box_mode == 1)
			emit signalInitCompressionRelaxation(NULL, box_random_orientation, box_wall_speed, box_stop_filling_factor, box_stop_dissipation_factor);
		else if(box_mode == 2)
			emit signalInitShockwave(NULL, box_random_orientation, box_wall_speed, box_stop_pressure, box_side_wall_mod);
		else if(box_mode == 3)
			emit signalInitDynamicCompression(NULL, box_random_orientation, box_wall_speed, 1e3 * box_dynamic_pressure, box_side_wall_mod);
		else
			emit signalInitOpenBox(NULL, box_random_orientation, box_side_walls, box_side_wall_mod);
	}
	else
	{
		if(box_mode == 0)
			emit signalInitCompression(filename, box_random_orientation, box_side_walls, box_moving_bottom_wall, box_wall_speed, box_stop_filling_factor, box_side_wall_mod);
		else if(box_mode == 1)
			emit signalInitCompressionRelaxation(filename, box_random_orientation, box_wall_speed, box_stop_filling_factor, box_stop_dissipation_factor);
		else if(box_mode == 2)
			emit signalInitShockwave(filename, box_random_orientation, box_wall_speed, box_stop_pressure, box_side_wall_mod);
		else if(box_mode == 3)
			emit signalInitDynamicCompression(filename, box_random_orientation, box_wall_speed, 1e3 * box_dynamic_pressure, box_side_wall_mod);
		else
			emit signalInitOpenBox(filename, box_random_orientation, box_side_walls, box_side_wall_mod);
	}
}

void InitBoxWidget::modeChanged(int value)
{
	box_mode = value;

	if(value == 0)
	{
		wallSpeedEdit->setEnabled(true);
		stopFillingFactorEdit->setEnabled(true);
		stopDissipationFactorEdit->setEnabled(false);
		stopPressureEdit->setEnabled(false);
		dynamicPressureEdit->setEnabled(false);
		sideWallsCheckBox->setEnabled(true);
		movingBottomWallCheckBox->setEnabled(true);
	}
	else if(value == 1)
	{
		wallSpeedEdit->setEnabled(true);
		stopFillingFactorEdit->setEnabled(true);
		stopDissipationFactorEdit->setEnabled(true);
		stopPressureEdit->setEnabled(false);
		dynamicPressureEdit->setEnabled(false);
		sideWallsCheckBox->setEnabled(true);
		movingBottomWallCheckBox->setEnabled(false);
	}
	else if(value == 2)
	{
		wallSpeedEdit->setEnabled(true);
		stopFillingFactorEdit->setEnabled(false);
		stopDissipationFactorEdit->setEnabled(false);
		stopPressureEdit->setEnabled(true);
		dynamicPressureEdit->setEnabled(false);
		sideWallsCheckBox->setEnabled(false);
		movingBottomWallCheckBox->setEnabled(false);
	}
	else if(value == 3)
	{
		wallSpeedEdit->setEnabled(true);
		stopFillingFactorEdit->setEnabled(false);
		stopDissipationFactorEdit->setEnabled(false);
		stopPressureEdit->setEnabled(false);
		dynamicPressureEdit->setEnabled(true);
		sideWallsCheckBox->setEnabled(true);
		movingBottomWallCheckBox->setEnabled(false);
	}
	else
	{
		wallSpeedEdit->setEnabled(false);
		stopFillingFactorEdit->setEnabled(false);
		stopDissipationFactorEdit->setEnabled(false);
		stopPressureEdit->setEnabled(false);
		dynamicPressureEdit->setEnabled(false);
		sideWallsCheckBox->setEnabled(true);
		movingBottomWallCheckBox->setEnabled(false);
	}
}
