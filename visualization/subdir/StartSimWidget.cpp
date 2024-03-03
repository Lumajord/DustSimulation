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

#include "StartSimWidget.h"
#include "Constants.h"

extern double sim_timestep;
extern double sim_duration;

extern int print_energy_interval;
extern int print_positions_interval;
extern int take_screenshot_interval;
extern Qt::CheckState follow_cms;
extern int sim_time_display_mode;
extern int sim_averaging_steps;
extern Qt::CheckState sim_use_gravity;
extern double gravity_modifier;
extern double wall_inertia;
extern double sim_azimuthal_acceleration;
extern Qt::CheckState sim_use_sim_azimuthal_acceleration;
extern double azimuthal_acceleration;

#ifdef ENABLE_GRAVITY
	extern bool gravity_enabled;
	extern double gravity_modifier;
#endif

extern bool damping_enabled;

StartSimWidget::StartSimWidget()
{
	setMinimumHeight(645);
	setMaximumHeight(645);
	setMinimumWidth(400);
	setMaximumWidth(400);

	setWindowTitle (tr("Start Simulation"));

	timelengthEdit = new QDoubleSpinBox(this);
	timelengthEdit->setGeometry(290, 10, 100, 30);
	timelengthEdit->setDecimals(0);
	timelengthEdit->setRange(0, 10000000);
	timelengthEdit->setValue(sim_duration);
	timelengthEdit->setSingleStep(100);
	connect(timelengthEdit, SIGNAL(valueChanged(double)), this, SLOT(setDuration(double)));

	timestepEdit = new QDoubleSpinBox(this);
	timestepEdit->setGeometry(290, 50, 100, 30);
	timestepEdit->setDecimals(4);
	timestepEdit->setRange(1e-4, 1000);
	timestepEdit->setValue(sim_timestep);
	timestepEdit->setSingleStep(1e-2);
	connect(timestepEdit, SIGNAL(valueChanged(double)), this, SLOT(setTimestep(double)));

	printEnergyIntervalEdit = new QSpinBox(this);
	printEnergyIntervalEdit->setGeometry(290, 90, 100, 30);
	printEnergyIntervalEdit->setRange(0, 10000);
	printEnergyIntervalEdit->setValue(print_energy_interval);
	printEnergyIntervalEdit->setSingleStep(10);
	printEnergyIntervalEdit->setToolTip(tr("After the specified number of update steps various kinds of energy are logged in an external file"));
	connect(printEnergyIntervalEdit, SIGNAL(valueChanged(int)), this, SLOT(setPrintEnergyInterval(int)));

	printPositionsIntervalEdit = new QSpinBox(this);
	printPositionsIntervalEdit->setGeometry(290, 130, 100, 30);
	printPositionsIntervalEdit->setRange(0, 10000);
	printPositionsIntervalEdit->setValue(print_positions_interval);
	printPositionsIntervalEdit->setSingleStep(10);
	printPositionsIntervalEdit->setToolTip(tr("After the specified number of update steps the positions of the particles are printed to an external file"));
	connect(printPositionsIntervalEdit, SIGNAL(valueChanged(int)), this, SLOT(setPrintPositionsInterval(int)));

	unitOfTimeDisplaySelectionBox = new QComboBox(this);
	unitOfTimeDisplaySelectionBox->setGeometry(290, 180, 100, 30);
	unitOfTimeDisplaySelectionBox->addItem("1 ns");
	unitOfTimeDisplaySelectionBox->addItem("10 ns");
	unitOfTimeDisplaySelectionBox->addItem("100 ns");
	unitOfTimeDisplaySelectionBox->addItem("1 µs");
	unitOfTimeDisplaySelectionBox->setCurrentIndex(sim_time_display_mode);
	unitOfTimeDisplaySelectionBox->setToolTip(tr("Accuracy of the time displayed in the info box (\"Display Options->Display Info\" must be enabled)"));
	connect(unitOfTimeDisplaySelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(unitOfTimeChanged(int)));

	averagingStepsEdit = new QSpinBox(this);
	averagingStepsEdit->setGeometry(290, 220, 100, 30);
	averagingStepsEdit->setRange(0, 10000);
	averagingStepsEdit->setValue(sim_averaging_steps);
	averagingStepsEdit->setSingleStep(10);
	averagingStepsEdit->setToolTip(tr("The forces used to calculate the pressure is averaged over the specified number of update steps (only affects compression mode)"));
	connect(averagingStepsEdit, SIGNAL(valueChanged(int)), this, SLOT(setAveragingSteps(int)));

	followCMSCheckBox = new QCheckBox("Camera follows center of mass", this);
	followCMSCheckBox->setGeometry(10, 255, 240, 30);
	followCMSCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	followCMSCheckBox->setCheckState(follow_cms);
	followCMSCheckBox->setToolTip(tr("If enabled camera will try follow the center of mass"));
	connect(followCMSCheckBox, SIGNAL(stateChanged(int)), this, SLOT(followCMSChanged(int)));

	videoModeCheckBox = new QCheckBox("Enable Video Mode", this);
	videoModeCheckBox->setGeometry(10, 295, 200, 30);
	videoModeCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	videoModeCheckBox->setCheckState(Qt::Unchecked);
	videoModeCheckBox->setToolTip(tr("If enabled snapshots of the current OpenGL window are stored in \"/images/\" subfolder"));
	connect(videoModeCheckBox, SIGNAL(stateChanged(int)), this, SLOT(checkBoxChanged(int)));

	takeScreenshotIntervalEdit = new QSpinBox(this);
	takeScreenshotIntervalEdit->setGeometry(290, 330, 100, 30);
	takeScreenshotIntervalEdit->setRange(10, 10000);
	takeScreenshotIntervalEdit->setValue(take_screenshot_interval);
	takeScreenshotIntervalEdit->setSingleStep(10);
	takeScreenshotIntervalEdit->setEnabled(false);
	takeScreenshotIntervalEdit->setToolTip(tr("Specifies after how many update steps a snapshot is taken"));
	connect(takeScreenshotIntervalEdit, SIGNAL(valueChanged(int)), this, SLOT(setTakeScreenshotInterval(int)));

	useGravityCheckBox = new QCheckBox("Enable Gravity", this);
	useGravityCheckBox->setGeometry(10, 370, 200, 30);
	useGravityCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	useGravityCheckBox->setCheckState(sim_use_gravity);
	useGravityCheckBox->setToolTip(tr("Activates a constant force acting upon all particles dragging them towards the bottom"));
	connect(useGravityCheckBox, SIGNAL(stateChanged(int)), this, SLOT(useGravityChanged(int)));

	gravityModifierEdit = new QDoubleSpinBox(this);
	gravityModifierEdit->setGeometry(290, 405, 100, 30);
	gravityModifierEdit->setRange(0, 10000000);
	gravityModifierEdit->setValue(gravity_modifier);
	gravityModifierEdit->setSingleStep(1);
	gravityModifierEdit->setToolTip(tr("Allows to modify the strength of the external force where a value of 1 is equialent to the gravitational force on earth"));
	connect(gravityModifierEdit, SIGNAL(valueChanged(double)), this, SLOT(gravityModifierChanged(double)));

	wallInertiaEdit = new QDoubleSpinBox(this);
	wallInertiaEdit->setGeometry(290, 445, 100, 30);
	wallInertiaEdit->setRange(0, 10000000);
	wallInertiaEdit->setValue(wall_inertia);
	wallInertiaEdit->setSingleStep(0.1);
	wallInertiaEdit->setToolTip(tr("Modifies the inertia of walls that use dynamic pressure mode"));
	connect(wallInertiaEdit, SIGNAL(valueChanged(double)), this, SLOT(wallInertiaChanged(double)));

	useAzimuthalAccelerationCheckBox = new QCheckBox("Enable Azimuthal Acceleration", this);
	useAzimuthalAccelerationCheckBox->setGeometry(10, 485, 250, 30);
	useAzimuthalAccelerationCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	useAzimuthalAccelerationCheckBox->setCheckState(sim_use_sim_azimuthal_acceleration);
	useAzimuthalAccelerationCheckBox->setToolTip(tr("Activates a constant acceleration of the particles in azimuthal direction (with respect to the z-Axis)"));
	connect(useAzimuthalAccelerationCheckBox, SIGNAL(stateChanged(int)), this, SLOT(useAzimuthalAccelerationChanged(int)));

	azimuthalAccelerationEdit = new QDoubleSpinBox(this);
	azimuthalAccelerationEdit->setGeometry(290, 520, 100, 30);
	azimuthalAccelerationEdit->setRange(0, 100000000);
	azimuthalAccelerationEdit->setValue(sim_azimuthal_acceleration);
	azimuthalAccelerationEdit->setSingleStep(1000);
	azimuthalAccelerationEdit->setToolTip(tr("Allows to modify the strength of the external force where a value of 1 is equivalent to the gravitational force on earth"));
	connect(azimuthalAccelerationEdit, SIGNAL(valueChanged(double)), this, SLOT(angularAccelerationChanged(double)));
	
	useDampingCheckBox = new QCheckBox("Enable Damping", this);
	useDampingCheckBox->setGeometry(10, 560, 250, 30);
	useDampingCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	useDampingCheckBox->setCheckState((Qt::CheckState)damping_enabled);
	useDampingCheckBox->setToolTip(tr("Activates damping of (angular) velocities"));
	connect(useDampingCheckBox, SIGNAL(stateChanged(int)), this, SLOT(enableDampingChanged(int)));

	startSimButton = new QPushButton(tr("Run Simulation"), this);
	startSimButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	startSimButton->setGeometry(70, 605, 260, 30);
	startSimButton->setFocus();
	connect(startSimButton, SIGNAL(clicked(bool)), this, SLOT(startSim()));

	useGravityChanged(sim_use_gravity);
}

StartSimWidget::~StartSimWidget(void)
{
}

void StartSimWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));
	painter.drawText(10, 30, tr("Time duration (in ns):"));
	painter.drawText(10, 70, tr("Update interval (in ns):"));
	painter.drawText(10, 110, tr("Energy print interval (in update steps):"));
	painter.drawText(10, 150, tr("Position print interval (in update steps):"));
	painter.drawText(10, 200, tr("Unit of time (for visualization):"));
	painter.drawText(10, 240, tr("Force averaging period (in update steps):"));
	painter.drawText(10, 350, tr("Screenshot interval (in update steps):"));
	painter.drawText(10, 425, tr("Gravity strength modifier:"));
	painter.drawText(10, 465, tr("Wall inertia modifier:"));
	painter.drawText(10, 540, tr("Azimuthal acceleration: (in 1e6/s^2)"));
}

void StartSimWidget::startSim()
{
	emit signalStartSim(1e-9 * sim_duration, 1e-9 * sim_timestep, print_energy_interval, print_positions_interval, sim_averaging_steps, take_screenshot_interval, videoModeCheckBox->isChecked(), sim_use_gravity);
}

void StartSimWidget::checkBoxChanged(int state)
{
	if(state == Qt::Checked)
		takeScreenshotIntervalEdit->setEnabled(true);
	else
		takeScreenshotIntervalEdit->setEnabled(false);
}

void StartSimWidget::followCMSChanged(int state)
{
	follow_cms = (Qt::CheckState)state;
}

void StartSimWidget::setPrintEnergyInterval(int value)
{
	print_energy_interval = value;
}

void StartSimWidget::setPrintPositionsInterval(int value)
{
	print_positions_interval = value;
}

void StartSimWidget::setTakeScreenshotInterval(int value)
{
	take_screenshot_interval = value;
}

void StartSimWidget::setTimestep(double value)
{
	sim_timestep = value;
}

void StartSimWidget::setDuration(double value)
{
	sim_duration = value;
}

void StartSimWidget::unitOfTimeChanged(int mode)
{
	sim_time_display_mode = mode;
}

void StartSimWidget::setAveragingSteps(int value)
{
	sim_averaging_steps = value;
}

void StartSimWidget::useGravityChanged(int state)
{
	sim_use_gravity = (Qt::CheckState)state;
	gravity_enabled = (bool)state;

	if(state == Qt::Checked)
		gravityModifierEdit->setEnabled(true);
	else
		gravityModifierEdit->setEnabled(false);

	emit signalUpdateView();
}

void StartSimWidget::gravityModifierChanged(double value)
{
	gravity_modifier = value;
	emit signalUpdateView();
}

void StartSimWidget::wallInertiaChanged(double value)
{
	wall_inertia = value;
}

void StartSimWidget::useAzimuthalAccelerationChanged(int state)
{
	sim_use_sim_azimuthal_acceleration = (Qt::CheckState)state;

	if(state == Qt::Checked)
	{
		azimuthalAccelerationEdit->setEnabled(true);
		azimuthal_acceleration = 1e6 * sim_azimuthal_acceleration;
	}
	else
	{
		azimuthalAccelerationEdit->setEnabled(false);
		azimuthal_acceleration = 0;
	}
}

void StartSimWidget::angularAccelerationChanged(double value)
{
	sim_azimuthal_acceleration = value;

	if(useAzimuthalAccelerationCheckBox->isChecked())
		azimuthal_acceleration = 1e6 * sim_azimuthal_acceleration;
}

void StartSimWidget::enableDampingChanged(int value)
{
	damping_enabled = (bool) value;
	emit signalUpdateView();
}
