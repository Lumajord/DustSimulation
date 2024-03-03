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

#include "InitCakeWidget.h"

extern double particle_radius;

extern unsigned int init_cake_x_particles;
extern unsigned int init_cake_y_particles;
extern unsigned int init_cake_z_particles;
extern double init_cake_filling_factor;
extern unsigned int init_cake_particles;
extern double init_cake_migration_probability1;
extern double init_cake_migration_probability2;
extern BAMSelectionMethod init_cake_bam_selection_method;
extern double init_cake_x_size;
extern double init_cake_y_size;
extern double init_cake_top_slice_factor;
extern double init_cake_bottom_slice_factor;

InitCakeWidget::InitCakeWidget()
{
	setMinimumHeight(695);
	setMaximumHeight(695);
	setMinimumWidth(260);
	setMaximumWidth(260);

	setWindowTitle (tr("Init Cake"));

	xParticlesEdit = new QSpinBox(this);
	xParticlesEdit->setGeometry(160, 15, 90, 30);
	xParticlesEdit->setRange(1, 100);
	xParticlesEdit->setValue(init_cake_x_particles);
	xParticlesEdit->setSingleStep(1);
	connect(xParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(xParticlesChanged(int)));

	yParticlesEdit = new QSpinBox(this);
	yParticlesEdit->setGeometry(160, 55, 90, 30);
	yParticlesEdit->setRange(1, 100);
	yParticlesEdit->setValue(init_cake_y_particles);
	yParticlesEdit->setSingleStep(1);
	connect(yParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(yParticlesChanged(int)));

	zParticlesEdit = new QSpinBox(this);
	zParticlesEdit->setGeometry(160, 95, 90, 30);
	zParticlesEdit->setRange(1, 100);
	zParticlesEdit->setValue(init_cake_z_particles);
	zParticlesEdit->setSingleStep(1);
	connect(zParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(zParticlesChanged(int)));

	fillingFactorEdit = new QDoubleSpinBox(this);
	fillingFactorEdit->setGeometry(160, 135, 90, 30);
	fillingFactorEdit->setRange(0.0, 1);
	fillingFactorEdit->setValue(init_cake_filling_factor);
	fillingFactorEdit->setToolTip(tr("Ratio of randomly selected particles that will be left out from grid"));
	fillingFactorEdit->setSingleStep(0.1);
	connect(fillingFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(fillingFactorChanged(double)));

	initSimpleCubicButton = new QPushButton(tr("Init Simple Cubic"), this);
	initSimpleCubicButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initSimpleCubicButton->setGeometry(10, 180, 240, 30);
	connect(initSimpleCubicButton, SIGNAL(clicked(bool)), this, SLOT(initSimpleCubicLattice()));

	initHexagonalButton = new QPushButton(tr("Init Hexagonal"), this);
	initHexagonalButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initHexagonalButton->setGeometry(10, 220, 240, 30);
	connect(initHexagonalButton, SIGNAL(clicked(bool)), this, SLOT(initHexagonalLattice()));

	numberOfParticlesEdit = new QSpinBox(this);
	numberOfParticlesEdit->setGeometry(160, 290, 90, 30);
    numberOfParticlesEdit->setRange(1, 10000000);
	numberOfParticlesEdit->setValue(init_cake_particles);
	numberOfParticlesEdit->setToolTip(tr("Number of particles that will be added from random directions"));
	numberOfParticlesEdit->setSingleStep(1000);
	connect(numberOfParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(particleNumberChanged(int)));

	migrationProbability1Edit = new QDoubleSpinBox(this);
	migrationProbability1Edit->setGeometry(160, 330, 90, 30);
	migrationProbability1Edit->setRange(0, 1);
	migrationProbability1Edit->setValue(init_cake_migration_probability1);
	migrationProbability1Edit->setToolTip(tr("After being deposited a monomer will migrate with this probability to a position where it establishes contacts with two existing particles"));
	migrationProbability1Edit->setSingleStep(0.1);
	connect(migrationProbability1Edit, SIGNAL(valueChanged(double)), this, SLOT(migrationProbability1Changed(double)));

	migrationProbability2Edit = new QDoubleSpinBox(this);
	migrationProbability2Edit->setGeometry(160, 370, 90, 30);
	migrationProbability2Edit->setRange(0, 1);
	migrationProbability2Edit->setValue(init_cake_migration_probability2);
	migrationProbability2Edit->setToolTip(tr("After being deposited a monomer will migrate with this probability to a position where it establishes contacts with three existing particles"));
	migrationProbability2Edit->setSingleStep(0.1);
	connect(migrationProbability2Edit, SIGNAL(valueChanged(double)), this, SLOT(migrationProbability2Changed(double)));

	bamMethodSelectionBox = new QComboBox(this);
	bamMethodSelectionBox->setGeometry(160, 410, 90, 30);
	bamMethodSelectionBox->addItem("Random");
	bamMethodSelectionBox->addItem("Shortest migration");
	bamMethodSelectionBox->addItem("Most interior");
	bamMethodSelectionBox->setCurrentIndex(init_cake_bam_selection_method);
	bamMethodSelectionBox->setToolTip(tr("Random: Random position\nShortest migration: Position closest to initial position\nMost interior: Position closest to center of mass"));
	connect(bamMethodSelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(bamSelectionMethodChanged(int)));

	xSizeEdit = new QDoubleSpinBox(this);
	xSizeEdit->setGeometry(160, 450, 90, 30);
	xSizeEdit->setRange(0.0, 10000);
	xSizeEdit->setValue(init_cake_x_size);
	xSizeEdit->setToolTip(tr("x-Size of the base plate in particle radii"));
	xSizeEdit->setSingleStep(1);
	connect(xSizeEdit, SIGNAL(valueChanged(double)), this, SLOT(xSizeChanged(double)));

	ySizeEdit = new QDoubleSpinBox(this);
	ySizeEdit->setGeometry(160, 490, 90, 30);
	ySizeEdit->setRange(0.0, 10000);
	ySizeEdit->setValue(init_cake_y_size);
	ySizeEdit->setToolTip(tr("y-Size of the base plate in particle radii"));
	ySizeEdit->setSingleStep(1);
	connect(ySizeEdit, SIGNAL(valueChanged(double)), this, SLOT(ySizeChanged(double)));

	topSliceFactorEdit = new QDoubleSpinBox(this);
	topSliceFactorEdit->setGeometry(160, 530, 90, 30);
	topSliceFactorEdit->setRange(0, 1);
	topSliceFactorEdit->setValue(init_cake_top_slice_factor);
	topSliceFactorEdit->setToolTip(tr("Specifies what amount of the RBD aggregate is sliced from the top (to prevent lower filling factors at the top)"));
	topSliceFactorEdit->setSingleStep(0.05);
	connect(topSliceFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(topSliceFactorChanged(double)));

	bottomSliceFactorEdit = new QDoubleSpinBox(this);
	bottomSliceFactorEdit->setGeometry(160, 570, 90, 30);
	bottomSliceFactorEdit->setRange(0, 1);
	bottomSliceFactorEdit->setValue(init_cake_bottom_slice_factor);
	bottomSliceFactorEdit->setToolTip(tr("Specifies what amount of the RBD aggregate is sliced from the bottom (to prevent higher filling factors at the bottom)"));
	bottomSliceFactorEdit->setSingleStep(0.05);
	connect(bottomSliceFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(bottomSliceFactorChanged(double)));

	initBAMCakeButton = new QPushButton(tr("Init BAM Cake"), this);
	initBAMCakeButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initBAMCakeButton->setGeometry(10, 615, 240, 30);
	initBAMCakeButton->setToolTip(tr("Random ballistic deposition with migration on the specified base plate"));
	connect(initBAMCakeButton, SIGNAL(clicked(bool)), this, SLOT(initBAMCake()));

	initRBDCakeButton = new QPushButton(tr("Init RBD Cake"), this);
	initRBDCakeButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initRBDCakeButton->setGeometry(10, 655, 240, 30);
	initRBDCakeButton->setToolTip(tr("Random ballistic deposition on the specified base plate"));
	connect(initRBDCakeButton, SIGNAL(clicked(bool)), this, SLOT(initRBDCake()));
}

InitCakeWidget::~InitCakeWidget(void)
{
}

void InitCakeWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 35, tr("x-Particles:"));
	painter.drawText(10, 75, tr("y-Particles:"));
	painter.drawText(10, 115, tr("z-Particles:"));
	painter.drawText(10, 155, tr("Filling Factor:"));

	painter.drawText(10, 310, tr("Number of particles:"));
	painter.drawText(10, 350, tr("Migration probability:"));
	painter.drawText(10, 390, tr("Migration probability:"));
	painter.drawText(10, 430, tr("Position selection:"));
	painter.drawText(10, 470, tr("x-Size (in µm):"));
	painter.drawText(10, 510, tr("y-Size (in µm):"));
	painter.drawText(10, 550, tr("Top slice factor:"));
	painter.drawText(10, 590, tr("Bottom slice factor:"));
}

void InitCakeWidget::xParticlesChanged(int value)
{
	init_cake_x_particles = (unsigned int)value;
}

void InitCakeWidget::yParticlesChanged(int value)
{
	init_cake_y_particles = (unsigned int)value;
}

void InitCakeWidget::zParticlesChanged(int value)
{
	init_cake_z_particles = (unsigned int)value;
}

void InitCakeWidget::fillingFactorChanged(double value)
{
	init_cake_filling_factor = value;
}

void InitCakeWidget::particleNumberChanged(int value)
{
	init_cake_particles = (unsigned int)value;
}

void InitCakeWidget::migrationProbability1Changed(double value)
{
	init_cake_migration_probability1 = value;
}

void InitCakeWidget::migrationProbability2Changed(double value)
{
	init_cake_migration_probability2 = value;
}

void InitCakeWidget::bamSelectionMethodChanged(int value)
{
	init_cake_bam_selection_method = (BAMSelectionMethod)value;
}

void InitCakeWidget::xSizeChanged(double value)
{
	init_cake_x_size = value;
}

void InitCakeWidget::ySizeChanged(double value)
{
	init_cake_y_size = value;
}

void InitCakeWidget::topSliceFactorChanged(double value)
{
	init_cake_top_slice_factor = value;
}

void InitCakeWidget::bottomSliceFactorChanged(double value)
{
	init_cake_bottom_slice_factor = value;
}

void InitCakeWidget::initSimpleCubicLattice()
{
	emit signalInitSimpleCubicLattice(init_cake_x_particles, init_cake_y_particles, init_cake_z_particles, init_cake_filling_factor);
}

void InitCakeWidget::initHexagonalLattice()
{
	emit signalInitHexagonalLattice(init_cake_x_particles, init_cake_y_particles, init_cake_z_particles, init_cake_filling_factor);
}

void InitCakeWidget::initBAMCake()
{
	if(init_cake_migration_probability1 > 0)
		emit signalInitBAMCake(init_cake_particles, init_cake_x_size * 1e-4, init_cake_y_size * 1e-4, init_cake_migration_probability1, init_cake_migration_probability2, init_cake_top_slice_factor, init_cake_bottom_slice_factor, init_cake_bam_selection_method);
	else
		emit signalInitRBDCake(init_cake_particles, init_cake_x_size * 1e-4, init_cake_y_size * 1e-4, init_cake_top_slice_factor, init_cake_bottom_slice_factor);
}

void InitCakeWidget::initRBDCake()
{
	emit signalInitRBDCake(init_cake_particles, init_cake_x_size * 1e-4, init_cake_y_size * 1e-4, init_cake_top_slice_factor, init_cake_bottom_slice_factor);
}
