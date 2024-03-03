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

#include "InitAggregateWidget.h"

extern double particle_radius;

extern unsigned int init_aggregate_particles;
extern double init_aggregate_migration_probability1;
extern double init_aggregate_migration_probability2;
extern BAMSelectionMethod init_aggregate_bam_selection_method;
extern double init_aggregate_chain_ratio;
extern unsigned int init_aggregate_min_chain_length;
extern unsigned int init_aggregate_max_chain_length;
extern double init_aggregate_contact_distribution[6];

InitAggregateWidget::InitAggregateWidget()
{
	setMinimumHeight(580);
	setMaximumHeight(580);
	setMinimumWidth(260);
	setMaximumWidth(260);

	setWindowTitle (tr("Init Aggregate"));

	numberOfParticlesEdit = new QSpinBox(this);
	numberOfParticlesEdit->setGeometry(160, 15, 90, 30);
	numberOfParticlesEdit->setRange(1, 1000000);
	numberOfParticlesEdit->setValue(init_aggregate_particles);
	numberOfParticlesEdit->setToolTip(tr("Number of particles that will be added"));
	numberOfParticlesEdit->setSingleStep(1000);
	connect(numberOfParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(numberOfParticlesChanged(int)));

	migrationProbability1Edit = new QDoubleSpinBox(this);
	migrationProbability1Edit->setGeometry(160, 55, 90, 30);
	migrationProbability1Edit->setRange(0, 1);
	migrationProbability1Edit->setValue(init_aggregate_migration_probability1);
	migrationProbability1Edit->setToolTip(tr("After being deposited a monomer will migrate with this probability to a position where it establishes contacts with two existing particles"));
	migrationProbability1Edit->setSingleStep(0.1);
	connect(migrationProbability1Edit, SIGNAL(valueChanged(double)), this, SLOT(migrationProbability1Changed(double)));

	migrationProbability2Edit = new QDoubleSpinBox(this);
	migrationProbability2Edit->setGeometry(160, 95, 90, 30);
	migrationProbability2Edit->setRange(0, 1);
	migrationProbability2Edit->setValue(init_aggregate_migration_probability2);
	migrationProbability2Edit->setToolTip(tr("After being deposited a monomer will migrate with this probability to a position where it establishes contacts with three existing particles"));
	migrationProbability2Edit->setSingleStep(0.1);
	connect(migrationProbability2Edit, SIGNAL(valueChanged(double)), this, SLOT(migrationProbability2Changed(double)));

	bamMethodSelectionBox = new QComboBox(this);
	bamMethodSelectionBox->setGeometry(160, 135, 90, 30);
	bamMethodSelectionBox->addItem("Random");
	bamMethodSelectionBox->addItem("Shortest migration");
	bamMethodSelectionBox->addItem("Most interior");
	bamMethodSelectionBox->setCurrentIndex(init_aggregate_bam_selection_method);
	bamMethodSelectionBox->setToolTip(tr("Random: Random position\nShortest migration: Position closest to initial position\nMost interior: Position closest to center of mass"));
	connect(bamMethodSelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(bamSelectionMethodChanged(int)));

	initBAMAggregateButton = new QPushButton(tr("Init BAM Aggregate"), this);
	initBAMAggregateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initBAMAggregateButton->setGeometry(10, 180, 240, 30);
	initBAMAggregateButton->setToolTip(tr("Particles will be added to existing aggregate from random directions"));
	connect(initBAMAggregateButton, SIGNAL(clicked(bool)), this, SLOT(initBAMAggregate()));

	chainRatioEdit = new QDoubleSpinBox(this);
	chainRatioEdit->setGeometry(160, 240, 90, 30);
	chainRatioEdit->setRange(0, 1);
	chainRatioEdit->setValue(init_aggregate_chain_ratio);
	chainRatioEdit->setSingleStep(0.1);
	connect(chainRatioEdit, SIGNAL(valueChanged(double)), this, SLOT(chainRatioChanged(double)));

	minChainLengthEdit = new QSpinBox(this);
	minChainLengthEdit->setGeometry(160, 280, 90, 30);
	minChainLengthEdit->setRange(1, 1000);
	minChainLengthEdit->setValue(init_aggregate_min_chain_length);
	minChainLengthEdit->setSingleStep(1);
	connect(minChainLengthEdit, SIGNAL(valueChanged(int)), this, SLOT(minChainLengthChanged(int)));

	maxChainLengthEdit = new QSpinBox(this);
	maxChainLengthEdit->setGeometry(160, 320, 90, 30);
	maxChainLengthEdit->setRange(1, 1000);
	maxChainLengthEdit->setValue(init_aggregate_max_chain_length);
	maxChainLengthEdit->setSingleStep(1);
	connect(maxChainLengthEdit, SIGNAL(valueChanged(int)), this, SLOT(maxChainLengthChanged(int)));

	initFractalAggregateButton = new QPushButton(tr("Init Fractal Aggregate"), this);
	initFractalAggregateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initFractalAggregateButton->setGeometry(10, 365, 240, 30);
	initFractalAggregateButton->setToolTip(tr("Particles will be added to existing aggregate from random directions"));
	connect(initFractalAggregateButton, SIGNAL(clicked(bool)), this, SLOT(initFractalAggregate()));

	treeContacts0Edit = new QDoubleSpinBox(this);
	treeContacts0Edit->setGeometry(10, 440, 48, 30);
	treeContacts0Edit->setRange(0, 1);
	treeContacts0Edit->setValue(init_aggregate_contact_distribution[0]);
	treeContacts0Edit->setSingleStep(0.1);
	connect(treeContacts0Edit, SIGNAL(valueChanged(double)), this, SLOT(treeContacts0Changed(double)));

	treeContacts1Edit = new QDoubleSpinBox(this);
	treeContacts1Edit->setGeometry(58, 440, 48, 30);
	treeContacts1Edit->setRange(0, 1);
	treeContacts1Edit->setValue(init_aggregate_contact_distribution[1]);
	treeContacts1Edit->setSingleStep(0.1);
	connect(treeContacts1Edit, SIGNAL(valueChanged(double)), this, SLOT(treeContacts1Changed(double)));

	treeContacts2Edit = new QDoubleSpinBox(this);
	treeContacts2Edit->setGeometry(106, 440, 48, 30);
	treeContacts2Edit->setRange(0, 1);
	treeContacts2Edit->setValue(init_aggregate_contact_distribution[2]);
	treeContacts2Edit->setSingleStep(0.1);
	connect(treeContacts2Edit, SIGNAL(valueChanged(double)), this, SLOT(treeContacts2Changed(double)));

	treeContacts3Edit = new QDoubleSpinBox(this);
	treeContacts3Edit->setGeometry(154, 440, 48, 30);
	treeContacts3Edit->setRange(0, 1);
	treeContacts3Edit->setValue(init_aggregate_contact_distribution[3]);
	treeContacts3Edit->setSingleStep(0.1);
	connect(treeContacts3Edit, SIGNAL(valueChanged(double)), this, SLOT(treeContacts3Changed(double)));

	treeContacts4Edit = new QDoubleSpinBox(this);
	treeContacts4Edit->setGeometry(202, 440, 48, 30);
	treeContacts4Edit->setRange(0, 1);
	treeContacts4Edit->setValue(init_aggregate_contact_distribution[4]);
	treeContacts4Edit->setSingleStep(0.1);
	connect(treeContacts4Edit, SIGNAL(valueChanged(double)), this, SLOT(treeContacts4Changed(double)));

	initTreeAggregateButton = new QPushButton(tr("Init Tree Aggregate"), this);
	initTreeAggregateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initTreeAggregateButton->setGeometry(10, 485, 240, 30);
	initTreeAggregateButton->setToolTip(tr("Particles will be added using a tree algorithm"));
	connect(initTreeAggregateButton, SIGNAL(clicked(bool)), this, SLOT(initTreeAggregate()));

	initNeedleButton = new QPushButton(tr("Init Needle Deposition"), this);
	initNeedleButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initNeedleButton->setGeometry(10, 540, 240, 30);
	initNeedleButton->setToolTip(tr("Random ballistic deposition on the head of a needle -> only one fractal chain"));
	connect(initNeedleButton, SIGNAL(clicked(bool)), this, SLOT(initNeedleDeposition()));
}

InitAggregateWidget::~InitAggregateWidget(void)
{
}

void InitAggregateWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 35, tr("Number of particles:"));
	painter.drawText(10, 75, tr("Migration probability:"));
	painter.drawText(10, 115, tr("Migration probability:"));
	painter.drawText(10, 155, tr("Position selection:"));

	painter.drawText(10, 260, tr("Chain ratio:"));
	painter.drawText(10, 300, tr("Min chain length:"));
	painter.drawText(10, 340, tr("Max chain length:"));

	painter.drawText(10, 430, tr("Contact distribution:"));

	painter.drawText(10, 565, tr("Number of particles:"));
}

void InitAggregateWidget::numberOfParticlesChanged(int value)
{
	init_aggregate_particles = (unsigned int)value;
}

void InitAggregateWidget::migrationProbability1Changed(double value)
{
	init_aggregate_migration_probability1 = value;
}

void InitAggregateWidget::migrationProbability2Changed(double value)
{
	init_aggregate_migration_probability2 = value;
}

void InitAggregateWidget::bamSelectionMethodChanged(int value)
{
	init_aggregate_bam_selection_method = (BAMSelectionMethod)value;
}

void InitAggregateWidget::chainRatioChanged(double value)
{
	init_aggregate_chain_ratio = value;
}

void InitAggregateWidget::minChainLengthChanged(int value)
{
	init_aggregate_min_chain_length = (unsigned int)value;

	if(init_aggregate_min_chain_length > init_aggregate_max_chain_length)
	{
		init_aggregate_max_chain_length = init_aggregate_min_chain_length;
		maxChainLengthEdit->setValue(init_aggregate_min_chain_length);
	}
}

void InitAggregateWidget::maxChainLengthChanged(int value)
{
	init_aggregate_max_chain_length = (unsigned int)value;
}

void InitAggregateWidget::treeContacts0Changed(double value)
{
	init_aggregate_contact_distribution[0] = value;
}

void InitAggregateWidget::treeContacts1Changed(double value)
{
	init_aggregate_contact_distribution[1] = value;
}

void InitAggregateWidget::treeContacts2Changed(double value)
{
	init_aggregate_contact_distribution[2] = value;
}

void InitAggregateWidget::treeContacts3Changed(double value)
{
	init_aggregate_contact_distribution[3] = value;
}

void InitAggregateWidget::treeContacts4Changed(double value)
{
	init_aggregate_contact_distribution[4] = value;
}

void InitAggregateWidget::initBAMAggregate()
{
	emit signalInitBAMAggregate(init_aggregate_particles, init_aggregate_migration_probability1, init_aggregate_migration_probability2, init_aggregate_bam_selection_method);
}

void InitAggregateWidget::initFractalAggregate()
{
	emit signalInitFractalAggregate(init_aggregate_particles, init_aggregate_migration_probability1, init_aggregate_migration_probability2, init_aggregate_chain_ratio, init_aggregate_min_chain_length, init_aggregate_max_chain_length);
}

void InitAggregateWidget::initTreeAggregate()
{
	emit signalInitTreeAggregate(init_aggregate_particles, init_aggregate_migration_probability1, init_aggregate_migration_probability2, init_aggregate_contact_distribution[0], init_aggregate_contact_distribution[1], init_aggregate_contact_distribution[2], init_aggregate_contact_distribution[3], init_aggregate_contact_distribution[4], init_aggregate_contact_distribution[5]);
}

void InitAggregateWidget::initNeedleDeposition()
{
	emit signalInitNeedleDeposition(init_aggregate_particles);
}