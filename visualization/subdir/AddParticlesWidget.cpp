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

#include "AddParticlesWidget.h"

extern double particle_radius;

extern unsigned int add_particles_number_of_particles;
extern double add_particles_migration_probability1;
extern double add_particles_migration_probability2;
extern BAMSelectionMethod add_particles_bam_selection_method;
extern unsigned int add_particles_number_of_chains;
extern unsigned int add_particles_min_chain_length;
extern unsigned int add_particles_max_chain_length;
extern QString file_path;

AddParticlesWidget::AddParticlesWidget()
{
	setMinimumHeight(475);
	setMaximumHeight(475);
	setMinimumWidth(260);
	setMaximumWidth(260);

	setWindowTitle (tr("Add Particles"));

	selectAgglomerateButton = new QPushButton(tr("Select Agglomerate"), this);
	selectAgglomerateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectAgglomerateButton->setGeometry(10, 15, 240, 30);
	connect(selectAgglomerateButton, SIGNAL(clicked(bool)), this, SLOT(selectAgglomerate()));

	useCurrentAgglomerateCheckBox = new QCheckBox("Use current Agglomerate", this);
	useCurrentAgglomerateCheckBox->setGeometry(10, 55, 240, 30);
	useCurrentAgglomerateCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	useCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	connect(useCurrentAgglomerateCheckBox, SIGNAL(stateChanged(int)), this, SLOT(useCurrentAgglomerateChanged(int)));

	numberOfParticlesEdit = new QSpinBox(this);
	numberOfParticlesEdit->setGeometry(160, 95, 90, 30);
	numberOfParticlesEdit->setRange(1, 1000000);
	numberOfParticlesEdit->setValue(add_particles_number_of_particles);
	numberOfParticlesEdit->setSingleStep(1000);
	numberOfParticlesEdit->setToolTip(tr("Number of particles that will be added to the selected aggregate from random directions"));
	connect(numberOfParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(particleNumberChanged(int)));

	migrationProbability1Edit = new QDoubleSpinBox(this);
	migrationProbability1Edit->setGeometry(160, 135, 90, 30);
	migrationProbability1Edit->setRange(0, 1);
	migrationProbability1Edit->setValue(add_particles_migration_probability1);
	migrationProbability1Edit->setSingleStep(0.1);
	migrationProbability1Edit->setToolTip(tr("After being deposited a monomer will migrate with this probability to a position where it establishes contacts with two existing particles"));
	connect(migrationProbability1Edit, SIGNAL(valueChanged(double)), this, SLOT(migrationProbability1Changed(double)));

	migrationProbability2Edit = new QDoubleSpinBox(this);
	migrationProbability2Edit->setGeometry(160, 175, 90, 30);
	migrationProbability2Edit->setRange(0, 1);
	migrationProbability2Edit->setValue(add_particles_migration_probability2);
	migrationProbability2Edit->setSingleStep(0.1);
	migrationProbability2Edit->setToolTip(tr("After being deposited a monomer will migrate with this probability to a position where it establishes contacts with three existing particles"));
	connect(migrationProbability2Edit, SIGNAL(valueChanged(double)), this, SLOT(migrationProbability2Changed(double)));

	bamMethodSelectionBox = new QComboBox(this);
	bamMethodSelectionBox->setGeometry(160, 215, 90, 30);
	bamMethodSelectionBox->addItem("Random");
	bamMethodSelectionBox->addItem("Shortest migration");
	bamMethodSelectionBox->addItem("Most interior");
	bamMethodSelectionBox->setCurrentIndex(add_particles_bam_selection_method);
	bamMethodSelectionBox->setToolTip(tr("Random: Random position\nShortest migration: Position closest to initial position\nMost interior: Position closest to center of mass"));
	connect(bamMethodSelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(bamSelectionMethodChanged(int)));
	
	addBAMButton = new QPushButton(tr("Add RBD Particles"), this);
	addBAMButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	addBAMButton->setGeometry(10, 265, 240, 30);
	connect(addBAMButton, SIGNAL(clicked(bool)), this, SLOT(addBAMParticles()));

	numberOfChainsEdit = new QSpinBox(this);
	numberOfChainsEdit->setGeometry(160, 310, 90, 30);
	numberOfChainsEdit->setRange(1, 100000);
	numberOfChainsEdit->setValue(add_particles_number_of_chains);
	numberOfChainsEdit->setSingleStep(10);
	numberOfChainsEdit->setToolTip(tr("Number of fractal chains that will be added to the selected aggregate from random directions"));
	connect(numberOfChainsEdit, SIGNAL(valueChanged(int)), this, SLOT(numberOfChainsChanged(int)));

	minChainLengthEdit = new QSpinBox(this);
	minChainLengthEdit->setGeometry(160, 350, 90, 30);
	minChainLengthEdit->setRange(1, 100000);
	minChainLengthEdit->setValue(add_particles_min_chain_length);
	minChainLengthEdit->setSingleStep(10);
	minChainLengthEdit->setToolTip(tr("Minimum number of particles of the fractal chains"));
	connect(minChainLengthEdit, SIGNAL(valueChanged(int)), this, SLOT(minChainLengthChanged(int)));

	maxChainLengthEdit = new QSpinBox(this);
	maxChainLengthEdit->setGeometry(160, 390, 90, 30);
	maxChainLengthEdit->setRange(1, 100000);
	maxChainLengthEdit->setValue(add_particles_max_chain_length);
	maxChainLengthEdit->setSingleStep(10);
	maxChainLengthEdit->setToolTip(tr("Maximum number of particles of the fractal chains"));
	connect(maxChainLengthEdit, SIGNAL(valueChanged(int)), this, SLOT(maxChainLengthChanged(int)));

	addFractalChainsButton = new QPushButton(tr("Add Fractal Chains"), this);
	addFractalChainsButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	addFractalChainsButton->setGeometry(10, 430, 240, 30);
	connect(addFractalChainsButton, SIGNAL(clicked(bool)), this, SLOT(addFractalChains()));

	strcpy(filename, "");
}

AddParticlesWidget::~AddParticlesWidget(void)
{
}

void AddParticlesWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 115, tr("Number of particles:"));
	painter.drawText(10, 155, tr("Migration probability:"));
	painter.drawText(10, 195, tr("Migration probability:"));
	painter.drawText(10, 235, tr("Position selection:"));
	painter.drawText(10, 330, tr("Number of chains:"));
	painter.drawText(10, 370, tr("Min chain length:"));
	painter.drawText(10, 410, tr("Max chain length:"));
}

void AddParticlesWidget::selectAgglomerate()
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

void AddParticlesWidget::useCurrentAgglomerateChanged(int value)
{
	if(value == Qt::Unchecked)
	{
		if( strcmp(filename, "") == 0)
			useCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	}
}

void AddParticlesWidget::particleNumberChanged(int value)
{
	add_particles_number_of_particles = (unsigned int)value;
}

void AddParticlesWidget::migrationProbability1Changed(double value)
{
	add_particles_migration_probability1 = value;
}

void AddParticlesWidget::migrationProbability2Changed(double value)
{
	add_particles_migration_probability2 = value;
}

void AddParticlesWidget::bamSelectionMethodChanged(int value)
{
	add_particles_bam_selection_method = (BAMSelectionMethod)value;
}

void AddParticlesWidget::numberOfChainsChanged(int value)
{
	add_particles_number_of_chains = (unsigned int)value;
}

void AddParticlesWidget::minChainLengthChanged(int value)
{
	add_particles_min_chain_length = (unsigned int)value;

	if(add_particles_max_chain_length < add_particles_min_chain_length)
	{
		add_particles_max_chain_length = add_particles_min_chain_length;
		maxChainLengthEdit->setValue(add_particles_min_chain_length);
	}
}

void AddParticlesWidget::maxChainLengthChanged(int value)
{
	add_particles_max_chain_length = (unsigned int)value;

	if(add_particles_max_chain_length < add_particles_min_chain_length)
	{
		add_particles_min_chain_length = add_particles_max_chain_length;
		minChainLengthEdit->setValue(add_particles_max_chain_length);
	}
}
	
void AddParticlesWidget::addBAMParticles()
{
	if(useCurrentAgglomerateCheckBox->checkState() == Qt::Unchecked)
		emit signalAddBAMParticles(add_particles_number_of_particles, add_particles_migration_probability1, add_particles_migration_probability2, filename, add_particles_bam_selection_method);
	else
		emit signalAddBAMParticles(add_particles_number_of_particles, add_particles_migration_probability1, add_particles_migration_probability2, NULL, add_particles_bam_selection_method);
}

void AddParticlesWidget::addFractalChains()
{
	if(useCurrentAgglomerateCheckBox->checkState() == Qt::Unchecked)
		emit signalAddFractalChains(filename, add_particles_number_of_chains, add_particles_min_chain_length, add_particles_max_chain_length);
	else
		emit signalAddFractalChains(NULL, add_particles_number_of_chains, add_particles_min_chain_length, add_particles_max_chain_length);
}
