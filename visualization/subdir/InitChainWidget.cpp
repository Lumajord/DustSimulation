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

#include "InitChainWidget.h"

extern int init_chain_particles;
extern int init_chain_impact_particles;
extern int init_chain_target_particle;
extern double init_chain_impact_speed;
extern double init_chain_angular_irregularity;
extern Qt::CheckState init_chain_impact_current_agglomerate;

InitChainWidget::InitChainWidget()
{
	resize(400, 280);
	setMinimumHeight(280);
	setMaximumHeight(280);
	setMinimumWidth(400);
	setMaximumWidth(400);

	setWindowTitle (tr("Init Chains"));

	particleNumberEdit = new QSpinBox(this);
	particleNumberEdit->setGeometry(310, 10, 80, 30); 
	particleNumberEdit->setRange(1, 1000000);
	particleNumberEdit->setValue(init_chain_particles);
	particleNumberEdit->setSingleStep(1);

	particleImpactNumberEdit = new QSpinBox(this);
	particleImpactNumberEdit->setGeometry(310, 50, 80, 30); 
	particleImpactNumberEdit->setRange(0, 10000);
	particleImpactNumberEdit->setValue(init_chain_impact_particles);
	particleImpactNumberEdit->setSingleStep(1);

	targetParticleEdit = new QSpinBox(this);
	targetParticleEdit->setGeometry(310, 90, 80, 30); 
	targetParticleEdit->setRange(1, particleNumberEdit->value());
	targetParticleEdit->setValue(init_chain_target_particle);
	targetParticleEdit->setSingleStep(1);

	impactSpeedEdit = new QDoubleSpinBox(this);
	impactSpeedEdit->setGeometry(310, 130, 80, 30); 
	impactSpeedEdit->setRange(0, 100000);
	impactSpeedEdit->setValue(init_chain_impact_speed);
	impactSpeedEdit->setSingleStep(10);

	angularIrregularityEdit = new QDoubleSpinBox(this);
	angularIrregularityEdit->setGeometry(310, 170, 80, 30); 
	angularIrregularityEdit->setRange(0, 1);
	angularIrregularityEdit->setValue(init_chain_angular_irregularity);
	angularIrregularityEdit->setSingleStep(0.1);

	impactOnCurrentAgglomerateCheckBox = new QCheckBox("Impact on current agglomerate", this);
	impactOnCurrentAgglomerateCheckBox->setGeometry(10, 200, 300, 30);
	impactOnCurrentAgglomerateCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	impactOnCurrentAgglomerateCheckBox->setCheckState(init_chain_impact_current_agglomerate);

	initChainButton = new QPushButton(tr("Init Chain"), this);
	initChainButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initChainButton->setGeometry(10, 240, 380, 30);

	connect(impactOnCurrentAgglomerateCheckBox, SIGNAL(stateChanged(int)), this, SLOT(checkBoxChanged(int)));
	connect(particleNumberEdit, SIGNAL(valueChanged(int)), this, SLOT(particleNumberChanged(int)));
	connect(particleImpactNumberEdit, SIGNAL(valueChanged(int)), this, SLOT(particleImpactNumberChanged(int)));
	connect(targetParticleEdit, SIGNAL(valueChanged(int)), this, SLOT(targetParticleChanged(int)));
	connect(impactSpeedEdit, SIGNAL(valueChanged(double)), this, SLOT(impactSpeedChanged(double)));
	connect(angularIrregularityEdit, SIGNAL(valueChanged(double)), this, SLOT(angularIrregularityChanged(double)));
}

InitChainWidget::~InitChainWidget(void)
{
}

void InitChainWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 30, tr("Number of particles of target:"));
	painter.drawText(10, 70, tr("Number of particles of impacting projectile:"));
	painter.drawText(10, 110, tr("Particle the projectile will impact on:"));
	painter.drawText(10, 150, tr("Impact speed (in cm/s):"));
	painter.drawText(10, 190, tr("Angular irregularity of chain:"));
}

void InitChainWidget::particleNumberChanged(int max_value)
{
	targetParticleEdit->setRange(1, max_value);

	init_chain_particles = max_value;
}

void InitChainWidget::checkBoxChanged(int state)
{
	init_chain_impact_current_agglomerate = (Qt::CheckState)state;

	if(state == Qt::Checked)
	{
		particleNumberEdit->setEnabled(false);
		targetParticleEdit->setEnabled(false);
	}
	else
	{
		particleNumberEdit->setEnabled(true);
		targetParticleEdit->setEnabled(true);
	}
}

void InitChainWidget::particleImpactNumberChanged(int value)
{
	init_chain_impact_particles = value;
}

void InitChainWidget::targetParticleChanged(int value)
{
	init_chain_target_particle = value;
}

void InitChainWidget::impactSpeedChanged(double value)
{
	init_chain_impact_speed = value;
}

void InitChainWidget::angularIrregularityChanged(double value)
{
	init_chain_angular_irregularity = value;
}
