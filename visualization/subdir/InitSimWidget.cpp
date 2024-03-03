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

#include "InitSimWidget.h"

InitSimWidget::InitSimWidget()
{
	resize(400, 320);
	setMinimumHeight(320);
	setMaximumHeight(320);
	setMinimumWidth(400);
	setMaximumWidth(400);

	setWindowTitle (tr("Init Simulation "));

	particleNumberEdit = new QSpinBox(this);
	particleNumberEdit->setGeometry(300, 10, 80, 30); 
	particleNumberEdit->setRange(1, 10000);
	particleNumberEdit->setValue(10);
	particleNumberEdit->setSingleStep(1);

	particleImpactNumberEdit = new QSpinBox(this);
	particleImpactNumberEdit->setGeometry(300, 50, 80, 30); 
	particleImpactNumberEdit->setRange(1, 1000);
	particleImpactNumberEdit->setValue(1);
	particleImpactNumberEdit->setSingleStep(1);

	targetParticleEdit = new QSpinBox(this);
	targetParticleEdit->setGeometry(300, 90, 80, 30); 
	targetParticleEdit->setRange(1, particleNumberEdit->value());
	targetParticleEdit->setValue(particleNumberEdit->value()/2);
	targetParticleEdit->setSingleStep(1);

	impactSpeedEdit = new QDoubleSpinBox(this);
	impactSpeedEdit->setGeometry(300, 130, 80, 30); 
	impactSpeedEdit->setRange(0.001, 100000);
	impactSpeedEdit->setValue(2000);
	impactSpeedEdit->setSingleStep(100);

	angularIrregularityEdit = new QDoubleSpinBox(this);
	angularIrregularityEdit->setGeometry(300, 170, 80, 30); 
	angularIrregularityEdit->setRange(0, 1);
	angularIrregularityEdit->setValue(0);
	angularIrregularityEdit->setSingleStep(0.1);

	fillingFactorEdit = new QDoubleSpinBox(this);
	fillingFactorEdit->setGeometry(300, 210, 80, 30); 
	fillingFactorEdit->setRange(0, 1);
	fillingFactorEdit->setValue(0);
	fillingFactorEdit->setSingleStep(0.05);

	initChainButton = new QPushButton(tr("Init Chain"), this);
	initChainButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initChainButton->setGeometry(20, 280, 170, 30);

	initClusterButton = new QPushButton(tr("Init Cluster"), this);
	initClusterButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initClusterButton->setGeometry(210, 280, 170, 30);

	connect(particleNumberEdit, SIGNAL(valueChanged(int)), this, SLOT(adjustTargetParticleEditRange(int)));
}

InitSimWidget::~InitSimWidget(void)
{
}

void InitSimWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 30, tr("Number of particles of target:"));
	painter.drawText(10, 70, tr("Number of particles of impacting projectile:"));
	painter.drawText(10, 110, tr("Particle the projectile will impact on:"));
	painter.drawText(10, 150, tr("Impact speed (in cm/s):"));
	painter.drawText(10, 190, tr("Angular irregularity of chain:"));
	painter.drawText(10, 230, tr("Filling factor of cluster aggregates:"));
}

void InitSimWidget::adjustTargetParticleEditRange(int max_value)
{
	targetParticleEdit->setRange(1, max_value);
}
