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

#include "ModifyWallWidget.h"
#include "SimulationThread.h"

ModifyWallWidget::ModifyWallWidget(SimulationThread *sim_thread)
{
	this->sim_thread = sim_thread;

	resize(300, 300);
	setMinimumWidth(300);
	setMaximumWidth(300);
	setMinimumHeight(380);
	setMaximumHeight(380);
	
	setWindowTitle (tr("Modify Walls"));

	wallSpeedBox = new QDoubleSpinBox(this);
	wallSpeedBox->setGeometry(210, 60, 80, 30); 
	wallSpeedBox->setRange(-100000, 100000);
	wallSpeedBox->setValue(0);
	wallSpeedBox->setSingleStep(10);
	wallSpeedBox->setToolTip(tr("Constant speed in direction of wall normal"));
	//connect(wallSpeedBox, SIGNAL(valueChanged(double)), this, SLOT(setSideWallMod(double)));

	dynamicPressureBox = new QDoubleSpinBox(this);
	dynamicPressureBox->setGeometry(210, 100, 80, 30); 
	dynamicPressureBox->setRange(0, 100000);
	dynamicPressureBox->setValue(0);
	dynamicPressureBox->setSingleStep(1);
	dynamicPressureBox->setToolTip(tr("Pressure where no force is exerted on the wall"));
	//connect(dynamicPressureBox, SIGNAL(valueChanged(double)), this, SLOT(setSideWallMod(double)));

	alphaBox = new QDoubleSpinBox(this);
	alphaBox->setGeometry(210, 140, 80, 30); 
	alphaBox->setRange(0, 1);
	alphaBox->setValue(1);
	alphaBox->setSingleStep(0.1);
	alphaBox->setToolTip(tr("Transparency of wall"));

	compressionModifierBox = new QDoubleSpinBox(this);
	compressionModifierBox->setGeometry(210, 180, 80, 30); 
	compressionModifierBox->setRange(0, 100000);
	compressionModifierBox->setDecimals(4);
	compressionModifierBox->setValue(1);
	compressionModifierBox->setSingleStep(0.1);
	compressionModifierBox->setToolTip(tr("Modifier for the particle<->wall normal interaction"));

	rollingModifierBox = new QDoubleSpinBox(this);
	rollingModifierBox->setGeometry(210, 220, 80, 30); 
	rollingModifierBox->setRange(0, 100000);
	rollingModifierBox->setDecimals(4);
	rollingModifierBox->setValue(1);
	rollingModifierBox->setSingleStep(0.1);
	rollingModifierBox->setToolTip(tr("Modifier for the particle<->wall rolling interaction"));

	slidingModifierBox = new QDoubleSpinBox(this);
	slidingModifierBox->setGeometry(210, 260, 80, 30); 
	slidingModifierBox->setRange(0, 100000);
	slidingModifierBox->setDecimals(4);
	slidingModifierBox->setValue(1);
	slidingModifierBox->setSingleStep(0.1);
	slidingModifierBox->setToolTip(tr("Modifier for the particle<->wall sliding interaction"));

	setWallPropertiesButton = new QPushButton(tr("Set Properties"), this);
	setWallPropertiesButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	setWallPropertiesButton->setGeometry(10, 340, 280, 30);
	setWallPropertiesButton->setEnabled(false);
	connect(setWallPropertiesButton, SIGNAL(clicked(bool)), this, SLOT(setWallProperties()));

	wallSelectionBox = new QComboBox(this);
	wallSelectionBox->setGeometry(160, 10, 130, 30);
	wallSelectionBox->setToolTip(tr("Static Compression: Top wall moves downward at constant speed until stop filling factor is reached.\nRelaxation: Similar to static compression but simulation is continued until the specified ammount of kinetic energy is dissipated\nDynamic Compression: Wall speed influenced by particle<->wall interaction\nOpen Box: No top wall"));
	connect(wallSelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(wallSelectionChanged(int)));
	setupWallSelectionBox(0);
}

ModifyWallWidget::~ModifyWallWidget(void)
{
}

void ModifyWallWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 30, tr("Select Wall:"));
	painter.drawText(10, 80, tr("Wall speed:"));
	painter.drawText(10, 120, tr("Dynamic pressure:"));
	painter.drawText(10, 160, tr("Alpha value:"));
	painter.drawText(10, 200, tr("Compression modifier:"));
	painter.drawText(10, 240, tr("Rolling modifier:"));
	painter.drawText(10, 280, tr("Sliding modifier:"));
}

void ModifyWallWidget::setupWallSelectionBox(int number_of_walls)
{
	int items = wallSelectionBox->count();

	for(int item = items-1; item >= 0; --item)
		wallSelectionBox->removeItem(item);

	if(number_of_walls > 0)
	{
		if(number_of_walls == 1)
			wallSelectionBox->addItem("Wall");
		else
		{
			wallSelectionBox->addItem("Bottom Wall");
			
			if(number_of_walls > 2)
			{
				wallSelectionBox->addItem("Left Wall");
				wallSelectionBox->addItem("Right Wall");
				wallSelectionBox->addItem("Front Wall");
				wallSelectionBox->addItem("Back Wall");

				if(number_of_walls == 6) 
					wallSelectionBox->addItem("Top Wall");
			}
			else
				wallSelectionBox->addItem("Top Wall");

			wallSelectionBox->addItem("All walls");
		}

		wallSelectionBox->setEnabled(true);
		wallSpeedBox->setEnabled(true);
		dynamicPressureBox->setEnabled(true);
		alphaBox->setEnabled(true);
		compressionModifierBox->setEnabled(true);
		rollingModifierBox->setEnabled(true);
		slidingModifierBox->setEnabled(true);
		setWallPropertiesButton->setEnabled(true);
	}
	else
	{
		wallSelectionBox->addItem("No wall available");
		wallSelectionBox->setEnabled(false);
		wallSpeedBox->setEnabled(false);
		dynamicPressureBox->setEnabled(false);
		alphaBox->setEnabled(false);
		compressionModifierBox->setEnabled(false);
		rollingModifierBox->setEnabled(false);
		slidingModifierBox->setEnabled(false);
		setWallPropertiesButton->setEnabled(false);
	}
}

void ModifyWallWidget::wallSelectionChanged(int value)
{
	if(value >= 0 && value < sim_thread->number_of_walls)
	{
		double v_wall = sim_thread->sim.walls[value].velocity[0] * sim_thread->sim.walls[value].normal[0] 
				+ sim_thread->sim.walls[value].velocity[1] * sim_thread->sim.walls[value].normal[1] 
				+ sim_thread->sim.walls[value].velocity[2] * sim_thread->sim.walls[value].normal[2];
		wallSpeedBox->setValue(v_wall);
		dynamicPressureBox->setValue(sim_thread->sim.walls[value].mass);
		alphaBox->setValue(sim_thread->sim.walls[value].alpha);
		compressionModifierBox->setValue(sim_thread->sim.walls[value].compression_modifier);
		rollingModifierBox->setValue(sim_thread->sim.walls[value].rolling_modifier);
		slidingModifierBox->setValue(sim_thread->sim.walls[value].sliding_modifier);
	}
		
}

void ModifyWallWidget::setWallProperties()
{
	emit signalSetWallProperties(wallSelectionBox->currentIndex(), wallSpeedBox->value(), dynamicPressureBox->value(), alphaBox->value(), compressionModifierBox->value(), rollingModifierBox->value(), slidingModifierBox->value());
}