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

#include "InitWallWidget.h"

extern QString file_path;
extern double wall_collision_impact_speed;
extern int wall_collision_impact_angle;
extern double wall_collision_impact_distance;
extern Qt::CheckState wall_collision_random_orientation;

WallCollisionWidget::WallCollisionWidget()
{
	resize(300, 300);
	setMinimumWidth(300);
	setMaximumWidth(300);
	setMinimumHeight(300);
	setMaximumHeight(300);

	setWindowTitle (tr("Setup Wall Collision"));

	impactSpeedEdit = new QDoubleSpinBox(this);
	impactSpeedEdit->setGeometry(210, 10, 80, 30); 
	impactSpeedEdit->setRange(0.001, 100000);
	impactSpeedEdit->setValue(wall_collision_impact_speed);
	impactSpeedEdit->setSingleStep(100);
	connect(impactSpeedEdit, SIGNAL(valueChanged(double)), this, SLOT(setImpactSpeed(double)));

	impactAngleEdit = new QSpinBox(this);
	impactAngleEdit->setGeometry(210, 50, 80, 30); 
	impactAngleEdit->setRange(-80, 80);
	impactAngleEdit->setValue(wall_collision_impact_angle);
	impactAngleEdit->setSingleStep(5);
	connect(impactAngleEdit, SIGNAL(valueChanged(int)), this, SLOT(setImpactAngle(int)));

	impactDistanceEdit = new QDoubleSpinBox(this);
	impactDistanceEdit->setGeometry(210, 90, 80, 30); 
	impactDistanceEdit->setRange(0, 100000);
	impactDistanceEdit->setValue(wall_collision_impact_distance);
	impactDistanceEdit->setToolTip(tr("Initial distance of the aggregate to the wall (in particle radii)"));
	impactDistanceEdit->setSingleStep(1);
	connect(impactDistanceEdit, SIGNAL(valueChanged(double)), this, SLOT(setImpactDistance(double)));

	selectAgglomerateButton = new QPushButton(tr("Select Agglomerate"), this);
	selectAgglomerateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectAgglomerateButton->setGeometry(10, 150, 280, 30);
	connect(selectAgglomerateButton, SIGNAL(clicked(bool)), this, SLOT(selectAgglomerate()));

	useCurrentAgglomerateCheckBox = new QCheckBox("Use current Agglomerate", this);
	useCurrentAgglomerateCheckBox->setGeometry(10, 190, 260, 30);
	useCurrentAgglomerateCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	useCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	connect(useCurrentAgglomerateCheckBox, SIGNAL(stateChanged(int)), this, SLOT(useCurrentAgglomerateChanged(int)));

	randomOrientationCheckBox = new QCheckBox("Random orientation of agglomerate", this);
	randomOrientationCheckBox->setGeometry(10, 220, 260, 30);
	randomOrientationCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomOrientationCheckBox->setCheckState(wall_collision_random_orientation);
	connect(randomOrientationCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRandomOrientation(int)));

	initButton = new QPushButton(tr("Init"), this);
	initButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initButton->setGeometry(10, 260, 280, 30);
	connect(initButton, SIGNAL(clicked(bool)), this, SLOT(initCollision()));

	strcpy(filename, "");
}

WallCollisionWidget::~WallCollisionWidget(void)
{
}

void WallCollisionWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 30, tr("Impact speed (in cm/s):"));
	painter.drawText(10, 70, tr("Impact angle (in degrees):"));
	painter.drawText(10, 110, tr("Impact distance:"));
}

void WallCollisionWidget::selectAgglomerate()
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

void WallCollisionWidget::useCurrentAgglomerateChanged(int value)
{
	if(value == Qt::Unchecked)
	{
		if( strcmp(filename, "") == 0)
			useCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	}
}

void WallCollisionWidget::setRandomOrientation(int value)
{
	wall_collision_random_orientation = (Qt::CheckState)value;
}

void WallCollisionWidget::setImpactAngle(int value)
{
	wall_collision_impact_angle = value;
}

void WallCollisionWidget::setImpactSpeed(double value)
{
	wall_collision_impact_speed = value;
}

void WallCollisionWidget::setImpactDistance(double value)
{
	wall_collision_impact_distance = value;
}

void WallCollisionWidget::initCollision()
{
	if(useCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
		emit signalCollideAgglomerateWithWall(NULL, wall_collision_impact_speed, wall_collision_impact_angle, wall_collision_impact_distance, wall_collision_random_orientation);
	else
		emit signalCollideAgglomerateWithWall(filename, wall_collision_impact_speed, wall_collision_impact_angle, wall_collision_impact_distance, wall_collision_random_orientation);
}
