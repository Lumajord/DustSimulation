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

#include "CollideAgglomeratesWidget.h"

extern double agglomerates_collision_impact_speed;
extern double agglomerates_collision_impact_parameter;
extern double agglomerates_collision_initial_separation;
extern Qt::CheckState agglomerates_collision_random_orientation1;
extern Qt::CheckState agglomerates_collision_random_orientation2;
extern Qt::CheckState agglomerates_collision_minimize_distance;
extern unsigned int agglomerates_collision_projectile_count;
extern unsigned int agglomerates_collision_projectile_size;
extern Qt::CheckState agglomerates_collision_use_PCA_projectiles;
extern Qt::CheckState agglomerates_collision_random_impact_parameter;
extern QString file_path;

CollideAgglomeratesWidget::CollideAgglomeratesWidget(void)
{
	setMinimumWidth(300);
	setMaximumWidth(300);
	setMinimumHeight(605);
	setMaximumHeight(605);

	setWindowTitle (tr("Setup Agglomerate Collision"));

	impactSpeedEdit = new QDoubleSpinBox(this);
	impactSpeedEdit->setGeometry(200, 10, 90, 30); 
	impactSpeedEdit->setRange(0.001, 100000);
	impactSpeedEdit->setValue(agglomerates_collision_impact_speed);
	impactSpeedEdit->setSingleStep(100);
	connect(impactSpeedEdit, SIGNAL(valueChanged(double)), this, SLOT(setImpactSpeed(double)));

	impactParameterEdit = new QDoubleSpinBox(this);
	impactParameterEdit->setGeometry(200, 50, 90, 30); 
	impactParameterEdit->setRange(-1, 1);
	impactParameterEdit->setValue(agglomerates_collision_impact_parameter);
	impactParameterEdit->setSingleStep(0.05);
	connect(impactParameterEdit, SIGNAL(valueChanged(double)), this, SLOT(setImpactParameter(double)));

	initialSeparationEdit = new QDoubleSpinBox(this);
	initialSeparationEdit->setGeometry(200, 90, 90, 30); 
	initialSeparationEdit->setRange(0.0, 10000);
	initialSeparationEdit->setValue(agglomerates_collision_initial_separation);
	initialSeparationEdit->setSingleStep(1);
	initialSeparationEdit->setToolTip(tr("Initial distance between agglomerates (in radii of a monomere)"));
	connect(initialSeparationEdit, SIGNAL(valueChanged(double)), this, SLOT(setInitialSeparation(double)));

	minimizeDistanceCheckBox = new QCheckBox("Minimize distance", this);
	minimizeDistanceCheckBox->setGeometry(10, 125, 280, 30);
	minimizeDistanceCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	minimizeDistanceCheckBox->setCheckState(agglomerates_collision_minimize_distance);
	minimizeDistanceCheckBox->setToolTip(tr("If enabled the initial separation of the agglomerates is as low as possible"));
	connect(minimizeDistanceCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setMinimizeDistance(int)));

	hitAndStickCheckBox = new QCheckBox("Hit and stick", this);
	hitAndStickCheckBox->setGeometry(10, 160, 280, 30);
	hitAndStickCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	hitAndStickCheckBox->setCheckState(Qt::Unchecked);
	hitAndStickCheckBox->setToolTip(tr("If enabled the selected aggregates will be attached to each other without any internal restructuring"));
	connect(hitAndStickCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setHitAndStick(int)));

	selectFirstAgglomerateButton = new QPushButton(tr("Select first agglomerate"), this);
	selectFirstAgglomerateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectFirstAgglomerateButton->setGeometry(10, 205, 280, 30);
	connect(selectFirstAgglomerateButton, SIGNAL(clicked(bool)), this, SLOT(selectFirstAgglomerate()));

	selectSecondAgglomerateButton = new QPushButton(tr("Select second agglomerate"), this);
	selectSecondAgglomerateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectSecondAgglomerateButton->setGeometry(10, 245, 280, 30);
	connect(selectSecondAgglomerateButton, SIGNAL(clicked(bool)), this, SLOT(selectSecondAgglomerate()));

	randomOrientation1CheckBox = new QCheckBox("Random orientation of 1st agglomerate", this);
	randomOrientation1CheckBox->setGeometry(10, 290, 280, 30);
	randomOrientation1CheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomOrientation1CheckBox->setCheckState(agglomerates_collision_random_orientation1);
	randomOrientation1CheckBox->setToolTip(tr("Toggles random rotation of first agglomerate"));
	connect(randomOrientation1CheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRandomOrientation1(int)));

	randomOrientation2CheckBox = new QCheckBox("Random orientation of 2nd agglomerate", this);
	randomOrientation2CheckBox->setGeometry(10, 320, 280, 30);
	randomOrientation2CheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomOrientation2CheckBox->setCheckState(agglomerates_collision_random_orientation2);
	randomOrientation2CheckBox->setToolTip(tr("Toggles random rotation of second agglomerate"));
	connect(randomOrientation2CheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRandomOrientation2(int)));

	initCollisionButton = new QPushButton(tr("Init Collision"), this);
	initCollisionButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initCollisionButton->setGeometry(10, 365, 280, 30);
	connect(initCollisionButton, SIGNAL(clicked(bool)), this, SLOT(initCollision()));

	projectileCountEdit = new QSpinBox(this);
	projectileCountEdit->setGeometry(200, 420, 90, 30); 
	projectileCountEdit->setRange(0, 10000);
	projectileCountEdit->setValue(agglomerates_collision_projectile_count);
	projectileCountEdit->setSingleStep(1);
	projectileCountEdit->setToolTip(tr("Number of projectiles"));
	connect(projectileCountEdit, SIGNAL(valueChanged(int)), this, SLOT(setProjectileCount(int)));

	projectileSizeEdit = new QSpinBox(this);
	projectileSizeEdit->setGeometry(200, 460, 90, 30); 
	projectileSizeEdit->setRange(0, 100000);
	projectileSizeEdit->setValue(agglomerates_collision_projectile_size);
	projectileSizeEdit->setSingleStep(1);
	projectileSizeEdit->setToolTip(tr("Number of monomers the PCA projectiles will be composed of"));
	connect(projectileSizeEdit, SIGNAL(valueChanged(int)), this, SLOT(setProjectileSize(int)));

	usePCAProjectilesCheckBox = new QCheckBox("Use PCA projectiles", this);
	usePCAProjectilesCheckBox->setGeometry(10, 495, 280, 30);
	usePCAProjectilesCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	usePCAProjectilesCheckBox->setCheckState(agglomerates_collision_use_PCA_projectiles);
	usePCAProjectilesCheckBox->setToolTip(tr("If checked, instead of the second aggregate, PCA aggregates of specified size will be used as projectiles"));
	connect(usePCAProjectilesCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setUsePCAProjectiles(int)));

	randomImpactParameterCheckBox = new QCheckBox("Random impact parameter", this);
	randomImpactParameterCheckBox->setGeometry(10, 525, 280, 30);
	randomImpactParameterCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomImpactParameterCheckBox->setCheckState(agglomerates_collision_random_impact_parameter);
	randomImpactParameterCheckBox->setToolTip(tr("If checked the impact parameter of incoming projectiles will be randomized"));
	connect(randomImpactParameterCheckBox, SIGNAL(stateChanged(int)), this, SLOT(setRandomImpactParameter(int)));

	initMultiCollisionButton = new QPushButton(tr("Init Multi Collision"), this);
	initMultiCollisionButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	initMultiCollisionButton->setGeometry(10, 565, 280, 30);
	connect(initMultiCollisionButton, SIGNAL(clicked(bool)), this, SLOT(initMultiCollision()));

	if(agglomerates_collision_use_PCA_projectiles == Qt::Checked)
		projectileSizeEdit->setEnabled(true);
	else
		projectileSizeEdit->setEnabled(false);

	strcpy(filename1, "");
	strcpy(filename2, "");
}

CollideAgglomeratesWidget::~CollideAgglomeratesWidget(void)
{
}

void CollideAgglomeratesWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);

	painter.setFont(QFont("Helvetica", 10, QFont::Bold));
	painter.drawText(10, 30, tr("Impact speed (in cm/s):"));
	painter.drawText(10, 70, tr("Impact parameter:"));
	painter.drawText(10, 110, tr("Initial separation:"));
	painter.drawText(10, 440, tr("Number of projectiles:"));
	painter.drawText(10, 480, tr("PCA projectile size:"));
}

void CollideAgglomeratesWidget::selectFirstAgglomerate()
{
	QString temp = QFileDialog::getOpenFileName(this, tr("Select First Agglomerate"), file_path, tr("Particle Files (*.dat)"));

	if(!temp.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(temp).path();
        strcpy(filename1, temp.toLatin1().data());

		selectFirstAgglomerateButton->setText(QFileInfo(temp).fileName());
	}
}

void CollideAgglomeratesWidget::selectSecondAgglomerate()
{
	QString temp = QFileDialog::getOpenFileName(this, tr("Select Second Agglomerate"), file_path, tr("Particle Files (*.dat)"));

	if(!temp.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(temp).path();
        strcpy(filename2, temp.toLatin1().data());

		selectSecondAgglomerateButton->setText(QFileInfo(temp).fileName());
	}
}

void CollideAgglomeratesWidget::setRandomOrientation1(int value)
{
	agglomerates_collision_random_orientation1 = (Qt::CheckState)value;
}

void CollideAgglomeratesWidget::setRandomOrientation2(int value)
{
	agglomerates_collision_random_orientation2 = (Qt::CheckState)value;
}

void CollideAgglomeratesWidget::setMinimizeDistance(int value)
{
	agglomerates_collision_minimize_distance = (Qt::CheckState)value;
}

void CollideAgglomeratesWidget::setHitAndStick(int value)
{
	if((Qt::CheckState)value == Qt::Checked)
	{
		impactSpeedEdit->setEnabled(false);
		initialSeparationEdit->setEnabled(false);
	}
	else
	{
		impactSpeedEdit->setEnabled(true);
		initialSeparationEdit->setEnabled(true);
	}
}

void CollideAgglomeratesWidget::setImpactParameter(double value)
{
	agglomerates_collision_impact_parameter = value;
}

void CollideAgglomeratesWidget::setImpactSpeed(double value)
{
	agglomerates_collision_impact_speed = value;
}

void CollideAgglomeratesWidget::setInitialSeparation(double value)
{
	agglomerates_collision_initial_separation = value;
}

void CollideAgglomeratesWidget::setProjectileCount(int value)
{
	agglomerates_collision_projectile_count = (unsigned int)value;
}

void CollideAgglomeratesWidget::setProjectileSize(int value)
{
	agglomerates_collision_projectile_size = (unsigned int)value;
}

void CollideAgglomeratesWidget::setUsePCAProjectiles(int value)
{
	if((Qt::CheckState)value == Qt::Checked)
		projectileSizeEdit->setEnabled(true);
	else
		projectileSizeEdit->setEnabled(false);

	agglomerates_collision_use_PCA_projectiles = (Qt::CheckState)value;
}

void CollideAgglomeratesWidget::setRandomImpactParameter(int value)
{
	agglomerates_collision_random_impact_parameter = (Qt::CheckState)value;
}

void CollideAgglomeratesWidget::initCollision()
{
	if(hitAndStickCheckBox->checkState() == Qt::Checked)
		emit signalHitAndStick(filename1, filename2, agglomerates_collision_random_orientation1, agglomerates_collision_random_orientation2, agglomerates_collision_impact_parameter);
	else
		emit signalCollideAgglomerates(filename1, filename2, agglomerates_collision_random_orientation1, agglomerates_collision_random_orientation2, agglomerates_collision_impact_speed, agglomerates_collision_impact_parameter, agglomerates_collision_initial_separation + 0.001, agglomerates_collision_minimize_distance);
}

void CollideAgglomeratesWidget::initMultiCollision()
{
	if(agglomerates_collision_use_PCA_projectiles)
		emit signalImpactMultiPCAProjectiles(filename1, agglomerates_collision_projectile_count, agglomerates_collision_projectile_size, agglomerates_collision_impact_speed, agglomerates_collision_impact_parameter, agglomerates_collision_random_impact_parameter, 1e-4 * (agglomerates_collision_initial_separation + 0.001), agglomerates_collision_minimize_distance);
	else
		emit signalImpactMultiProjectiles(filename1, filename2, agglomerates_collision_projectile_count, agglomerates_collision_impact_speed, agglomerates_collision_impact_parameter, agglomerates_collision_random_impact_parameter, 1e-4 * (agglomerates_collision_initial_separation + 0.001), agglomerates_collision_minimize_distance);
}
