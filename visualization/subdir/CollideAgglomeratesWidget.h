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

#ifndef COLLIDEAGGLOMERATESWIDGET_H
#define COLLIDEAGGLOMERATESWIDGET_H

#include <QtGui>
#include <QWidget>
#include <QPushButton>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QString>

class CollideAgglomeratesWidget :	public QWidget
{
	Q_OBJECT

public:
	CollideAgglomeratesWidget(void);
	~CollideAgglomeratesWidget(void);

	QPushButton *selectFirstAgglomerateButton;
	QPushButton *selectSecondAgglomerateButton;
	QPushButton *initCollisionButton;
	QPushButton *initMultiCollisionButton;

	QDoubleSpinBox *impactSpeedEdit;
	QDoubleSpinBox *impactParameterEdit;
	QDoubleSpinBox *initialSeparationEdit;
	QCheckBox *randomOrientation1CheckBox;
	QCheckBox *randomOrientation2CheckBox;
	QCheckBox *minimizeDistanceCheckBox;
	QCheckBox *hitAndStickCheckBox;
	QSpinBox *projectileCountEdit;
	QSpinBox *projectileSizeEdit;
	QCheckBox *randomImpactParameterCheckBox;
	QCheckBox *usePCAProjectilesCheckBox;

	char filename1[300];
	char filename2[300];

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void selectFirstAgglomerate();
	void selectSecondAgglomerate();	
	void setRandomOrientation1(int value);
	void setRandomOrientation2(int value);
	void setMinimizeDistance(int value);
	void setHitAndStick(int value);
	void setImpactParameter(double value);
	void setImpactSpeed(double value);
	void setInitialSeparation(double value);
	void setProjectileCount(int value);
	void setProjectileSize(int value);
	void setUsePCAProjectiles(int value);
	void setRandomImpactParameter(int value);
	void initCollision();
	void initMultiCollision();

signals:
	void signalCollideAgglomerates(const char*, const char*, bool, bool, double, double, double, bool);
	void signalImpactMultiProjectiles(const char*, const char*, unsigned int, double, double, bool, double, bool);
	void signalImpactMultiPCAProjectiles(const char*, unsigned int, unsigned int, double, double, bool, double, bool);
	void signalHitAndStick(const char*, const char*, bool, bool, double);
};

#endif
