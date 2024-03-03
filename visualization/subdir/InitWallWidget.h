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

#ifndef INITWALLWIDGET_H
#define INITWALLWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QDoubleSpinBox>
#include <QCheckBox>
#include <QPushButton>
#include <QFileDialog>
#include <QString>

class WallCollisionWidget : public QWidget
{
	Q_OBJECT

public:
	WallCollisionWidget();
	~WallCollisionWidget(void);

	QDoubleSpinBox *impactSpeedEdit;
	QSpinBox *impactAngleEdit;
	QDoubleSpinBox *impactDistanceEdit;
	QCheckBox *randomOrientationCheckBox;
	QCheckBox *useCurrentAgglomerateCheckBox;

	QPushButton *selectAgglomerateButton;
	QPushButton *initButton;

	char filename[300];

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void selectAgglomerate();
	void setRandomOrientation(int value);
	void useCurrentAgglomerateChanged(int value);
	void setImpactAngle(int value);
	void setImpactSpeed(double value);
	void setImpactDistance(double value);
	void initCollision();

signals:
	void signalCollideAgglomerateWithWall(const char*, double, int, double, bool);
};

#endif

