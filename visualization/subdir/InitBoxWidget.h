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

#ifndef INITBOXWIDGET_H
#define INITBOXWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QComboBox>
#include <QFileDialog>
#include <QString>

class InitBoxWidget : public QWidget
{
	Q_OBJECT

public:
	InitBoxWidget();
	~InitBoxWidget(void);

	char filename[300];

	QPushButton *selectAgglomerateButton;
	QCheckBox* useCurrentAgglomerateCheckBox;
	QCheckBox* randomOrientationCheckBox;
	QCheckBox* movingBottomWallCheckBox;
	QCheckBox* sideWallsCheckBox;

	QDoubleSpinBox* wallSpeedEdit;
	QDoubleSpinBox* dynamicPressureEdit;
	QDoubleSpinBox* stopFillingFactorEdit;
	QDoubleSpinBox* stopDissipationFactorEdit;
	QDoubleSpinBox* stopPressureEdit;
	QDoubleSpinBox* sideWallModEdit;

	QComboBox* modeSelectionBox;

	QPushButton *initButton;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void selectAgglomerate();
	void setWallSpeed(double value);
	void setDynamicPressure(double value);
	void useCurrentAgglomerateCheckBoxChanged(int value);
	void setRandomOrientation(int value);
	void setMovingBottomWall(int value);
	void setSideWalls(int value);
	void modeChanged(int value);
	void init();
	void setStopFillingFactor(double value);
	void setStopDissipationFactor(double value);
	void setPenetrationDepth(double value);
	void setSideWallMod(double value);

signals:
	void signalInitCompression(const char*, bool, bool, bool, double, double, double);
	void signalInitCompressionRelaxation(const char*, bool, double, double, double);
	void signalInitShockwave(const char*, bool, double, double, double);
	void signalInitDynamicCompression(const char*, bool, double, double, double);
	void signalInitOpenBox(const char*, bool, bool, double);
};

#endif

