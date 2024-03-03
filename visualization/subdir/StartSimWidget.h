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

#ifndef STARTSIMWIDGET_H
#define STARTSIMWIDGET_H


#include <QtGui>
#include <QPushButton>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>

class StartSimWidget : public QWidget
{
	Q_OBJECT

public:

	StartSimWidget();
	~StartSimWidget(void);

	QPushButton *startSimButton;

	QDoubleSpinBox *timelengthEdit;
	QDoubleSpinBox *timestepEdit;

	QSpinBox *printEnergyIntervalEdit;
	QSpinBox *printPositionsIntervalEdit;

	QComboBox* unitOfTimeDisplaySelectionBox;
	QSpinBox *averagingStepsEdit;
	QCheckBox *followCMSCheckBox;

	QCheckBox *videoModeCheckBox;
	QSpinBox *takeScreenshotIntervalEdit;

	QCheckBox *useGravityCheckBox;
	QDoubleSpinBox *gravityModifierEdit;
	QDoubleSpinBox *wallInertiaEdit;

	QCheckBox *useAzimuthalAccelerationCheckBox;
	QDoubleSpinBox *azimuthalAccelerationEdit;

	QCheckBox *useDampingCheckBox;

public slots:
	void checkBoxChanged(int state);
	void followCMSChanged(int state);
	void setPrintEnergyInterval(int value);
	void setPrintPositionsInterval(int value);
	void setTakeScreenshotInterval(int value);
	void setTimestep(double value);
	void setDuration(double value);
	void unitOfTimeChanged(int);
	void setAveragingSteps(int);
	void useGravityChanged(int);
	void gravityModifierChanged(double);
	void wallInertiaChanged(double);
	void startSim();
	void useAzimuthalAccelerationChanged(int);
	void angularAccelerationChanged(double);
	void enableDampingChanged(int);

protected:
	void paintEvent(QPaintEvent *event);

signals:
	void signalStartSim(double, double, int, int, int, int, bool, bool);
	void signalUpdateView();
};

#endif
