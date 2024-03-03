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

#ifndef MODIFYWALLWIDGET_H
#define MODIFYWALLWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QPushButton>

class SimulationThread;

class ModifyWallWidget : public QWidget
{
	Q_OBJECT

public:
	ModifyWallWidget(SimulationThread *sim_thread);
	~ModifyWallWidget(void);

	QComboBox* wallSelectionBox;
	QDoubleSpinBox *wallSpeedBox;
	QDoubleSpinBox *dynamicPressureBox;
	QDoubleSpinBox *alphaBox;
	QDoubleSpinBox *compressionModifierBox;
	QDoubleSpinBox *rollingModifierBox;
	QDoubleSpinBox *slidingModifierBox;

	QPushButton *setWallPropertiesButton;

	SimulationThread *sim_thread;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void setupWallSelectionBox(int);
	void wallSelectionChanged(int);
	void setWallProperties();

signals:
	void signalSetWallProperties(int, double, double, double, double, double, double);
};

#endif

