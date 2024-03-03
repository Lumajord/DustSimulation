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

#ifndef INITCAKEWIDGET_H
#define INITCAKEWIDGET_H


#include "Constants.h"

#include <QtGui>
#include <QWidget>
#include <QSpinBox>
#include <QPushButton>
#include <QComboBox>

class InitCakeWidget : public QWidget
{
	Q_OBJECT

public:
	InitCakeWidget();
	~InitCakeWidget(void);

	QSpinBox* xParticlesEdit;
	QSpinBox* yParticlesEdit;
	QSpinBox* zParticlesEdit;
	QDoubleSpinBox* fillingFactorEdit;
	QPushButton *initSimpleCubicButton;
	QPushButton *initHexagonalButton;

	QSpinBox* numberOfParticlesEdit;
	QDoubleSpinBox* migrationProbability1Edit;
	QDoubleSpinBox* migrationProbability2Edit;
	QComboBox* bamMethodSelectionBox;
	QDoubleSpinBox* xSizeEdit;
	QDoubleSpinBox* ySizeEdit;
	QDoubleSpinBox* topSliceFactorEdit;
	QDoubleSpinBox* bottomSliceFactorEdit;
	QPushButton *initBAMCakeButton;
	QPushButton *initRBDCakeButton;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void xParticlesChanged(int value);
	void yParticlesChanged(int value);
	void zParticlesChanged(int value);
	void fillingFactorChanged(double value);

	void particleNumberChanged(int value);
	void migrationProbability1Changed(double value);
	void migrationProbability2Changed(double value);
	void bamSelectionMethodChanged(int value);
	void xSizeChanged(double value);
	void ySizeChanged(double value);
	void topSliceFactorChanged(double value);
	void bottomSliceFactorChanged(double value);

	void initSimpleCubicLattice();
	void initHexagonalLattice();

	void initBAMCake();
	void initRBDCake();

signals:
	void signalInitSimpleCubicLattice(unsigned int, unsigned int, unsigned int, double);
	void signalInitHexagonalLattice(unsigned int, unsigned int, unsigned int, double);
	void signalInitBAMCake(unsigned int, double, double, double, double, double, double, BAMSelectionMethod);
	void signalInitRBDCake(unsigned int, double, double, double, double);
};

#endif

