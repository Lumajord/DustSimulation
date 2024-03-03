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

#ifndef INITAGGREGATEWIDGET_H
#define INITAGGREGATEWIDGET_H

#include "Constants.h"

#include <QtGui>
#include <QWidget>
#include <QSpinBox>
#include <QPushButton>
#include <QComboBox>

class InitAggregateWidget : public QWidget
{
	Q_OBJECT

public:
    InitAggregateWidget();
    ~InitAggregateWidget(void);

	QSpinBox* numberOfParticlesEdit;

	QPushButton *initNeedleButton;

	QDoubleSpinBox* migrationProbability1Edit;
	QDoubleSpinBox* migrationProbability2Edit;
	QComboBox* bamMethodSelectionBox;
	QPushButton *initBAMAggregateButton;

	QDoubleSpinBox* chainRatioEdit;
	QSpinBox *minChainLengthEdit;
	QSpinBox *maxChainLengthEdit;
	QPushButton *initFractalAggregateButton;

	QDoubleSpinBox* treeContacts0Edit;
	QDoubleSpinBox* treeContacts1Edit;
	QDoubleSpinBox* treeContacts2Edit;
	QDoubleSpinBox* treeContacts3Edit;
	QDoubleSpinBox* treeContacts4Edit;
	QPushButton *initTreeAggregateButton;


protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void numberOfParticlesChanged(int value);
	void migrationProbability1Changed(double value);
	void migrationProbability2Changed(double value);
	void bamSelectionMethodChanged(int value);
	
	void chainRatioChanged(double value);
	void minChainLengthChanged(int value);
	void maxChainLengthChanged(int value);
	void treeContacts0Changed(double value);
	void treeContacts1Changed(double value);
	void treeContacts2Changed(double value);
	void treeContacts3Changed(double value);
	void treeContacts4Changed(double value);

	void initNeedleDeposition();
	void initBAMAggregate();
	void initFractalAggregate();
	void initTreeAggregate();


signals:
	void signalInitNeedleDeposition(unsigned int);
	void signalInitBAMAggregate(unsigned int, double, double, BAMSelectionMethod);
	void signalInitFractalAggregate(unsigned int, double, double, double, unsigned int, unsigned int);
	void signalInitTreeAggregate(unsigned int, double, double, double, double, double, double, double, double);	
};

#endif

