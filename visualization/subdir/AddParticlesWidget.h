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

#ifndef ADDPARTICLESWIDGET_H
#define ADDPARTICLESWIDGET_H

#include "Constants.h"
#include <QtGui>
#include <QWidget>
#include <QSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>
#include <QFileDialog>

class AddParticlesWidget : public QWidget
{
	Q_OBJECT

public:
	AddParticlesWidget();
	~AddParticlesWidget(void);

	char filename[300];

	QSpinBox* numberOfParticlesEdit;
	QDoubleSpinBox* migrationProbability1Edit;
	QDoubleSpinBox* migrationProbability2Edit;
	QComboBox* bamMethodSelectionBox;
	QPushButton *selectAgglomerateButton;
	QCheckBox *useCurrentAgglomerateCheckBox;
	QPushButton *addBAMButton;
	QSpinBox* numberOfChainsEdit;
	QSpinBox* minChainLengthEdit;
	QSpinBox* maxChainLengthEdit;
	QPushButton *addFractalChainsButton;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void particleNumberChanged(int);
	void migrationProbability1Changed(double);
	void migrationProbability2Changed(double);
	void bamSelectionMethodChanged(int value);
	void selectAgglomerate();
	void useCurrentAgglomerateChanged(int);
	void numberOfChainsChanged(int);
	void minChainLengthChanged(int);
	void maxChainLengthChanged(int);
	void addBAMParticles();
	void addFractalChains();


signals:
	void signalAddBAMParticles(unsigned int, double, double, const char*, BAMSelectionMethod);
	void signalAddFractalChains(const char*, unsigned int, unsigned int, unsigned int);
};

#endif

