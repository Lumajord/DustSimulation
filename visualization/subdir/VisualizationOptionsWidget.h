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

#ifndef VISUALIZATIONOPTIONSWIDGET_H
#define VISUALIZATIONOPTIONSWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QComboBox>

class VisualizationOptionsWidget :
	public QWidget
{
	Q_OBJECT

public:
	VisualizationOptionsWidget(void);
	~VisualizationOptionsWidget(void);

	QComboBox *displayModeSelectionBox;

	QCheckBox *trackFragmentsCheckBox;
	QCheckBox *displayDepthCheckBox;

	QSpinBox *neighbourSearchDistEdit;
	QSpinBox *neighbourMaxParticlesEdit;
	QSpinBox *neighbourMinOffsetEdit;
	QSpinBox *neighbourMaxOffsetEdit;
	QDoubleSpinBox *densityMaxFillingFactorEdit;
	QDoubleSpinBox *dislocationMinValueEdit;
	QDoubleSpinBox *dislocationMaxValueEdit;

public slots:
	void displayModeChanged(int mode);
	void displayDepthChanged(int state);
	void trackFragmentsChanged(int state);
	void neighbourSearchDistChanged(int value);
	void neighbourMaxParticlesChanged(int value);
	void neighbourMinOffsetChanged(int value);
	void neighbourMaxOffsetChanged(int value);
	void densityMaxFillingFactorChanged(double value);
	void dislocationMinValueChanged(double value);
	void dislocationMaxValueChanged(double value);

signals:
	void signalUpdateColors();
	
protected:
	void paintEvent(QPaintEvent *event);
};

#endif
