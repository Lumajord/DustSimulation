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

#ifndef SLICEWIDGET_H
#define SLICEWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QCheckBox>
#include <QFileDialog>
#include <QString>

class SliceWidget : public QWidget
{
	Q_OBJECT

public:
	SliceWidget();
	~SliceWidget(void);

	QPushButton *sliceBoxButton;
	QPushButton *sliceSphereButton;
	QPushButton *sliceCylinderButton;
	QPushButton *selectAgglomerateButton;
	QPushButton *sliceTopButton;
	QPushButton *sliceBottomButton;

	char filename[300];

	QDoubleSpinBox *boxXSizeEdit;
	QDoubleSpinBox *boxYSizeEdit;
	QDoubleSpinBox *boxZSizeEdit;
	QDoubleSpinBox *sphereRadiusModifierEdit;
	QDoubleSpinBox *topSliceFactorEdit;

	QCheckBox *randomOrientationCheckBox;
	QCheckBox *sliceCurrentAgglomerateCheckBox;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void selectAgglomerate();
	void boxXSizeChanged(double value);
	void boxYSizeChanged(double value);
	void boxZSizeChanged(double value);
	void topSliceFactorChanged(double value);
	void sphereRadiusModifierChanged(double value);
	void sliceCurrentAgglomerateChanged(int value);
	void randomOrientationChanged(int value);

	void sliceBox();
	void sliceSphere();
	void sliceCylinder();
	void sliceTop();
	void sliceBottom();

signals:
	void signalSliceBox(const char*, double, double, double);
	void signalSliceSphere(const char*, double, bool);
	void signalSliceCylinder(const char*, double);
	void signalSliceTop(const char*, double);
	void signalSliceBottom(const char*, double);
};

#endif

