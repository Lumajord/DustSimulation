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

#ifndef ROTATIONWIDGET_H
#define ROTATIONWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QPushButton>
#include <QDoubleSpinBox>

class RotationWidget : public QWidget
{
	Q_OBJECT

public:
	RotationWidget();
	~RotationWidget(void);

	QPushButton *rotateXAxisButton;
	QPushButton *rotateYAxisButton;
	QPushButton *rotateZAxisButton;
	QPushButton *rotateAggregateButton;

	QDoubleSpinBox *rotationAxisXEdit;
	QDoubleSpinBox *rotationAxisYEdit;
	QDoubleSpinBox *rotationAxisZEdit;
	QDoubleSpinBox *rotationAngleEdit;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void rotationAxisXChanged(double value);
	void rotationAxisYChanged(double value);
	void rotationAxisZChanged(double value);
	void rotationAngleChanged(double value);
	void rotateXAxis();
	void rotateYAxis();
	void rotateZAxis();
	void rotateAggregate();

signals:
	void signalRotateAggregate(double, double, double, double);
};

#endif

