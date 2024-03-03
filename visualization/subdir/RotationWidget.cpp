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

#include "RotationWidget.h"

extern double rotation_axis_x;
extern double rotation_axis_y;
extern double rotation_axis_z;
extern double rotation_angle;

RotationWidget::RotationWidget()
{
	resize(300, 310);
	setMinimumWidth(300);
	setMaximumWidth(300);
	setMinimumHeight(200);
	setMaximumHeight(200);

	setWindowTitle(tr("Rotate Aggregate"));

	rotationAngleEdit = new QDoubleSpinBox(this);
	rotationAngleEdit->setGeometry(200, 20, 90, 30);
	rotationAngleEdit->setRange(-360, 360);
	rotationAngleEdit->setValue(rotation_angle);
	rotationAngleEdit->setSingleStep(45);
	connect(rotationAngleEdit, SIGNAL(valueChanged(double)), this, SLOT(rotationAngleChanged(double)));

	rotationAxisXEdit = new QDoubleSpinBox(this);
	rotationAxisXEdit->setGeometry(130, 60, 50, 30);
	rotationAxisXEdit->setRange(-1, 1);
	rotationAxisXEdit->setValue(rotation_axis_x);
	rotationAxisXEdit->setSingleStep(0.1);
	connect(rotationAxisXEdit, SIGNAL(valueChanged(double)), this, SLOT(rotationAxisXChanged(double)));

	rotationAxisYEdit = new QDoubleSpinBox(this);
	rotationAxisYEdit->setGeometry(185, 60, 50, 30);
	rotationAxisYEdit->setRange(-1, 1);
	rotationAxisYEdit->setValue(rotation_axis_y);
	rotationAxisYEdit->setSingleStep(0.1);
	connect(rotationAxisYEdit, SIGNAL(valueChanged(double)), this, SLOT(rotationAxisYChanged(double)));

	rotationAxisZEdit = new QDoubleSpinBox(this);
	rotationAxisZEdit->setGeometry(240, 60, 50, 30);
	rotationAxisZEdit->setRange(-1, 1);
	rotationAxisZEdit->setValue(rotation_axis_z);
	rotationAxisZEdit->setSingleStep(0.1);
	connect(rotationAxisZEdit, SIGNAL(valueChanged(double)), this, SLOT(rotationAxisZChanged(double)));

	rotateXAxisButton = new QPushButton(tr("X Axis"), this);
	rotateXAxisButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	rotateXAxisButton->setGeometry(20, 100, 80, 30);
	connect(rotateXAxisButton, SIGNAL(clicked(bool)), this, SLOT(rotateXAxis()));

	rotateYAxisButton = new QPushButton(tr("Y Axis"), this);
	rotateYAxisButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	rotateYAxisButton->setGeometry(110, 100, 80, 30);
	connect(rotateYAxisButton, SIGNAL(clicked(bool)), this, SLOT(rotateYAxis()));

	rotateZAxisButton = new QPushButton(tr("Z Axis"), this);
	rotateZAxisButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	rotateZAxisButton->setGeometry(200, 100, 80, 30);
	connect(rotateZAxisButton, SIGNAL(clicked(bool)), this, SLOT(rotateZAxis()));

	rotateAggregateButton = new QPushButton(tr("Rotate Aggregate"), this);
	rotateAggregateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	rotateAggregateButton->setGeometry(10, 160, 280, 30);
	connect(rotateAggregateButton, SIGNAL(clicked(bool)), this, SLOT(rotateAggregate()));
}

RotationWidget::~RotationWidget(void)
{
}

void RotationWidget::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 40, tr("Angle (in degree):"));
	painter.drawText(10, 80, tr("Axis (x,y,z):"));
}

void RotationWidget::rotationAxisXChanged(double value)
{
	rotation_axis_x = value;
}

void RotationWidget::rotationAxisYChanged(double value)
{
	rotation_axis_y = value;
}

void RotationWidget::rotationAxisZChanged(double value)
{
	rotation_axis_z = value;
}

void RotationWidget::rotationAngleChanged(double value)
{
	rotation_angle = value;
}

void RotationWidget::rotateXAxis()
{
	rotationAxisXEdit->setValue(1.0);
	rotationAxisYEdit->setValue(0);
	rotationAxisZEdit->setValue(0);
}

void RotationWidget::rotateYAxis()
{
	rotationAxisXEdit->setValue(0);
	rotationAxisYEdit->setValue(1.0);
	rotationAxisZEdit->setValue(0);
}

void RotationWidget::rotateZAxis()
{
	rotationAxisXEdit->setValue(0);
	rotationAxisYEdit->setValue(0);
	rotationAxisZEdit->setValue(1.0);
}

void RotationWidget::rotateAggregate()
{
	emit signalRotateAggregate(rotation_axis_x, rotation_axis_y, rotation_axis_z, rotation_angle);
}
