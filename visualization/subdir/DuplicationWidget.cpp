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

#include "DuplicationWidget.h"

extern int duplicate_x_duplications;
extern int duplicate_y_duplications;
extern int duplicate_z_duplications;
extern Qt::CheckState duplicate_random_orientation;
extern Qt::CheckState duplicate_mirror;
extern QString file_path;

DuplicationWidget::DuplicationWidget()
{
	resize(260, 260);
	setMinimumHeight(270);
	setMaximumHeight(270);
	setMinimumWidth(260);
	setMaximumWidth(260);

	setWindowTitle (tr("Duplicate Aggregate"));

	xDuplicationsEdit = new QSpinBox(this);
	xDuplicationsEdit->setGeometry(160, 10, 90, 30); 
	xDuplicationsEdit->setRange(0, 1000);
	xDuplicationsEdit->setValue(duplicate_x_duplications);
	xDuplicationsEdit->setSingleStep(1);
	connect(xDuplicationsEdit, SIGNAL(valueChanged(int)), this, SLOT(xDuplicationsChanged(int)));

	yDuplicationsEdit = new QSpinBox(this);
	yDuplicationsEdit->setGeometry(160, 50, 90, 30); 
	yDuplicationsEdit->setRange(0, 1000);
	yDuplicationsEdit->setValue(duplicate_y_duplications);
	yDuplicationsEdit->setSingleStep(1);
	connect(yDuplicationsEdit, SIGNAL(valueChanged(int)), this, SLOT(yDuplicationsChanged(int)));

	zDuplicationsEdit = new QSpinBox(this);
	zDuplicationsEdit->setGeometry(160, 90, 90, 30); 
	zDuplicationsEdit->setRange(0, 1000);
	zDuplicationsEdit->setValue(duplicate_z_duplications);
	zDuplicationsEdit->setSingleStep(1);
	connect(zDuplicationsEdit, SIGNAL(valueChanged(int)), this, SLOT(zDuplicationsChanged(int)));

	randomOrientationCheckBox = new QCheckBox("Random Orientation", this);
	randomOrientationCheckBox->setGeometry(10, 160, 280, 30);
	randomOrientationCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomOrientationCheckBox->setCheckState(duplicate_random_orientation);
	connect(randomOrientationCheckBox, SIGNAL(stateChanged(int)), this, SLOT(randomOrientationChanged(int)));

	mirrorCheckBox = new QCheckBox("Mirror aggregates", this);
	mirrorCheckBox->setGeometry(10, 190, 280, 30);
	mirrorCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	mirrorCheckBox->setCheckState(duplicate_mirror);
	connect(mirrorCheckBox, SIGNAL(stateChanged(int)), this, SLOT(mirrorChanged(int)));

	duplicateButton = new QPushButton(tr("Duplicate"), this);
	duplicateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	duplicateButton->setGeometry(10, 230, 240, 30);
	connect(duplicateButton, SIGNAL(clicked(bool)), this, SLOT(duplicate()));
}

DuplicationWidget::~DuplicationWidget(void)
{
}

void DuplicationWidget::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 30, tr("x-Duplications:"));
	painter.drawText(10, 70, tr("y-Duplications:"));
	painter.drawText(10, 110, tr("z-Duplications:"));
}

void DuplicationWidget::xDuplicationsChanged(int value)
{
	duplicate_x_duplications = value;
}

void DuplicationWidget::yDuplicationsChanged(int value)
{
	duplicate_y_duplications = value;
}

void DuplicationWidget::zDuplicationsChanged(int value)
{
	duplicate_z_duplications = value;
}

void DuplicationWidget::randomOrientationChanged(int value)
{
	duplicate_random_orientation = (Qt::CheckState)value;
}

void DuplicationWidget::mirrorChanged(int value)
{
	duplicate_mirror = (Qt::CheckState)value;
}

void DuplicationWidget::duplicate()
{
	emit signalDuplicate(duplicate_x_duplications, duplicate_y_duplications, duplicate_z_duplications, duplicate_mirror, duplicate_random_orientation);
}
