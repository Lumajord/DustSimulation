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

#ifndef DUPLICATIONWIDGET_H
#define DUPLICATIONWIDGET_H

#include <QtGui>
#include <QWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QCheckBox>

class DuplicationWidget : public QWidget
{
	Q_OBJECT

public:
	DuplicationWidget();
	~DuplicationWidget(void);

	QPushButton *duplicateButton;

	char file_name[200];

	QSpinBox *xDuplicationsEdit;
	QSpinBox *yDuplicationsEdit;
	QSpinBox *zDuplicationsEdit;

	QCheckBox *randomOrientationCheckBox;
	QCheckBox *mirrorCheckBox;

protected:
	void paintEvent(QPaintEvent *event);

public slots:
	void xDuplicationsChanged(int value);
	void yDuplicationsChanged(int value);
	void zDuplicationsChanged(int value);
	void randomOrientationChanged(int value);
	void mirrorChanged(int value);
	void duplicate();

signals:
	void signalDuplicate(int, int, int, bool, bool);
};

#endif

