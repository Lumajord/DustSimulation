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

#ifndef INITSIMWIDGET_H
#define INITSIMWIDGET_H


#include <QtGui>
#include <QWidget>
#include <QPushButton>
#include <QDoubleSpinBox>

class InitSimWidget : public QWidget
{
	Q_OBJECT

public:

	InitSimWidget();
	~InitSimWidget(void);

	QPushButton *initChainButton;
	QPushButton *initClusterButton;

	QSpinBox *particleNumberEdit;
	QSpinBox *particleImpactNumberEdit;
	QSpinBox *targetParticleEdit;
	QDoubleSpinBox *impactSpeedEdit;
	QDoubleSpinBox *angularIrregularityEdit;
	QDoubleSpinBox *fillingFactorEdit;

protected:
	void paintEvent(QPaintEvent *event);

private slots:
		void adjustTargetParticleEditRange(int max_value);
};

#endif

