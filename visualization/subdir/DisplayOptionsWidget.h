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

#ifndef DISPLAYOPTIONSWIDGET_H
#define DISPLAYOPTIONSWIDGET_H

#include <QtGui>
#include <QWidget>
#include <QPushButton>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QColorDialog>

class DisplayOptionsWidget :
	public QWidget
{
	Q_OBJECT

public:
	DisplayOptionsWidget(void);
	~DisplayOptionsWidget(void);

	QPushButton *selectBackgroundColorButton;
	QPushButton *selectBorderColorButton;
	QPushButton *selectWallColorButton;
	QPushButton *selectTextColorButton;

	QCheckBox *displayInfoCheckBox;
	QCheckBox *displayPressureCheckBox;
	QCheckBox *displayForceCheckBox;
	QCheckBox *displayChangedContactsCheckBox;
	QCheckBox *drawParticlesCheckBox;
	QCheckBox *drawBordersCheckBox;
	QCheckBox *drawKeyCheckBox;

	QDoubleSpinBox *borderDistanceFactorEdit;
	QDoubleSpinBox *borderMinDepthEdit;

	QColor q_background_color;
	QColor q_border_color;
	QColor q_wall_color;
	QColor q_text_color;

public slots:
	void getBackgroundColor();
	void getBorderColor();
	void getWallColor();
	void getTextColor();
	void borderMinDepthChanged(double value);
	void borderDistanceFactorChanged(double value);
	void drawBordersChanged(int state);
	void drawParticlesChanged(int state);
	void displayInfoChanged(int state);
	void displayPressureChanged(int state);
	void displayForceChanged(int state);
	void displayChangedContactsChanged(int state);
	void drawKeyChanged(int state);

signals:
	void signalUpdateView();
	void signalUpdateShaders();
	void signalUpdateColors();
	
protected:
	void paintEvent(QPaintEvent *event);
};

#endif
