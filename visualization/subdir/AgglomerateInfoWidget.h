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

#ifndef AGGLOMERATEINFOWIDGET_H
#define AGGLOMERATEINFOWIDGET_H

#include <QWidget>
#include <QCheckBox>
#include <QPushButton>
#include <QtGui>

class AgglomerateInfoWidget :	public QWidget
{
	Q_OBJECT

public:
	AgglomerateInfoWidget(void);
	~AgglomerateInfoWidget(void);

	// used to determine text width/height
	QFont *font;
	QFontMetrics *fontMetrics;

	QCheckBox *onlyInnerParticlesCheckBox;
	QPushButton *refreshButton;

protected:
	void paintEvent(QPaintEvent *event);


public slots:
	void onlyInnerParticlesChanged(int);

signals:
	void signalUpdateAgglomerateInfo();
	

};

#endif
