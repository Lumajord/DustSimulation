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

#include "AgglomerateInfoWidget.h"
#include "SimulationVisualization.h"

extern DisplayInfo sim_info;
extern AgglomerateInfo agg_info;
extern Qt::CheckState agg_info_only_inner_particles;
extern double mass;

AgglomerateInfoWidget::AgglomerateInfoWidget(void)
{
	setMinimumWidth(320);
	setMaximumWidth(320);
	setMinimumHeight(690);
	setMaximumHeight(690);

	setWindowTitle (tr("Agglomerate Properties"));

	onlyInnerParticlesCheckBox = new QCheckBox("Use only inner particles", this);
	onlyInnerParticlesCheckBox->setGeometry(10, 610, 200, 30);
	onlyInnerParticlesCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	onlyInnerParticlesCheckBox->setCheckState(agg_info_only_inner_particles);
	onlyInnerParticlesCheckBox->setToolTip(tr("If checked only the particles within the box are taken into account for the contact histogram"));
	connect(onlyInnerParticlesCheckBox, SIGNAL(stateChanged(int)), this, SLOT(onlyInnerParticlesChanged(int)));

	refreshButton = new QPushButton(tr("Refresh"), this);
	refreshButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	refreshButton->setGeometry(70, 650, 160, 30);

	font = new QFont("Helvetica", 10, QFont::Bold);
	fontMetrics = new QFontMetrics(*font);
}

AgglomerateInfoWidget::~AgglomerateInfoWidget(void)
{
	delete fontMetrics;
	delete font;
}

void AgglomerateInfoWidget::paintEvent(QPaintEvent *)
{
	char buffer[50];
	int h_pos;
	int width = this->width();
	int edge_dist = 35;
	int edge_dist_units = 30;

	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.setFont(QFont("Helvetica", 10, QFont::Bold));
	painter.drawText(10, 30, tr("Particles:"));
	painter.drawText(10, 60, tr("Total Mass:"));
	painter.drawText(10, 90, tr("Fragments:"));
	painter.drawText(10, 120, tr("Gyration Radius:"));
	painter.drawText(10, 150, tr("Outer Radius:"));
	painter.drawText(10, 180, tr("Density (in g/cm^3):"));
	painter.drawText(10, 210, tr("Box height:"));
	painter.drawText(10, 240, tr("Box base:"));
	painter.drawText(10, 270, tr("Enclosing box filling factor:"));
	painter.drawText(10, 300, tr("Enclosing sphere filling factor:"));

	sprintf(buffer, "%i", sim_info.particles);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 30, buffer);

	if(sim_info.particles > 0)
		sprintf(buffer, "%.3le", sim_info.particles * mass);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 60, buffer);
	painter.drawText(width - edge_dist_units, 60, tr("g"));

	sprintf(buffer, "%i", agg_info.fragments);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 90, buffer);

	sprintf(buffer, "%i", agg_info.fragments);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 90, buffer);

	if(agg_info.gyration_radius > 0)
		sprintf(buffer, "%.3le", 0.01 * agg_info.gyration_radius);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 120, buffer);
	painter.drawText(width - edge_dist_units, 120, tr("m"));

	if(agg_info.outer_radius > 0)
		sprintf(buffer, "%.3le", 0.01 * agg_info.outer_radius);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 150, buffer);
	painter.drawText(width - edge_dist_units, 150, tr("m"));

	if(agg_info.gyration_radius > 0)
	{
		// characteristic radius (Mukai et al., 1992)
		double r = sqrt(5.0 / 3.0) * agg_info.gyration_radius;
		double rho = 3.0 * (double)sim_info.particles * mass / (4.0 * M_PI * r*r*r);
		sprintf(buffer, "%.3le", rho);
	}
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 180, buffer);

	if(agg_info.box_height > 0)
		sprintf(buffer, "%.3le", 0.01 * agg_info.box_height);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 210, buffer);
	painter.drawText(width - edge_dist_units, 210, tr("m"));

	if(agg_info.box_base > 0)
		sprintf(buffer, "%.3le", 0.0001 * agg_info.box_base);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 240, buffer);
	painter.drawText(width - edge_dist_units, 240, tr("m^2"));

	if(agg_info.filling_factor_box > 0)
		sprintf(buffer, "%.3lf", agg_info.filling_factor_box);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 270, buffer);

	if(agg_info.filling_factor_sphere > 0)
		sprintf(buffer, "%.3lf", agg_info.filling_factor_sphere);
	else
		sprintf(buffer, "-");
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 300, buffer);

	// contact histogram
	int max_particles = 0;

	for(int i = 0; i < 13; ++i)
	{
		if(agg_info.contact_histogram[i] > max_particles)
			max_particles = agg_info.contact_histogram[i];
	}

	if(max_particles > 0)
	{
		h_pos = 10;
		
                for(int i = 0; i < 13; ++i)
		{
			int height = (int)(250.0 * (double)agg_info.contact_histogram[i] / (double)max_particles);
            int v_pos = 330 + (250 - height);

			if(i%2==1)
				painter.fillRect(h_pos, v_pos, 25, height, QColor(255, 0, 0));
			else
				painter.fillRect(h_pos, v_pos, 25, height, QColor(0, 150, 0));

			h_pos += 25;
		}
	}

	// draw frame
	painter.drawRect(10, 320, 300, 260 );
        painter.drawText(18, 605, tr("0"));
        painter.drawText(68, 605, tr("2"));
        painter.drawText(118, 605, tr("4"));
        painter.drawText(168, 605, tr("6"));
        painter.drawText(218, 605, tr("8"));
        painter.drawText(263, 605, tr("10"));
        painter.drawText(313, 605, tr("12"));
}

void AgglomerateInfoWidget::onlyInnerParticlesChanged(int state)
{
	agg_info_only_inner_particles = (Qt::CheckState)state;
	emit signalUpdateAgglomerateInfo();
}
