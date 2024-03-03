#include "MaterialInfoWidget.h"

extern char material_name[200];

extern double particle_radius;
extern double reduced_radius;
extern double density;
extern double mass;
extern double moment_of_inertia;

extern double surface_energy;
extern double nu;
extern double young_mod;
extern double shear_mod;
extern double bulk_mod;

extern double equilibrium_contact_radius;
extern double delta_0;
extern double delta_c;
extern double F_c;
extern double t_c;
extern double T_normal;
extern double v_c;
extern double ENERGY_UNIT;
extern double bulk_sound_speed;

extern double T_vis;
extern double crit_rolling_displacement;
extern double osc_damping_factor;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;

extern double mass;
extern double k_s;

MaterialInfoWidget::MaterialInfoWidget(void)
{
	setMinimumHeight(560);
	setMaximumHeight(560);
	setMinimumWidth(320);
	setMaximumWidth(320);
	resize(320, 410);

	setWindowTitle (tr("Material Info"));

	font = new QFont("Helvetica", 10, QFont::Bold);
	fontMetrics = new QFontMetrics(*font);
}

MaterialInfoWidget::~MaterialInfoWidget(void)
{
	delete fontMetrics;
	delete font;
}

void MaterialInfoWidget::paintEvent(QPaintEvent *)
{
	char buffer[50]; 
	int h_pos;
	int width = this->width();
	int edge_dist = 75;
	int edge_dist_units = 70;
	
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 20, tr("Material:"));

	painter.drawText(10, 50, tr("Particle Radius:"));
	painter.drawText(10, 70, tr("Density:"));
	painter.drawText(10, 90, tr("Mass:"));
	painter.drawText(10, 110, tr("Moment of Inertia:"));

	painter.drawText(10, 140, tr("Surface Energy:"));
	painter.drawText(10, 160, tr("Poisson's Number:"));
	painter.drawText(10, 180, tr("Young's Modulus:"));
	painter.drawText(10, 200, tr("Shear Modulus:"));
	painter.drawText(10, 220, tr("Bulk Modulus:"));

	painter.drawText(10, 250, tr("Viscous Constant:"));
	painter.drawText(10, 270, tr("Normal Damping Factor:"));
	painter.drawText(10, 290, tr("Rolling Modifier:"));
	painter.drawText(10, 310, tr("Sliding Modifier:"));
	painter.drawText(10, 330, tr("Twisting Modifier:"));

	painter.drawText(10, 360, tr("Breaking Force F_c:"));
	painter.drawText(10, 380, tr("Breaking Distance:"));
	painter.drawText(10, 400, tr("Equilibrium Distance:"));
	painter.drawText(10, 420, tr("Sticking Velocity v_c:"));
	painter.drawText(10, 440, tr("Characteristic Time t_c:"));
	painter.drawText(10, 460, tr("Bulk Sound Speed:"));
	painter.drawText(10, 480, tr("Normal Osc. Period:"));
	painter.drawText(10, 500, tr("Sliding Osc. Period:"));

	painter.drawText(10, 530, tr("Energy Unit:"));
	painter.drawText(10, 550, tr("E_roll:"));
	painter.drawText(10, 570, tr("E_break:"));

	h_pos = width - fontMetrics->width(material_name) - edge_dist;
	painter.drawText(h_pos, 20, material_name);

	sprintf(buffer, "%g", 0.01 * particle_radius);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 50, buffer);
	painter.drawText(width - edge_dist_units, 50, tr("m"));

	sprintf(buffer, "%.2lf", density);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 70, buffer);
	painter.drawText(width - edge_dist_units, 70, tr("g/cm^3"));

	sprintf(buffer, "%.3g", mass);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 90, buffer);
	painter.drawText(width - edge_dist_units, 90, tr("g"));

	sprintf(buffer, "%.3g", moment_of_inertia);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 110, buffer);
	painter.drawText(width - edge_dist_units, 110, tr("g cm^2"));

	sprintf(buffer, "%g", surface_energy);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 140, buffer);
	painter.drawText(width - edge_dist_units, 140, tr("mJ/m^2"));

	sprintf(buffer, "%.2lf", nu);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 160, buffer);

	sprintf(buffer, "%.2lf", young_mod * 1e-10);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 180, buffer);
	painter.drawText(width - edge_dist_units, 180, tr("GPa"));

	sprintf(buffer, "%.2lf", shear_mod * 1e-10);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 200, buffer);
	painter.drawText(width - edge_dist_units, 200, tr("GPa"));

	sprintf(buffer, "%.2lf", bulk_mod * 1e-10);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 220, buffer);
	painter.drawText(width - edge_dist_units, 220, tr("GPa"));
	

	sprintf(buffer, "%.3g", T_vis);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 250, buffer);

	sprintf(buffer, "%.3g", osc_damping_factor);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 270, buffer);

	sprintf(buffer, "%.3g", rolling_modifier);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 290, buffer);

	sprintf(buffer, "%.3g", sliding_modifier);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 310, buffer);

	sprintf(buffer, "%.3g", twisting_modifier);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 330, buffer);


	sprintf(buffer, "%.3g", 1e-5*F_c);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 360, buffer);
	painter.drawText(width - edge_dist_units, 360, tr("N"));

	sprintf(buffer, "%.3g", 0.01 * delta_c);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 380, buffer);
	painter.drawText(width - edge_dist_units, 380, tr("m"));

	sprintf(buffer, "%.3g", 0.01 * delta_0);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 400, buffer);
	painter.drawText(width - edge_dist_units, 400, tr("m"));

	sprintf(buffer, "%.3g", 0.01 * v_c);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 420, buffer);
	painter.drawText(width - edge_dist_units, 420, tr("m/s"));

	sprintf(buffer, "%.3g", t_c);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 440, buffer);
	painter.drawText(width - edge_dist_units, 440, tr("s"));

	sprintf(buffer, "%.3g", bulk_sound_speed * 0.01);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 460, buffer);
	painter.drawText(width - edge_dist_units, 460, tr("m/s"));

	sprintf(buffer, "%.3g", T_normal);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 480, buffer);
	painter.drawText(width - edge_dist_units, 440, tr("s"));

	sprintf(buffer, "%.3g", 2.0 * M_PI * sqrt(mass / k_s));
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 500, buffer);
	painter.drawText(width - edge_dist_units, 480, tr("s"));

	sprintf(buffer, "%.3g", F_c * delta_c * 1e-7);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 530, buffer);
	painter.drawText(width - edge_dist_units, 510, tr("J"));

	sprintf(buffer, "%.3g", 6.0 * M_PI*M_PI * surface_energy * reduced_radius * crit_rolling_displacement * 1e-7);
    h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 550, buffer);
	painter.drawText(width - edge_dist_units, 530, tr("J"));

	sprintf(buffer, "%.3g", 1.8 * F_c * delta_c * 1e-7);
	h_pos = width - fontMetrics->width(buffer) - edge_dist;
	painter.drawText(h_pos, 570, buffer);
	painter.drawText(width - edge_dist_units, 550, tr("J"));
}
