#include "VisualizationOptionsWidget.h"
#include "Constants.h"

extern DisplayMode vis_display_mode;
extern int vis_neighbour_search_dist;
extern int vis_neighbour_max_particles;
extern int vis_neighbour_min_offset;
extern int vis_neighbour_max_offset;
extern double vis_density_max_filling_factor;
extern double vis_dislocation_min_value;
extern double vis_dislocation_max_value;
extern Qt::CheckState display_depth;
extern Qt::CheckState track_fragments;

VisualizationOptionsWidget::VisualizationOptionsWidget(void)
{
	setMinimumHeight(440);
	setMaximumHeight(440);
	setMinimumWidth(260);
	setMaximumWidth(260);
	resize(240, 440);
	setWindowTitle (tr("Visualization Options"));

	int h_pos = width() - 70;

	displayModeSelectionBox = new QComboBox(this);
	displayModeSelectionBox->setGeometry(70, 40, 170, 30);
	displayModeSelectionBox->setFont(QFont("Helvetica", 8, QFont::Bold));
	displayModeSelectionBox->addItem(" Nothing");
	displayModeSelectionBox->addItem(" Relative Particle Speed");
	displayModeSelectionBox->addItem(" Dissipated Energy");
	displayModeSelectionBox->addItem(" Volume Filling Factor");
	displayModeSelectionBox->addItem(" Coordination Number");
	displayModeSelectionBox->addItem(" Particle-Wall Contacts");
	displayModeSelectionBox->addItem(" Velocity-Wall Angle");
	displayModeSelectionBox->addItem(" Fragments");
	displayModeSelectionBox->addItem(" Traveled Distance");
        displayModeSelectionBox->addItem(" 50/50 for collision");
	displayModeSelectionBox->setCurrentIndex(vis_display_mode);
	connect(displayModeSelectionBox, SIGNAL(currentIndexChanged(int)), this, SLOT(displayModeChanged(int)));

	trackFragmentsCheckBox = new QCheckBox("Track Fragments", this);
	trackFragmentsCheckBox->setGeometry(10, 90, 240, 30);
	trackFragmentsCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	trackFragmentsCheckBox->setCheckState(track_fragments);
	trackFragmentsCheckBox->setToolTip(tr("If enabled, fragments will be updated after every integration step.\nAttention: This slows down the simulation significantly!"));
	connect(trackFragmentsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(trackFragmentsChanged(int)));

	displayDepthCheckBox = new QCheckBox("Depth", this);
	displayDepthCheckBox->setGeometry(10, 120, 200, 30);
	displayDepthCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	displayDepthCheckBox->setCheckState(display_depth);
	displayDepthCheckBox->setToolTip(tr("If enabled, particles in areas with higher density appear darker"));
	connect(displayDepthCheckBox, SIGNAL(stateChanged(int)), this, SLOT(displayDepthChanged(int)));

	neighbourSearchDistEdit = new QSpinBox(this);
	neighbourSearchDistEdit->setGeometry(h_pos, 160, 60, 30); 
	neighbourSearchDistEdit->setRange(0, 5);
	neighbourSearchDistEdit->setValue(vis_neighbour_search_dist);
	neighbourSearchDistEdit->setSingleStep(1);
	connect(neighbourSearchDistEdit, SIGNAL(valueChanged(int)), this, SLOT(neighbourSearchDistChanged(int)));

	neighbourMaxParticlesEdit = new QSpinBox(this);
	neighbourMaxParticlesEdit->setGeometry(h_pos, 200, 60, 30); 
	neighbourMaxParticlesEdit->setRange(0, 100000);
	neighbourMaxParticlesEdit->setValue(vis_neighbour_max_particles);
	neighbourMaxParticlesEdit->setSingleStep(10);
	connect(neighbourMaxParticlesEdit, SIGNAL(valueChanged(int)), this, SLOT(neighbourMaxParticlesChanged(int)));

	neighbourMinOffsetEdit = new QSpinBox(this);
	neighbourMinOffsetEdit->setGeometry(h_pos, 240, 60, 30); 
	neighbourMinOffsetEdit->setRange(0, 10000);
	neighbourMinOffsetEdit->setValue(vis_neighbour_min_offset);
	neighbourMinOffsetEdit->setSingleStep(10);
	connect(neighbourMinOffsetEdit, SIGNAL(valueChanged(int)), this, SLOT(neighbourMinOffsetChanged(int)));

	neighbourMaxOffsetEdit = new QSpinBox(this);
	neighbourMaxOffsetEdit->setGeometry(h_pos, 280, 60, 30); 
	neighbourMaxOffsetEdit->setRange(0, 10000);
	neighbourMaxOffsetEdit->setValue(vis_neighbour_max_offset);
	neighbourMaxOffsetEdit->setSingleStep(10);
	connect(neighbourMaxOffsetEdit, SIGNAL(valueChanged(int)), this, SLOT(neighbourMaxOffsetChanged(int)));

	densityMaxFillingFactorEdit = new QDoubleSpinBox(this);
	densityMaxFillingFactorEdit->setGeometry(h_pos, 320, 60, 30); 
	densityMaxFillingFactorEdit->setRange(0.01, 0.75);
	densityMaxFillingFactorEdit->setValue(vis_density_max_filling_factor);
	densityMaxFillingFactorEdit->setSingleStep(0.05);
	connect(densityMaxFillingFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(densityMaxFillingFactorChanged(double)));

	dislocationMinValueEdit = new QDoubleSpinBox(this);
	dislocationMinValueEdit->setGeometry(h_pos, 360, 60, 30); 
	dislocationMinValueEdit->setRange(0.0, 100000);
	dislocationMinValueEdit->setValue(vis_dislocation_min_value);
	dislocationMinValueEdit->setSingleStep(1);
	connect(dislocationMinValueEdit, SIGNAL(valueChanged(double)), this, SLOT(dislocationMinValueChanged(double)));

	dislocationMaxValueEdit = new QDoubleSpinBox(this);
	dislocationMaxValueEdit->setGeometry(h_pos, 400, 60, 30); 
	dislocationMaxValueEdit->setRange(0.0, 100000);
	dislocationMaxValueEdit->setValue(vis_dislocation_max_value);
	dislocationMaxValueEdit->setSingleStep(1);
	connect(dislocationMaxValueEdit, SIGNAL(valueChanged(double)), this, SLOT(dislocationMaxValueChanged(double)));

	displayDepthChanged(display_depth);
}

VisualizationOptionsWidget::~VisualizationOptionsWidget(void)
{
}

void VisualizationOptionsWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 30, tr("Select visualization mode:"));
	painter.drawText(10, 180, tr("Neighbour search dist:"));
	painter.drawText(10, 220, tr("Cap neighbours  at:"));
	painter.drawText(10, 260, tr("Min. neighbour offset:"));
	painter.drawText(10, 300, tr("Max. neighbour offset:"));
	painter.drawText(10, 340, tr("Max. filling factor:"));
	painter.drawText(10, 380, tr("Min. dislocation (in delta_0):"));
	painter.drawText(10, 420, tr("Max. dislocation (in r_p):"));
}

void VisualizationOptionsWidget::displayModeChanged(int mode)
{
	vis_display_mode = (DisplayMode)mode;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::displayDepthChanged(int state)
{
	display_depth = (Qt::CheckState)state;

	if(state == Qt::Checked)
	{
		neighbourSearchDistEdit->setEnabled(true);
		neighbourMaxParticlesEdit->setEnabled(true);
		neighbourMinOffsetEdit->setEnabled(true);
		neighbourMaxOffsetEdit->setEnabled(true);
	}
	else
	{
		neighbourSearchDistEdit->setEnabled(false);
		neighbourMaxParticlesEdit->setEnabled(false);
		neighbourMinOffsetEdit->setEnabled(false);
		neighbourMaxOffsetEdit->setEnabled(false);
	}

	emit signalUpdateColors();
}

void VisualizationOptionsWidget::trackFragmentsChanged(int state)
{
	track_fragments = (Qt::CheckState)state;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::neighbourSearchDistChanged(int value)
{
	vis_neighbour_search_dist = value;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::neighbourMaxParticlesChanged(int value)
{
	vis_neighbour_max_particles = value;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::neighbourMinOffsetChanged(int value)
{
	vis_neighbour_min_offset = value;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::neighbourMaxOffsetChanged(int value)
{
	vis_neighbour_max_offset = value;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::densityMaxFillingFactorChanged(double value)
{
	vis_density_max_filling_factor = (float)value;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::dislocationMinValueChanged(double value)
{
	vis_dislocation_min_value = value;
	emit signalUpdateColors();
}

void VisualizationOptionsWidget::dislocationMaxValueChanged(double value)
{
	vis_dislocation_max_value = value;
	emit signalUpdateColors();
}
