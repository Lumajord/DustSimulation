#include "DisplayOptionsWidget.h"
#include "SimulationVisualization.h"

extern float clear_color[3];
extern float border_color[3];
extern float wall_color[3];
extern float text_color[3];
extern float border_distance_factor;
extern float border_min_depth;
extern Qt::CheckState display_text_info;
extern Qt::CheckState display_pressure;
extern Qt::CheckState display_force;
extern Qt::CheckState display_changed_contacts;
extern Qt::CheckState display_key;
extern Qt::CheckState draw_borders;
extern Qt::CheckState draw_particles;

extern DisplayInfo sim_info;

DisplayOptionsWidget::DisplayOptionsWidget(void)
{
	setMinimumHeight(550);
	setMaximumHeight(550);
	setMinimumWidth(240);
	setMaximumWidth(240);
	resize(240, 460);
	setWindowTitle (tr("Display Options"));

	selectBackgroundColorButton = new QPushButton(tr("Select"), this);
	selectBackgroundColorButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectBackgroundColorButton->setGeometry(130, 25, 100, 30);
	connect(selectBackgroundColorButton, SIGNAL(clicked(bool)), this, SLOT(getBackgroundColor()));

	selectBorderColorButton = new QPushButton(tr("Select"), this);
	selectBorderColorButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectBorderColorButton->setGeometry(130, 85, 100, 30);
	selectBorderColorButton->setEnabled(draw_borders);
	connect(selectBorderColorButton, SIGNAL(clicked(bool)), this, SLOT(getBorderColor()));

	selectWallColorButton = new QPushButton(tr("Select"), this);
	selectWallColorButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectWallColorButton->setGeometry(130, 145, 100, 30);
	connect(selectWallColorButton, SIGNAL(clicked(bool)), this, SLOT(getWallColor()));

	selectTextColorButton = new QPushButton(tr("Select"), this);
	selectTextColorButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectTextColorButton->setGeometry(130, 205, 100, 30);
	selectTextColorButton->setEnabled(display_text_info);
	connect(selectTextColorButton, SIGNAL(clicked(bool)), this, SLOT(getTextColor()));

	displayInfoCheckBox = new QCheckBox("Display Info", this);
	displayInfoCheckBox->setGeometry(10, 260, 200, 30);
	displayInfoCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	displayInfoCheckBox->setCheckState(display_text_info);
	connect(displayInfoCheckBox, SIGNAL(stateChanged(int)), this, SLOT(displayInfoChanged(int)));
	
	displayPressureCheckBox = new QCheckBox("Display Pressure", this);
	displayPressureCheckBox->setGeometry(10, 290, 200, 30);
	displayPressureCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	displayPressureCheckBox->setCheckState(display_pressure);
	displayPressureCheckBox->setToolTip(tr("If checked pressure on top wall is displayed"));
	connect(displayPressureCheckBox, SIGNAL(stateChanged(int)), this, SLOT(displayPressureChanged(int)));

	displayForceCheckBox = new QCheckBox("Display Force", this);
	displayForceCheckBox->setGeometry(10, 320, 200, 30);
	displayForceCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	displayForceCheckBox->setCheckState(display_force);
	displayForceCheckBox->setToolTip(tr("If checked force on top wall is displayed"));
	connect(displayForceCheckBox, SIGNAL(stateChanged(int)), this, SLOT(displayForceChanged(int)));

	displayChangedContactsCheckBox = new QCheckBox("Display Contacts", this);
	displayChangedContactsCheckBox->setGeometry(10, 350, 200, 30);
	displayChangedContactsCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	displayChangedContactsCheckBox->setCheckState(display_changed_contacts);
	displayChangedContactsCheckBox->setToolTip(tr("If checked the number of broken/established contacts is displayed"));
	connect(displayChangedContactsCheckBox, SIGNAL(stateChanged(int)), this, SLOT(displayChangedContactsChanged(int)));

	drawKeyCheckBox = new QCheckBox("Display Key", this);
	drawKeyCheckBox->setGeometry(10, 380, 200, 30);
	drawKeyCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	drawKeyCheckBox->setCheckState(display_key);
	connect(drawKeyCheckBox, SIGNAL(stateChanged(int)), this, SLOT(drawKeyChanged(int)));

	drawParticlesCheckBox = new QCheckBox("Draw Particles", this);
	drawParticlesCheckBox->setGeometry(10, 410, 200, 30);
	drawParticlesCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	drawParticlesCheckBox->setCheckState(draw_particles);
	connect(drawParticlesCheckBox, SIGNAL(stateChanged(int)), this, SLOT(drawParticlesChanged(int)));

	drawBordersCheckBox = new QCheckBox("Draw Borders", this);
	drawBordersCheckBox->setGeometry(10, 440, 200, 30);
	drawBordersCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	drawBordersCheckBox->setCheckState(draw_borders);
	connect(drawBordersCheckBox, SIGNAL(stateChanged(int)), this, SLOT(drawBordersChanged(int)));

	borderMinDepthEdit = new QDoubleSpinBox(this);
	borderMinDepthEdit->setGeometry(170, 470, 60, 30); 
	borderMinDepthEdit->setRange(0, 10000.0f);
	borderMinDepthEdit->setValue(border_min_depth);
	borderMinDepthEdit->setSingleStep(1.0f);
	connect(borderMinDepthEdit, SIGNAL(valueChanged(double)), this, SLOT(borderMinDepthChanged(double)));

	borderDistanceFactorEdit = new QDoubleSpinBox(this);
	borderDistanceFactorEdit->setGeometry(170, 510, 60, 30); 
	borderDistanceFactorEdit->setRange(0, 10000.0f);
	borderDistanceFactorEdit->setValue(border_distance_factor);
	borderDistanceFactorEdit->setSingleStep(1.0f);
	connect(borderDistanceFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(borderDistanceFactorChanged(double)));

	q_background_color.setRedF(clear_color[0]);
	q_background_color.setGreenF(clear_color[1]);
	q_background_color.setBlueF(clear_color[2]);
	
	q_border_color.setRedF(border_color[0]);
	q_border_color.setGreenF(border_color[1]);
	q_border_color.setBlueF(border_color[2]);

	q_wall_color.setRedF(wall_color[0]);
	q_wall_color.setGreenF(wall_color[1]);
	q_wall_color.setBlueF(wall_color[2]);

	q_text_color.setRedF(text_color[0]);
	q_text_color.setGreenF(text_color[1]);
	q_text_color.setBlueF(text_color[2]);

	drawBordersChanged(draw_borders);
}

DisplayOptionsWidget::~DisplayOptionsWidget(void)
{
}

void DisplayOptionsWidget::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 20, tr("Background color:"));
	painter.drawText(10, 80, tr("Border color:"));
	painter.drawText(10, 140, tr("Wall color:"));
	painter.drawText(10, 200, tr("Text color:"));
	painter.fillRect (10, 25, 110, 30, q_background_color);
	painter.fillRect (10, 85, 110, 30, q_border_color);
	painter.fillRect (10, 145, 110, 30, q_wall_color);
	painter.fillRect (10, 205, 110, 30, q_text_color);

	painter.drawText(10, 490, tr("Minimum border depth:"));
	painter.drawText(10, 530, tr("Border thickness factor:"));
}

void DisplayOptionsWidget::getBackgroundColor()
{
	q_background_color = QColorDialog::getColor(q_background_color, this, "Select Background Color");

	clear_color[0] = q_background_color.redF();
	clear_color[1] = q_background_color.greenF();
	clear_color[2] = q_background_color.blueF();

	emit signalUpdateView();
}

void DisplayOptionsWidget::getBorderColor()
{
	q_border_color = QColorDialog::getColor(q_border_color, this, "Select Border Color");

	border_color[0] = q_border_color.redF();
	border_color[1] = q_border_color.greenF();
	border_color[2] = q_border_color.blueF();

	emit signalUpdateView();
}

void DisplayOptionsWidget::getWallColor()
{
	q_wall_color = QColorDialog::getColor(q_wall_color, this, "Select Wall Color");

	wall_color[0] = q_wall_color.redF();
	wall_color[1] = q_wall_color.greenF();
	wall_color[2] = q_wall_color.blueF();

	emit signalUpdateView();
}

void DisplayOptionsWidget::getTextColor()
{
	q_text_color = QColorDialog::getColor(q_text_color, this, "Select Text Color");

	text_color[0] = q_text_color.redF();
	text_color[1] = q_text_color.greenF();
	text_color[2] = q_text_color.blueF();

	emit signalUpdateView();
}

void DisplayOptionsWidget::displayInfoChanged(int state)
{
	display_text_info = (Qt::CheckState)state;

	if(state == Qt::Checked)
		selectTextColorButton->setEnabled(true);
	else
		selectTextColorButton->setEnabled(false);

	emit signalUpdateView();
}

void DisplayOptionsWidget::displayPressureChanged(int state)
{
	display_pressure = (Qt::CheckState)state;

	if(state == Qt::Unchecked)
		sim_info.pressure = -1.0;

	emit signalUpdateView();
}

void DisplayOptionsWidget::displayForceChanged(int state)
{
	display_force = (Qt::CheckState)state;

	if(state == Qt::Unchecked)
		sim_info.force = -1.0;

	emit signalUpdateView();
}

void DisplayOptionsWidget::displayChangedContactsChanged(int state)
{
	display_changed_contacts = (Qt::CheckState)state;
	emit signalUpdateView();
}

void DisplayOptionsWidget::borderMinDepthChanged(double value)
{
	border_min_depth = value;
	emit signalUpdateView();
}

void DisplayOptionsWidget::borderDistanceFactorChanged(double value)
{
	border_distance_factor = value;
	emit signalUpdateView();
}

void DisplayOptionsWidget::drawBordersChanged(int state)
{
	draw_borders = (Qt::CheckState)state;

	if(state == Qt::Checked)
	{
		selectBorderColorButton->setEnabled(true);
		borderMinDepthEdit->setEnabled(true);
		borderDistanceFactorEdit->setEnabled(true);
	}
	else
	{
		selectBorderColorButton->setEnabled(false);
		borderMinDepthEdit->setEnabled(false);
		borderDistanceFactorEdit->setEnabled(false);
	}

	emit signalUpdateShaders();
	emit signalUpdateView();
}

void DisplayOptionsWidget::drawParticlesChanged(int state)
{
	draw_particles = (Qt::CheckState)state;
	emit signalUpdateShaders();
	emit signalUpdateView();
}

void DisplayOptionsWidget::drawKeyChanged(int state)
{
	display_key = (Qt::CheckState)state;
	emit signalUpdateView();
}
