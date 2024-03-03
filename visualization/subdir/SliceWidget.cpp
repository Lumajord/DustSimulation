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

#include "SliceWidget.h"

extern Qt::CheckState slice_random_orientation;
extern double slice_box_x_size;
extern double slice_box_y_size;
extern double slice_box_z_size;
extern double slice_sphere_radius;
extern double slice_top_slice_factor;
extern QString file_path;

SliceWidget::SliceWidget()
{
	resize(300, 310);
	setMinimumWidth(300);
	setMaximumWidth(300);
	setMinimumHeight(560);
	setMaximumHeight(560);

	setWindowTitle(tr("Slice Aggregate"));

	selectAgglomerateButton = new QPushButton(tr("Select Agglomerate"), this);
	selectAgglomerateButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	selectAgglomerateButton->setGeometry(10, 20, 280, 30);
	connect(selectAgglomerateButton, SIGNAL(clicked(bool)), this, SLOT(selectAgglomerate()));

	sliceCurrentAgglomerateCheckBox = new QCheckBox("Slice current Agglomerate", this);
	sliceCurrentAgglomerateCheckBox->setGeometry(10, 55, 280, 30);
	sliceCurrentAgglomerateCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	sliceCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	connect(sliceCurrentAgglomerateCheckBox, SIGNAL(stateChanged(int)), this, SLOT(sliceCurrentAgglomerateChanged(int)));

	randomOrientationCheckBox = new QCheckBox("Random Orientation", this);
	randomOrientationCheckBox->setGeometry(10, 85, 280, 30);
	randomOrientationCheckBox->setFont(QFont("Helvetica", 10, QFont::Bold));
	randomOrientationCheckBox->setCheckState(slice_random_orientation);
	randomOrientationCheckBox->setToolTip(tr("Rotate aggregate randomly before slicing (only affects sphere)"));
	connect(randomOrientationCheckBox, SIGNAL(stateChanged(int)), this, SLOT(randomOrientationChanged(int)));

	boxXSizeEdit = new QDoubleSpinBox(this);
	boxXSizeEdit->setGeometry(200, 130, 90, 30); 
	boxXSizeEdit->setRange(0, 1000);
	boxXSizeEdit->setValue(slice_box_x_size);
	boxXSizeEdit->setSingleStep(0.1);
	boxXSizeEdit->setToolTip(tr("x-size (in µm) of the slice box"));
	connect(boxXSizeEdit, SIGNAL(valueChanged(double)), this, SLOT(boxXSizeChanged(double)));

	boxYSizeEdit = new QDoubleSpinBox(this);
	boxYSizeEdit->setGeometry(200, 170, 90, 30); 
	boxYSizeEdit->setRange(0, 1000);
	boxYSizeEdit->setValue(slice_box_y_size);
	boxYSizeEdit->setSingleStep(0.1);
	boxYSizeEdit->setToolTip(tr("x-size (in µm) of the slice box"));
	connect(boxYSizeEdit, SIGNAL(valueChanged(double)), this, SLOT(boxYSizeChanged(double)));

	boxZSizeEdit = new QDoubleSpinBox(this);
	boxZSizeEdit->setGeometry(200, 210, 90, 30); 
	boxZSizeEdit->setRange(0, 1000);
	boxZSizeEdit->setValue(slice_box_z_size);
	boxZSizeEdit->setSingleStep(0.1);
	boxZSizeEdit->setToolTip(tr("x-size (in µm) of the slice box"));
	connect(boxZSizeEdit, SIGNAL(valueChanged(double)), this, SLOT(boxZSizeChanged(double)));

	sliceBoxButton = new QPushButton(tr("Slice Box"), this);
	sliceBoxButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	sliceBoxButton->setGeometry(10, 250, 280, 30);
	connect(sliceBoxButton, SIGNAL(clicked(bool)), this, SLOT(sliceBox()));

	sphereRadiusModifierEdit = new QDoubleSpinBox(this);
	sphereRadiusModifierEdit->setGeometry(200, 305, 90, 30); 
	sphereRadiusModifierEdit->setRange(0, 1000);
	sphereRadiusModifierEdit->setValue(slice_sphere_radius);
	sphereRadiusModifierEdit->setSingleStep(1);
	sphereRadiusModifierEdit->setToolTip(tr("Modifies the slice radius of the sphere / cylinder"));
	connect(sphereRadiusModifierEdit, SIGNAL(valueChanged(double)), this, SLOT(sphereRadiusModifierChanged(double)));

	sliceSphereButton = new QPushButton(tr("Slice Sphere"), this);
	sliceSphereButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	sliceSphereButton->setGeometry(10, 345, 280, 30);
	connect(sliceSphereButton, SIGNAL(clicked(bool)), this, SLOT(sliceSphere()));

	sliceCylinderButton = new QPushButton(tr("Slice Cylinder"), this);
	sliceCylinderButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	sliceCylinderButton->setGeometry(10, 385, 280, 30);
	sliceCylinderButton->setToolTip(tr("Slices a cylinder around the Y-Axis"));
	connect(sliceCylinderButton, SIGNAL(clicked(bool)), this, SLOT(sliceCylinder()));

	topSliceFactorEdit = new QDoubleSpinBox(this);
	topSliceFactorEdit->setGeometry(200, 440, 90, 30); 
	topSliceFactorEdit->setRange(0, 1);
	topSliceFactorEdit->setValue(slice_top_slice_factor);
	topSliceFactorEdit->setSingleStep(0.05);
	topSliceFactorEdit->setToolTip(tr("Determines what percentage of the aggregate is sliced from the top/bottom"));
	connect(topSliceFactorEdit, SIGNAL(valueChanged(double)), this, SLOT(topSliceFactorChanged(double)));

	sliceTopButton = new QPushButton(tr("Slice Top"), this);
	sliceTopButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	sliceTopButton->setGeometry(10, 480, 280, 30);
	sliceTopButton->setToolTip(tr("Slices the top of the aggregte"));
	connect(sliceTopButton, SIGNAL(clicked(bool)), this, SLOT(sliceTop()));

	sliceBottomButton = new QPushButton(tr("Slice Bottom"), this);
	sliceBottomButton->setFont(QFont("Helvetica", 10, QFont::Bold));
	sliceBottomButton->setGeometry(10, 520, 280, 30);
	sliceBottomButton->setToolTip(tr("Slices the bottom of the aggregate"));
	connect(sliceBottomButton, SIGNAL(clicked(bool)), this, SLOT(sliceBottom()));

	strcpy(filename, "");
}

SliceWidget::~SliceWidget(void)
{
}

void SliceWidget::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	painter.setFont(QFont("Helvetica", 10, QFont::Bold));

	painter.drawText(10, 150, tr("Box x-size (in µm):"));
	painter.drawText(10, 190, tr("Box y-size (in µm):"));
	painter.drawText(10, 230, tr("Box z-size (in µm):"));
	painter.drawText(10, 325, tr("Slice radius (in µm):"));
	painter.drawText(10, 460, tr("Top/Bottom slice factor:"));
}

void SliceWidget::selectAgglomerate()
{
	QString temp = QFileDialog::getOpenFileName(this, tr("Select Agglomerate"), file_path, tr("Particle Files (*.dat)"));

	if(!temp.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(temp).path();
        strcpy(filename, temp.toLatin1().data());

		selectAgglomerateButton->setText(QFileInfo(temp).fileName());
		sliceCurrentAgglomerateCheckBox->setCheckState(Qt::Unchecked);
	}
}

void SliceWidget::sliceCurrentAgglomerateChanged(int value)
{
	if(value == Qt::Unchecked)
	{
		if( strcmp(filename, "") == 0)
			sliceCurrentAgglomerateCheckBox->setCheckState(Qt::Checked);
	}
}

void SliceWidget::randomOrientationChanged(int value)
{
	slice_random_orientation = (Qt::CheckState)value;
}

void SliceWidget::boxXSizeChanged(double value)
{
	slice_box_x_size = value;
}

void SliceWidget::boxYSizeChanged(double value)
{
	slice_box_y_size = value;
}

void SliceWidget::boxZSizeChanged(double value)
{
	slice_box_z_size = value;
}

void SliceWidget::sphereRadiusModifierChanged(double value)
{
	slice_sphere_radius = value;
}

void SliceWidget::topSliceFactorChanged(double value)
{
	slice_top_slice_factor = value;
}

void SliceWidget::sliceBox()
{
	if(sliceCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
		emit signalSliceBox(NULL, slice_box_x_size*1e-4, slice_box_z_size*1e-4, slice_box_y_size*1e-4);
	else
		emit signalSliceBox(filename, slice_box_x_size, slice_box_y_size, slice_box_z_size);
}

void SliceWidget::sliceSphere()
{
	if(sliceCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
		emit signalSliceSphere(NULL, slice_sphere_radius * 1e-4, slice_random_orientation);
	else
		emit signalSliceSphere(filename, slice_sphere_radius * 1e-4, slice_random_orientation);
}

void SliceWidget::sliceCylinder()
{
	if(sliceCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
		emit signalSliceCylinder(NULL, slice_sphere_radius * 1e-4);
	else
		emit signalSliceCylinder(filename, slice_sphere_radius * 1e-4);
}

void SliceWidget::sliceTop()
{
	if(sliceCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
		emit signalSliceTop(NULL, slice_top_slice_factor);
	else
		emit signalSliceTop(filename, slice_top_slice_factor);
}

void SliceWidget::sliceBottom()
{
	if(sliceCurrentAgglomerateCheckBox->checkState() == Qt::Checked)
		emit signalSliceBottom(NULL, slice_top_slice_factor);
	else
		emit signalSliceBottom(filename, slice_top_slice_factor);
}

