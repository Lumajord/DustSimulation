#include "OpenGLNavigationWidget.h"

OpenGLNavigationWidget::OpenGLNavigationWidget(QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f)
: QGLWidget(parent, shareWidget, f)
{
	initialize();
}

OpenGLNavigationWidget::OpenGLNavigationWidget(QGLContext* context, QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f)
: QGLWidget(context, parent, shareWidget, f)
{
	initialize();
}

OpenGLNavigationWidget::OpenGLNavigationWidget(const QGLFormat& format, QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f)
: QGLWidget(format, parent, shareWidget, f)
{
	initialize();
}

OpenGLNavigationWidget::~OpenGLNavigationWidget()
{
}

/**
	initialize everything (called from constructors)
*/
void OpenGLNavigationWidget::initialize()
{
    makeCurrent();
    initializeOpenGLFunctions();

	cameraDefaultPosition(0) = 0.0;
	cameraDefaultPosition(1) = 0.0;
	cameraDefaultPosition(2) = -1.0;

	cameraDefaultLookAt(0) = 0.0;
	cameraDefaultLookAt(1) = 0.0;
	cameraDefaultLookAt(2) = 0.0;

	cameraDefaultUp(0) = 0.0;
	cameraDefaultUp(1) = 1.0;
	cameraDefaultUp(2) = 0.0;
	
	normalMoveFactor = 1.0;
	fastMoveFactor = 5.0;
	normalRotateFactor = 1.0;
	fastRotateFactor = 5.0;
	normalZoomFactor = 1.0;
	fastZoomFactor = 5.0;

	resetCamera();
}

/**
	set camera default position
*/
void OpenGLNavigationWidget::setCameraDefaultPosition(const Vector<GLdouble, 3> &v)
{
	cameraDefaultPosition = v;
}

/**
	set camera default target
*/
void OpenGLNavigationWidget::setCameraDefaultLookAt(const Vector<GLdouble, 3> &v)
{
	cameraDefaultLookAt = v;
}

/**
	set camera default up vector
*/
void OpenGLNavigationWidget::setCameraDefaultUp(const Vector<GLdouble, 3> &v)
{
	cameraDefaultUp = v;
}

/**
	set camera position
*/
void OpenGLNavigationWidget::setCameraPosition(const Vector<GLdouble, 3> &v)
{
	cameraPosition = v;
	updateCameraBasis();
}

/**
	set camera target
*/
void OpenGLNavigationWidget::setCameraLookAt(const Vector<GLdouble, 3> &v)
{
	cameraLookAt = v;
	updateCameraBasis();
}

/**
	set camera up vector
*/
void OpenGLNavigationWidget::setCameraUp(const Vector<GLdouble, 3> &v)
{
	cameraUp = v;
	cameraUp.normalize();
	updateCameraBasis();
}


/**
	updates the camera basis (camera forward, camera side, camera up)
*/
void OpenGLNavigationWidget::updateCameraBasis()
{
	cameraForward = cameraLookAt - cameraPosition;
	cameraForward.normalize();

	cameraSide = cameraForward^cameraUp;
	cameraSide.normalize();

	cameraUp = cameraSide^cameraForward;
	cameraUp.normalize();

	emit cameraUpdated();
}

/**
	reset camera
*/
void OpenGLNavigationWidget::resetCamera()
{
	cameraPosition = cameraDefaultPosition;
	cameraLookAt = cameraDefaultLookAt;
	cameraUp = cameraDefaultUp;

	updateCameraBasis();
	update();
}

/**
	setup camera in OpenGL
*/
void OpenGLNavigationWidget::setupCamera()
{
	glLoadIdentity();

	Matrix<GLdouble, 4, 4> matrix;

	matrix(0,0) = cameraSide(0);
	matrix(1,0) = cameraSide(1);
	matrix(2,0) = cameraSide(2);
	matrix(3,0) = 0.0;

	matrix(0,1) = cameraUp(0);
	matrix(1,1) = cameraUp(1);
	matrix(2,1) = cameraUp(2);
	matrix(3,1) = 0.0;

	matrix(0,2) = -(cameraForward(0));
	matrix(1,2) = -(cameraForward(1));
	matrix(2,2) = -(cameraForward(2));
	matrix(3,2) = 0.0;

	matrix(0,3) = 0.0;
	matrix(1,3) = 0.0;
	matrix(2,3) = 0.0;
	matrix(3,3) = 1.0;

	glMultMatrixd(matrix);
	glTranslated(-cameraPosition(0), -cameraPosition(1), -cameraPosition(2));
}

void OpenGLNavigationWidget::rotateCamera()
{
	Matrix<GLdouble, 4, 4> matrix;

	matrix(0,0) = cameraSide(0);
	matrix(1,0) = cameraSide(1);
	matrix(2,0) = cameraSide(2);
	matrix(3,0) = 0.0;

	matrix(0,1) = cameraUp(0);
	matrix(1,1) = cameraUp(1);
	matrix(2,1) = cameraUp(2);
	matrix(3,1) = 0.0;

	matrix(0,2) = -(cameraForward(0));
	matrix(1,2) = -(cameraForward(1));
	matrix(2,2) = -(cameraForward(2));
	matrix(3,2) = 0.0;

	matrix(0,3) = 0.0;
	matrix(1,3) = 0.0;
	matrix(2,3) = 0.0;
	matrix(3,3) = 1.0;

	glMultMatrixd(matrix);
}

/**
	handle mouse press
*/
void OpenGLNavigationWidget::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton) {
		lastDragPosition = event->pos();
	}

	if (event->button() == Qt::RightButton) {
		lastDragPosition = event->pos();
	}
}

/**
	handle mouse move
*/
void OpenGLNavigationWidget::mouseMoveEvent(QMouseEvent *event)
{
	// rotation with right mouse button
	if (event->buttons() & Qt::RightButton) {
		
		Matrix<GLdouble, 3, 3> r;
		r = rotationMatrix(cameraUp, -(double)(event->pos().x() - lastDragPosition.x()) / 180.0 * (event->modifiers() & Qt::ControlModifier ? fastRotateFactor : normalRotateFactor));
		cameraPosition = r*cameraPosition;

		r = rotationMatrix(cameraSide, -(double)(event->pos().y() - lastDragPosition.y()) / 180.0 * (event->modifiers() & Qt::ControlModifier ? fastRotateFactor : normalRotateFactor));
		cameraPosition = r*cameraPosition;

		updateCameraBasis();
		lastDragPosition = event->pos();
		update();
	}

	// translation with left mouse button
	if (event->buttons() & Qt::LeftButton) {
		cameraLookAt -= cameraSide * (GLdouble)(event->pos().x() - lastDragPosition.x()) * (event->modifiers() & Qt::ControlModifier ? fastMoveFactor : normalMoveFactor);
		cameraPosition -= cameraSide * (GLdouble)(event->pos().x() - lastDragPosition.x()) * (event->modifiers() & Qt::ControlModifier ? fastMoveFactor : normalMoveFactor);

		cameraLookAt += cameraUp * (GLdouble)(event->pos().y() - lastDragPosition.y()) * (event->modifiers() & Qt::ControlModifier ? fastMoveFactor : normalMoveFactor);
		cameraPosition += cameraUp * (GLdouble)(event->pos().y() - lastDragPosition.y()) * (event->modifiers() & Qt::ControlModifier ? fastMoveFactor : normalMoveFactor);

		updateCameraBasis();
		lastDragPosition = event->pos();
		update();
	}
}

/**
	handle mouse wheel
*/
void OpenGLNavigationWidget::wheelEvent(QWheelEvent *event)
{
	double scrollFactor = 1.0 - (GLdouble)(event->delta()/120)*0.01*(event->modifiers() & Qt::ControlModifier ? fastZoomFactor : normalZoomFactor);

	cameraPosition = cameraLookAt - scrollFactor * (cameraLookAt - cameraPosition);

	updateCameraBasis();
	update();
}
