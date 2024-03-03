/***************************************************************************
 *   OpenGL Navigation Widget                                              *
 *   Copyright (C) 2010 by Tobias Müller                                   *
 *   Tobias_Mueller@twam.info                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _OPENGLNAVIGATIONWIDGET_H_
#define _OPENGLNAVIGATIONWIDGET_H_


#include <QGLWidget>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QOpenGLFunctions>
#include "Vector.h"
#include "Matrix.h"

class OpenGLNavigationWidget : public QGLWidget, protected QOpenGLFunctions
{
	Q_OBJECT

	public:
		OpenGLNavigationWidget(QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags f = 0);
		OpenGLNavigationWidget(QGLContext* context, QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags f = 0);
		OpenGLNavigationWidget(const QGLFormat& format, QWidget* parent = 0, const QGLWidget* shareWidget = 0, Qt::WindowFlags f = 0);
		~OpenGLNavigationWidget();

		// setter
		void setCameraDefaultPosition(const Vector<GLdouble, 3> &);
		void setCameraDefaultLookAt(const Vector<GLdouble, 3> &);
		void setCameraDefaultUp(const Vector<GLdouble, 3> &);
		void setCameraPosition(const Vector<GLdouble, 3> &);
		void setCameraLookAt(const Vector<GLdouble, 3> &);
		void setCameraUp(const Vector<GLdouble, 3> &);
		inline void setNormalMoveFactor(double c) { normalMoveFactor = c; }
		inline void setFastMoveFactor(double c) { fastMoveFactor = c; }
		inline void setNormalRotationFactor(double c) { normalRotateFactor = c; }
		inline void setFastRotationFactor(double c) { fastRotateFactor = c; }
		inline void setNormalZoomFactor(double c) { normalZoomFactor = c; }
		inline void setFastZoomFactor(double c) { fastZoomFactor = c; }

		// getter
		inline const Vector<GLdouble ,3>& getCameraDefaultPosition() const { return cameraDefaultPosition; }
		inline const Vector<GLdouble ,3>& getCameraDefaultLookAt() const { return cameraDefaultLookAt; }
		inline const Vector<GLdouble ,3>& getCameraDefaultUp() const { return cameraDefaultUp; }
		inline const Vector<GLdouble ,3>& getCameraPosition() const { return cameraPosition; }
		inline const Vector<GLdouble ,3>& getCameraLookAt() const { return cameraLookAt; }
		inline const Vector<GLdouble ,3>& getCameraUp() const { return cameraUp; }
		inline double getNormalMoveFactor() const { return normalMoveFactor; }
		inline double getFastMoveFactor() const { return fastMoveFactor; }
		inline double getNormalRotationFactor() const { return normalRotateFactor; }
		inline double getFastRotationFactor() const { return fastRotateFactor; }
		inline double getNormalZoomFactor() const { return normalZoomFactor; }
		inline double getFastZoomFactor() const { return fastZoomFactor; }

		void updateCameraBasis();
		void resetCamera();
		void setupCamera();
		void rotateCamera();

	protected:
		/// position of last mouse drag
		QPoint lastDragPosition;

		/// camera default position
		Vector<GLdouble, 3> cameraDefaultPosition;
		/// camera default target
		Vector<GLdouble, 3> cameraDefaultLookAt;
		/// camera default up vector
		Vector<GLdouble, 3> cameraDefaultUp;

		/// camera position
		Vector<GLdouble, 3> cameraPosition;
		/// camera target
		Vector<GLdouble, 3> cameraLookAt;

		/// camera up vector (direction which will be upwards on the screen), must be normalized
		Vector<GLdouble, 3> cameraUp;
		/// camera side vector (direction which will be sidewards on the screen), must be normalized
		Vector<GLdouble, 3> cameraSide;
		/// camera forward vector (direction which will be in camera direction), must be normalized
		Vector<GLdouble, 3> cameraForward;

		/// factor for normal moving
		double normalMoveFactor;
		/// factor for fast moving (ctrl pressed)
		double fastMoveFactor;
		/// factor for normal rotating
		double normalRotateFactor;
		/// factor for fast rotating
		double fastRotateFactor;
		/// factor for normal zooming
		double normalZoomFactor;
		/// factor for fast zooming (ctrl pressed)
		double fastZoomFactor;

		virtual void mousePressEvent(QMouseEvent *event);
		virtual void mouseMoveEvent(QMouseEvent *event);
		virtual void wheelEvent(QWheelEvent *event);

	private:
		void initialize();

	signals:
		void cameraUpdated();
};

#endif
