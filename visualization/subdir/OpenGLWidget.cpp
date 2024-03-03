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

#include "OpenGLWidget.h"
#include "Constants.h"
#include "SimulationThread.h"
#include <QtGui>

extern float clear_color[3];
extern float border_color[3];
extern float wall_color[3];
extern float text_color[3];
extern float border_distance_factor;
extern float border_min_depth;
extern double particle_radius;
extern Qt::CheckState display_key;
extern Qt::CheckState display_text_info;
extern Qt::CheckState display_changed_contacts;
extern Qt::CheckState draw_borders;
extern Qt::CheckState draw_particles;
extern int font_size;
extern bool follow_cms;

extern char material_name[200];
extern DisplayInfo sim_info;
extern QString info_text;
extern int sim_time_display_mode;

extern DisplayMode vis_display_mode;
extern double vis_density_max_filling_factor;
extern double vis_dislocation_min_value;
extern double vis_dislocation_max_value;
extern double delta_0;
extern double particle_radius;

#ifdef ENABLE_GRAVITY
    extern double gravity_strength;
    extern vec3 gravity_direction;
    extern bool gravity_enabled;
    extern double gravity_modifier;
#endif

extern bool damping_enabled;

#ifdef DRAW_CONTACT_POINTERS
OpenGLWidget::OpenGLWidget(QWidget *parent, multimap<int, vec3> **contact_pointers)
    : OpenGLNavigationWidget(parent)
{
    this->walls = walls;
    this->contact_pointers = contact_pointers;

    quadric = gluNewQuadric();
    gluQuadricNormals(quadric, GLU_SMOOTH);

    vbPosistion = NULL;
    vbColor = NULL;

    cam_pos_x = 0.0f;
    cam_pos_y = 0.0f;
    cam_pos_z = -150.0f * (float)particle_radius;

    cam_rot_x = 0;
    cam_rot_y = 0;
    cam_rot_z = 0;

    array_size = 0;
}
#else
OpenGLWidget::OpenGLWidget(QWidget *parent, SimulationThread *sim_thread)
    : OpenGLNavigationWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    this->sim_thread = sim_thread;

    quadric = gluNewQuadric();
    gluQuadricNormals(quadric, GLU_SMOOTH);

    vbPosition = NULL;
    vbColor = NULL;

    setNormalMoveFactor(0.05 * particle_radius);
    setFastMoveFactor(0.25 * particle_radius);
    setNormalRotationFactor(1.0);
    setFastRotationFactor(5.0);
    setNormalZoomFactor(3.0);
    setFastZoomFactor(15.0);

    setCameraDefaultPosition(Vector<GLdouble, 3>(3, 0.0, 0.0, 150.0 * particle_radius));
    setCameraDefaultLookAt(Vector<GLdouble, 3>(3, 0.0, 0.0, 0.0));
    setCameraDefaultUp(Vector<GLdouble, 3>(3, 0.0, 1.0, 0.0));

    resetCamera();

    particle_array_size = 0;

    font = QFont("Helvetica", font_size, QFont::Bold);
    font_script = QFont("Helvetica", font_size*3.0/4.0, QFont::Bold);
    font_metrics = new QFontMetrics(font);

    point_size_scale_factor = 0.5f;

    noise_texture = 0;
    noise_texture_ptr = NULL;
}
#endif


OpenGLWidget::~OpenGLWidget()
{
    makeCurrent();

    if(particle_array_size > 0)
    {
        glDeleteBuffers(1, (const GLuint*)&vbColor);
        glDeleteBuffers(1, (const GLuint*)&vbPosition);
    }

    delete font_metrics;

}

void OpenGLWidget::initializeGL()
{
    glEnable(GL_DEPTH_TEST);
    glClearColor(clear_color[0], clear_color[1], clear_color[2], clear_color[3]);

    glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glEnable(GL_LIGHTING);

    GLfloat global_ambient[] = {0.6f, 0.6f, 0.6f, 1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

    // lighting for wall
    GLfloat diffuse_wall[] = {0.5f, 0.5f, 0.5f , 1.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse_wall);

    GLfloat position_wall[] = {float(-10.0 * particle_radius), float(32.0 * particle_radius), float(-32.0 * particle_radius), 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, position_wall);

    // lighting for gravity vector
    GLfloat diffuse_arrow[] = {0.9f, 0.9f, 0.9f , 1.0f};
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse_arrow);

    GLfloat position_arrow[] = {80.0f, 60.0f, 80.0f, 1.0f};
    glLightfv(GL_LIGHT1, GL_POSITION, position_arrow);

    glShadeModel(GL_SMOOTH);

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glDepthFunc(GL_LEQUAL);

    glewInit();

    // make sure our graphic cards meets the requirements
    if(!glewIsSupported("GL_VERSION_2_0 GL_VERSION_1_5 GL_ARB_multitexture GL_ARB_vertex_buffer_object"))
    {
        QMessageBox::about(this, tr("Error"), tr("Required OpenGL extensions missing"));
        initialized = false;
    }
    else
        initialized = true;

#if defined (_WIN32)
    if (wglewIsSupported("WGL_EXT_swap_control")) {
        // disable vertical sync
        wglSwapIntervalEXT(0);
    }
#endif

    // default = sphereShader
    selectShader();

    // enable AA
    glEnable(GL_MULTISAMPLE_ARB);
}

void OpenGLWidget::selectShader()
{
    if(draw_borders == Qt::Checked)
    {
        if(draw_particles == Qt::Checked)
            program = compileProgram(vertexShader, sphereBoundariesPixelShader);
        else
            program = compileProgram(vertexShader, boundariesPixelShader);

        point_size_scale_factor = 0.25f;
    }
    else
    {
        program = compileProgram(vertexShader, spherePixelShader);
        point_size_scale_factor = 0.5f;
    }
}

void OpenGLWidget::resizeGL(int width, int height)
{
    // prevent division by 0
    if(height == 0)
        height = 1;

    GLfloat z_near = (GLfloat)particle_radius;
    GLfloat z_far = (GLfloat)(5000.0 * particle_radius);

    glViewport(0,0,width,height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(FIELD_OF_VIEW,(GLfloat)width/(GLfloat)height, z_near, z_far);

    glMatrixMode(GL_MODELVIEW);
}

void OpenGLWidget::paintGL()
{
//    makeCurrent();

    glEnable(GL_MULTISAMPLE);

    // clear view
    glClearColor(clear_color[0], clear_color[1], clear_color[2], clear_color[3]);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //////////////////////////////////////////////////////////////////////////////////
    // prepare drawing
    //////////////////////////////////////////////////////////////////////////////////

    glEnable(GL_POINT_SPRITE_ARB);
    glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);

    setupCamera();

    //////////////////////////////////////////////////////////////////////////////////
    // draw particles
    //////////////////////////////////////////////////////////////////////////////////

    glColor3f(1, 1, 1);

    // prepare shaders
    glUseProgram(program);
    glUniform1f( glGetUniformLocation(program, "pointScale"), this->height() / tanf(FIELD_OF_VIEW * point_size_scale_factor * (float)M_PI/180.0f) );
    glUniform1f( glGetUniformLocation(program, "pointRadius"), particle_radius );
    glUniform4f( glGetUniformLocation(program, "clear_color"), clear_color[0], clear_color[1], clear_color[2], 1.0f);
    glUniform4f( glGetUniformLocation(program, "border_color"), border_color[0], border_color[1], border_color[2], 1.0f);
    glUniform1f( glGetUniformLocation(program, "distance_factor"), particle_radius * border_distance_factor);
    glUniform1f( glGetUniformLocation(program, "min_border_depth"), particle_radius * border_min_depth);

    // activate noise texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, noise_texture);

    // bind buffers
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbPosition);
    glVertexPointer(3, GL_DOUBLE, 0, 0);
    glEnableClientState(GL_VERTEX_ARRAY);

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbColor);
    glColorPointer(3, GL_FLOAT, 0, 0);
    glEnableClientState(GL_COLOR_ARRAY);

    glDrawArrays(GL_POINTS, 0, particle_array_size);

    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);

    glUseProgram(0);
    glDisable(GL_POINT_SPRITE_ARB);

    //////////////////////////////////////////////////////////////////////////////////
    // draw walls
    //////////////////////////////////////////////////////////////////////////////////

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_LIGHT0);

    for(int w = 0; w < sim_thread->number_of_walls; ++w)
    {
        if(sim_thread->alpha_values[w] != 0)
        {
            glColor4f(wall_color[0], wall_color[1], wall_color[2], sim_thread->alpha_values[w]);

            // uppper
            glBegin(GL_QUADS);
            glNormal3d(sim_thread->wall_normals[3*w], sim_thread->wall_normals[3*w+1], sim_thread->wall_normals[3*w+2]);
            glVertex3d(sim_thread->wall_positions[24*w], sim_thread->wall_positions[24*w+1], sim_thread->wall_positions[24*w+2]);
            glVertex3d(sim_thread->wall_positions[24*w+3], sim_thread->wall_positions[24*w+4], sim_thread->wall_positions[24*w+5]);
            glVertex3d(sim_thread->wall_positions[24*w+6], sim_thread->wall_positions[24*w+7], sim_thread->wall_positions[24*w+8]);
            glVertex3d(sim_thread->wall_positions[24*w+9], sim_thread->wall_positions[24*w+10], sim_thread->wall_positions[24*w+11]);
            glEnd();

            // lower
            //glBegin(GL_QUADS);
            //glVertex3f(walls[w].pos1[0] - walls[w].normal[0], walls[w].pos1[1] - walls[w].normal[1], walls[w].pos1[2] - walls[w].normal[2]);
            //glVertex3f(walls[w].pos4[0] - walls[w].normal[0], walls[w].pos4[1] - walls[w].normal[1], walls[w].pos4[2] - walls[w].normal[2]);
            //glVertex3f(walls[w].pos3[0] - walls[w].normal[0], walls[w].pos3[1] - walls[w].normal[1], walls[w].pos3[2] - walls[w].normal[2]);
            //glVertex3f(walls[w].pos2[0] - walls[w].normal[0], walls[w].pos2[1] - walls[w].normal[1], walls[w].pos2[2] - walls[w].normal[2]);
            //glEnd();
        }
    }

    glDisable(GL_BLEND);

    //////////////////////////////////////////////////////////////////////////////////
    // draw gyration radii
    //////////////////////////////////////////////////////////////////////////////////

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for(int i = 0; i < sim_thread->gyration_radii.size(); ++i)
    {
        glColor4fv((GLfloat*) (&sim_thread->gyration_radii[i].color) );

        glPushMatrix();
        glTranslatef(sim_thread->gyration_radii[i].pos[0], sim_thread->gyration_radii[i].pos[1], sim_thread->gyration_radii[i].pos[2]);
        gluSphere(quadric, (float) sim_thread->gyration_radii[i].radius, 24, 24);
        glPopMatrix();
    }

    glDisable(GL_BLEND);
    glDisable(GL_LIGHT0);

    //////////////////////////////////////////////////////////////////////////////////
    // display gravity vector
    //////////////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_GRAVITY
    if(gravity_enabled)
    {
        GLfloat global_ambient_arrow[] = {0.2f, 0.2f, 0.2f, 1.0f};
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient_arrow);
        glEnable(GL_LIGHT1);

        double length = 30;
        double length2 = length * (double)width() / (double)height();
        const double x_center = 50;
        const double y_center = 60;
        const double z_center = 0;

        const double arrow_shaft_width = 0.15;
        const double arrow_head_width = 0.4;
        const double arrow_head_length = 0.3;

        double x_start = - length * gravity_direction[0];
        double y_start = - length * gravity_direction[1];
        double z_start = - length * gravity_direction[2];
        double x_middle = (0.5 - arrow_head_length) * length * gravity_direction[0];
        double y_middle = (0.5 - arrow_head_length) * length * gravity_direction[1];
        double z_middle = (0.5 - arrow_head_length) * length * gravity_direction[2];
        double x_end = length * gravity_direction[0];
        double y_end = length * gravity_direction[1];
        double z_end = length * gravity_direction[2];

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, width(), 0, height(), -150, 150);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glTranslated(x_center, y_center, z_center);
        rotateCamera();

        glColor3f(0.7f, 0.7f, 0.7f);

        // arrow shaft cap
        glBegin(GL_TRIANGLE_FAN);
        glNormal3d(-gravity_direction[0], -gravity_direction[1], -gravity_direction[2]);
        glVertex3d(x_start, y_start, z_start);
        for(int i = 0; i < LOD; ++i)
        {
            glNormal3d(-gravity_direction[0], -gravity_direction[1], -gravity_direction[2]);
            double x = x_start - cos_table[i] * arrow_shaft_width * length * gravity_direction[1] - sin_table[i] * arrow_shaft_width * length2 * gravity_direction[2];
            double y = y_start + cos_table[i] * arrow_shaft_width * length * gravity_direction[0];
            double z = z_start + sin_table[i] * arrow_shaft_width * length2 * gravity_direction[0];
            glVertex3d(x,y,z);
        }
        glEnd();

        // arrow shaft
        glBegin(GL_QUAD_STRIP);
        for(int i = 0; i < LOD; ++i)
        {
            glNormal3d(cos_table[i] * gravity_direction[1] - sin_table[i] * gravity_direction[2], cos_table[i] * gravity_direction[0], sin_table[i] * gravity_direction[0]);
            double x = x_start - cos_table[i] * arrow_shaft_width * length * gravity_direction[1] - sin_table[i] * arrow_shaft_width * length2 * gravity_direction[2];
            double y = y_start + cos_table[i] * arrow_shaft_width * length * gravity_direction[0];
            double z = z_start + sin_table[i] * arrow_shaft_width * length2 * gravity_direction[0];
            glVertex3d(x,y,z);

            glNormal3d(cos_table[i] * gravity_direction[1] - sin_table[i] * gravity_direction[2], cos_table[i] * gravity_direction[0], sin_table[i] * gravity_direction[0]);
            x = x_middle - cos_table[i] * arrow_shaft_width * length * gravity_direction[1] - sin_table[i] * arrow_shaft_width * length2 * gravity_direction[2];
            y = y_middle + cos_table[i] * arrow_shaft_width * length * gravity_direction[0];
            z = z_middle + sin_table[i] * arrow_shaft_width * length2 * gravity_direction[0];
            glVertex3d(x,y,z);
        }
        glEnd();

        // arrow head cap
        glBegin(GL_TRIANGLE_FAN);
        glNormal3d(-gravity_direction[0], -gravity_direction[1], -gravity_direction[2]);
        glVertex3d(x_middle, y_middle, z_middle);
        for(int i = 0; i < LOD; ++i)
        {
            glNormal3d(-gravity_direction[0], -gravity_direction[1], -gravity_direction[2]);
            double x = x_middle - cos_table[i] * arrow_head_width * length * gravity_direction[1] - sin_table[i] * arrow_head_width * length2 * gravity_direction[2];
            double y = y_middle + cos_table[i] * arrow_head_width * length * gravity_direction[0];
            double z = z_middle + sin_table[i] * arrow_head_width * length2 * gravity_direction[0];
            glVertex3d(x,y,z);
        }
        glEnd();

        // arrow head
        glBegin(GL_TRIANGLE_FAN);
        glNormal3d(gravity_direction[0], gravity_direction[1], gravity_direction[2]);
        glVertex3d(x_end, y_end, z_end);
        for(int i = 0; i < LOD; ++i)
        {
            glNormal3d(cos_table[i] * gravity_direction[1] - sin_table[i] * gravity_direction[2], cos_table[i] * gravity_direction[0], sin_table[i] * gravity_direction[0]);
            double x = x_middle - cos_table[i] * arrow_head_width * length * gravity_direction[1] - sin_table[i] * arrow_head_width * length2 * gravity_direction[2];
            double y = y_middle + cos_table[i] * arrow_head_width * length * gravity_direction[0];
            double z = z_middle + sin_table[i] * arrow_head_width * length2 * gravity_direction[0];
            glVertex3d(x,y,z);
        }
        glEnd();

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();

        glColor3f(text_color[0], text_color[1], text_color[2]);
        renderText(80, height() - 70, tr("Gravity enabled"), QFont("Helvetica", font_size-2, QFont::Bold));
        if(gravity_modifier > 1.001)
        {
            char buffer[50];
            sprintf(buffer, "Modifier: %.0f", gravity_modifier);
            renderText(80, height() - 50, buffer, QFont("Helvetica", font_size-2, QFont::Bold));
        }
        else if(gravity_modifier < 0.999)
        {
            char buffer[50];
            sprintf(buffer, "Modifier: %.2f", gravity_modifier);
            renderText(80, height() - 50, buffer, QFont("Helvetica", font_size-2, QFont::Bold));
        }

        glDisable(GL_LIGHT1);
        GLfloat global_ambient_wall[] = {0.6f, 0.6f, 0.6f, 1.0f};
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient_wall);
    }
#endif

    //////////////////////////////////////////////////////////////////////////////////
    // display additional information
    //////////////////////////////////////////////////////////////////////////////////

    glDisable(GL_LIGHTING);
    glColor3f(text_color[0], text_color[1], text_color[2]);
    renderText(30, 30, info_text, QFont("Helvetica", font_size, QFont::Bold));

    if(damping_enabled)
        renderText(width() - 200, height() - 70, tr("Damping enabled"), QFont("Helvetica", font_size-2, QFont::Bold));

    if(display_text_info)
    {
        char buffer[50];
        int h_pos2;
        int h_pos3 = 55;
        int h_pos1 = font_metrics->width("Coord. Number:") + font_metrics->width("100000") + h_pos3;
        int v_pos = font_metrics->width("kPa") + 10;
        int v_space = font_metrics->height();

        if(sim_info.particles > 0)
        {
            renderText(width()-h_pos1, v_pos, "Material:", font);
            h_pos2 = h_pos3 + font_metrics->width(material_name) + 5;
            renderText(width()-h_pos2, v_pos, material_name, font);
            v_pos += v_space;

            renderText(width()-h_pos1, v_pos, "Particles:", font);
            sprintf(buffer, "%i", sim_info.particles);
            h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
            renderText(width()-h_pos2, v_pos, buffer, font);
            v_pos += v_space;

            if(display_changed_contacts)
            {
                if(sim_info.created_contacts > 9999)
                    renderText(width()-h_pos1, v_pos, "Created cont.:", font);
                else
                    renderText(width()-h_pos1, v_pos, "Created contacts:", font);

                sprintf(buffer, "%i", sim_info.created_contacts);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                v_pos += v_space;

                if(sim_info.broken_contacts > 9999)
                    renderText(width()-h_pos1, v_pos, "Broken cont.:", font);
                else
                    renderText(width()-h_pos1, v_pos, "Broken contacts:", font);
                sprintf(buffer, "%i", sim_info.broken_contacts);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                v_pos += v_space;
            }
        }
        if(sim_info.collision_speed != 0)
        {
            renderText(width()-h_pos1, v_pos, "Collision speed:", font);

            if(sim_info.collision_speed >= 100)
            {
                sprintf(buffer, "%.2lf", sim_info.collision_speed * 0.01);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font );
                renderText(width()-h_pos3, v_pos, "m/s", font);
            }
            else
            {
                sprintf(buffer, "%.0lf", sim_info.collision_speed);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font );
                renderText(width()-h_pos3, v_pos, "cm/s", font);
            }

            v_pos += v_space;
        }

        v_pos += 10;

        if(sim_info.time != 0)
        {
            // determine unit of time display
            if(sim_time_display_mode == 0)
            {
                sprintf(buffer, "%.0lf", sim_info.time*1e9);
                renderText(width()-h_pos3, v_pos, "ns", font);
            }
            else if(sim_time_display_mode == 1)
            {
                sprintf(buffer, "%.2lf", sim_info.time*1e6);
                renderText(width()-h_pos3, v_pos, "\u03BCs", font);
            }
            else if(sim_time_display_mode == 2)
            {
                sprintf(buffer, "%.1lf", sim_info.time*1e6);
                renderText(width()-h_pos3, v_pos, "\u03BCs", font);
            }
            else if(sim_time_display_mode == 3)
            {
                sprintf(buffer, "%.0lf", sim_info.time*1e6);
                renderText(width()-h_pos3, v_pos, "\u03BCs", font);
            }

            renderText(width()-h_pos1, v_pos, "Time:", font);
            h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
            renderText(width()-h_pos2, v_pos, buffer, font);

            v_pos += v_space;
        }

        if(sim_info.wall_speed != 0)
        {
            renderText(width()-h_pos1, v_pos, "Wall speed:", font);

            if(sim_info.wall_speed >= 100.0)
            {
                sprintf(buffer, "%.2lf", sim_info.wall_speed * 0.01);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "m/s", font);
            }
            else
            {
                sprintf(buffer, "%.0lf", sim_info.wall_speed);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "cm/s", font);
            }

            v_pos += v_space;
        }

        if(sim_info.filling_factor != 0)
        {
            renderText(width()-h_pos1, v_pos, "Filling factor:", font);
            sprintf(buffer, "%.3lf", sim_info.filling_factor);
            h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
            renderText(width()-h_pos2, v_pos, buffer, font);
            v_pos += v_space;
        }

        if(sim_info.coordination_number != 0)
        {
            renderText(width()-h_pos1, v_pos, "Coord. Number:", font);
            sprintf(buffer, "%.2lf", sim_info.coordination_number);
            h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
            renderText(width()-h_pos2, v_pos, buffer, font);
            v_pos += v_space;
        }

        if(sim_info.pressure >= 0)
        {
            renderText(width()-h_pos1, v_pos, "Pressure:", font);

            if(sim_info.pressure >= 1e6)
            {
                sprintf(buffer, "%.2lf", 1e-6 * sim_info.pressure);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "MPa", font);
            }
            else if(sim_info.pressure >= 1e3)
            {
                sprintf(buffer, "%.2lf", 1e-3 * sim_info.pressure);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "kPa", font);
            }
            else
            {
                sprintf(buffer, "%.2lf", sim_info.pressure);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "Pa", font);
            }

            v_pos += v_space;
        }

        if(sim_info.force >= 0)
        {
            renderText(width()-h_pos1, v_pos, "Force:", font);

            if(sim_info.force >= 1e3)
            {
                sprintf(buffer, "%.2lf", 1e-3 * sim_info.force);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "kN", font);
            }
            else if(sim_info.force >= 1)
            {
                sprintf(buffer, "%.2lf", sim_info.force);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "N", font);
            }
            else if(sim_info.force >= 1e-3)
            {
                sprintf(buffer, "%.2lf", 1e3 * sim_info.force);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "mN", font);
            }
            else
            {
                sprintf(buffer, "%.2lf", 1e6 * sim_info.force);
                h_pos2 = h_pos3 + font_metrics->width(buffer) + 5;
                renderText(width()-h_pos2, v_pos, buffer, font);
                renderText(width()-h_pos3, v_pos, "\u03BCN", font);
            }

            v_pos += v_space;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // display key
    //////////////////////////////////////////////////////////////////////////////////

    if(display_key && sim_info.particles > 0)
        renderKey();

    glEnable(GL_LIGHTING);

//    doneCurrent();
}

void OpenGLWidget::renderKey()
{
    // do not render key if not initialized
    if(key_tic_labels.size() == 0)
        return;

    GLfloat key_margin_right = 10;
    GLfloat key_margin_top = 260;
    GLfloat key_margin_bottom = 50;
    GLfloat key_width = 30.0;
    GLfloat key_height = height() - key_margin_top - key_margin_bottom;

    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, width(), 0, height(), -150, 150);
    glMatrixMode(GL_MODELVIEW);

    //////////////////////////////////////////////////////////////////////////////////
    // render color stripe
    //////////////////////////////////////////////////////////////////////////////////

    glBegin(GL_QUAD_STRIP);
    float color[3];

    for(int i = 0; i < 11; ++i)
    {
        double k = 0.1 * (double)i;

        // determine color
        if(vis_display_mode == DISPLAY_DENSITY)
            SimulationVis::getColorDensity(color, k);
        else if(vis_display_mode == DISPLAY_PARTICLE_CONTACTS)
            SimulationVis::getColorCoordinationNumber(color, k);
        else if(vis_display_mode == DISPLAY_DISLOCATION)
            SimulationVis::getColorDislocation(color, k+0.02);

        glColor3f(color[0], color[1], color[2]);

        GLfloat cellHeight = key_height * (1.0 - k);
        glVertex2f(width()-key_margin_right-key_width, height()-key_margin_top-cellHeight);
        glVertex2f(width()-key_margin_right, height()-key_margin_top-cellHeight);
    }

    glEnd();

    // draw frame
    glLineWidth(1.0f);
    glColor3f(text_color[0], text_color[1], text_color[2]);

    glBegin(GL_LINE_LOOP);
    glVertex2f(width()-key_margin_right-key_width, height()-key_margin_top-key_height);
    glVertex2f(width()-key_margin_right, height()-key_margin_top-key_height);
    glVertex2f(width()-key_margin_right, height()-key_margin_top+1);
    glVertex2f(width()-key_margin_right-key_width, height()-key_margin_top+1);
    glEnd();

    //////////////////////////////////////////////////////////////////////////////////
    // display tics
    //////////////////////////////////////////////////////////////////////////////////

    glLineWidth(2.0f);
    {
        double h_pos;
        double v_pos;
        glColor3f(text_color[0], text_color[1], text_color[2]);

        for (unsigned int tic = 0; tic < key_tics; ++tic)
        {
            h_pos = width() - key_margin_right - key_width;
            v_pos = key_margin_top + key_height * (GLfloat)tic/(GLfloat)(key_tics-1);

            glBegin(GL_LINE_STRIP);
            glVertex2f(h_pos - 5.0, height() - key_margin_top - (key_height * (GLfloat)tic / (GLfloat)(key_tics-1) ));
            glVertex2f(h_pos, height() - key_margin_top - (key_height * (GLfloat)tic / (GLfloat)(key_tics-1) ));
            glEnd();

            renderText(h_pos - 8.0 - font_metrics->width(key_tic_labels[tic].data), v_pos + 0.5*(double)font_size, key_tic_labels[tic].data, font);
        }
    }

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

void OpenGLWidget::updateKey()
{

    switch(vis_display_mode)
    {
        case DISPLAY_DENSITY:
        {
            key_log_scale = false;
            key_min_value = 0;
            key_max_value = vis_density_max_filling_factor;
            key_tics = 10;

            key_tic_labels.resize(key_tics);

            for(int i = 0; i < key_tics; ++i)
            {
                sprintf(key_tic_labels[key_tics-1-i].data, "%.2lf", key_min_value + (key_max_value - key_min_value) * (double)i / (double)(key_tics-1) );
            }

            break;
        }
        case DISPLAY_PARTICLE_CONTACTS:
        {
            key_log_scale = false;
            key_min_value = 0;
            key_max_value = 10;
            key_tics = 11;
            key_tic_labels.resize(key_tics);

            for(int i = 0; i < key_tics; ++i)
            {
                sprintf(key_tic_labels[key_tics-1-i].data, "%i", i);
            }

            break;
        }
        case DISPLAY_DISLOCATION:
        {
            key_log_scale = false;
            key_min_value = vis_dislocation_min_value * delta_0;
            key_max_value = vis_dislocation_max_value * particle_radius;
            key_tics = 9;
            key_tic_labels.resize(key_tics);

            for(int i = 0; i < key_tics; ++i)
            {
                double value = key_min_value + (key_max_value - key_min_value) * 0.125 * (double)i;

                if(value < 1e-4)
                    sprintf(key_tic_labels[key_tics-1-i].data, "%5.0fnm", 1e7*value);
                else
                    sprintf(key_tic_labels[key_tics-1-i].data, "%5.2f\u03BCm", 1e4*value);
            }

            break;
        }
        default:
        {
            key_log_scale = false;
            key_min_value = 0;
            key_max_value = 0;
            key_tics = 0;
            key_tic_labels.clear();
        }
    }
    update();
}

GLuint OpenGLWidget::compileProgram(const char *vsource, const char *fsource)
{
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    glShaderSource(vertexShader, 1, &vsource, 0);
    glShaderSource(fragmentShader, 1, &fsource, 0);

    glCompileShader(vertexShader);
    glCompileShader(fragmentShader);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    glLinkProgram(program);

    // check if program linked
    GLint success = 0;
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success)
    {
        char temp[256];
        glGetProgramInfoLog(program, 256, 0, temp);
        QMessageBox::about(this, tr("Error"), tr("Preparing shader failed"));
        glDeleteProgram(program);
        program = 0;
    }

    return program;
}

void OpenGLWidget::setNumberOfParticles(int number_of_particles)
{
    // number of particles changed -> resize vertex buffer
    if(particle_array_size != number_of_particles)
    {
        // delete old vertex buffer object
        if(particle_array_size > 0)
        {
            glDeleteBuffers(1, (const GLuint*)&vbColor);
            glDeleteBuffers(1, (const GLuint*)&vbPosition);
        }

        particle_array_size = number_of_particles;

        // setup new vertex buffer object
        if(number_of_particles > 0)
        {
            unsigned int pos_size = number_of_particles * 3 * sizeof(double);
            unsigned int color_size = number_of_particles * 3 * sizeof(float);

            double* pos_buffer = new double[3 * number_of_particles];
            memset(pos_buffer, 1.0, pos_size);

            float* color_buffer = new float[3 * number_of_particles];
            memset(color_buffer, 1.0f, color_size);

            // color vbo
            glGenBuffers(1, &vbColor);
            glBindBuffer(GL_ARRAY_BUFFER, vbColor);
            glBufferData(GL_ARRAY_BUFFER, color_size, color_buffer, GL_DYNAMIC_DRAW);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            // position vbo
            glGenBuffers(1, &vbPosition);
            glBindBuffer(GL_ARRAY_BUFFER, vbPosition);
            glBufferData(GL_ARRAY_BUFFER, pos_size, pos_buffer, GL_DYNAMIC_DRAW);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            delete [] pos_buffer;
            delete [] color_buffer;
        }
    }
}

void OpenGLWidget::updateVertexBuffers()
{
    // wait until drawing is allowed
    //while(!sim_thread->can_draw){}

    // reset particle/color vbos if necessary
    if(particle_array_size != sim_thread->number_of_particles)
        setNumberOfParticles(sim_thread->number_of_particles);

    if(particle_array_size > 0)
    {
        // position buffer
        glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbPosition);
        double* p_ptr = (double*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);

        // if the pointer is valid(mapped), update VBO
        if(p_ptr)
        {
            size_t size = particle_array_size * 3 * sizeof(double);
            memcpy(p_ptr, sim_thread->particle_positions, size);

            // unmap it after use
            glUnmapBufferARB(GL_ARRAY_BUFFER_ARB);
        }

        // color buffer
        glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbColor);
        float* c_ptr = (float*)glMapBufferARB(GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);

        // if the pointer is valid(mapped), update VBO
        if(c_ptr)
        {
            size_t size = particle_array_size * 3 * sizeof(float);
            memcpy(c_ptr, sim_thread->particle_colors, size);
            glUnmapBufferARB(GL_ARRAY_BUFFER_ARB); // unmap it after use
        }
    }

    updateGL();
}
