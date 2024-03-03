 /***************************************************************************
 *   Simulation of particle collisions/agglomeration                        *
 *   Copyright (C) 2009 Alexander Seizinger                                 *
 *   alexander.seizinger[...]gmx.net                                        *
 *                                                                          *
 *   This program is free software; you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation; either version 2 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANYs WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with this program; if not, write to the                          *
 *   Free Software Foundation, Inc.                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              *
 ***************************************************************************/
//#include <QDebug>

#include "SimulationLib.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <list>

#include <vector>



extern double particle_radius;
extern double reduced_radius;
extern double density;
extern double mass;
extern double mass_inv;
extern double moment_of_inertia;

extern double surface_energy;
extern double nu;
extern double young_mod;
extern double crit_rolling_displacement;
extern double crit_rolling_displacement_squared;
extern double crit_sliding_displacement_squared;
extern double osc_damping_factor;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;
extern double crit_sliding_displacement_modifier;
extern double crit_wall_sliding_displacement_modifier;

extern double k_r;
extern double k_s;
extern double k_t;
extern double wall_k_s;
extern double wall_k_r;
extern double wall_k_t;
extern double ENERGY_UNIT;
extern double delta_0;
extern double delta_c;
extern double F_c;
extern double equilibrium_distance;

extern bool damping_enabled;

#include <clocale>


ErrorCode SimLib::disturbParticles(Simulation *sim, double max_dist)
{

    // get number of contacts per particle
    std::vector< int > contacts(sim->number_of_particles);
    ContactListEntry *cl_entry = NULL;

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        cl_entry = sim->contact_list[p];

        while(cl_entry)
        {
            // ignore contacts with walls
            if(cl_entry->id >= 0)
            {
                ++contacts[p];
                ++contacts[cl_entry->id];
            }

            cl_entry = cl_entry->next;
        }
    }



    sim->deleteContacts();

    for(int p = 0; p < sim->number_of_particles; ++p)
    {

        if(contacts[p] <= 3) // dont stirr particles that are not elastically charged
            continue;

        // determine random direction
        double phi = sim->get_random_zero_twoPi();
        double cos_theta = sim->get_random_cos_theta();
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        vec3 dir;
        dir[0] = sin_theta * cos(phi);
        dir[1] = sin_theta * sin(phi);
        dir[2] = cos_theta;

        double dist = max_dist * sim->get_random_zero_one();

        sim->pos_old[X_COORD(p)] += dir[0] * dist;
        sim->pos_old[Y_COORD(p)] += dir[1] * dist;
        sim->pos_old[Z_COORD(p)] += dir[2] * dist;
    }

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    return EC_OK;
}

ErrorCode SimLib::duplicateSlices(Simulation *sim, int x_duplications, int y_duplications, int z_duplications, bool mirror, bool random_orientation)
{
    if(sim->number_of_particles < 1)
        return EC_NO_PARTICLES;

    vec3 center_of_mass;
    SimLib::getCenterOfMass(&center_of_mass, *sim, 0, sim->number_of_particles-1);

    // translate agglomerates to make sure cms = 0
    SimLib::centerCMS(sim);

    if(random_orientation)
        SimLib::rotateSimRandomly(sim);

    // determine new size
    double x_min = sim->pos_old[0];
    double x_max = sim->pos_old[0];
    double y_min = sim->pos_old[1];
    double y_max = sim->pos_old[1];
    double z_min = sim->pos_old[2];
    double z_max = sim->pos_old[2];

    for(int p = 1; p < sim->number_of_particles; ++p)
    {
        if(sim->pos_old[X_COORD(p)] > x_max)
            x_max = sim->pos_old[X_COORD(p)];

        if(sim->pos_old[X_COORD(p)] < x_min)
            x_min = sim->pos_old[X_COORD(p)];

        if(sim->pos_old[Y_COORD(p)] > y_max)
            y_max = sim->pos_old[Y_COORD(p)];

        if(sim->pos_old[Y_COORD(p)] < y_min)
            y_min = sim->pos_old[Y_COORD(p)];

        if(sim->pos_old[Z_COORD(p)] > z_max)
            z_max = sim->pos_old[Z_COORD(p)];

        if(sim->pos_old[Z_COORD(p)] < z_min)
            z_min = sim->pos_old[Z_COORD(p)];
    }

    x_min -= particle_radius;
    x_max += particle_radius;
    y_min -= particle_radius;
    y_max += particle_radius;
    z_min -= particle_radius;
    z_max += particle_radius;


    // translate -> particles are located between 0 and x/y/z_max-x/y/z_min
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        sim->pos_old[X_COORD(p)] -= x_min;
        sim->pos_old[Y_COORD(p)] -= y_min;
        sim->pos_old[Z_COORD(p)] -= z_min;
    }

    sim->saveToFile("temp_agglomerate.dat");

    Simulation sim2;
    sim2.loadFromFile("temp_agglomerate.dat");

    if(mirror && sim2.contact_list)
    {
        ContactListEntry *cl_entry, *last_entry;

        for(int p = 0; p < sim2.number_of_particles; ++p)
        {
            cl_entry = sim2.contact_list[p];

            while(cl_entry)
            {
                last_entry = cl_entry;
                cl_entry = cl_entry->next;

                delete last_entry->contact;
                delete last_entry;
            }

            sim2.contact_list[p] = NULL;
        }
    }

    // for mirror
    double x_center = (x_max - x_min)/2.0;
    double y_center = (y_max - y_min)/2.0;
    double z_center = (z_max - z_min)/2.0;

    //double particle_dist = 2.0 * particle_radius; //particle_radius;
    int x_duplications_initial = x_duplications;

    // duplicate in x direction
    for(int z = 0; z < z_duplications+1; ++z)
    {
        for(int y = 0; y < y_duplications+1; ++y)
        {
            for(int x = 0; x < x_duplications; ++x)
            {
                for(int p = 0; p < sim->number_of_particles; ++p)
                    sim->pos_old[X_COORD(p)] -= (x_max - x_min);

                if(mirror)
                {
                    for(int p = 0; p < sim2.number_of_particles; ++p)
                        sim2.pos_old[X_COORD(p)] = x_center - (sim2.pos_old[X_COORD(p)] - x_center);
                }

                sim->addParticlesFromSim(sim2);
            }

            for(int p = 0; p < sim->number_of_particles; ++p)
            {
                sim->pos_old[X_COORD(p)] += (double)(x_duplications_initial+1) * (x_max - x_min);   // reset x_pos
                sim->pos_old[Y_COORD(p)] -= (y_max - y_min);
            }

            if(mirror)
            {
                for(int p = 0; p < sim2.number_of_particles; ++p)
                    sim2.pos_old[Y_COORD(p)] = y_center - (sim2.pos_old[Y_COORD(p)] - y_center);
            }

            if(y == 0 && z == 0)
                ++x_duplications;
        }

        for(int p = 0; p < sim->number_of_particles; ++p)
        {
            sim->pos_old[Y_COORD(p)] += (double)(y_duplications+1) * (y_max - y_min);   // reset y pos
            sim->pos_old[Z_COORD(p)] -= (z_max - z_min);
        }

        if(mirror)
        {
            for(int p = 0; p < sim2.number_of_particles; ++p)
                sim2.pos_old[Z_COORD(p)] = z_center - (sim2.pos_old[Z_COORD(p)] - z_center);
        }

    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

ErrorCode SimLib::particleClusterAggregation(Simulation *sim, int final_cluster_size, double fractal_prefactor, double fractal_dimension)
{
    if(final_cluster_size < 2)
        return EC_INVALID_PARAMETER;

    // prepare sim
    sim->resizeArrays(final_cluster_size, 0);
    sim->number_of_walls = 0;
    sim->walls.clear();

    // prepare grid
    sim->grid.init(256, 2.1 * particle_radius, final_cluster_size);
    sim->grid.resetGrid();

    vec3 new_pos, cms, delta_r, m;
    double desired_dist, desired_dist_squared;
    double dist_squared, min_dist_squared, max_dist_squared;
    double k;
    bool particle_added;
    std::list<int> particles;
    double particle_dist = 2.0 * particle_radius - delta_0;

    // place first particle
    sim->pos_old[0] = 0; sim->pos_old[1] = 0; sim->pos_old[2] = 0;
    new_pos[0] = 0; new_pos[1] = 0; new_pos[2] = 0;
    sim->grid.addParticle(new_pos, 0);

    sim->pos_old[X_COORD(1)] = particle_dist; sim->pos_old[Y_COORD(1)] = 0; sim->pos_old[Z_COORD(1)] = 0;
    new_pos[0] = particle_dist; new_pos[1] = 0; new_pos[2] = 0;
    sim->grid.addParticle(new_pos, 1);

    for(int particle_count = 3; particle_count <= final_cluster_size; ++particle_count)
    {
        particle_added = false;

        // get center of mass
        SimLib::getCenterOfMass(&cms, *sim, 0, particle_count-2);

        // determine required distance of new particle from cms
        desired_dist_squared = (double)(particle_count*particle_count) * particle_radius*particle_radius / (double)(particle_count-1) * pow((double)particle_count / fractal_prefactor, 2.0 / fractal_dimension)
                                - (double)(particle_count) * particle_radius*particle_radius / (double)(particle_count-1)
                                - (double)(particle_count) * particle_radius*particle_radius * pow((double)(particle_count-1)/fractal_prefactor, 2.0 / fractal_dimension);

        desired_dist = sqrt(desired_dist_squared);
        min_dist_squared = (desired_dist - particle_dist) * (desired_dist - particle_dist);
        max_dist_squared = (desired_dist + particle_dist) * (desired_dist + particle_dist);

        // geneate particle list
        particles.clear();

        if( particle_count%3 == 0 && particle_count > 20)
        {
            for(int p = 0; p < particle_count-1; ++p)
                particles.push_back(p);
        }
        else
        {
            for(int p = particle_count-2; p >= 0; --p)
                particles.push_back(p);
        }

        // find suitable particles where a new particle can be added
        //for(int p = 0; p < particle_count-1; ++p)
        for(std::list<int>::iterator particle = particles.begin(); particle != particles.end(); ++particle)
        {
            int p = *particle;

            delta_r[0] = sim->pos_old[X_COORD(p)] - cms[0];
            delta_r[1] = sim->pos_old[Y_COORD(p)] - cms[1];
            delta_r[2] = sim->pos_old[Z_COORD(p)] - cms[2];

             dist_squared = norm_squared(delta_r);

            // add particle to list if its center is close enough to desired distance
            if(dist_squared > min_dist_squared && dist_squared < max_dist_squared)
            {
                // determine center of intersection circle
                double temp = sim->pos_old[X_COORD(p)]*sim->pos_old[X_COORD(p)] + sim->pos_old[Y_COORD(p)]*sim->pos_old[Y_COORD(p)] + sim->pos_old[Z_COORD(p)]*sim->pos_old[Z_COORD(p)];
                k = ( 0.5 * (desired_dist_squared - particle_dist*particle_dist + temp - norm_squared(cms)) - dot_product(cms, delta_r) ) / dist_squared;

                m[0] = cms[0] + k * delta_r[0];
                m[1] = cms[1] + k * delta_r[1];
                m[2] = cms[2] + k * delta_r[2];

                // determine radius of intersection circle
                double radius = sqrt(particle_dist*particle_dist - (sim->pos_old[X_COORD(p)]-m[0])*(sim->pos_old[X_COORD(p)]-m[0]) - (sim->pos_old[Y_COORD(p)]-m[1])*(sim->pos_old[Y_COORD(p)]-m[1]) - (sim->pos_old[Z_COORD(p)]-m[2])*(sim->pos_old[Z_COORD(p)]-m[2]) );

                // get random point on intersection circle and try to place new particle there
                int counter = 0;
                while(counter < 10)
                {
                    ++counter;

                    vec3 n;
                    getOrthoVector(sim, &n, delta_r);

                    new_pos[0] = m[0] + radius * n[0];
                    new_pos[1] = m[1] + radius * n[1];
                    new_pos[2] = m[2] + radius * n[2];

                    // check if other particles are too close
                    if( sim->grid.canAddParticleAt(new_pos, sim->pos_old) )
                    {
                        sim->grid.addParticle(new_pos, particle_count-1);
                        sim->pos_old[X_COORD(particle_count-1)] = new_pos[0];
                        sim->pos_old[Y_COORD(particle_count-1)] = new_pos[1];
                        sim->pos_old[Z_COORD(particle_count-1)] = new_pos[2];

                        particle_added = true;
                        break;
                    }
                }
            }

            if(particle_added)
                break;
        }

        // abort if unable to add particle
        if(!particle_added)
            break;
    }

    // debug
    /*double r_g = 0;

    cms = SimLib::getCenterOfMass(sim, 0, sim->number_of_particles-1);

    for(int p = 0; p < sim->number_of_particles; ++p)
    r_g += (sim->pos_old[X_COORD(p)]-cms.x)*(sim->pos_old[X_COORD(p)]-cms.x) + (sim->pos_old[Y_COORD(p)]-cms.y)*(sim->pos_old[Y_COORD(p)]-cms.y) + (sim->pos_old[Z_COORD(p)]-cms.z)*(sim->pos_old[Z_COORD(p)]-cms.z);

    r_g /= (double)sim->number_of_particles;

    r_g = sqrt(r_g);

    double N = fractal_prefactor * pow(r_g / particle_radius, fractal_dimension);*/

    sim->initSimData(SIM_TYPE_GENERAL);
    return EC_OK;
}

void invert(double *a)
{
    double a11 = a[0];
    double a12 = a[1];
    double a13 = a[2];
    double a21 = a[3];
    double a22 = a[4];
    double a23 = a[5];
    double a31 = a[6];
    double a32 = a[7];
    double a33 = a[8];

    double det = a11 * (a22 * a33 - a23 * a32) - a12 * (a33 * a21 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
    double det_inv = 1.0 / det;

    if(det != 0)
    {
        a[0] = det_inv * (a22*a33 - a23*a32);
        a[1] = det_inv * (a13*a32 - a12*a33);
        a[2] = det_inv * (a12*a23 - a13*a22);
        a[3] = det_inv * (a23*a31 - a21*a33);
        a[4] = det_inv * (a11*a33 - a13*a31);
        a[5] = det_inv * (a13*a21 - a11*a23);
        a[6] = det_inv * (a21*a32 - a22*a31);
        a[7] = det_inv * (a31*a12 - a11*a32);
        a[8] = det_inv * (a11*a22 - a12*a21);
    }
}

ErrorCode SimLib::initBAMAggregate(Simulation *sim, const char *filename, unsigned int number_of_particles, double migration_probability1, double migration_probability2, BAMSelectionMethod bam_selection_method)
{
    if(number_of_particles < 1)
        return EC_INVALID_PARAMETER;

    if(migration_probability1 < 0 || migration_probability1 > 1 || migration_probability2 < 0 || migration_probability2 > 1)
        return EC_INVALID_PARAMETER;

    double eq_dist = 2.0 * particle_radius - delta_0;
    double outer_radius = 0;

    if(filename)
    {
        sim->loadFromFile(filename);
        sim->removeWalls();
        SimLib::centerCMS(sim);
    }

    SimLib::getSize(*sim, NULL, &outer_radius);

    int first_add_particle_id = sim->number_of_particles;

    sim->addParticles(number_of_particles);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // place first two particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    vec3 pos = {0,0,0};

    if(first_add_particle_id == 0)
    {
        sim->pos_old[X_COORD(first_add_particle_id)] = pos[0];
        sim->pos_old[Y_COORD(first_add_particle_id)] = pos[1];
        sim->pos_old[Z_COORD(first_add_particle_id)] = pos[2];
        sim->grid.addParticle(pos, first_add_particle_id);
        ++first_add_particle_id;

        outer_radius = particle_radius;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // place other particles
    /////////////////////////////////////////////////////////////////////////////////////////////


    double done = 0.0;


    for(int p = first_add_particle_id; p < sim->number_of_particles; ++p)
    {

        if(p - first_add_particle_id > int(done*(sim->number_of_particles - first_add_particle_id))-1)
        {
            printf("Generating agglomerate %d percent done %d/%d.\n", (int)(done*100), p, sim->number_of_particles);
            done += 0.1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // shoot particle from random direction at aggregate
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        int p1 = -1;

        // determine random direction
        double phi = sim->get_random_zero_twoPi();
        double cos_theta = sim->get_random_cos_theta();
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        vec3 dir;
        dir[0] = sin_theta * cos(phi);
        dir[1] = sin_theta * sin(phi);
        dir[2] = cos_theta;

        // select random direction
        pos[0] = - (outer_radius + 4.0 * particle_radius) * dir[0];
        pos[1] = - (outer_radius + 4.0 * particle_radius) * dir[1];
        pos[2] = - (outer_radius + 4.0 * particle_radius) * dir[2];

        // determine next collision with existing particles
        double min_value = 1.e100;


        for(int j = 0; j < p; ++j)
        {
            double a = (sim->pos_old[X_COORD(j)]-pos[0]) * dir[0] + (sim->pos_old[Y_COORD(j)]-pos[1]) * dir[1] + (sim->pos_old[Z_COORD(j)]-pos[2]) * dir[2];
            double b = (sim->pos_old[X_COORD(j)]-pos[0])*(sim->pos_old[X_COORD(j)]-pos[0]) + (sim->pos_old[Y_COORD(j)]-pos[1])*(sim->pos_old[Y_COORD(j)]-pos[1]) + (sim->pos_old[Z_COORD(j)]-pos[2])*(sim->pos_old[Z_COORD(j)]-pos[2]);
            double c = a*a-b+eq_dist*eq_dist;

            if(c >= 0)
            {
                double t = a - sqrt(c);
                if(t < min_value)
                {
                    min_value = t;
                    p1 = j;
                }
            }
        }

        if(p1 < 0)
        {
            sim->cleanUp();
            return EC_INTERNAL_ERROR;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check migration
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        bool migration_successful = false;

        if(migration_probability1 >= sim->get_random_zero_one_incl() && p > 1)
        {
            vec3 default_pos;
            default_pos[0] = (min_value - (outer_radius + 4.0 * particle_radius)) * dir[0];
            default_pos[1] = (min_value - (outer_radius + 4.0 * particle_radius)) * dir[1];
            default_pos[2] = (min_value - (outer_radius + 4.0 * particle_radius)) * dir[2];

            if(migration_probability2 >= sim->get_random_zero_one_incl() && p > 3)
                migration_successful = SimLib::getMigrationPos2(&pos, *sim, p1, default_pos, bam_selection_method);
            else
                migration_successful = SimLib::getMigrationPos1(&pos, *sim, p1, default_pos, bam_selection_method);
        }

        // no migration - add particle to default pos
        if(!migration_successful)
        {
            pos[0] = (min_value - (outer_radius + 4.0 * particle_radius)) * dir[0];
            pos[1] = (min_value - (outer_radius + 4.0 * particle_radius)) * dir[1];
            pos[2] = (min_value - (outer_radius + 4.0 * particle_radius)) * dir[2];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // add particle to sim
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        sim->pos_old[X_COORD(p)] = pos[0];
        sim->pos_old[Y_COORD(p)] = pos[1];
        sim->pos_old[Z_COORD(p)] = pos[2];
        sim->grid.addParticle(pos, p);

        double radius = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];

        if(radius > outer_radius*outer_radius)
            outer_radius = sqrt(radius);
    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // remove misplaced particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    int removed_particles = 0;
    bool *remove_list = new bool[sim->number_of_particles];

    double max_deviation = 0;
    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        remove_list[p] = false;

        ContactListEntry *next_cl_entry = sim->contact_list[p];

        while(next_cl_entry)
        {
            double dist_squared = (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)]) * (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)]) * (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]) * (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]);


            int deviation = fabs( (dist_squared - eq_dist_squared)  / eq_dist_squared ) > 0.01;

            if(deviation > max_deviation)
                max_deviation = deviation;

            if(deviation > 0.01 )
            {
                ++removed_particles;
                remove_list[p] = true;

                next_cl_entry = NULL;
            }
            else
                next_cl_entry = next_cl_entry->next;
        }
    }

    sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();
    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    printf("100 percent .\n");

    return EC_OK;
}

bool SimLib::getNeighbourPos(vec3 *pos, Simulation &sim, int particle_id)
{
    for(int i = 0; i < 64; ++i)
    {
        // determine random direction
        double phi = sim.get_random_zero_twoPi();
        double cos_theta = sim.get_random_cos_theta();
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        vec3 dir;
        dir[0] = sin_theta * cos(phi);
        dir[1] = sin_theta * sin(phi);
        dir[2] = cos_theta;

        (*pos)[0] = sim.pos_old[X_COORD(particle_id)] - equilibrium_distance * dir[0];
        (*pos)[1] = sim.pos_old[Y_COORD(particle_id)] - equilibrium_distance * dir[1];
        (*pos)[2] = sim.pos_old[Z_COORD(particle_id)] - equilibrium_distance * dir[2];

        if(sim.grid.canAddParticleAt(*pos, sim.pos_old))
            return true;

        (*pos)[0] = sim.pos_old[X_COORD(particle_id)] + equilibrium_distance * dir[0];
        (*pos)[1] = sim.pos_old[Y_COORD(particle_id)] + equilibrium_distance * dir[1];
        (*pos)[2] = sim.pos_old[Z_COORD(particle_id)] + equilibrium_distance * dir[2];

        if(sim.grid.canAddParticleAt(*pos, sim.pos_old))
            return true;
    }

    return false;
}

bool SimLib::getMigrationPos1(vec3 *migration_pos, Simulation &sim, int particle_id, vec3 &default_pos, BAMSelectionMethod bam_selection_method)
{
    // get neighbouring particles
    std::list<int> neighbours;
    sim.grid.getNeighbours(&neighbours, sim.pos_old, particle_id, 2.0 * equilibrium_distance);

    if( sim.get_random_zero_one() > 0.5)
        neighbours.reverse();

    double min_dist = 100.0;
    vec3 best_pos;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // check migration
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(std::list<int>::iterator p2 = neighbours.begin(); p2 != neighbours.end(); ++p2)
    {
        // determine contact normal/distance between p1 and p2
        vec3 n_c;
        n_c[0] = sim.pos_old[X_COORD(*p2)] - sim.pos_old[X_COORD(particle_id)];
        n_c[1] = sim.pos_old[Y_COORD(*p2)] - sim.pos_old[Y_COORD(particle_id)];
        n_c[2] = sim.pos_old[Z_COORD(*p2)] - sim.pos_old[Z_COORD(particle_id)];

        double d = sqrt(n_c[0]*n_c[0] + n_c[1]*n_c[1] +n_c[2]*n_c[2]);

        n_c[0] /=d;
        n_c[1] /=d;
        n_c[2] /=d;

        vec3 center;
        center[0] = sim.pos_old[X_COORD(particle_id)] + 0.5 * d * n_c[0];
        center[1] = sim.pos_old[Y_COORD(particle_id)] + 0.5 * d * n_c[1];
        center[2] = sim.pos_old[Z_COORD(particle_id)] + 0.5 * d * n_c[2];

        vec3 k1, k2;
        SimLib::getOrthoVector(&sim, &k1, n_c);
        k2[0] = n_c[1] * k1[2] - n_c[2] * k1[1];
        k2[1] = n_c[2] * k1[0] - n_c[0] * k1[2];
        k2[2] = n_c[0] * k1[1] - n_c[1] * k1[0];

        double a = sqrt(equilibrium_distance*equilibrium_distance - 0.25 * d*d);
        double phi_0 = 0; //2.0 * M_PI * (double)rand() / (double)(RAND_MAX+1);

        for(int i = 0; i < 64; ++i)
        {
            double phi = phi_0 + 2.0 * M_PI * (double)i / 64.0;
            double a1 = a * cos(phi);
            double a2 = a * sin(phi);

            (*migration_pos)[0] = center[0] + a1 * k1[0] + a2 * k2[0];
            (*migration_pos)[1] = center[1] + a1 * k1[1] + a2 * k2[1];
            (*migration_pos)[2] = center[2] + a1 * k1[2] + a2 * k2[2];

            if(sim.grid.canAddParticleAt(*migration_pos, sim.pos_old))
            {
                if(bam_selection_method == BAM_SELECT_RANDOM)
                    return true;
                else
                {
                    double dist;

                    if(bam_selection_method == BAM_SELECT_CLOSEST)
                        dist = ((*migration_pos)[0] - default_pos[0])*((*migration_pos)[0] - default_pos[0]) + ((*migration_pos)[1] - default_pos[1])*((*migration_pos)[1] - default_pos[1]) + ((*migration_pos)[2] - default_pos[2])*((*migration_pos)[2] - default_pos[2]);
                    else
                        dist = (*migration_pos)[0]*(*migration_pos)[0] + (*migration_pos)[1]*(*migration_pos)[1] + (*migration_pos)[2]*(*migration_pos)[2];

                    if(dist < min_dist)
                    {
                        min_dist = dist;
                        best_pos[0] = (*migration_pos)[0];
                        best_pos[1] = (*migration_pos)[1];
                        best_pos[2] = (*migration_pos)[2];
                    }
                }
            }
        }
    }

    if(min_dist < 99.0)
    {
        (*migration_pos)[0] = best_pos[0];
        (*migration_pos)[1] = best_pos[1];
        (*migration_pos)[2] = best_pos[2];
        return true;
    }
    else
        return false;
}

bool SimLib::getMigrationPos2(vec3 *migration_pos, Simulation &sim, int particle_id, vec3 &default_pos, BAMSelectionMethod bam_selection_method)
{
    // get neighbouring particles
    std::list<int> neighbours;
    sim.grid.getNeighbours(&neighbours, sim.pos_old, particle_id, 2.0 * equilibrium_distance);

    if( sim.get_random_zero_one() > 0.5)
        neighbours.reverse();

    double min_dist = 100.0;
    vec3 best_pos;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // do double migration
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(std::list<int>::iterator p2 = neighbours.begin(); p2 != neighbours.end(); ++p2)
    {
        // determine contact normal/distance between p1 and p2
        vec3 n_c;
        n_c[0] = sim.pos_old[X_COORD(*p2)] - sim.pos_old[X_COORD(particle_id)];
        n_c[1] = sim.pos_old[Y_COORD(*p2)] - sim.pos_old[Y_COORD(particle_id)];
        n_c[2] = sim.pos_old[Z_COORD(*p2)] - sim.pos_old[Z_COORD(particle_id)];

        double d12 = sqrt(n_c[0]*n_c[0] + n_c[1]*n_c[1] + n_c[2]*n_c[2]);

        n_c[0] /= d12;
        n_c[1] /= d12;
        n_c[2] /= d12;

        vec3 center1;
        center1[0] = sim.pos_old[X_COORD(particle_id)] + 0.5 * d12 * n_c[0];
        center1[1] = sim.pos_old[Y_COORD(particle_id)] + 0.5 * d12 * n_c[1];
        center1[2] = sim.pos_old[Z_COORD(particle_id)] + 0.5 * d12 * n_c[2];

        double a1 = sqrt(equilibrium_distance*equilibrium_distance - 0.25 * d12*d12);

        for(std::list<int>::iterator p3 = neighbours.begin(); p3 != neighbours.end(); ++p3)
        {
            if(*p3 != *p2)
            {
                double k = (center1[0] - sim.pos_old[X_COORD(*p3)]) * n_c[0] + (center1[1] - sim.pos_old[Y_COORD(*p3)]) * n_c[1] + (center1[2] - sim.pos_old[Z_COORD(*p3)]) * n_c[2];

                vec3 center2;
                center2[0] = sim.pos_old[X_COORD(*p3)] + k * n_c[0];
                center2[1] = sim.pos_old[Y_COORD(*p3)] + k * n_c[1];
                center2[2] = sim.pos_old[Z_COORD(*p3)] + k * n_c[2];

                double a2 = equilibrium_distance*equilibrium_distance - (center2[0] - sim.pos_old[X_COORD(*p3)])*(center2[0] - sim.pos_old[X_COORD(*p3)])
                            - (center2[1] - sim.pos_old[Y_COORD(*p3)])*(center2[1] - sim.pos_old[Y_COORD(*p3)]) - (center2[2] - sim.pos_old[Z_COORD(*p3)])*(center2[2] - sim.pos_old[Z_COORD(*p3)]);

                if(a2 >= 0)
                {
                    a2 = sqrt(a2);
                    k = a1*a1-a2*a2;

                    vec3 delta_center;
                    delta_center[0] = center2[0] - center1[0];
                    delta_center[1] = center2[1] - center1[1];
                    delta_center[2] = center2[2] - center1[2];

                    double b_squared = delta_center[0]*delta_center[0] + delta_center[1]*delta_center[1] + delta_center[2]*delta_center[2];
                    double c = 0.5 + 0.5 * k / b_squared;
                    double a3 = a1*a1 - 0.25*b_squared - 0.5 * k - 0.25 * k*k / b_squared;

                    if(a3 >= 0)
                    {
                        a3 = sqrt(a3);

                        vec3 n;
                        n[0] = n_c[1] * delta_center[2] - n_c[2] * delta_center[1];
                        n[1] = n_c[2] * delta_center[0] - n_c[0] * delta_center[2];
                        n[2] = n_c[0] * delta_center[1] - n_c[1] * delta_center[0];
                        normalize(&n);

                        // try to add particle at first possible position
                        (*migration_pos)[0] = center1[0] + c * delta_center[0] + a3 * n[0];
                        (*migration_pos)[1] = center1[1] + c * delta_center[1] + a3 * n[1];
                        (*migration_pos)[2] = center1[2] + c * delta_center[2] + a3 * n[2];

                        if(sim.grid.canAddParticleAt(*migration_pos, sim.pos_old))
                        {
                            if(bam_selection_method == BAM_SELECT_RANDOM)
                                return true;
                            else
                            {
                                // valid migration pos found -> compare with other possible positions when center or shortest migration is selected
                                double dist;

                                if(bam_selection_method == BAM_SELECT_CLOSEST)
                                    dist = ((*migration_pos)[0] - default_pos[0])*((*migration_pos)[0] - default_pos[0]) + ((*migration_pos)[1] - default_pos[1])*((*migration_pos)[1] - default_pos[1]) + ((*migration_pos)[2] - default_pos[2])*((*migration_pos)[2] - default_pos[2]);
                                else
                                    dist = (*migration_pos)[0]*(*migration_pos)[0] + (*migration_pos)[1]*(*migration_pos)[1] + (*migration_pos)[2]*(*migration_pos)[2];

                                if(dist < min_dist)
                                {
                                    min_dist = dist;
                                    best_pos[0] = (*migration_pos)[0];
                                    best_pos[1] = (*migration_pos)[1];
                                    best_pos[2] = (*migration_pos)[2];
                                }
                            }
                        }

                        // try to add particle at second possible position
                        (*migration_pos)[0] = center1[0] + c * delta_center[0] - a3 * n[0];
                        (*migration_pos)[1] = center1[1] + c * delta_center[1] - a3 * n[1];
                        (*migration_pos)[2] = center1[2] + c * delta_center[2] - a3 * n[2];

                        if(sim.grid.canAddParticleAt(*migration_pos, sim.pos_old))
                        {
                            if(bam_selection_method == BAM_SELECT_RANDOM)
                                return true;
                            else
                            {
                                // valid migration pos found -> compare with other possible positions when center or shortest migration is selected
                                double dist;

                                if(bam_selection_method == BAM_SELECT_CLOSEST)
                                    dist = ((*migration_pos)[0] - default_pos[0])*((*migration_pos)[0] - default_pos[0]) + ((*migration_pos)[1] - default_pos[1])*((*migration_pos)[1] - default_pos[1]) + ((*migration_pos)[2] - default_pos[2])*((*migration_pos)[2] - default_pos[2]);
                                else
                                    dist = (*migration_pos)[0]*(*migration_pos)[0] + (*migration_pos)[1]*(*migration_pos)[1] + (*migration_pos)[2]*(*migration_pos)[2];

                                if(dist < min_dist)
                                {
                                    min_dist = dist;
                                    best_pos[0] = (*migration_pos)[0];
                                    best_pos[1] = (*migration_pos)[1];
                                    best_pos[2] = (*migration_pos)[2];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(min_dist < 99.0)
    {
        (*migration_pos)[0] = best_pos[0];
        (*migration_pos)[1] = best_pos[1];
        (*migration_pos)[2] = best_pos[2];
        return true;
    }
    else
        return false;
}

int getUppermostParticle(Simulation *sim, double pos_x, double pos_z, double y_max, double y_min)
{
    // determine start cell
    int center_cell_id = sim->grid.getCellID(pos_x, y_max, pos_z);
    int max_y_cell_id = sim->grid.getYID(y_max);
    int min_y_cell_id = sim->grid.getYID(y_min);
    double current_y_pos = y_min;
    int closest_id = -1;

    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);
    bool stop_search = false;

    // check cells around the current cell for particles
    // if none are found, decrease y and continue search in next level
    for(int y = max_y_cell_id; y >= min_y_cell_id; --y)
    {
        for(int x = -1; x <= 1; ++x)
        {
            for(int z = -1; z <= 1; ++z)
            {
                // check particles in this cell
                GridEntry *grid_entry = sim->grid.particles_in_cell[center_cell_id + x + (y-max_y_cell_id) * sim->grid.x_cells + z * sim->grid.x_cells * sim->grid.y_cells];

                while(grid_entry)
                {
                    int id = grid_entry->particle;

                    // calc dist to "path" of falling particle
                    double dist_squared = (pos_x - sim->pos_old[X_COORD(id)])*(pos_x - sim->pos_old[X_COORD(id)]) + (pos_z - sim->pos_old[Z_COORD(id)])*(pos_z - sim->pos_old[Z_COORD(id)]);

                    // only consider particles that are close enough to establish a contact
                    if(dist_squared < eq_dist_squared)
                    {
                        // calc pos of new particle
                        double k = sim->pos_old[Y_COORD(id)]*sim->pos_old[Y_COORD(id)] + (sim->pos_old[X_COORD(id)]-pos_x)*(sim->pos_old[X_COORD(id)]-pos_x) + (sim->pos_old[Z_COORD(id)]-pos_z)*(sim->pos_old[Z_COORD(id)]-pos_z) - eq_dist_squared;

                        double y_pos = sim->pos_old[Y_COORD(id)] + sqrt(sim->pos_old[Y_COORD(id)]*sim->pos_old[Y_COORD(id)] - k);

                        if(y_pos > current_y_pos)
                        {
                            current_y_pos = y_pos;
                            closest_id = id;
                        }
                    }

                    grid_entry = grid_entry->next;
                }
            }
        }

        if(stop_search)
            return closest_id;

        // if particle has been found, end search after checking the subsequent slice
        // (particles in slices with lower y id will be too far away)
        if(closest_id >= 0)
            stop_search = true;
    }

    if(closest_id >= 0)
        return closest_id;
    else
        return -1;
}


ErrorCode SimLib::initBAMCake(Simulation *sim, unsigned int number_of_particles, double x_size, double y_size, double migration_probability1, double migration_probability2)
{
    if(number_of_particles < 2)
        return EC_INVALID_PARAMETER;

    if(migration_probability1 < 0 || migration_probability2 < 0)
        return EC_INVALID_PARAMETER;

    double eq_dist = 2.0 * particle_radius - delta_0;
    double current_height = 0;

    sim->resizeArrays(number_of_particles, 0);

    double x_min = particle_radius;
    double x_max = x_size - particle_radius;
    double y_min = particle_radius;
    double y_max = y_size - particle_radius;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // place particles
    /////////////////////////////////////////////////////////////////////////////////////////////


    double done = 0.0;

    for(int p = 0; p < sim->number_of_particles; ++p)
    {


        if(p > int(done*sim->number_of_particles))
        {
            printf("%d percent done.\n", (int)(done*100));
            done += 0.1;
        }

        // determine random position
        vec3 pos, init_pos;
        init_pos[0] = x_min + (x_max - x_min) * sim->get_random_zero_one_incl();
        init_pos[1] = current_height + 1.001 * eq_dist;
        init_pos[2] = y_min + (y_max - y_min) * sim->get_random_zero_one_incl();

        int p1 = getUppermostParticle(sim, init_pos[0], init_pos[2], init_pos[1], 0);

        // no particle on the way to bottom found
        if(p1 >= 0)
        {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check migration
        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(migration_probability1 == 1.0 || migration_probability1 >= sim->get_random_zero_one())
        {
            // get particles that are close to p1
            std::list<int> neighbours;
            sim->grid.getNeighbours(&neighbours, sim->pos_old, p1, 2.0 * eq_dist);

            if( sim->get_random_zero_one() > 0.5)
                neighbours.reverse();

            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            // do double migration
            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            if(migration_probability2 == 1.0 || migration_probability2 >= sim->get_random_zero_one())
            {
                for(std::list<int>::iterator p2 = neighbours.begin(); p2 != neighbours.end(); ++p2)
                {
                    // determine contact normal/distance between p1 and p2
                    vec3 n_c;
                    n_c[0] = sim->pos_old[X_COORD(*p2)] - sim->pos_old[X_COORD(p1)];
                    n_c[1] = sim->pos_old[Y_COORD(*p2)] - sim->pos_old[Y_COORD(p1)];
                    n_c[2] = sim->pos_old[Z_COORD(*p2)] - sim->pos_old[Z_COORD(p1)];

                    double d12 = sqrt(n_c[0]*n_c[0] + n_c[1]*n_c[1] + n_c[2]*n_c[2]);

                    n_c[0] /= d12;
                    n_c[1] /= d12;
                    n_c[2] /= d12;

                    vec3 center1;
                    center1[0] = sim->pos_old[X_COORD(p1)] + 0.5 * d12 * n_c[0];
                    center1[1] = sim->pos_old[Y_COORD(p1)] + 0.5 * d12 * n_c[1];
                    center1[2] = sim->pos_old[Z_COORD(p1)] + 0.5 * d12 * n_c[2];

                    double a1 = sqrt(eq_dist*eq_dist - 0.25 * d12*d12);

                    for(std::list<int>::iterator p3 = neighbours.begin(); p3 != neighbours.end(); ++p3)
                    {
                        if(*p3 != *p2)
                        {
                            double k = (center1[0] - sim->pos_old[X_COORD(*p3)]) * n_c[0] + (center1[1] - sim->pos_old[Y_COORD(*p3)]) * n_c[1] + (center1[2] - sim->pos_old[Z_COORD(*p3)]) * n_c[2];

                            vec3 center2;
                            center2[0] = sim->pos_old[X_COORD(*p3)] + k * n_c[0];
                            center2[1] = sim->pos_old[Y_COORD(*p3)] + k * n_c[1];
                            center2[2] = sim->pos_old[Z_COORD(*p3)] + k * n_c[2];

                            double a2 = eq_dist*eq_dist - (center2[0] - sim->pos_old[X_COORD(*p3)])*(center2[0] - sim->pos_old[X_COORD(*p3)])
                                                        - (center2[1] - sim->pos_old[Y_COORD(*p3)])*(center2[1] - sim->pos_old[Y_COORD(*p3)])
                                                        - (center2[2] - sim->pos_old[Z_COORD(*p3)])*(center2[2] - sim->pos_old[Z_COORD(*p3)]);

                            if(a2 >= 0)
                            {
                                a2 = sqrt(a2);
                                k = a1*a1-a2*a2;

                                vec3 delta_center;
                                delta_center[0] = center2[0] - center1[0];
                                delta_center[1] = center2[1] - center1[1];
                                delta_center[2] = center2[2] - center1[2];

                                double b_squared = delta_center[0]*delta_center[0] + delta_center[1]*delta_center[1] + delta_center[2]*delta_center[2];
                                double c = 0.5 + 0.5 * k / b_squared;
                                double a3 = a1*a1 - 0.25*b_squared - 0.5 * k - 0.25 * k*k / b_squared;

                                if(a3 >= 0)
                                {
                                    a3 = sqrt(a3);

                                    vec3 n;
                                    n[0] = n_c[1] * delta_center[2] - n_c[2] * delta_center[1];
                                    n[1] = n_c[2] * delta_center[0] - n_c[0] * delta_center[2];
                                    n[2] = n_c[0] * delta_center[1] - n_c[1] * delta_center[0];
                                    normalize(&n);

                                    pos[0] = center1[0] + c * delta_center[0] + a3 * n[0];
                                    pos[1] = center1[1] + c * delta_center[1] + a3 * n[1];
                                    pos[2] = center1[2] + c * delta_center[2] + a3 * n[2];

                                    if(sim->grid.canAddParticleAt(pos, sim->pos_old) && pos[0] > x_min && pos[0] < x_max && pos[2] > y_min && pos[2] < y_max)
                                        goto migration_successful;

                                    pos[0] = center1[0] + c * delta_center[0] - a3 * n[0];
                                    pos[1] = center1[1] + c * delta_center[1] - a3 * n[1];
                                    pos[2] = center1[2] + c * delta_center[2] - a3 * n[2];

                                    if(sim->grid.canAddParticleAt(pos, sim->pos_old) && pos[0] > x_min && pos[0] < x_max && pos[2] > y_min && pos[2] < y_max)
                                        goto migration_successful;
                                }
                            }
                        }
                    }
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            // only one migration
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            for(std::list<int>::iterator p2 = neighbours.begin(); p2 != neighbours.end(); ++p2)
            {
                // determine contact normal/distance between p1 and p2
                vec3 n_c;
                n_c[0] = sim->pos_old[X_COORD(*p2)] - sim->pos_old[X_COORD(p1)];
                n_c[1] = sim->pos_old[Y_COORD(*p2)] - sim->pos_old[Y_COORD(p1)];
                n_c[2] = sim->pos_old[Z_COORD(*p2)] - sim->pos_old[Z_COORD(p1)];

                double d = sqrt(n_c[0]*n_c[0] + n_c[1]*n_c[1] +n_c[2]*n_c[2]);

                n_c[0] /=d;
                n_c[1] /=d;
                n_c[2] /=d;

                vec3 center;
                center[0] = sim->pos_old[X_COORD(p1)] + 0.5 * d * n_c[0];
                center[1] = sim->pos_old[Y_COORD(p1)] + 0.5 * d * n_c[1];
                center[2] = sim->pos_old[Z_COORD(p1)] + 0.5 * d * n_c[2];

                vec3 k1, k2;
                SimLib::getOrthoVector(sim, &k1, n_c);
                k2[0] = n_c[1] * k1[2] - n_c[2] * k1[1];
                k2[1] = n_c[2] * k1[0] - n_c[0] * k1[2];
                k2[2] = n_c[0] * k1[1] - n_c[1] * k1[0];

                double a = sqrt(eq_dist*eq_dist - 0.25 * d*d);
                double phi_0 = sim->get_random_zero_twoPi();

                for(int i = 0; i < 64; ++i)
                {
                    double phi = phi_0 + 2.0 * M_PI * (double)i / 64.0;
                    double a1 = a * cos(phi);
                    double a2 = a * sin(phi);

                    pos[0] = center[0] + a1 * k1[0] + a2 * k2[0];
                    pos[1] = center[1] + a1 * k1[1] + a2 * k2[1];
                    pos[2] = center[2] + a1 * k1[2] + a2 * k2[2];

                    if(sim->grid.canAddParticleAt(pos, sim->pos_old) && pos[0] > x_min && pos[0] < x_max && pos[2] > y_min && pos[2] < y_max )
                        goto migration_successful;
                }
            }
        }
        }

        // no migraton occurred
        pos[0] = init_pos[0];
        pos[2] = init_pos[2];

        if(p1 >= 0)
        {
            double k = sim->pos_old[Y_COORD(p1)]*sim->pos_old[Y_COORD(p1)] + (sim->pos_old[X_COORD(p1)]-init_pos[0])*(sim->pos_old[X_COORD(p1)]-init_pos[0]) + (sim->pos_old[Z_COORD(p1)]-init_pos[2])*(sim->pos_old[Z_COORD(p1)]-init_pos[2]) - eq_dist*eq_dist;
            pos[1] = sim->pos_old[Y_COORD(p1)] + sqrt(sim->pos_old[Y_COORD(p1)]*sim->pos_old[Y_COORD(p1)] - k);
        }
        else
            pos[1] = 0;

migration_successful:

        sim->pos_old[X_COORD(p)] = pos[0];
        sim->pos_old[Y_COORD(p)] = pos[1];
        sim->pos_old[Z_COORD(p)] = pos[2];
        sim->grid.addParticle(pos, p);

        if(pos[1] > current_height)
            current_height = pos[1];
    }


    printf("100 percent done.\n");

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // remove misplaced particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    int removed_particles = 0;
    bool *remove_list = new bool[sim->number_of_particles];

    double max_deviation = 0;
    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        remove_list[p] = false;

        ContactListEntry *next_cl_entry = sim->contact_list[p];

        while(next_cl_entry)
        {
            double dist_squared = (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)]) * (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)]) * (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]) * (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]);


            int deviation = fabs( (dist_squared - eq_dist_squared)  / eq_dist_squared ) > 0.01;

            if(deviation > max_deviation)
                max_deviation = deviation;

            if(deviation > 0.01 )
            {
                ++removed_particles;
                remove_list[p] = true;

                next_cl_entry = NULL;
            }
            else
                next_cl_entry = next_cl_entry->next;
        }
    }

    sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();
    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

ErrorCode SimLib::initFractalAggregate(Simulation *sim, unsigned int number_of_particles, double migration_probability1, double migration_probability2, double chain_ratio, unsigned int min_chain_length, unsigned int max_chain_length)
{
    if(number_of_particles < 1)
        return EC_INVALID_PARAMETER;

    if(migration_probability1 < 0 || migration_probability2 < 0 || chain_ratio < 0)
        return EC_INVALID_PARAMETER;

    double eq_dist = 2.0 * particle_radius - delta_0;

    sim->resizeArrays(number_of_particles, 0);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // build aggregate
    /////////////////////////////////////////////////////////////////////////////////////////////

    // place first particle
    vec3 pos = {0,0,0};
    sim->pos_old[X_COORD(0)] = pos[0];
    sim->pos_old[Y_COORD(0)] = pos[1];
    sim->pos_old[Z_COORD(0)] = pos[2];
    sim->grid.addParticle(pos, 0);

    double outer_radius = particle_radius;
    unsigned int current_particles = 1;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // add chains first
    //////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int total_chain_paticles = (unsigned int)((double)number_of_particles * chain_ratio);

    while(current_particles < total_chain_paticles)
    {
        std::uniform_int_distribution<int> uniform_dist(0, current_particles-1);
        unsigned int particle_id = uniform_dist(sim->rand_generator);
        unsigned int chain_particles = min_chain_length;

        if (min_chain_length < max_chain_length)
        {
            std::uniform_int_distribution<int> uniform_dist2(0, max_chain_length - min_chain_length - 1);
            chain_particles += uniform_dist2(sim->rand_generator);
        }

        if(chain_particles >= number_of_particles - current_particles)
            chain_particles = number_of_particles - current_particles;

        // determine direction of chain
        /*pos[0] = sim->pos_old[X_COORD(particle_id)];
        pos[1] = sim->pos_old[Y_COORD(particle_id)];
        pos[2] = sim->pos_old[Z_COORD(particle_id)];

        // prevent error in atan2
        if(pos[0] == 0 && pos[1] == 0)
            pos[0] = 2.0 * particle_radius;

        normalize(&pos);

        double phi0 = atan2(pos[0], pos[1]) + 1.0 * M_PI * (-0.5 + (double)rand() / ((double)RAND_MAX));
        double theta0 = acos(pos[2]) + 0.5 * M_PI * (-0.5 + (double)rand() / ((double)RAND_MAX));*/

        // determine random direction
        double phi0 = sim->get_random_zero_twoPi();
        double theta0 = sim->get_random_cos_theta();

        for(unsigned int p = 0; p < chain_particles; ++p)
        {
            for(int i = 0; i < 32; ++i)
            {
                double phi = phi0 + M_PI * (-0.5 + sim->get_random_zero_one_incl());

                double theta = theta0 + ( 0.5 * sim->get_random_zero_one_incl() - 0.25);
                if(theta < -1.0)
                    theta = -1.0;
                else if(theta  > 1.0)
                    theta = 1.0;
                theta = asin(theta) + M_PI_2;

                vec3 dir;
                dir[0] = sin(theta) * cos(phi);
                dir[1] = sin(theta) * sin(phi);
                dir[2] = cos(theta);

                pos[0] = sim->pos_old[X_COORD(particle_id)] + equilibrium_distance * dir[0];
                pos[1] = sim->pos_old[Y_COORD(particle_id)] + equilibrium_distance * dir[1];
                pos[2] = sim->pos_old[Z_COORD(particle_id)] + equilibrium_distance * dir[2];

                if(sim->grid.canAddParticleAt(pos, sim->pos_old))
                {
                    sim->pos_old[X_COORD(current_particles)] = pos[0];
                    sim->pos_old[Y_COORD(current_particles)] = pos[1];
                    sim->pos_old[Z_COORD(current_particles)] = pos[2];
                    sim->grid.addParticle(pos, current_particles);
                    particle_id = current_particles;
                    ++current_particles;
                    break;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // add single particles at random locations
    //////////////////////////////////////////////////////////////////////////////////////////////

    while(current_particles < number_of_particles)
    {
        std::uniform_int_distribution<int> uniform_dist(0, current_particles - 1);
        unsigned int particle_id = uniform_dist(sim->rand_generator);

        if(SimLib::getNeighbourPos(&pos, *sim, particle_id))
        {
            bool migration_successful = false;
            vec3 migration_pos;

            if(migration_probability1 >= sim->get_random_zero_one_incl())
            {
                if(migration_probability2 >= sim->get_random_zero_one_incl() && current_particles > 3)
                    migration_successful = SimLib::getMigrationPos2(&migration_pos, *sim, particle_id, pos, BAM_SELECT_RANDOM);
                else
                    migration_successful = SimLib::getMigrationPos1(&migration_pos, *sim, particle_id, pos, BAM_SELECT_RANDOM);
            }

            if(migration_successful)
            {
                sim->pos_old[X_COORD(current_particles)] = migration_pos[0];
                sim->pos_old[Y_COORD(current_particles)] = migration_pos[1];
                sim->pos_old[Z_COORD(current_particles)] = migration_pos[2];
                sim->grid.addParticle(migration_pos, current_particles);
            }
            else
            {
                sim->pos_old[X_COORD(current_particles)] = pos[0];
                sim->pos_old[Y_COORD(current_particles)] = pos[1];
                sim->pos_old[Z_COORD(current_particles)] = pos[2];
                sim->grid.addParticle(pos, current_particles);
            }

            ++current_particles;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();
    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;
    SimLib::centerCMS(sim);

    return EC_OK;
}

ErrorCode SimLib::initTreeAggregate(Simulation *sim, unsigned int number_of_particles, double contact_distribution[6], double migration_rate1, double migration_rate2)
{
    if(number_of_particles < 1)
        return EC_INVALID_PARAMETER;

    // normalize contact distribution
    double sum = contact_distribution[0] + contact_distribution[1] + contact_distribution[2] + contact_distribution[3] + contact_distribution[4] + contact_distribution[5];

    contact_distribution[0] /= sum;

    for(int c = 1; c < 6; ++c)
    {
        contact_distribution[c] /= sum;
        contact_distribution[c] += contact_distribution[c-1];
    }

    double eq_dist = 2.0 * particle_radius - delta_0;

    sim->resizeArrays(number_of_particles, 0);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // add particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    // place first particle
    vec3 pos = {0,0,0};
    sim->pos_old[X_COORD(0)] = pos[0];
    sim->pos_old[Y_COORD(0)] = pos[1];
    sim->pos_old[Z_COORD(0)] = pos[2];
    sim->grid.addParticle(pos, 0);

    double outer_radius = particle_radius;
    unsigned int current_particles = 1;

    int first_seed_id = 0;
    int last_seed_id = 0;

    while(current_particles < number_of_particles)
    {
        for(int p = first_seed_id; p <= last_seed_id; ++p)
        {
            // determine number of contacts
            unsigned int coordination_number = 5;
            double value = sim->get_random_zero_one();

            int c = 0;
            if(current_particles < 25 && contact_distribution[1] > 0)
                c = 1;

            for(; c < 5; ++c)
            {
                if(value < contact_distribution[c])
                {
                    coordination_number = c;
                    break;
                }
            }

            // ensure that at least one particle is added
            if(coordination_number == 0 && p == first_seed_id)
                coordination_number = 1;

            for(c = 0; c < coordination_number; ++c)
            {
                bool valid_pos = false;

                if(sim->get_random_zero_one() < migration_rate1 && c < coordination_number-1)
                {
                    valid_pos = SimLib::getMigrationPos1(&pos, *sim, p, pos, BAM_SELECT_RANDOM);

                    if(valid_pos)
                        c += 1;
                }
                else if(sim->get_random_zero_one() < migration_rate2 && c < coordination_number-2)
                {
                    valid_pos = SimLib::getMigrationPos2(&pos, *sim, p, pos, BAM_SELECT_RANDOM);

                    if(valid_pos)
                        c += 2;
                }

                if(!valid_pos)
                {
                    for(int i = 0; i < 20; ++i)
                    {
                        valid_pos = SimLib::getNeighbourPos(&pos, *sim, p);

                        if(valid_pos)
                            break;
                    }
                }

                if(valid_pos)
                {
                    sim->pos_old[X_COORD(current_particles)] = pos[0];
                    sim->pos_old[Y_COORD(current_particles)] = pos[1];
                    sim->pos_old[Z_COORD(current_particles)] = pos[2];
                    sim->grid.addParticle(pos, current_particles);
                    ++current_particles;

                    // stop aggregation if maximum number of particles has been reached
                    if(current_particles == number_of_particles)
                        goto tree_aggregation_finished;
                }
            }
        }

        // prevent getting stuck in an infinite loop if no particles have been added
        if(current_particles == last_seed_id+1)
            goto tree_aggregation_failed;

        // neighbours to all particles have been added -> next iteration
        first_seed_id = last_seed_id+1;
        last_seed_id = current_particles-1;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

tree_aggregation_finished:

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();
    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    SimLib::centerCMS(sim);
    return EC_OK;

tree_aggregation_failed:
    sim->cleanUp();
    return EC_INTERNAL_ERROR;
}

ErrorCode SimLib::addFractalChainsToAggregate(Simulation *sim, const char *filename, unsigned int number_of_chains, unsigned int min_chain_length, unsigned int max_chain_length, double migration_probability)
{
    if(migration_probability < 0)
        return EC_INVALID_PARAMETER;

    if(filename)
    {
        sim->loadFromFile(filename);
        sim->removeWalls();
    }
    else
    {
        if(sim->number_of_particles == 0)
            return EC_NO_PARTICLES;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // place chains
    /////////////////////////////////////////////////////////////////////////////////////////////

    Simulation sim2;
    sim2.setMaterialConstants(particle_radius, density, surface_energy, nu, young_mod, crit_rolling_displacement, osc_damping_factor, rolling_modifier, sliding_modifier, twisting_modifier, crit_sliding_displacement_modifier, crit_wall_sliding_displacement_modifier);

    std::uniform_int_distribution<int> uniform_dist(min_chain_length, max_chain_length-1);
    for(unsigned int chain = 0; chain < number_of_chains; ++chain)
    {
        // generate chain
        unsigned int chain_length = (unsigned int)uniform_dist(sim->rand_generator);
        SimLib::initNeedleDeposition(&sim2, chain_length);
        vec3 axis = {0,0, 1.0};
        SimLib::rotateParticles(&sim2, axis, - M_PI_2, 0, sim2.number_of_particles-1);

        ErrorCode error_code = SimLib::hitAndStick(sim, &sim2, 0, true, false);

        if(error_code != EC_OK)
            return error_code;
    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();
    sim->initSimData(SIM_TYPE_GENERAL);

    return EC_OK;
}

ErrorCode SimLib::buildSimpleCubicGrid(Simulation *sim, double x_size, double y_size, double z_size, double filling_factor)
{
    int particles = 0;
    int estimated_particles = 0;

    // calculate estimated number of particles
    estimated_particles = (int) (x_size*y_size*z_size  / (8.0 * particle_radius*particle_radius*particle_radius) );

    double *positions = new double[3 * (int)(1.25 * (double)estimated_particles)];
    double dist = 2.0 * particle_radius - delta_0;
    double x_pos = 0.5 * dist;
    double y_pos = 0.5 * dist;
    double z_pos = 0.5 * dist;
    bool stop = false;

    while(!stop)
    {
        if( filling_factor > sim->get_random_zero_one_incl())
        {
            positions[X_COORD(particles)] = x_pos;
            positions[Y_COORD(particles)] = z_pos;
            positions[Z_COORD(particles)] = y_pos;

            ++particles;
        }

        x_pos += dist;

        if(x_pos > x_size - 0.5 * dist)
        {
            x_pos = 0.5 * dist;
            y_pos += dist;
        }

        if(y_pos > y_size - 0.5 * dist)
        {
            y_pos = 0.5 * dist;
            z_pos += dist;
        }

        if(z_pos > z_size - 0.5 * dist)
            stop = true;
    }

    // copy to sim & clean up
    sim->resizeArrays(particles, 0);
    memcpy(sim->pos_old, positions, 3 * particles * sizeof(double) );
    SimLib::centerCMS(sim);
    delete [] positions;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

ErrorCode SimLib::buildHexagonalGrid(Simulation *sim, int x_particles, int y_particles, int z_particles, double filling_factor)
{
    int particles = 0;

    // temporary array (entire memory will only be needed if all grid positions will be filled with a particle)
    double *positions = new double[3 * x_particles * y_particles * z_particles];

    double dist = 2.0 * particle_radius - delta_0;

    // base vectors describing the grid
    vec3 x_base_vec, y_base_vec, z_base_vec;

    x_base_vec[0] = dist;
    x_base_vec[1] = 0;
    x_base_vec[2] = 0;

    y_base_vec[0] = 0.5 * dist;
    y_base_vec[1] = 0.5 * sqrt(3.0) * dist;
    y_base_vec[2] = 0;

    z_base_vec[0] = 0.5 * dist;
    z_base_vec[1] = sqrt(3.0) / 6.0 * dist;
    z_base_vec[2] = sqrt(6.0) / 3.0 * dist;

    for(int x = 0; x < x_particles; ++x)
    {
        for(int y = 0; y < y_particles; ++y)
        {
            for(int z = 0; z < z_particles; ++z)
            {
                if( filling_factor > sim->get_random_zero_one_incl())
                {
                    positions[X_COORD(particles)] = (double)(x-y/2-z/3) * x_base_vec[0] + (double)(y-z/3) * y_base_vec[0] + (double)z * z_base_vec[0];
                    positions[Y_COORD(particles)] = (double)(x-y/2-z/3) * x_base_vec[2] + (double)(y-z/3) * y_base_vec[2] + (double)z * z_base_vec[2];
                    positions[Z_COORD(particles)] = (double)(x-y/2-z/3) * x_base_vec[1] + (double)(y-z/3) * y_base_vec[1] + (double)z * z_base_vec[1];

                    ++particles;
                }
            }
        }
    }

    // copy to sim & clean up
    sim->resizeArrays(particles, 0);
    memcpy(sim->pos_old, positions, 3 * particles * sizeof(double) );
    SimLib::centerCMS(sim);
    delete [] positions;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

double SimLib::getLowerXPos(Simulation *sim, double pos_y, double pos_z, double x_max, double x_min)
{
    // determine start cell
    int center_cell_id = sim->grid.getCellID(x_max, pos_y, pos_z);
    int max_x_cell_id = sim->grid.getXID(x_max);
    int min_x_cell_id = sim->grid.getXID(x_min);
    double current_x_pos = x_min;
    int closest_id = -1;

    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);
    bool stop_search = false;

    // check cells around the current cell for particles
    // if none are found, decrease x and continue search in next level
    for(int x = max_x_cell_id; x >= min_x_cell_id; --x)
    {
        for(int y = -1; y <= 1; ++y)
        {
            for(int z = -1; z <= 1; ++z)
            {
                // check particles in this cell
                GridEntry *grid_entry = sim->grid.particles_in_cell[center_cell_id + (x-max_x_cell_id) + y * sim->grid.x_cells + z * sim->grid.x_cells * sim->grid.y_cells];

                while(grid_entry)
                {
                    int id = grid_entry->particle;

                    // calc dist to "path" of particle moving along x-axis on (x, pos_y, pos_z)
                    double dist_squared = (pos_y - sim->pos_old[Y_COORD(id)])*(pos_y - sim->pos_old[Y_COORD(id)]) + (pos_z - sim->pos_old[Z_COORD(id)])*(pos_z - sim->pos_old[Z_COORD(id)]);

                    // only consider particles that are close enough to establish a contact
                    if(dist_squared < eq_dist_squared)
                    {
                        // calc pos of new particle
                        double k = sim->pos_old[X_COORD(id)]*sim->pos_old[X_COORD(id)] + (sim->pos_old[Y_COORD(id)]-pos_y)*(sim->pos_old[Y_COORD(id)]-pos_y) + (sim->pos_old[Z_COORD(id)]-pos_z)*(sim->pos_old[Z_COORD(id)]-pos_z) - eq_dist_squared;

                        double x_pos = sim->pos_old[X_COORD(id)] + sqrt(sim->pos_old[X_COORD(id)]*sim->pos_old[X_COORD(id)] - k);

                        if(x_pos > current_x_pos)
                        {
                            current_x_pos = x_pos;
                            closest_id = id;
                        }
                    }

                    grid_entry = grid_entry->next;
                }
            }
        }

        if(stop_search)
            return current_x_pos;

        // if particle has been found, end search after checking the subsequent slice
        // (particles in slices with lower x id will be too far away)
        if(closest_id >= 0)
            stop_search = true;
    }

    if(closest_id >= 0)
        return current_x_pos;
    else
        return x_min;
}

double SimLib::getYPos(Simulation *sim, double pos_x, double pos_z, double y_max, double y_min)
{
    // determine start cell
    int center_cell_id = sim->grid.getCellID(pos_x, y_max, pos_z);
    int max_y_cell_id = sim->grid.getYID(y_max);
    int min_y_cell_id = sim->grid.getYID(y_min);
    double current_y_pos = y_min;
    int closest_id = -1;

    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);
    bool stop_search = false;

    // check cells around the current cell for particles
    // if none are found, decrease y and continue search in next level
    for(int y = max_y_cell_id; y >= min_y_cell_id; --y)
    {
        for(int x = -1; x <= 1; ++x)
        {
            for(int z = -1; z <= 1; ++z)
            {
                // check particles in this cell
                GridEntry *grid_entry = sim->grid.particles_in_cell[center_cell_id + x + (y-max_y_cell_id) * sim->grid.x_cells + z * sim->grid.x_cells * sim->grid.y_cells];

                while(grid_entry)
                {
                    int id = grid_entry->particle;

                    // calc dist to "path" of falling particle
                    double dist_squared = (pos_x - sim->pos_old[X_COORD(id)])*(pos_x - sim->pos_old[X_COORD(id)]) + (pos_z - sim->pos_old[Z_COORD(id)])*(pos_z - sim->pos_old[Z_COORD(id)]);

                    // only consider particles that are close enough to establish a contact
                    if(dist_squared < eq_dist_squared)
                    {
                        // calc pos of new particle
                        double y_pos = sim->pos_old[Y_COORD(id)] + sqrt(eq_dist_squared - dist_squared);

                        if(y_pos > current_y_pos)
                        {
                            current_y_pos = y_pos;
                            closest_id = id;
                        }
                    }

                    grid_entry = grid_entry->next;
                }
            }
        }

        if(stop_search)
            return current_y_pos;

        // if particle has been found, end search after checking the subsequent slice
        // (particles in slices with lower y id will be too far away)
        if(closest_id >= 0)
            stop_search = true;
    }

    if(closest_id >= 0)
        return current_y_pos;
    else
        return y_min;
}



ErrorCode SimLib::initRandomBallisticDeposition(Simulation *sim, int number_of_particles, double x_size, double y_size, double dest_filling_factor)
{
    if(number_of_particles <= 0)
        return EC_NO_PARTICLES;

    // init sim
    sim->resizeArrays(number_of_particles, 0);


    double x_min = particle_radius;
    double x_max = x_size - particle_radius;
    double y_min = particle_radius;
    double y_max = y_size - particle_radius;

    // place first particle
    vec3 pos;
    pos[0] = x_min + (x_max - x_min) * sim->get_random_zero_one_incl(); //(x_size - 2.0*particle_radius) * sim->get_random_zero_one_incl();
    pos[2] = y_min + (y_max - y_min) * sim->get_random_zero_one_incl(); //(y_size - 2.0*particle_radius) * sim->get_random_zero_one_incl();
    pos[1] = particle_radius;

    sim->pos_old[X_COORD(0)] = pos[0];
    sim->pos_old[Y_COORD(0)] = pos[1];
    sim->pos_old[Z_COORD(0)] = pos[2];

    sim->grid.addParticle(pos, 0);

    double current_height = 2.0 * particle_radius;

    // place other particles
    for(int p = 1; p < number_of_particles; ++p)
    {
        // drop particle from current height
        pos[0] = x_min + (x_max - x_min) * sim->get_random_zero_one_incl(); //x_size * sim->get_random_zero_one_incl();
        pos[2] = y_min + (y_max - y_min) * sim->get_random_zero_one_incl(); //y_size * sim->get_random_zero_one_incl();
        pos[1] = current_height + 1.1 * particle_radius;

        // determine y_pos of next particle
        pos[1] = SimLib::getYPos(sim, pos[0], pos[2], pos[1], particle_radius);

        // add particle
        sim->pos_old[X_COORD(p)] = pos[0];
        sim->pos_old[Y_COORD(p)] = pos[1];
        sim->pos_old[Z_COORD(p)] = pos[2];
        sim->grid.addParticle(pos, p);

        if(pos[1] + particle_radius > current_height)
            current_height = pos[1] + particle_radius;
    }

    // artificially reduce filling factor
    if(dest_filling_factor > 0)
    {

        vec3 lower_pos, upper_pos;
        sim->getEnclosingBox(&lower_pos, &upper_pos);

        double volume = (upper_pos[0] - lower_pos[0]) * (upper_pos[1] - lower_pos[1]) * (upper_pos[2] - lower_pos[2]);

        //adjust size of box in which the filling factor is determined
        lower_pos[0] += 3.0 * particle_radius;
        lower_pos[1] += 3.0 * particle_radius;
        lower_pos[2] += 3.0 * particle_radius;
        upper_pos[0] -= 3.0 * particle_radius;
        upper_pos[1] -= 3.0 * particle_radius;
        upper_pos[2] -= 3.0 * particle_radius;

        double current_filling_factor = SimLib::getFillingFactorOfBox(*sim, lower_pos, upper_pos);

        if(current_filling_factor > dest_filling_factor)
        {
            int removed_particles = (current_filling_factor - dest_filling_factor) * 3.0/(4.0*M_PI) * volume / (particle_radius*particle_radius*particle_radius);

            bool *remove_list = new bool[sim->number_of_particles];
            for(int p = 0; p < sim->number_of_particles; ++p)
                remove_list[p] = false;

            std::set<int> selected_particles;

            // select particles randomly
            std::uniform_int_distribution<int> distribution(0, number_of_particles-1);
            for(int p = 0; p < removed_particles; ++p)
            {
                // select a particle that is not already on the remove list
                int id = distribution(sim->rand_generator);
                while(selected_particles.find(id) != selected_particles.end())
                    id = distribution(sim->rand_generator);

                remove_list[id] = true;
                selected_particles.insert(id);
            }

            sim->removeParticles(remove_list, removed_particles);

            delete [] remove_list;
        }

    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    /////////////////////////////////////////////////////////////////////////////////////////////
    // check for misplaced particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    int removed_particles = 0;
    bool *remove_list = new bool[sim->number_of_particles];

    double max_deviation = 0;
    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        remove_list[p] = false;

        ContactListEntry *next_cl_entry = sim->contact_list[p];

        while(next_cl_entry)
        {
            double dist_squared = (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)]) * (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)]) * (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]) * (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]);


            int deviation = fabs( (dist_squared - eq_dist_squared)  / eq_dist_squared ) > 0.01;

            if(deviation > max_deviation)
                max_deviation = deviation;

            if(deviation > 0.01 )
            {
                ++removed_particles;
                remove_list[p] = true;

                next_cl_entry = NULL;
            }
            else
                next_cl_entry = next_cl_entry->next;
        }
    }

    sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}



ErrorCode SimLib::reduceBoxFillingFactor(Simulation *sim, double dest_filling_factor, bool side_walls)
{
    if(sim->number_of_particles <= 0)
        return EC_NO_PARTICLES;

    const int number_of_particles = sim->number_of_particles;

    // artificially reduce filling factor
    if(dest_filling_factor <= 0.0)
        return EC_OK;

    {
        double N_remove;
        if(side_walls)
        {
            vec3 lower_pos, upper_pos;
            sim->getEnclosingBox(&lower_pos, &upper_pos);

            //adjust size of box in which the filling factor is determined
            lower_pos[0] += 3.0 * particle_radius;
            lower_pos[1] += 3.0 * particle_radius;
            lower_pos[2] += 3.0 * particle_radius;
            upper_pos[0] -= 3.0 * particle_radius;
            upper_pos[1] -= 3.0 * particle_radius;
            upper_pos[2] -= 3.0 * particle_radius;

            double volume_i = (upper_pos[0] - lower_pos[0]) * (upper_pos[1] - lower_pos[1]) * (upper_pos[2] - lower_pos[2]);

            double filling_factor_i = SimLib::getFillingFactorOfBox(*sim, lower_pos, upper_pos);

            double volume_p =  (4.0*M_PI) / 3.0 * (particle_radius*particle_radius*particle_radius);

            double N_i = filling_factor_i * volume_i / volume_p;


            N_remove = (filling_factor_i - dest_filling_factor) * volume_i / volume_p;

            N_remove = N_remove / N_i * sim->number_of_particles;
        }
        else
        {

            vec3 lower_pos, upper_pos;
            sim->getEnclosingBox(&lower_pos, &upper_pos);

            double volume = (upper_pos[0] - lower_pos[0]) * (upper_pos[1] - lower_pos[1]) * (upper_pos[2] - lower_pos[2]);


            double current_filling_factor = (double)sim->number_of_particles * 4.0/3.0 * M_PI * particle_radius*particle_radius*particle_radius / volume;

            N_remove = (current_filling_factor - dest_filling_factor) * 3.0/(4.0*M_PI) * volume / (particle_radius*particle_radius*particle_radius);

        }


        if(N_remove >= 1.0)
        {
            int removed_particles = (int)(N_remove+0.5);



            bool *remove_list = new bool[sim->number_of_particles];
            for(int p = 0; p < sim->number_of_particles; ++p)
                remove_list[p] = false;

            std::set<int> selected_particles;

            // select particles randomly
            std::uniform_int_distribution<int> distribution(0, number_of_particles-1);
            for(int p = 0; p < removed_particles; ++p)
            {
                // select a particle that is not already on the remove list
                int id = distribution(sim->rand_generator);
                while(selected_particles.find(id) != selected_particles.end())
                {
                    id = distribution(sim->rand_generator);
                }
                remove_list[id] = true;
                selected_particles.insert(id);
            }

            sim->removeParticles(remove_list, removed_particles);

            delete [] remove_list;
        }

    }


    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    /////////////////////////////////////////////////////////////////////////////////////////////
    // check for misplaced particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    int removed_particles = 0;
    bool *remove_list = new bool[sim->number_of_particles];

    double max_deviation = 0;
    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        remove_list[p] = false;

        ContactListEntry *next_cl_entry = sim->contact_list[p];

        while(next_cl_entry)
        {
            double dist_squared = (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)]) * (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)]) * (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]) * (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]);


            int deviation = fabs( (dist_squared - eq_dist_squared)  / eq_dist_squared ) > 0.01;

            if(deviation > max_deviation)
                max_deviation = deviation;

            if(deviation > 0.01 )
            {
                ++removed_particles;
                remove_list[p] = true;

                next_cl_entry = NULL;
            }
            else
                next_cl_entry = next_cl_entry->next;
        }
    }

    sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}


/*
ErrorCode SimLib::initRandomBallisticDeposition_old(Simulation *sim, int number_of_particles, double x_size, double y_size, double dest_filling_factor)
{
    if(number_of_particles <= 0)
        return EC_NO_PARTICLES;

    // init sim
    sim->resizeArrays(number_of_particles, 0);

    // place first particle
    vec3 pos;
    pos[0] = (x_size - 2.0*particle_radius) * sim->get_random_zero_one_incl();
    pos[2] = (y_size - 2.0*particle_radius) * sim->get_random_zero_one_incl();
    pos[1] = particle_radius;

    sim->pos_old[X_COORD(0)] = pos[0];
    sim->pos_old[Y_COORD(0)] = pos[1];
    sim->pos_old[Z_COORD(0)] = pos[2];

    sim->grid.addParticle(pos, 0);

    double current_height = 2.0 * particle_radius;

    // place other particles
    for(int p = 1; p < number_of_particles; ++p)
    {
        // drop particle from current height
        pos[0] = x_size * sim->get_random_zero_one_incl();
        pos[2] = y_size * sim->get_random_zero_one_incl();
        pos[1] = current_height + 1.1 * particle_radius;

        // determine y_pos of next particle
        pos[1] = SimLib::getYPos(sim, pos[0], pos[2], pos[1], particle_radius);

        // add particle
        sim->pos_old[X_COORD(p)] = pos[0];
        sim->pos_old[Y_COORD(p)] = pos[1];
        sim->pos_old[Z_COORD(p)] = pos[2];
        sim->grid.addParticle(pos, p);

        if(pos[1] + particle_radius > current_height)
            current_height = pos[1] + particle_radius;
    }

    // artificially reduce filling factor
    if(dest_filling_factor > 0)
    {
        vec3 lower_pos, upper_pos;
        sim->getEnclosingBox(&lower_pos, &upper_pos);

        double volume = (upper_pos[0] - lower_pos[0]) * (upper_pos[1] - lower_pos[1]) * (upper_pos[2] - lower_pos[2]);
        double current_filling_factor = (double)sim->number_of_particles * 4.0/3.0 * M_PI * particle_radius*particle_radius*particle_radius / volume;

        if(current_filling_factor > dest_filling_factor)
        {
            int removed_particles = (current_filling_factor - dest_filling_factor) * 3.0/(4.0*M_PI) * volume / (particle_radius*particle_radius*particle_radius);

            bool *remove_list = new bool[sim->number_of_particles];
            for(int p = 0; p < sim->number_of_particles; ++p)
                remove_list[p] = false;

            std::set<int> selected_particles;

            // select particles randomly
            std::uniform_int_distribution<int> distribution(0, number_of_particles-1);
            for(int p = 0; p < removed_particles; ++p)
            {
                // select a particle that is not already on the remove list
                int id = distribution(sim->rand_generator);
                while(selected_particles.find(id) != selected_particles.end())
                    id = distribution(sim->rand_generator);

                remove_list[id] = true;
                selected_particles.insert(id);
            }

            sim->removeParticles(remove_list, removed_particles);

            delete [] remove_list;
        }
    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    /////////////////////////////////////////////////////////////////////////////////////////////
    // check for misplaced particles
    /////////////////////////////////////////////////////////////////////////////////////////////

    int removed_particles = 0;
    bool *remove_list = new bool[sim->number_of_particles];

    double max_deviation = 0;
    double eq_dist_squared = (2.0 * particle_radius - delta_0) * (2.0 * particle_radius - delta_0);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        remove_list[p] = false;

        ContactListEntry *next_cl_entry = sim->contact_list[p];

        while(next_cl_entry)
        {
            double dist_squared = (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)]) * (sim->pos_new[X_COORD(p)] - sim->pos_new[X_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)]) * (sim->pos_new[Y_COORD(p)] - sim->pos_new[Y_COORD(next_cl_entry->id)])
                         + (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]) * (sim->pos_new[Z_COORD(p)] - sim->pos_new[Z_COORD(next_cl_entry->id)]);


            int deviation = fabs( (dist_squared - eq_dist_squared)  / eq_dist_squared ) > 0.01;

            if(deviation > max_deviation)
                max_deviation = deviation;

            if(deviation > 0.01 )
            {
                ++removed_particles;
                remove_list[p] = true;

                next_cl_entry = NULL;
            }
            else
                next_cl_entry = next_cl_entry->next;
        }
    }

    sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}
*/

ErrorCode SimLib::initNeedleDeposition(Simulation *sim, int number_of_particles)
{
    if(number_of_particles <= 0)
        return EC_NO_PARTICLES;

    sim->resizeArrays(number_of_particles, 0);

    // place first particle
    vec3 pos;
    pos[0] = 0;
    pos[2] = 0;
    pos[1] = particle_radius;

    sim->pos_old[X_COORD(0)] = pos[0];
    sim->pos_old[Y_COORD(0)] = pos[1];
    sim->pos_old[Z_COORD(0)] = pos[2];

    sim->grid.addParticle(pos, 0);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // // drop particles on top of the needle
    /////////////////////////////////////////////////////////////////////////////////////////////

    double dist;
    double outer_radius = 2.0 * particle_radius;
    double current_height = 2.0 * particle_radius;

    for(int p = 1; p < number_of_particles; ++p)
    {
POUR_DOWN_NEW_PARTICLE:

        do
        {
            pos[0] = outer_radius * (1.0 - 2.0 * sim->get_random_zero_one_incl());
            pos[2] = outer_radius * (1.0 - 2.0 * sim->get_random_zero_one_incl());
            dist = pos[0]*pos[0] + pos[2]*pos[2];
        }while(dist > outer_radius*outer_radius);

        pos[1] = current_height + 2.1 * particle_radius;

        // determine y_pos of next particle
        pos[1] = SimLib::getYPos(sim, pos[0], pos[2], pos[1], particle_radius);

        if(pos[1] < 1.01 * particle_radius)
            goto POUR_DOWN_NEW_PARTICLE;

        // add particle
        sim->pos_old[X_COORD(p)] = pos[0];
        sim->pos_old[Y_COORD(p)] = pos[1];
        sim->pos_old[Z_COORD(p)] = pos[2];
        sim->grid.addParticle(pos, p);

        // check outer radius
        double dist = 2.0 * particle_radius + sqrt(pos[0]*pos[0] + pos[2]*pos[2]);

        if(dist > outer_radius)
            outer_radius = dist;

        if(pos[1] > current_height)
            current_height = pos[1];
    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;
    return EC_OK;
}

ErrorCode SimLib::initCylindricalRBD(Simulation *sim, int number_of_particles, double radius, double slice_factor)
{
    if(number_of_particles <= 0)
        return EC_NO_PARTICLES;

    if(slice_factor < 0 || radius < 0)
        return EC_INVALID_PARAMETER;

    // init sim
    sim->resizeArrays(number_of_particles, 0);

    vec3 pos;
    double phi;
    double r;

    // place first particle
    phi = sim->get_random_zero_twoPi();
    r = sqrt(sim->get_random_zero_one_incl()) * radius;

    pos[0] = r * sin(phi);
    pos[2] = r * cos(phi);
    pos[1] = particle_radius;

    sim->pos_old[0] = pos[0];
    sim->pos_old[1] = pos[1];
    sim->pos_old[2] = pos[2];

    sim->grid.addParticle(pos, 0);

    double current_height = 2.0 * particle_radius;

    double percentage;
    printf("0.00%% done");

    // place other particles
    for(int p = 1; p < number_of_particles; ++p)
    {
        // drop particle from current height
        phi = sim->get_random_zero_twoPi();
        r = sqrt(sim->get_random_zero_one_incl()) * radius;
        pos[0] = r * sin(phi);
        pos[2] = r * cos(phi);
        pos[1] = current_height + 1.1 * particle_radius;

        // determine y_pos of next particle
        pos[1] = SimLib::getYPos(sim, pos[0], pos[2], pos[1], particle_radius);

        // add particle
        sim->pos_old[X_COORD(p)] = pos[0];
        sim->pos_old[Y_COORD(p)] = pos[1];
        sim->pos_old[Z_COORD(p)] = pos[2];
        sim->grid.addParticle(pos, p);

        if(pos[1] + particle_radius > current_height)
            current_height = pos[1] + particle_radius;

        printf("\b\b\b\b\b\b\b\b\b\b\b");
        percentage = 100.0 * (double)p/(double)number_of_particles;
        printf("%3.2lf%% done", percentage);
    }

    if(slice_factor > 0)
    {
        // determine cake height
        vec3 lower_pos, upper_pos;
        sim->getEnclosingBox(&lower_pos, &upper_pos);

        upper_pos[1] -= (upper_pos[1] - lower_pos[1]) * slice_factor;
        SimLib::sliceBox(sim, lower_pos, upper_pos);
    }

    SimLib::centerCMS(sim);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

ErrorCode SimLib::initChain(Simulation *sim, int number_of_chain_particles, int number_of_impact_particles, int target_id, double angular_irregularity, double impact_speed)
{
    angular_irregularity *= 0.5 * M_PI;
    sim->resizeArrays(number_of_chain_particles + number_of_impact_particles, 0);

    double phi, theta;
    double chain_dist = 2.0 * particle_radius - delta_0;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up first chain
    /////////////////////////////////////////////////////////////////////////////////////////////

    sim->pos_old[X_COORD(0)] = 0;
    sim->pos_old[Y_COORD(0)] = 0.5 * number_of_chain_particles * chain_dist;
    sim->pos_old[Z_COORD(0)] = 0;

    for(int p = 1; p < number_of_chain_particles; ++p)
    {
        phi = angular_irregularity - 2.0 * sim->get_random_zero_one_incl() * angular_irregularity;
        theta = angular_irregularity - 2.0 * sim->get_random_zero_one_incl() * angular_irregularity;

        sim->pos_old[X_COORD(p)] = sim->pos_old[X_COORD(p-1)] - chain_dist * sin(phi) * cos(theta);
        sim->pos_old[Y_COORD(p)] = sim->pos_old[Y_COORD(p-1)] - chain_dist * cos(phi) * cos(theta);
        sim->pos_old[Z_COORD(p)] = sim->pos_old[Z_COORD(p-1)] + chain_dist * sin(theta);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up second chain
    /////////////////////////////////////////////////////////////////////////////////////////////

    if(number_of_impact_particles > 0)
    {
        // locate center of chain
        sim->pos_old[X_COORD(number_of_chain_particles)] = sim->pos_old[X_COORD(target_id)] + 1.5 * chain_dist;
        sim->pos_old[Y_COORD(number_of_chain_particles)] = sim->pos_old[Y_COORD(target_id)];
        sim->pos_old[Z_COORD(number_of_chain_particles)] = sim->pos_old[Z_COORD(target_id)];

        sim->vel[X_COORD(number_of_chain_particles)] = -impact_speed;

        for(int p = number_of_chain_particles+1; p < sim->number_of_particles; ++p)
        {
            phi = angular_irregularity - 2.0 * sim->get_random_zero_one_incl() * angular_irregularity;
            theta = angular_irregularity - 2.0 * sim->get_random_zero_one_incl() * angular_irregularity;

            sim->pos_old[X_COORD(p)] = sim->pos_old[X_COORD(p-1)] + chain_dist * cos(phi) * cos(theta);
            sim->pos_old[Y_COORD(p)] = sim->pos_old[Y_COORD(p-1)] - chain_dist * sin(phi) * cos(theta);
            sim->pos_old[Z_COORD(p)] = sim->pos_old[Z_COORD(p-1)] + chain_dist * sin(theta);

            sim->vel[X_COORD(p)] = -impact_speed;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_COLLISION, impact_speed, 0.0, 1.0, 0.0, 0.0);

    return EC_OK;
}

ErrorCode SimLib::initCluster(Simulation *sim, int number_of_cluster_particles, int number_of_impact_particles, double filling_factor, double angular_irregularity, double impact_speed)
{
    if(filling_factor <= 0)
        return EC_INVALID_PARAMETER;

    number_of_impact_particles = 0;

    // prepare sim
    sim->resizeArrays(number_of_cluster_particles + number_of_impact_particles, 0);

    double particle_dist = 2.0 * particle_radius - delta_0;

    int particle_counter = 0;
    int x_pos = 0;
    int y_pos = 0;

    int particles_per_side = (int)pow((double)number_of_cluster_particles / filling_factor, 1.0/3.0);

    vec3 temp;
    temp[0] = - 0.5 * (double)particles_per_side * particle_dist;
    temp[1] = - 0.5 * (double)particles_per_side * particle_dist;
    temp[2] = - 0.5 * (double)particles_per_side * particle_dist;

    // fluffy target

    while(particle_counter < number_of_cluster_particles)
    {
        if( sim->get_random_zero_one_incl() <= filling_factor)
        {
            sim->pos_old[X_COORD(particle_counter)] = temp[0];
            sim->pos_old[Y_COORD(particle_counter)] = temp[1];
            sim->pos_old[Z_COORD(particle_counter)] = temp[2];

            ++particle_counter;
        }

        temp[0] += particle_dist;
        ++x_pos;

        if(x_pos >= particles_per_side)
        {
            x_pos = 0;
            ++y_pos;

            temp[0] = - 0.5 * (double)particles_per_side * particle_dist;
            temp[1] += particle_dist;

            if(y_pos >= particles_per_side)
            {
                y_pos = 0;
                temp[1] = - 0.5 * (double)particles_per_side * particle_dist;
                temp[2] += particle_dist;
            }
        }
    }

    // impacting projectile
    /*sim->pos_old[4*number_of_cluster_particles] = (0.5 * particles_per_side + 1) * particle_dist;
    sim->pos_old[4*number_of_cluster_particles+1] = 0.5 * particle_dist * sin(M_PI / 3.0);
    sim->pos_old[4*number_of_cluster_particles+2] = 0;

    sim->vel[4*number_of_cluster_particles] = -impact_speed;


    double phi;

    for(int p = number_of_cluster_particles+1; p < sim->number_of_particles; ++p)
    {
        phi = angular_irregularity - (double)(rand()%100) * 0.02 * angular_irregularity;

        sim->pos_old[X_COORD(p)] = sim->pos_old[4*(p-1)] + particle_dist * cos(phi);
        sim->pos_old[Y_COORD(p)] = sim->pos_old[4*(p-1)+1] - particle_dist * sin(phi);
        sim->pos_old[Z_COORD(p)] = 0;

        sim->vel[X_COORD(p)] = -impact_speed;
    }*/

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

ErrorCode SimLib::initChainBox(Simulation *sim, int number_of_particles, double filling_factor, double angular_irregularity)
{
    if(filling_factor <= 0)
        return EC_INVALID_PARAMETER;

    angular_irregularity *= M_PI;

    // prepare sim
    sim->resizeArrays(number_of_particles, 0);

    // determine sizes
    double particle_dist = 2.0 * particle_radius - delta_0;

    int chain_length = pow((double)number_of_particles * (0.5 / filling_factor)*(0.5 / filling_factor), 1.0/3.0);
    int x_seeds = 1 + (double)chain_length*filling_factor * 2.0;
    int y_seeds = 1 + (double)chain_length*filling_factor * 2.0;

    double x_seed_dist = 2.0 * particle_radius / (2.0 * filling_factor);
    double y_seed_dist = 2.0 * particle_radius / (2.0 * filling_factor);

    int particle_counter = 0;
    int pos_counter;
    bool valid_pos = true;
    int x_pos = 0;
    int y_pos = 0;
    int z_pos = 0;

    double phi, theta;

    vec3 last_pos, next_pos;
    next_pos[0] = - 0.5 * (double)x_seeds * particle_dist;
    next_pos[1] = - 0.5 * (double)y_seeds * particle_dist;
    next_pos[2] = - 0.5 * (double)chain_length * particle_dist;

    // fluffy target
    while(particle_counter < number_of_particles)
    {
        // add particle
        sim->grid.addParticle(next_pos, particle_counter);
        sim->pos_old[X_COORD(particle_counter)] = next_pos[0];
        sim->pos_old[Y_COORD(particle_counter)] = next_pos[2];
        sim->pos_old[Z_COORD(particle_counter)] = next_pos[1];
        memcpy(last_pos, next_pos, sizeof(vec3));
        ++particle_counter;

        // select pos of next particle
        pos_counter = 0;
        valid_pos = false;

        while(pos_counter < 10)
        {
            ++pos_counter;

            phi = ( sim->get_random_zero_one_incl() - 0.5) * angular_irregularity;
            theta = (sim->get_random_zero_one_incl() - 0.5) * angular_irregularity; // TODO: check if this cosine distribution is intentional

            next_pos[0] = last_pos[0] + particle_dist * sin(theta);
            next_pos[1] = last_pos[1] + particle_dist * sin(phi) * cos(theta);
            next_pos[2] = last_pos[2] + particle_dist * cos(phi) * cos(theta);

            // check if not too close to existing particles
            if(sim->grid.canAddParticleAt(next_pos, sim->pos_old))
            {
                valid_pos = true;
                break;
            }
        }

        ++z_pos;

        if(z_pos >= chain_length || !valid_pos)
        {
            z_pos = 0;
            ++x_pos;

            if(x_pos >= x_seeds)
            {
                x_pos = 0;
                ++y_pos;
            }

            next_pos[0] = - 0.5 * (double)x_seeds * particle_dist + (double)x_pos * x_seed_dist;
            next_pos[1] = - 0.5 * (double)y_seeds * particle_dist + (double)y_pos * y_seed_dist;
            next_pos[2] = - 0.5 * (double)chain_length * particle_dist;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    sim->initSimData(SIM_TYPE_GENERAL);
    sim->sim_info.info_storage[5] = 2.0*(double)sim->getNumberOfContacts()/(double)sim->number_of_particles;

    return EC_OK;
}

ErrorCode SimLib::initImpactChainOnAgglomerate(Simulation *sim, int number_of_chain_particles, double angular_irregularity, double impact_speed, bool random_orientation)
{
    if(sim->number_of_particles < 1)
        return EC_NO_PARTICLES;

    angular_irregularity *= 0.5 * M_PI;

    int existing_particles = sim->number_of_particles;

    //copy positions of existing particles
    double* buffer = new double[3*existing_particles];

    memcpy(buffer, sim->pos_old, 3 * existing_particles * sizeof(double));

    // set up arrays for new number of particles
    sim->resizeArrays(existing_particles + number_of_chain_particles, 0);

    // paste old positions
    memcpy(sim->pos_old, buffer, 3 * existing_particles * sizeof(double));
    delete [] buffer;

    // set pos of new particles to zero
    memset(&(sim->pos_old[X_COORD(existing_particles)]), 0, 3 * number_of_chain_particles * sizeof(double));

    /////////////////////////////////////////////////////////////////////////////////////////////
    // get spatial data of current particle agglomerate
    /////////////////////////////////////////////////////////////////////////////////////////////

    SimLib::centerCMS(sim);

    if(random_orientation)
        SimLib::rotateSimRandomly(sim);

    double x_max = sim->pos_old[X_COORD(0)];

    for(int p = 1; p < existing_particles; ++p)
    {
        if(sim->pos_old[X_COORD(p)] > x_max)
            x_max = sim->pos_old[X_COORD(p)];
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up impacting chain
    /////////////////////////////////////////////////////////////////////////////////////////////

    double phi, theta;
    double chain_dist = 2.0 * particle_radius - delta_0;

    sim->pos_old[X_COORD(existing_particles)] = 6.0 * particle_radius + x_max;
    sim->pos_old[Y_COORD(existing_particles)] = 0.0;
    sim->pos_old[Z_COORD(existing_particles)] = 0.0;
    sim->vel[X_COORD(existing_particles)] = -impact_speed;
    sim->vel[Y_COORD(existing_particles)] = 0.0;
    sim->vel[Z_COORD(existing_particles)] = 0.0;

    for(int p = existing_particles+1; p < sim->number_of_particles; ++p)
    {
        phi = angular_irregularity - 2.0 * sim->get_random_zero_one_incl() * angular_irregularity;
        theta = angular_irregularity - 2.0 * sim->get_random_zero_one_incl() * angular_irregularity;

        sim->pos_old[X_COORD(p)] = sim->pos_old[X_COORD(p-1)] + chain_dist * cos(phi) * cos(theta);
        sim->pos_old[Y_COORD(p)] = sim->pos_old[Y_COORD(p-1)] - chain_dist * sin(theta);
        sim->pos_old[Z_COORD(p)] = sim->pos_old[Z_COORD(p-1)] - chain_dist * sin(phi) * cos(theta);

        sim->vel[X_COORD(p)] = -impact_speed;
        sim->vel[Y_COORD(p)] = 0.0;
        sim->vel[Z_COORD(p)] = 0.0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    sim->initSimData(SIM_TYPE_COLLISION, impact_speed, 0.0, 0.0, 0.0, 0.0);

    return EC_OK;
}

ErrorCode SimLib::collideAgglomerateWithWall(Simulation *sim, const char *filename, double impact_speed, int impact_angle, double impact_distance, bool random_orientation, double wall_rolling_modifier, double wall_sliding_modifier, double wall_twisting_modifier)
{
    if(filename)
    {
        ErrorCode error_code = sim->loadFromFile(filename);

        if(error_code != EC_OK)
            return error_code;
    }

    if(sim->number_of_particles < 1)
        return EC_NO_PARTICLES;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // put particles in proper locations
    /////////////////////////////////////////////////////////////////////////////////////////////

    vec3 center_of_mass;
    SimLib::getCenterOfMass(&center_of_mass, *sim, 0, sim->number_of_particles-1);

    // translate agglomerates to make sure cms = 0
    SimLib::centerCMS(sim);

    if(random_orientation)
        SimLib::rotateSimRandomly(sim);

    // determine spatial extension
    vec3 lower, upper;
    sim->getEnclosingBox(&lower, &upper);

    // set impact speed
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        sim->vel[X_COORD(p)] = impact_speed * sin( (double)impact_angle * M_PI / 180.0 );
        sim->vel[Y_COORD(p)] = - impact_speed * cos( (double)impact_angle * M_PI / 180.0 );
        sim->vel[Z_COORD(p)] = 0.0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up the wall
    /////////////////////////////////////////////////////////////////////////////////////////////

    sim->walls.clear();
    sim->number_of_walls = 1;

    // determine offset
    double v_dist = impact_distance * particle_radius - lower[1];
    double h_dist = v_dist * tan( (double)impact_angle * M_PI / 180.0 );

    Wall wall;
    wall.pos[0] = -128.0 * particle_radius + h_dist;
    wall.pos[1] = -v_dist;
    wall.pos[2] = -128.0 * particle_radius;

    wall.normal[0] = 0.0;
    wall.normal[1] = 1.0;
    wall.normal[2] = 0.0;
    wall.alpha = 1.0;

    vec3 x_dir = {1.0, 0.0, 0.0};
    vec3 y_dir = {0.0, 0.0, 1.0};
    wall.setSize(x_dir, y_dir, 256.0 * particle_radius, 256.0 * particle_radius);

    wall.rolling_modifier = wall_rolling_modifier;
    wall.sliding_modifier = wall_sliding_modifier;
    wall.twisting_modifier = wall_twisting_modifier;

    sim->walls.push_back(wall);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    sim->initSimData(SIM_TYPE_WALL_COLLISION, impact_speed, 0.0, 0.0, 0.0, 0.0);

    return EC_OK;
}

ErrorCode SimLib::initBox(Simulation *sim, const char *filename, bool side_walls, double side_wall_height_modifier, bool top_wall, double top_wall_speed, double dynamic_pressure)
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // load agglomerate & determine spatial extension
    /////////////////////////////////////////////////////////////////////////////////////////////

    ErrorCode error_code;

    if(filename)
    {
        error_code = sim->loadFromFile(filename);

        if(error_code != EC_OK)
            return error_code;
    }

    if(sim->number_of_particles < 1)
        return EC_NO_PARTICLES;

    sim->removeWalls();

    vec3 center_of_mass;
    SimLib::getCenterOfMass(&center_of_mass, *sim, 0, sim->number_of_particles-1);

    // translate agglomerates to make sure cms = 0
    SimLib::centerCMS(sim);

    // get enclosing box
    vec3 min_pos, max_pos;

    error_code = sim->getEnclosingBox(&min_pos, &max_pos);

    if(error_code != EC_OK)
        return error_code;

    // reset speed
    //memset(sim->vel, 0, 3*sim->number_of_particles * sizeof(double));
    //memset(sim->vel_angular, 0, 3*sim->number_of_particles * sizeof(double));

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up the bounding box
    /////////////////////////////////////////////////////////////////////////////////////////////

    double wall_x_min = min_pos[0];
    double wall_z_min = min_pos[2];
    double wall_x_max = max_pos[0];
    double wall_z_max = max_pos[2];

    /*if(!side_walls)
    {
        wall_x_min -= 0.75 * (max_pos[0]-min_pos[0]);
        wall_x_max += 0.75 * (max_pos[0]-min_pos[0]);
        wall_z_min -= 0.75 * (max_pos[2]-min_pos[2]);
        wall_z_max += 0.75 * (max_pos[2]-min_pos[2]);
    }*/

    Wall wall;
    vec3 x_dir;
    vec3 y_dir;

    // bottom
    wall.pos[0] = wall_x_min;
    wall.pos[1] = min_pos[1] - 1.0001 * particle_radius;
    wall.pos[2] = wall_z_min;

    wall.normal[0] = 0.0;
    wall.normal[1] = 1.0;
    wall.normal[2] = 0.0;

    x_dir[0] = 1.0;
    x_dir[1] = 0.0;
    x_dir[2] = 0.0;

    y_dir[0] = 0.0;
    y_dir[1] = 0.0;
    y_dir[2] = 1.0;

    wall.setSize(x_dir, y_dir, wall_x_max - wall_x_min, wall_z_max - wall_z_min);
    wall.velocity[0] = 0; wall.velocity[1] = 0; wall.velocity[2] = 0;
    wall.mass = 0;

    wall.alpha = 1.0;
    wall.track_force = true;

    sim->walls.push_back(wall);

    if(side_walls)
    {
        wall.track_force = false;

        // left
        wall.pos[0] = min_pos[0] - 1.0001 * particle_radius;
        wall.pos[1] = min_pos[1];
        wall.pos[2] = min_pos[2];
        wall.normal[0] = 1.0;
        wall.normal[1] = 0.0;
        wall.normal[2] = 0.0;
        x_dir[0] = 0.0;
        x_dir[1] = 0.0;
        x_dir[2] = 1.0;
        y_dir[0] = 0.0;
        y_dir[1] = 1.0;
        y_dir[2] = 0.0;
        wall.setSize(x_dir, y_dir, max_pos[2] - min_pos[2], (max_pos[1] - min_pos[1]) * side_wall_height_modifier);
        wall.alpha = 0.65f;
        sim->walls.push_back(wall);

        //right
        wall.pos[0] = max_pos[0] + 1.0001 * particle_radius;
        wall.pos[1] = min_pos[1];
        wall.pos[2] = min_pos[2];

        wall.normal[0] = -1.0;
        wall.normal[1] = 0.0;
        wall.normal[2] = 0.0;

        x_dir[0] = 0.0;
        x_dir[1] = 1.0;
        x_dir[2] = 0.0;

        y_dir[0] = 0.0;
        y_dir[1] = 0.0;
        y_dir[2] = 1.0;
        wall.setSize(x_dir, y_dir,  (max_pos[1] - min_pos[1]) * side_wall_height_modifier, max_pos[2] - min_pos[2]);
        wall.alpha = 1.0;
        sim->walls.push_back(wall);

        // front
        wall.pos[0] = min_pos[0];
        wall.pos[1] = min_pos[1];
        wall.pos[2] = max_pos[2] + 1.0001 * particle_radius;
        wall.normal[0] = 0.0;
        wall.normal[1] = 0.0;
        wall.normal[2] = -1.0;

        x_dir[0] = 1.0;
        x_dir[1] = 0.0;
        x_dir[2] = 0.0;

        y_dir[0] = 0.0;
        y_dir[1] = 1.0;
        y_dir[2] = 0.0;
        wall.setSize(x_dir, y_dir, max_pos[0] - min_pos[0], (max_pos[1] - min_pos[1]) * side_wall_height_modifier);
        wall.alpha = 0.0f;
        sim->walls.push_back(wall);

        // back
        wall.pos[0] = min_pos[0];
        wall.pos[1] = min_pos[1];
        wall.pos[2] = min_pos[2] - 1.0001 * particle_radius;

        wall.normal[0] = 0.0;
        wall.normal[1] = 0.0;
        wall.normal[2] = 1.0;

        x_dir[0] = 0.0;
        x_dir[1] = 1.0;
        x_dir[2] = 0.0;

        y_dir[0] = 1.0;
        y_dir[1] = 0.0;
        y_dir[2] = 0.0;
        wall.setSize(x_dir, y_dir,  (max_pos[1] - min_pos[1]) * side_wall_height_modifier, max_pos[0] - min_pos[0]);
        wall.alpha = 1.0;
        sim->walls.push_back(wall);
    }

    // top
    if(top_wall)
    {
        if(side_walls)
            wall.alpha = 1.0f;
        else
            wall.alpha = 0.7f;

        if(dynamic_pressure > 0)
            //wall.mass = density * top_wall_thickness * particle_radius * (max_pos[0] - min_pos[0]) * (max_pos[2] - min_pos[2]);
            wall.mass = 10.0 * dynamic_pressure * (max_pos[0] - min_pos[0]) * (max_pos[2] - min_pos[2]);

        wall.pos[0] = wall_x_min;
        wall.pos[1] = min_pos[1] + (max_pos[1] - min_pos[1]) * side_wall_height_modifier + 1.0001 * particle_radius; //max_pos[1] + 1.01 * particle_radius * side_wall_height_modifier;
        wall.pos[2] = wall_z_min;
        wall.normal[0] = 0.0;
        wall.normal[1] = -1.0;
        wall.normal[2] = 0.0;
        wall.velocity[0] = 0.0;
        wall.velocity[1] = -top_wall_speed;
        wall.velocity[2] = 0.0;

        x_dir[0] = 0.0;
        x_dir[1] = 0.0;
        x_dir[2] = 1.0;

        y_dir[0] = 1.0;
        y_dir[1] = 0.0;
        y_dir[2] = 0.0;
        wall.setSize(x_dir, y_dir, wall_z_max - wall_z_min, wall_x_max - wall_x_min);

        wall.track_force = true;
        sim->walls.push_back(wall);
    }

    sim->number_of_walls = (int)sim->walls.size();

    /////////////////////////////////////////////////////////////////////////////////////////////
    // init box
    /////////////////////////////////////////////////////////////////////////////////////////////

    if(top_wall)
    {
        if(side_walls)
            sim->initBox(max_pos[1] - min_pos[1], max_pos[0] - min_pos[0],  max_pos[2] - min_pos[2], min_pos[0], min_pos[1], min_pos[2], 0, 5);
        else
            sim->initBox(max_pos[1] - min_pos[1], max_pos[0] - min_pos[0],  max_pos[2] - min_pos[2], min_pos[0], min_pos[1], min_pos[2], 0, 1);
    }
    else
        sim->initBox(max_pos[1] - min_pos[1], max_pos[0] - min_pos[0],  max_pos[2] - min_pos[2], min_pos[0], min_pos[1], min_pos[2], 0, -1);

    return EC_OK;
}

ErrorCode SimLib::setWallInteractionModifiers(Simulation *sim, double side_wall_compression_modifier, double side_wall_rolling_modifier, double side_wall_sliding_modifier, double top_wall_compression_modifier, double top_wall_rolling_modifier, double top_wall_sliding_modifier)
{
    if(!sim->box)
        return EC_NO_BOX;

    if(sim->box->bottom_wall_id >= 0)
    {
        sim->walls[sim->box->bottom_wall_id].compression_modifier = top_wall_compression_modifier;
        sim->walls[sim->box->bottom_wall_id].rolling_modifier = top_wall_rolling_modifier;
        sim->walls[sim->box->bottom_wall_id].sliding_modifier = top_wall_sliding_modifier;
    }

    if(sim->box->top_wall_id >= 0)
    {
        sim->walls[sim->box->top_wall_id].compression_modifier = top_wall_compression_modifier;
        sim->walls[sim->box->top_wall_id].rolling_modifier = top_wall_rolling_modifier;
        sim->walls[sim->box->top_wall_id].sliding_modifier = top_wall_sliding_modifier;
    }

    if(sim->walls.size() == 5 || sim->walls.size() == 6)
    {
        sim->walls[1].compression_modifier = side_wall_compression_modifier;
        sim->walls[1].rolling_modifier = side_wall_rolling_modifier;
        sim->walls[1].sliding_modifier = side_wall_sliding_modifier;

        sim->walls[2].compression_modifier = side_wall_compression_modifier;
        sim->walls[2].rolling_modifier = side_wall_rolling_modifier;
        sim->walls[2].sliding_modifier = side_wall_sliding_modifier;

        sim->walls[3].compression_modifier = top_wall_compression_modifier;
        sim->walls[3].rolling_modifier = side_wall_rolling_modifier;
        sim->walls[3].sliding_modifier = side_wall_sliding_modifier;

        sim->walls[4].compression_modifier = side_wall_compression_modifier;
        sim->walls[4].rolling_modifier = side_wall_rolling_modifier;
        sim->walls[4].sliding_modifier = side_wall_sliding_modifier;
    }

    return EC_OK;
}

ErrorCode SimLib::initCompressionBox(Simulation *sim, const char *filename, bool side_walls, bool moving_bottom_wall, double wall_speed, double stop_filling_factor, double side_wall_compression_modifier, double side_wall_rolling_modifier, double side_wall_sliding_modifier, double top_wall_compression_modifier, double top_wall_rolling_modifier, double top_wall_sliding_modifier)
{

    ErrorCode error_code = initBox(sim, filename, side_walls, 1.0, true, wall_speed, 0);

    if(error_code == EC_OK)
    {
        SimLib::setWallInteractionModifiers(sim, side_wall_compression_modifier, side_wall_rolling_modifier, side_wall_sliding_modifier, top_wall_compression_modifier, top_wall_rolling_modifier, top_wall_sliding_modifier);

        if(moving_bottom_wall)
            sim->walls[sim->box->bottom_wall_id].velocity[1] = wall_speed;

        // set up data for sim
        if(side_walls)
            error_code = sim->initSimData(SIM_TYPE_COMPRESSION_WITH_SIDE_WALLS, 0.0, stop_filling_factor, 0.0, 0.0, 0.0);
        else
            error_code = sim->initSimData(SIM_TYPE_COMPRESSION_NO_SIDE_WALLS, 0.0, stop_filling_factor, 0.0, 0.0, 0.0);


        return error_code;
    }
    else
        return error_code;
}

ErrorCode SimLib::initCompressionRelaxationBox(Simulation *sim, const char *filename, double wall_speed, double stop_filling_factor, double stop_dissipation_factor, double side_wall_compression_modifier, double side_wall_rolling_modifier, double side_wall_sliding_modifier)
{
    ErrorCode error_code = initBox(sim, filename, true, 1.0, true, wall_speed, 0.0);

    if(error_code == EC_OK)
    {
        SimLib::setWallInteractionModifiers(sim, side_wall_compression_modifier, side_wall_rolling_modifier, side_wall_sliding_modifier, 1.0, 1.0, 1.0);

        // set up data for sim
        return  sim->initSimData(SIM_TYPE_COMPRESSION_RELAXATION, 0.0, stop_filling_factor, stop_dissipation_factor, 0.0, 0.0);
    }
    else
        return error_code;
}

ErrorCode SimLib::initShockwaveBox(Simulation *sim, const char *filename, double compression_speed, double perturbation_layer_thickness, double stop_pressure, double side_wall_rolling_modifier, double side_wall_sliding_modifier)
{
    ErrorCode error_code = initBox(sim, filename, true, 1.0, false, 0.0, 0.0);

    if(error_code == EC_OK)
    {
        SimLib::setWallInteractionModifiers(sim, 1.0, side_wall_rolling_modifier, side_wall_sliding_modifier, 1.0, 1.0, 1.0);

        // get particles on top
        vec3 lower, upper;
        sim->getEnclosingBox(&lower, &upper);

        for(int p = 0; p < sim->number_of_particles; ++p)
        {
            if(sim->pos_old[Y_COORD(p)] > upper[1] - perturbation_layer_thickness)
                sim->vel[Y_COORD(p)] = -compression_speed;
        }

        // set up data for sim
        return sim->initSimData(SIM_TYPE_SHOCKWAVE, 0.0, 0.0, 0.0, 0.0, stop_pressure);
    }
    else
        return error_code;
}

ErrorCode SimLib::initDynamicCompressionBox(Simulation *sim, const char *filename, double wall_speed, double dynamic_pressure, double side_wall_compression_modifier, double side_wall_rolling_modifier, double side_wall_sliding_modifier)
{
    ErrorCode error_code = initBox(sim, filename, true, 1.0, true, wall_speed, dynamic_pressure);

    if(error_code == EC_OK)
    {
        SimLib::setWallInteractionModifiers(sim, side_wall_compression_modifier, side_wall_rolling_modifier, side_wall_sliding_modifier, 1.0, 1.0, 1.0);

        return sim->initSimData(SIM_TYPE_DYNAMIC_COMPRESSION, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    else
        return error_code;
}

ErrorCode SimLib::initOpenBox(Simulation *sim, const char* filename, bool side_walls, double side_wall_height_modifier, double side_wall_compression_modifier, double side_wall_rolling_modifier, double side_wall_sliding_modifier)
{
    ErrorCode error_code = initBox(sim, filename, side_walls, side_wall_height_modifier, false, 0.0, 0.0);

    if(error_code == EC_OK)
    {
        SimLib::setWallInteractionModifiers(sim, side_wall_compression_modifier, side_wall_rolling_modifier, side_wall_sliding_modifier, 1.0, 1.0, 1.0);

        return sim->initSimData(SIM_TYPE_OPEN_BOX, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    else
        return error_code;
}

ErrorCode SimLib::initCompactionBox(Simulation *sim, const char *filename, double wall_speed, double stop_filling_factor)
{
    ErrorCode error_code = initBox(sim, filename, true, 1.0, true, 1.0, 0.0);

    if(error_code == EC_OK)
    {
        SimLib::setWallInteractionModifiers(sim, 1.0, 0.001, 0.001, 1.0, 0.001, 0.001);

        for(size_t w = 0; w < sim->walls.size(); ++w)
        {
            sim->walls[w].velocity[0] = wall_speed * sim->walls[w].normal[0];
            sim->walls[w].velocity[1] = wall_speed * sim->walls[w].normal[1];
            sim->walls[w].velocity[2] = wall_speed * sim->walls[w].normal[2];
        }
        return sim->initSimData(SIM_TYPE_COMPACTION, 0.0, stop_filling_factor, 0.0, 0.0, 0.0);
    }
    else
        return error_code;
}

ErrorCode SimLib::initPullStrengthTest(Simulation *sim, double pull_speed)
{
    if(sim->walls.size() != 6)
        return EC_INVALID_BOX_TYPE;

    // increase stickiness of particles on the walls
    for(int w = 0; w < 6; ++w)
    {
        sim->walls[w].rolling_modifier = 0.0;
        sim->walls[w].sliding_modifier = 0.0;
        sim->walls[w].compression_modifier = 1.0;
        sim->walls[w].alpha = 0;
    }

#ifdef ENABLE_WALL_GLUE
    sim->walls[sim->box->top_wall_id].compression_modifier = wall_glue_strength;
    sim->walls[sim->box->bottom_wall_id].compression_modifier = wall_glue_strength;
#else
    sim->walls[sim->box->top_wall_id].compression_modifier = 1.0;
    sim->walls[sim->box->bottom_wall_id].compression_modifier = 1.0;
#endif

    sim->walls[sim->box->top_wall_id].rolling_modifier = 1.0;
    sim->walls[sim->box->top_wall_id].sliding_modifier = 1.0;
    sim->walls[sim->box->top_wall_id].alpha = 0.7f;

    sim->walls[sim->box->bottom_wall_id].rolling_modifier = 1.0;
    sim->walls[sim->box->bottom_wall_id].sliding_modifier = 1.0;
    sim->walls[sim->box->bottom_wall_id].alpha = 1.0f;

    // change speed of top wall
    sim->walls[sim->box->top_wall_id].velocity[1] = pull_speed;
    sim->walls[sim->box->bottom_wall_id].velocity[1] = 0.0;

    sim->initSimData(SIM_TYPE_PULL_STRENGTH_TEST, 0.0, 0.0, 0.0, 0.0, 0.0);

    return EC_OK;
}

ErrorCode SimLib::initShearStrengthTest(Simulation *sim, double pull_speed)
{
    if(sim->walls.size() != 2)
        return EC_INVALID_BOX_TYPE;

#ifdef ENABLE_WALL_GLUE
    sim->walls[sim->box->top_wall_id].compression_modifier = wall_glue_strength;
    sim->walls[sim->box->bottom_wall_id].compression_modifier = wall_glue_strength;
#else
    sim->walls[sim->box->top_wall_id].compression_modifier = 1.0;
    sim->walls[sim->box->bottom_wall_id].compression_modifier = 1.0;
#endif

    sim->walls[sim->box->top_wall_id].rolling_modifier = 1.0;
    sim->walls[sim->box->top_wall_id].sliding_modifier = 1.0;
    sim->walls[sim->box->top_wall_id].alpha = 0.7f;

    sim->walls[sim->box->bottom_wall_id].rolling_modifier = 1.0;
    sim->walls[sim->box->bottom_wall_id].sliding_modifier = 1.0;
    sim->walls[sim->box->bottom_wall_id].alpha = 1.0f;

    // change speed of top wall
    sim->walls[sim->box->top_wall_id].velocity[0] = pull_speed;
    sim->walls[sim->box->top_wall_id].velocity[1] = 0.0;
    sim->walls[sim->box->bottom_wall_id].velocity[1] = 0.0;

    sim->initSimData(SIM_TYPE_SHEAR_STRENGTH_TEST, 0.0, 0.0, 0.0, 0.0, 0.0);

    return EC_OK;
}

ErrorCode SimLib::collideTwoAgglomerates(Simulation *sim, const char *filename1, const char *filename2, double impact_speed, double impact_parameter, double impact_distance, bool random_orientation1, bool random_orientation2, bool smart_distance, bool cms, int seed)
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // load both simulations
    /////////////////////////////////////////////////////////////////////////////////////////////

    ErrorCode error_code = sim->loadFromFile(filename1);

    if(error_code != EC_OK)
        return error_code;

    Simulation sim2;
    memcpy(&(sim2.sim_info), &(sim->sim_info), sizeof(SimInfo));

    /* not working in windows
    if(seed == -1)
        sim2.rand_generator.seed(time(NULL));
    else
    */
    sim2.rand_generator.seed(seed);

    error_code = sim2.loadFromFile(filename2);

    if(error_code != EC_OK)
        return error_code;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // determine spatial extension
    /////////////////////////////////////////////////////////////////////////////////////////////

    // translate agglomerates to make sure cms = 0
    SimLib::centerCMS(sim);
    SimLib::centerCMS(&sim2);

    if(random_orientation1)
        SimLib::rotateSimRandomly(sim);

    if(random_orientation2)
        SimLib::rotateSimRandomly(&sim2);

    // determine spatial extension
    vec3 lower1, lower2, upper1, upper2;

    sim->getEnclosingBox(&lower1, &upper1);
    sim2.getEnclosingBox(&lower2, &upper2);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // determine shifts
    /////////////////////////////////////////////////////////////////////////////////////////////

    double x1_shift = 0.0;
    double x2_shift = 0.0;
    double y_shift = 0.0;
    double v1 = impact_speed;
    double v2 = 0.0;

    if(impact_parameter > 0.0)
        y_shift = (upper1[1] - lower2[1]) * impact_parameter;
    else if (impact_parameter < 0.0)
        y_shift =  (upper2[1] - lower1[1]) * impact_parameter;

    if(smart_distance)
    {
        // start with distance where agglomerates cannot overlap
        double current_x_pos = upper1[0] - lower2[0] + particle_radius;

        vec3 pos;
        bool stop = false;

        // try if placing agg2 at a smaller distance to agg1 where agglomerates do not overlap
        while(!stop && current_x_pos > lower1[0])
        {
            current_x_pos -= 0.25 * particle_radius;

            for(int p = 0; p < sim2.number_of_particles; ++p)
            {
                pos[0] = sim2.pos_old[X_COORD(p)] + current_x_pos;
                pos[1] = sim2.pos_old[Y_COORD(p)] + y_shift;
                pos[2] = sim2.pos_old[Z_COORD(p)];

                if( ! sim->grid.canAddParticleAt(pos, sim->pos_old ) )
                {
                    current_x_pos += 0.25 * particle_radius;
                    stop = true;
                    break;
                }
            }
        }

        current_x_pos += impact_distance * particle_radius;

        if(cms)
        {
            x1_shift -= 0.5 * current_x_pos;
            x2_shift += 0.5 * current_x_pos;
        }
        else
            x1_shift -= current_x_pos;
    }
    else
    {
        x1_shift = -upper1[0];
        x2_shift = -lower2[0];

        if(cms)
        {
            x1_shift -= 0.5 * impact_distance * particle_radius;
            x2_shift += 0.5 * impact_distance * particle_radius;

            v1 = impact_speed / (1.0 + (double)sim->number_of_particles / (double)sim2.number_of_particles);
            v2 = -impact_speed / (1.0 + (double)sim2.number_of_particles / (double)sim->number_of_particles);
        }
        else
            x1_shift -= impact_distance * particle_radius;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // setup aggregates
    /////////////////////////////////////////////////////////////////////////////////////////////

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        sim->pos_old[X_COORD(p)] += x1_shift;

        sim->vel[X_COORD(p)] = v1;
        sim->vel[Y_COORD(p)] = 0.0;
        sim->vel[Z_COORD(p)] = 0.0;
    }

    for(int p = 0; p < sim2.number_of_particles; ++p)
    {
        sim2.pos_old[X_COORD(p)] += x2_shift;
        sim2.pos_old[Y_COORD(p)] += y_shift;

        sim2.vel[X_COORD(p)] = v2;
        sim2.vel[Y_COORD(p)] = 0.0;
        sim2.vel[Z_COORD(p)] = 0.0;
    }

    sim->addParticlesFromSim(sim2);
    sim->initSimData(SIM_TYPE_COLLISION, impact_speed, 0.0, 0.0, 0.0, 0.0);
    sim->sim_info.info_storage[6] = sim2.sim_info.info_storage[5];

    return EC_OK;
}

ErrorCode SimLib::collideTwoAgglomerates(
        Simulation *sim,
        const char *material_name,
        const char *file1,
        const char *file2,
        const unsigned int seed1,
        const unsigned int seed2,
        const double impact_speed,
        const double impact_parameter,
        const double impact_distance,
        const bool smart_distance,
        const bool cms
        )
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // load both agglomerates
    /////////////////////////////////////////////////////////////////////////////////////////////

    ErrorCode error_code = sim->loadMaterial(material_name);

    if(error_code != EC_OK)
        return error_code;

    sim->rand_generator.seed(seed1);

    printf("%s\n", file1);

    sim->loadFromFile(file1);


    printf("Sphere1 number of particles/contacts = %d %d\n", sim->number_of_particles, sim->number_of_contacts);

    Simulation sim2;

    error_code = sim2.loadMaterial(material_name);

    if(error_code != EC_OK)
        return error_code;

    sim2.rand_generator.seed(seed2);

    sim2.loadFromFile(file2);
    //printf("num = %d %d\n", sim2.number_of_particles, sim2.number_of_contacts);



    /////////////////////////////////////////////////////////////////////////////////////////////
    // determine spatial extension
    /////////////////////////////////////////////////////////////////////////////////////////////

    // translate agglomerates to make sure cms = 0
    SimLib::centerCMS(sim);
    SimLib::centerCMS(&sim2);



    SimLib::rotateSimRandomly(sim);
    SimLib::rotateSimRandomly(&sim2);


    double cross_section1;
    double cross_section2;
    double sigma_cross_section;

    SimLib::getCrossSection(*sim, 0.1*particle_radius, 12, &cross_section1, &sigma_cross_section);

    SimLib::getCrossSection(sim2, 0.1*particle_radius, 12, &cross_section2, &sigma_cross_section);


    double radius1 = sqrt(cross_section1 / M_PI);
    double radius2 = sqrt(cross_section2 / M_PI);

    /*
    double outer_radius;
    double gyration_radius;

    SimLib::getSize(*sim, &gyration_radius, &outer_radius);
    double radius1 = sqrt(5.0/3.0)*gyration_radius;


    SimLib::getSize(sim2, &gyration_radius, &outer_radius);
    double radius2 = sqrt(5.0/3.0)*gyration_radius;
    */

    // determine spatial extension
    vec3 lower1, lower2, upper1, upper2;

    sim->getEnclosingBox(&lower1, &upper1);
    sim2.getEnclosingBox(&lower2, &upper2);






    /////////////////////////////////////////////////////////////////////////////////////////////
    // determine shifts
    /////////////////////////////////////////////////////////////////////////////////////////////

    double x1_shift = 0.0;
    double x2_shift = 0.0;
    double y_shift = 0.0;
    double v1 = impact_speed;
    double v2 = 0.0;


    /*
    if(impact_parameter > 0.0)
        y_shift = (upper1[1] - lower2[1]) * impact_parameter;
    else if (impact_parameter < 0.0)
        y_shift =  (upper2[1] - lower1[1]) * impact_parameter;
    */

    y_shift = (radius1 + radius2)*impact_parameter;


    if(smart_distance)
    {
        // start with distance where agglomerates cannot overlap
        double current_x_pos = upper1[0] - lower2[0] + particle_radius;

        vec3 pos;
        bool stop = false;

        // try if placing agg2 at a smaller distance to agg1 where agglomerates do not overlap
        while(!stop && current_x_pos > lower1[0])
        {
            current_x_pos -= 0.25 * particle_radius;

            for(int p = 0; p < sim2.number_of_particles; ++p)
            {
                pos[0] = sim2.pos_old[X_COORD(p)] + current_x_pos;
                pos[1] = sim2.pos_old[Y_COORD(p)] + y_shift;
                pos[2] = sim2.pos_old[Z_COORD(p)];

                if( ! sim->grid.canAddParticleAt(pos, sim->pos_old ) )
                {
                    current_x_pos += 0.25 * particle_radius;
                    stop = true;
                    break;
                }
            }
        }

        current_x_pos += impact_distance;

        if(cms)
        {
            x1_shift -= 0.5 * current_x_pos;
            x2_shift += 0.5 * current_x_pos;

            v1 = impact_speed / (1.0 + (double)sim->number_of_particles / (double)sim2.number_of_particles);
            v2 = -impact_speed / (1.0 + (double)sim2.number_of_particles / (double)sim->number_of_particles);
        }
        else
            x1_shift -= current_x_pos;
    }
    else
    {
        x1_shift = -upper1[0];
        x2_shift = -lower2[0];

        if(cms)
        {
            x1_shift -= 0.5 * impact_distance * particle_radius;
            x2_shift += 0.5 * impact_distance * particle_radius;

            v1 = impact_speed / (1.0 + (double)sim->number_of_particles / (double)sim2.number_of_particles);
            v2 = -impact_speed / (1.0 + (double)sim2.number_of_particles / (double)sim->number_of_particles);
        }
        else
            x1_shift -= impact_distance * particle_radius;
    }


    printf("velocities = %e %e  num_particles = %d  %d	impact_dist = %e\n", v1, v2, sim->number_of_particles, sim2.number_of_particles, x1_shift/particle_radius);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // setup aggregates
    /////////////////////////////////////////////////////////////////////////////////////////////

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        sim->pos_old[X_COORD(p)] += x1_shift;

        sim->vel[X_COORD(p)] = v1;
        sim->vel[Y_COORD(p)] = 0.0;
        sim->vel[Z_COORD(p)] = 0.0;
    }

    for(int p = 0; p < sim2.number_of_particles; ++p)
    {
        sim2.pos_old[X_COORD(p)] += x2_shift;
        sim2.pos_old[Y_COORD(p)] += y_shift;

        sim2.vel[X_COORD(p)] = v2;
        sim2.vel[Y_COORD(p)] = 0.0;
        sim2.vel[Z_COORD(p)] = 0.0;
    }


    sim->addParticlesFromSim(sim2);

    sim->initSimData(SIM_TYPE_COLLISION, impact_speed, 0.0, 0.0, 0.0, 0.0);
    sim->sim_info.info_storage[6] = sim2.sim_info.info_storage[5];

    return EC_OK;
}

ErrorCode SimLib::generateRBDSphere(
        Simulation *sim,
        const char *material_name,
        const double agg_mass,
        const unsigned int seed1
        )
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // generate both agglomerates
    /////////////////////////////////////////////////////////////////////////////////////////////

    printf("Generating RBD sphere\n");
    ErrorCode error_code = sim->loadMaterial(material_name);

    if(error_code != EC_OK)
        return error_code;

    sim->rand_generator.seed(seed1);



    double init_fill_factor = 0.14; // for RBD


    const int number_of_particles = int(0.5 + agg_mass / (4.0/3.0 * 3.1415 * density * particle_radius*particle_radius*particle_radius));

    double agg_size = pow(3.0 / (4.0 * M_PI) * agg_mass / (density * init_fill_factor), 1.0/3.0);

    double init_size = 2.0 * agg_size + 10.0e-4; // make init_size bigger to guarantee that there are enough particles to get correct agglomerate size

    const int init_num = (init_size + 10.0e-4) * init_size*init_size*init_fill_factor / (4.0/3.0 * 3.1415 * particle_radius*particle_radius*particle_radius);

    printf("Size = %e   Mass =  %e\n", agg_size, agg_mass);
    printf("init number of particles:   %d  %d\n",  init_num, number_of_particles);
    error_code = SimLib::initBAMCake(sim, init_num, init_size, init_size, 0.0, 0.0);


    vec3 center;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;

 slice:
    agg_size -= 1.e-5;
    SimLib::sliceSphere(sim,  center, agg_size);

    if(sim->number_of_particles > number_of_particles)
      goto slice;

    /*
    double rr, gr;

    SimLib::getSize(*sim, &gr, &rr);



    double cross_section1;
    double sigma_cross_section;
    SimLib::getCrossSection(*sim, 0.1*particle_radius, 12, &cross_section1, &sigma_cross_section);


    double radius1 = sqrt(cross_section1 / M_PI);


    gr *= sqrt(5.0/3.0);
    printf("radius = %e %e  %e\n", gr, rr, radius1);
    */

    error_code = SimLib::initBAMAggregate(sim, NULL, number_of_particles - sim->number_of_particles, 0.0, 0.0, BAM_SELECT_INTERIOR);

    if(error_code != EC_OK)
        return error_code;

    return EC_OK;
}



ErrorCode SimLib::generateBAMSphere(
        Simulation *sim,
        const char *material_name,
        const double agg_mass,
        const unsigned int seed1
        )
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // generate both agglomerates
    /////////////////////////////////////////////////////////////////////////////////////////////

    ErrorCode error_code = sim->loadMaterial(material_name);

    if(error_code != EC_OK)
        return error_code;

    sim->rand_generator.seed(seed1);



    double init_fill_factor = 0.40; // for RBD random migration


    const int number_of_particles = int(0.5 + agg_mass / (4.0/3.0 * 3.1415 * density * particle_radius*particle_radius*particle_radius));
    double agg_size = pow(3.0 / (4.0 * M_PI) * agg_mass / (density * init_fill_factor), 1.0/3.0);

    double init_size = 2.0 * agg_size + 10.0e-4; // make init_size bigger to guarantee that there are enough particles to get correct agglomerate size

    init_fill_factor = 0.45; // for RBD random migration

    const int init_num = (init_size + 10.0e-4) * init_size*init_size*init_fill_factor / (4.0/3.0 * 3.1415 * particle_radius*particle_radius*particle_radius);

    printf("generating BAM sphere   %e  %e\n", agg_size, agg_mass);
    printf("init number of particles:   %d  %d\n",  init_num, number_of_particles);
    error_code = SimLib::initBAMCake(sim, init_num, init_size, init_size, 1.0, 1.0);


    vec3 center;
    center[0] = 0.0;
    center[1] = 0.0;
    center[2] = 0.0;

 slice:
    if(sim->number_of_particles > 5000 + number_of_particles)
        agg_size -= 1.e-5;
    else if (sim->number_of_particles > 500 + number_of_particles)
        agg_size -= 1.e-6;
    else
        agg_size -= 1.e-7;

    SimLib::sliceSphere(sim,  center, agg_size);
    printf("Slicing off, current number of particles = %d\n", sim->number_of_particles);

    if(sim->number_of_particles - number_of_particles > 10)
      goto slice;

    /*
    double rr, gr;

    SimLib::getSize(*sim, &gr, &rr);



    double cross_section1;
    double sigma_cross_section;
    SimLib::getCrossSection(*sim, 0.1*particle_radius, 12, &cross_section1, &sigma_cross_section);


    double radius1 = sqrt(cross_section1 / M_PI);


    gr *= sqrt(5.0/3.0);
    printf("radius = %e %e  %e\n", gr, rr, radius1);
    */

    printf("Particles: %d / %d\n", sim->number_of_particles, number_of_particles);
    //error_code = SimLib::initBAMAggregate(sim, NULL, number_of_particles - sim->number_of_particles, 1.0, 1.0, BAM_SELECT_RANDOM);

    if(error_code != EC_OK)
        return error_code;

    return EC_OK;
}




ErrorCode SimLib::generateBAMsphereAgglomerate(
        Simulation *sim,
        const char *material_name,
        const char *file1,
        double agg_size,
        const double agg_mass,
        const unsigned int seed1
        )
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // generate both agglomerates
    /////////////////////////////////////////////////////////////////////////////////////////////

    ErrorCode error_code = sim->loadMaterial(material_name);

    if(error_code != EC_OK)
        return error_code;

    sim->rand_generator.seed(seed1);


    double rim = 0.0 * particle_radius; // width of compacted ring build with BAM most interior
    double init_fill_factor = 0.4; // arbitrary chosen upper limit for the bam2 filling factor


    const double init_size = agg_size + 20.0e-4 - rim; // make init_size bigger to guarantee that there are enough particles to get correct agglomerate size

    int number_of_particles = 4.0/3.0 * 3.1415 * init_size*init_size*init_size*init_fill_factor / (4.0/3.0 * 3.1415 * particle_radius*particle_radius*particle_radius);

    if(agg_mass > 0.0)
    {
            number_of_particles = int(0.5 + agg_mass / (4.0/3.0 * 3.1415 * density * particle_radius*particle_radius*particle_radius));
            agg_size = pow(3.0 / (4.0 * M_PI) * agg_mass / density, 1.0/3.0);
    }

    const double extra_fill_factor = 0.6; // arbitrary chosen upper limit for the bam2 filling factor
    int extra_particles = int( 4.0/3.0 * 3.1415 * (agg_size*agg_size*agg_size*extra_fill_factor - (agg_size-rim)*(agg_size-rim)*(agg_size-rim)*extra_fill_factor)
            / (4.0/3.0 * 3.1415 * particle_radius*particle_radius*particle_radius));

    printf("init number of particles:   %d  extra particles in rim: %d\n", number_of_particles, extra_particles);


    printf("generating BAM sphere\n");
    error_code = SimLib::initBAMAggregate(sim, NULL, number_of_particles - extra_particles, 1.0, 1.0, BAM_SELECT_CLOSEST);

    if(error_code != EC_OK)
        return error_code;


    //  do not slice aggregate to create a rough surface
    //vec3 center_of_mass;
    //SimLib::getCenterOfMass(&center_of_mass, *sim, 0, sim->number_of_particles-1);
    //SimLib::sliceSphere(sim, center_of_mass, agg_size - rim);
    

    if(extra_particles > 0)
    {
        error_code = SimLib::initBAMAggregate(sim, NULL, extra_particles, 1.0, 1.0, BAM_SELECT_INTERIOR);
        if(error_code != EC_OK)
            return error_code;
    }


    sim->saveToFile(file1, true);

    return EC_OK;
}



ErrorCode SimLib::impactMultiProjectiles(Simulation *sim, const char *target_filename, bool random_target_orientation, const char *projectile_filename, unsigned int projectile_count, unsigned int projectile_size, double impact_velocity, double impact_parameter, bool random_impact_parameter, double min_distance, double max_distance, bool smart_distance)
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // load both simulations
    /////////////////////////////////////////////////////////////////////////////////////////////

    double outer_radius1, outer_radius2;

    ErrorCode error_code = sim->loadFromFile(target_filename);

    if(error_code != EC_OK)
        return error_code;

    SimLib::centerCMS(sim);
    SimLib::getSize(*sim, NULL, &outer_radius1);

    if(random_target_orientation)
        SimLib::rotateSimRandomly(sim);

    Simulation sim2;
    memcpy(&(sim2.sim_info), &(sim->sim_info), sizeof(SimInfo));

    if(projectile_filename)
    {
        error_code = sim2.loadFromFile(projectile_filename);

        if(error_code != EC_OK)
            return error_code;

        SimLib::centerCMS(&sim2);
        SimLib::getSize(sim2, NULL, &outer_radius2);

        projectile_size = sim2.number_of_particles;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // add projectiles
    /////////////////////////////////////////////////////////////////////////////////////////////

    int target_size = sim->number_of_particles;
    int next_particle_id = sim->number_of_particles;
    sim->addParticles(projectile_count * projectile_size);
    double eq_dist = 2.0 * particle_radius - delta_0;

    for(int projectile = 0; projectile < projectile_count; ++projectile)
    {
        if(projectile_filename)
            SimLib::rotateSimRandomly(&sim2);
        else
        {
            sim2.cleanUp();
            SimLib::initBAMAggregate(&sim2, NULL, projectile_size, 0.0, 0.0, BAM_SELECT_RANDOM);
            SimLib::centerCMS(&sim2);
            SimLib::getSize(sim2, NULL, &outer_radius2);
        }

        vec3 dir;
        vec3 pos, particle_pos;
        vec3 n;

        bool projectile_added = false;

        while(!projectile_added)
        {
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // find trajectory that hits the target
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

find_trajectory:

            // determine random direction
            double phi = sim->get_random_zero_twoPi();
            double cos_theta = sim->get_random_cos_theta();
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            dir[0] = sin_theta * cos(phi);
            dir[1] = sin_theta * sin(phi);
            dir[2] = cos_theta;

            SimLib::getOrthoVector(sim, &n, dir);

            if(random_impact_parameter)
            {
                double x = sim->get_random_zero_one();
                impact_parameter = exp( log(x * (pow(0.9,3)-pow(0.0,3)) + pow(0.0,3)) /(3.0) ); // TODO: pow(0.0,3) ???
            }

            pos[0] = -(outer_radius1+outer_radius2) * dir[0] + impact_parameter * outer_radius1 * n[0];
            pos[1] = -(outer_radius1+outer_radius2) * dir[1] + impact_parameter * outer_radius1 * n[1];
            pos[2] = -(outer_radius1+outer_radius2) * dir[2] + impact_parameter * outer_radius1 * n[2];

            double min_value = 1000.0;
            for(int j = 0; j < target_size; ++j)
            {
                double a = (sim->pos_old[X_COORD(j)]-pos[0]) * dir[0] + (sim->pos_old[Y_COORD(j)]-pos[1]) * dir[1] + (sim->pos_old[Z_COORD(j)]-pos[2]) * dir[2];
                double b = (sim->pos_old[X_COORD(j)]-pos[0])*(sim->pos_old[X_COORD(j)]-pos[0]) + (sim->pos_old[Y_COORD(j)]-pos[1])*(sim->pos_old[Y_COORD(j)]-pos[1]) + (sim->pos_old[Z_COORD(j)]-pos[2])*(sim->pos_old[Z_COORD(j)]-pos[2]);
                double c = a*a-b+eq_dist*eq_dist;

                if(c >= 0)
                {
                    double t = a - sqrt(c);
                    if(t < min_value)
                        goto find_distance;
                }
            }

            // trajetory does not hit the target -> retry
            goto find_trajectory;

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // determine distance
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

find_distance:

            double initial_dist = outer_radius1 + outer_radius2; // - outer_radius1 * (1.0 - sqrt(1.0 - impact_parameter*impact_parameter));
            double dist = initial_dist;

            if(smart_distance)
            {
                bool target_hit = false;

                while(!target_hit)
                {
                    pos[0] = dist * dir[0] + impact_parameter * outer_radius1 * n[0];
                    pos[1] = dist * dir[1] + impact_parameter * outer_radius1 * n[1];
                    pos[2] = dist * dir[2] + impact_parameter * outer_radius1 * n[2];

                    for(int p = 0; p < sim2.number_of_particles; ++p)
                    {
                        particle_pos[0] = pos[0] + sim2.pos_old[X_COORD(p)];
                        particle_pos[1] = pos[1] + sim2.pos_old[Y_COORD(p)];
                        particle_pos[2] = pos[2] + sim2.pos_old[Z_COORD(p)];

                        if(!sim->grid.canAddParticleAt(particle_pos, sim->pos_old))
                        {
                            target_hit = true;
                            dist += 0.1 * particle_radius;
                            break;
                        }
                    }

                    dist -= 0.05 * particle_radius;

                    if(dist < -initial_dist)
                    {
                        goto find_trajectory;
                        //dist = initial_dist;
                        //target_hit = true;
                    }
                }
            }

            dist += min_distance + (max_distance - min_distance) * sim->get_random_zero_one();
            //dist =  exp( log( dist * (pow( (outer_radius1+outer_radius2 + 1e-4 * impact_distance), 3 ) - pow(outer_radius1+outer_radius2 ,3)) + pow(outer_radius1+outer_radius2,3)) / 3.0 );

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // determine final pos
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            pos[0] = dist * dir[0] + impact_parameter * outer_radius1 * n[0];
            pos[1] = dist * dir[1] + impact_parameter * outer_radius1 * n[1];
            pos[2] = dist * dir[2] + impact_parameter * outer_radius1 * n[2];

            bool can_add_projectile = true;

            for(int p = 0; p < sim2.number_of_particles; ++p)
            {
                particle_pos[0] = pos[0] + sim2.pos_old[X_COORD(p)];
                particle_pos[1] = pos[1] + sim2.pos_old[Y_COORD(p)];
                particle_pos[2] = pos[2] + sim2.pos_old[Z_COORD(p)];

                if(!sim->grid.canAddParticleAt(particle_pos, sim->pos_old))
                {
                    can_add_projectile = false;
                    break;
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // add aggregate
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(can_add_projectile)
            {
                int id_offset = target_size + projectile * sim2.number_of_particles;

                for(int p = 0; p < sim2.number_of_particles; ++p)
                {
                    // add particles
                    particle_pos[0] = pos[0] + sim2.pos_old[X_COORD(p)];
                    particle_pos[1] = pos[1] + sim2.pos_old[Y_COORD(p)];
                    particle_pos[2] = pos[2] + sim2.pos_old[Z_COORD(p)];

                    sim->grid.addParticle(particle_pos, next_particle_id);
                    sim->pos_old[X_COORD(next_particle_id)] = particle_pos[0];
                    sim->pos_old[Y_COORD(next_particle_id)] = particle_pos[1];
                    sim->pos_old[Z_COORD(next_particle_id)] = particle_pos[2];
                    sim->vel[X_COORD(next_particle_id)] = - impact_velocity * dir[0];
                    sim->vel[Y_COORD(next_particle_id)] = - impact_velocity * dir[1];
                    sim->vel[Z_COORD(next_particle_id)] = - impact_velocity * dir[2];

                    // add contacts
                    sim->contact_list[next_particle_id] = NULL;

                    ContactListEntry *new_cl_entry = NULL;
                    ContactListEntry *cl_entry = sim2.contact_list[p];

                    while(cl_entry)
                    {
                        Contact *new_contact = new Contact();
                        *new_contact = *(cl_entry->contact);

                        new_contact->id1 += id_offset;
                        new_contact->id2 += id_offset;

                        if(sim->contact_list[p + id_offset] == NULL) // particle has no other contacts yet
                        {
                            // create contact list entry
                            sim->contact_list[p + id_offset] = new ContactListEntry;
                            sim->contact_list[p + id_offset]->next = NULL;
                            sim->contact_list[p + id_offset]->id = cl_entry->id + id_offset;
                            sim->contact_list[p + id_offset]->contact = new_contact;
                            new_cl_entry = sim->contact_list[p + id_offset];
                        }
                        else
                        {
                            new_cl_entry->next = new ContactListEntry;
                            new_cl_entry->next->next = NULL;
                            new_cl_entry->next->id = cl_entry->id + id_offset;
                            new_cl_entry->next->contact = new_contact;
                            new_cl_entry = new_cl_entry->next;
                        }

                        cl_entry = cl_entry->next;
                    }

                    ++next_particle_id;
                }

                projectile_added = true;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // setup aggregates
    /////////////////////////////////////////////////////////////////////////////////////////////

    sim->initSimData(SIM_TYPE_COLLISION, impact_velocity, 0.0, 0.0, 0.0, 0.0);

    return EC_OK;
}

ErrorCode SimLib::hitAndStick(Simulation *sim, const char *filename1, const char *filename2, double impact_parameter, bool random_orientation1, bool random_orientation2)
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // load both simulations
    /////////////////////////////////////////////////////////////////////////////////////////////

    ErrorCode error_code;

    if(filename1)
    {
        error_code = sim->loadFromFile(filename1);

        if(error_code != EC_OK)
            return error_code;
    }

    Simulation sim2;
    sim2.setMaterialConstants(particle_radius, density, surface_energy, nu, young_mod, crit_rolling_displacement, osc_damping_factor, rolling_modifier, sliding_modifier, twisting_modifier, crit_sliding_displacement_modifier, crit_wall_sliding_displacement_modifier);
    error_code = sim2.loadFromFile(filename2);

    if(error_code != EC_OK)
        return error_code;

    return SimLib::hitAndStick(sim, &sim2, impact_parameter, random_orientation1, random_orientation2);
}

ErrorCode SimLib::hitAndStick(Simulation *sim, Simulation *sim2, double impact_parameter, bool random_orientation1, bool random_orientation2)
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    // determine spatial extension
    /////////////////////////////////////////////////////////////////////////////////////////////

    memset(sim->vel, 0, sim->number_of_particles * 3 * sizeof(double));
    memset(sim2->vel, 0, sim2->number_of_particles * 3 * sizeof(double));

    // translate agglomerates to make sure cms = 0
    SimLib::centerCMS(sim);
    SimLib::centerCMS(sim2);

    if(random_orientation1)
        SimLib::rotateSimRandomly(sim);

    if(random_orientation2)
        SimLib::rotateSimRandomly(sim2);

    // determine spatial extension
    vec3 lower1, lower2, upper1, upper2;

    sim->getEnclosingBox(&lower1, &upper1);
    sim2->getEnclosingBox(&lower2, &upper2);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // determine pos for aggregate 2
    /////////////////////////////////////////////////////////////////////////////////////////////

    double x_shift;
    double y_shift = impact_parameter * (upper2[1] - lower1[1]);
    double max_x_shift = upper1[0] - lower1[0];

    for(int p = 0; p < sim2->number_of_particles; ++p)
    {
        double new_x_pos = getLowerXPos(sim, sim2->pos_old[Y_COORD(p)] - y_shift, sim2->pos_old[Z_COORD(p)], upper1[0] + 2.0 * particle_radius, lower1[0]);

        x_shift = sim2->pos_old[X_COORD(p)] + upper1[0] - lower2[0] - new_x_pos;

        if(x_shift < max_x_shift)
            max_x_shift = x_shift;
    }

    // shift agglomerate 2
    for(int p = 0; p < sim2->number_of_particles; ++p)
    {
        sim2->pos_old[X_COORD(p)] += upper1[0] - lower2[0] - max_x_shift;
        sim2->pos_old[Y_COORD(p)] -= y_shift;
    }

    sim->addParticlesFromSim(*sim2);
    sim->initSimData(SIM_TYPE_COLLISION, 0.0, 0.0, 0.0, 0.0, 0.0);

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();

    SimLib::centerCMS(sim);

    return EC_OK;
}

ErrorCode SimLib::sandblastAggregate(Simulation *sim, const char *filename, int number_of_projectiles, double impact_velocity, double impact_distance)
{
    double eq_dist = 2.0 * particle_radius - delta_0;

    if(filename)
    {
        ErrorCode error_code = sim->loadFromFile(filename);

        if(error_code != EC_OK)
            return error_code;
    }

    if(sim->number_of_particles <= 0)
        return EC_NO_PARTICLES;

    SimLib::centerCMS(sim);
    double outer_radius;
    SimLib::getSize(*sim, NULL, &outer_radius);

    int target_size = sim->number_of_particles;
    int next_particle_id = sim->number_of_particles;
    sim->addParticles(number_of_projectiles);

    while(next_particle_id < sim->number_of_particles)
    {
        vec3 dir;
        vec3 pos;
        vec3 n;
        double impact_parameter;
        bool impact_on_target = false;

        while(!impact_on_target)
        {
            // determine random direction
            double phi = sim->get_random_zero_twoPi();
            double cos_theta = sim->get_random_cos_theta();
            double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

            dir[0] = sin_theta * cos(phi);
            dir[1] = sin_theta * sin(phi);
            dir[2] = cos_theta;

            // determine impact parameter
            double x = sim->get_random_zero_one();
            impact_parameter = exp( log(x * (pow(1.0,3)-pow(0.0,3)) + pow(0.0,3)) /(3.0) ); // TODO pow(1.0, 3), pow(0.0, 3) ??????
            SimLib::getOrthoVector(sim, &n, dir);

            // check if trajectory hits the target
            pos[0] = -outer_radius * dir[0] + impact_parameter * outer_radius * n[0];
            pos[1] = -outer_radius * dir[1] + impact_parameter * outer_radius * n[1];
            pos[2] = -outer_radius * dir[2] + impact_parameter * outer_radius * n[2];

            double min_value = 1000.0;
            for(int j = 0; j < target_size; ++j)
            {
                double a = (sim->pos_old[X_COORD(j)]-pos[0]) * dir[0] + (sim->pos_old[Y_COORD(j)]-pos[1]) * dir[1] + (sim->pos_old[Z_COORD(j)]-pos[2]) * dir[2];
                double b = (sim->pos_old[X_COORD(j)]-pos[0])*(sim->pos_old[X_COORD(j)]-pos[0]) + (sim->pos_old[Y_COORD(j)]-pos[1])*(sim->pos_old[Y_COORD(j)]-pos[1]) + (sim->pos_old[Z_COORD(j)]-pos[2])*(sim->pos_old[Z_COORD(j)]-pos[2]);
                double c = a*a-b+eq_dist*eq_dist;

                if(c >= 0)
                {
                    double t = a - sqrt(c);
                    if(t < min_value)
                    {
                        impact_on_target = true;
                        break;
                    }
                }
            }
        }

        // determine distance
        double dist = sim->get_random_zero_one();
        dist =  exp( log( dist * (pow( (outer_radius + 1e-4 * impact_distance), 3 ) - pow(outer_radius,3)) + pow(outer_radius,3)) / 3.0 );

        // set final pos
        pos[0] = dist * dir[0] + impact_parameter * outer_radius * n[0];
        pos[1] = dist * dir[1] + impact_parameter * outer_radius * n[1];
        pos[2] = dist * dir[2] + impact_parameter * outer_radius * n[2];

        if(sim->grid.canAddParticleAt(pos, sim->pos_old))
        {
            sim->grid.addParticle(pos, next_particle_id);
            sim->pos_old[X_COORD(next_particle_id)] = pos[0];
            sim->pos_old[Y_COORD(next_particle_id)] = pos[1];
            sim->pos_old[Z_COORD(next_particle_id)] = pos[2];
            sim->vel[X_COORD(next_particle_id)] = - impact_velocity * dir[0];
            sim->vel[Y_COORD(next_particle_id)] = - impact_velocity * dir[1];
            sim->vel[Z_COORD(next_particle_id)] = - impact_velocity * dir[2];
            ++next_particle_id;
        }
    }

    sim->initSimData(SIM_TYPE_GENERAL, impact_velocity, 0.0, 0.0, 0.0, 0.0);

    return EC_OK;
}

void SimLib::detectFragments(Simulation &sim, std::vector<int> *fragment_id, std::vector<int> *size_of_fragment, std::vector< std::list<int> > *particles_of_fragment)
{
    if(sim.number_of_particles < 1)
        return;

    fragment_id->resize(sim.number_of_particles, -1);

    int *todo_list = new int[sim.number_of_particles];
    int current_todo_pos = 0;                   // current pos in the todo list
    int next_free_todo_pos = 0;                 // next free pos in the todo list
    int current_agglomerate = -1;               // id of the current agglomerate
    ContactListEntry *cl_entry;

    // build contact list in which all contacts with other particles (not only those with lower id) are stored
    std::vector< std::list<int> > contacts(sim.number_of_particles);

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        cl_entry = sim.contact_list[p];

        while(cl_entry)
        {
            // ignore contacts with walls
            if(cl_entry->id >= 0)
            {
                contacts[p].push_back(cl_entry->id);
                contacts[cl_entry->id].push_back(p);
            }

            cl_entry = cl_entry->next;
        }
    }

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        // particle does not belong to any agglomerate yet
        if((*fragment_id)[p] == -1)
        {
            // put it on the todo list
            todo_list[0] = p;

            // create a new agglomerate
            ++current_agglomerate;

            // the new agglomerate consist of one particle
            size_of_fragment->push_back(1);

            if(particles_of_fragment)
                particles_of_fragment->push_back( std::list<int>(1, p) );

            // the particle now belongs to the new agglomerate
            (*fragment_id)[p] = current_agglomerate;

            // set the pointers to the todo list
            current_todo_pos = 0;
            next_free_todo_pos = 1;

            // repeat until there is no particle left on the todo list
            while(current_todo_pos < next_free_todo_pos)
            {
                for(std::list<int>::iterator contact = contacts[todo_list[current_todo_pos]].begin(); contact != contacts[todo_list[current_todo_pos]].end(); ++contact)
                {
                    // check if particle has not been added to an agglomerate yet
                    if((*fragment_id)[*contact] == -1)
                    {
                        // add particle to the current agglomerate
                        (*fragment_id)[*contact] = current_agglomerate;

                        // add particle id to the list of ids of the current agglomerate
                        if(particles_of_fragment)
                            (*particles_of_fragment)[current_agglomerate].push_back(*contact);

                         // increase number of particles of that agglomerate
                        (*size_of_fragment)[current_agglomerate] += 1;

                        // add particle to the todo list
                        todo_list[next_free_todo_pos] = *contact;
                        ++next_free_todo_pos;
                    }
                }

                // work on this todo list entry is done
                ++current_todo_pos;
            }
        }
    }

    delete [] todo_list;
}

int SimLib::detectFragments(Simulation &sim, int *fragment_ids)
{
    if(fragment_ids == NULL || sim.number_of_particles < 1)
        return -1;

    int *todo_list = new int[sim.number_of_particles];
    int current_todo_pos = 0;                   // current pos in the todo list
    int next_free_todo_pos = 0;                 // next free pos in the todo list
    int current_agglomerate = -1;               // id of the current agglomerate
    int biggest_fragment = 0;
    int current_size = 0;
    int size_of_biggest_fragment = 0;
    ContactListEntry *cl_entry;

    // build contact list in which all contacts with other particles (not only those with lower id) are stored
    std::vector< std::list<int> > contacts(sim.number_of_particles);

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        fragment_ids[p] = -1;

        cl_entry = sim.contact_list[p];

        while(cl_entry)
        {
            // ignore contacts with walls
            if(cl_entry->id >= 0)
            {
                contacts[p].push_back(cl_entry->id);
                contacts[cl_entry->id].push_back(p);
            }

            cl_entry = cl_entry->next;
        }
    }

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        // particle does not belong to any agglomerate yet
        if(fragment_ids[p] == -1)
        {
            // put it on the todo list
            todo_list[0] = p;

            // create a new agglomerate
            ++current_agglomerate;
            current_size = 1;

            // the particle now belongs to the new agglomerate
            fragment_ids[p] = current_agglomerate;

            // set the pointers to the todo list
            current_todo_pos = 0;
            next_free_todo_pos = 1;

            // repeat until there is no particle left on the todo list
            while(current_todo_pos < next_free_todo_pos)
            {
                for(std::list<int>::iterator contact = contacts[todo_list[current_todo_pos]].begin(); contact != contacts[todo_list[current_todo_pos]].end(); ++contact)
                {
                    // check if particle has not been added to an agglomerate yet
                    if(fragment_ids[*contact] == -1)
                    {
                        // add particle to the current agglomerate
                        fragment_ids[*contact] = current_agglomerate;
                        ++current_size;

                        // add particle to the todo list
                        todo_list[next_free_todo_pos] = *contact;
                        ++next_free_todo_pos;
                    }
                }

                // work on this todo list entry is done
                ++current_todo_pos;
            }

            if(current_size > size_of_biggest_fragment)
            {
                biggest_fragment = current_agglomerate;
                size_of_biggest_fragment = current_size;
            }
        }
    }

    delete [] todo_list;

    return biggest_fragment;
}

void SimLib::filterFragments(Simulation *sim)
{
    if(sim->number_of_particles < 1)
        return;

    std::vector<int> fragment_ids;                          // array storing the fragment id of every particle
    std::vector<int> size_of_fragment;                      // number of particles of the fragment
    std::vector< std::list<int> > particles_of_fragment;    // ids of the particles of a specific fragment
    int largest_fragment = 0;                               // id of particles of the biggest fragment that has been detected so far

    SimLib::detectFragments(*sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

    // determine biggest agglomerate
    for(unsigned int agg = 1; agg < size_of_fragment.size(); ++agg)
    {
        if(size_of_fragment[agg] > size_of_fragment[largest_fragment])
            largest_fragment = (int)agg;
    }

    removeFragments(sim, largest_fragment, fragment_ids, size_of_fragment);
}

void SimLib::removeFragments(Simulation *sim, int remaining_fragment_id, std::vector<int> &fragment_ids, std::vector<int> &size_of_fragment)
{
    sim->removeWalls();

    // all particles but those of the biggest agglomerate will be removed
    int removed_particles = sim->number_of_particles - size_of_fragment[remaining_fragment_id];

    if(removed_particles > 0)
    {
        bool *remove_list = new bool[sim->number_of_particles];

        for(int p = 0; p < sim->number_of_particles; ++p)
        {
            // copy particle
            if(fragment_ids[p] == remaining_fragment_id)
                remove_list[p] = false;
            else
                remove_list[p] = true;
        }

        sim->removeParticles(remove_list, removed_particles);

        delete [] remove_list;
    }
}

void SimLib::sliceBox(Simulation *sim, vec3 &pos_lower, vec3 &pos_upper)
{
    sim->removeWalls();

    bool *remove_list = new bool[sim->number_of_particles];
    int removed_particles = 0;

    vec3 tmp_lower_pos;
    tmp_lower_pos[0] = pos_lower[0] + particle_radius;
    tmp_lower_pos[1] = pos_lower[1] + particle_radius;
    tmp_lower_pos[2] = pos_lower[2] + particle_radius;

    vec3 tmp_upper_pos;
    tmp_upper_pos[0] = pos_upper[0] - particle_radius;
    tmp_upper_pos[1] = pos_upper[1] - particle_radius;
    tmp_upper_pos[2] = pos_upper[2] - particle_radius;

    // determine which particles will be removed
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        if( sim->pos_old[X_COORD(p)] < tmp_lower_pos[0] || sim->pos_old[Y_COORD(p)] < tmp_lower_pos[1] || sim->pos_old[Z_COORD(p)] < tmp_lower_pos[2] ||
            sim->pos_old[X_COORD(p)] > tmp_upper_pos[0] || sim->pos_old[Y_COORD(p)] > tmp_upper_pos[1] || sim->pos_old[Z_COORD(p)] > tmp_upper_pos[2] )
        {
            remove_list[p] = true;
            ++removed_particles;
        }
        else
            remove_list[p] = false;
    }

    // remove particles from sim
    if(removed_particles > 0)
        sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;
}

void SimLib::sliceSphere(Simulation *sim, vec3 &center, double slice_radius)
{
    sim->removeWalls();

    if(sim->number_of_particles < 1)
        return;

    bool *remove_list = new bool[sim->number_of_particles];
    int removed_particles = 0;
    double dist;

    double slice_radius_sq = (slice_radius - particle_radius) * (slice_radius - particle_radius);

    // determine which particles will be removed
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        dist = (sim->pos_old[X_COORD(p)] - center[0])*(sim->pos_old[X_COORD(p)] - center[0])
            + (sim->pos_old[Y_COORD(p)] - center[1])*(sim->pos_old[Y_COORD(p)] - center[1])
            + (sim->pos_old[Z_COORD(p)] - center[2])*(sim->pos_old[Z_COORD(p)] - center[2]);

        if( dist > slice_radius_sq) //slice_radius*slice_radius )
        {
            remove_list[p] = true;
            ++removed_particles;
        }
        else
            remove_list[p] = false;
    }

    // remove particles from sim
    if(removed_particles > 0)
        sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;
}

void SimLib::sliceCylinder(Simulation *sim, vec3 &center, double slice_radius)
{
    sim->removeWalls();

    if(sim->number_of_particles < 1)
        return;

    bool *remove_list = new bool[sim->number_of_particles];
    int removed_particles = 0;
    double dist;

    double slice_radius_sq = (slice_radius - particle_radius) * (slice_radius - particle_radius);

    // determine which particles will be removed
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        dist = (sim->pos_old[X_COORD(p)] - center[0])*(sim->pos_old[X_COORD(p)] - center[0])
            + (sim->pos_old[Z_COORD(p)] - center[2])*(sim->pos_old[Z_COORD(p)] - center[2]);

        if( dist > slice_radius_sq) // slice_radius*slice_radius)
        {
            remove_list[p] = true;
            ++removed_particles;
        }
        else
            remove_list[p] = false;
    }

    // remove particles from sim
    if(removed_particles > 0)
        sim->removeParticles(remove_list, removed_particles);

    delete [] remove_list;
}

ErrorCode SimLib::sliceTop(Simulation *sim, double top_slice_factor)
{
    if(top_slice_factor < 0)
        return EC_INVALID_PARAMETER;

    // determine spatial extension
    vec3 lower_pos, upper_pos;
    sim->getEnclosingBox(&lower_pos, &upper_pos);

    upper_pos[1] -= (upper_pos[1] - lower_pos[1]) * top_slice_factor;
    SimLib::sliceBox(sim, lower_pos, upper_pos);

    return EC_OK;
}

ErrorCode SimLib::sliceBottom(Simulation *sim, double bottom_slice_factor)
{
    if(bottom_slice_factor < 0)
        return EC_INVALID_PARAMETER;

    // determine spatial extension
    vec3 lower_pos, upper_pos;
    sim->getEnclosingBox(&lower_pos, &upper_pos);

    lower_pos[1] += (upper_pos[1] - lower_pos[1]) * bottom_slice_factor;
    SimLib::sliceBox(sim, lower_pos, upper_pos);

    return EC_OK;
}

void SimLib::getCenterOfMass(vec3 *cms, const Simulation &sim, int first_particle, int last_particle)
{
    (*cms)[0] = 0;
    (*cms)[1] = 0;
    (*cms)[2] = 0;

    if(first_particle <= last_particle && last_particle < sim.number_of_particles)
    {
        for(int p = first_particle; p <= last_particle; ++p)
        {
            (*cms)[0] += sim.pos_old[X_COORD(p)];
            (*cms)[1] += sim.pos_old[Y_COORD(p)];
            (*cms)[2] += sim.pos_old[Z_COORD(p)];
        }

        double particles_inv = 1.0 / (double)(last_particle - first_particle + 1);

        (*cms)[0] *= particles_inv;
        (*cms)[1] *= particles_inv;
        (*cms)[2] *= particles_inv;
    }
}

void SimLib::rotateSimRandomly(Simulation *sim)
{
    // determine two random angles to get axis
    double phi = sim->get_random_zero_twoPi();
    double cos_theta = sim->get_random_cos_theta();
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

    vec3 axis;
    axis[0] = sin_theta * cos(phi);
    axis[1] = sin_theta * sin(phi);
    axis[2] = cos_theta;

    // random angle of rotation
    double angle = sim->get_random_zero_twoPi();

    // rotate around center of mass
    SimLib::rotateParticles(sim, axis, angle, 0, sim->number_of_particles-1);
}

void SimLib::rotateParticles(Simulation *sim, vec3 &axis, double angle, int first_particle, int last_particle)
{
    double cos_a = cos(angle);
    double sin_a = sin(angle);

    // rotation matrix
    double a11 = cos_a + axis[0] * axis[0] * (1.0 - cos_a);
    double a12 = axis[0] * axis[1] * (1.0 - cos_a) - axis[2] * sin_a;
    double a13 = axis[0] * axis[2] * (1.0 - cos_a) + axis[1] * sin_a;
    double a21 = axis[1] * axis[0] * (1.0 - cos_a) + axis[2] * sin_a;
    double a22 = cos_a + axis[1] * axis[1] * (1.0 - cos_a);
    double a23 = axis[1] * axis[2] * (1.0 - cos_a) - axis[0] * sin_a;
    double a31 = axis[2] * axis[0] * (1.0 - cos_a) - axis[1] * sin_a;
    double a32 = axis[2] * axis[1] * (1.0 - cos_a) + axis[0] * sin_a;
    double a33 = cos_a + axis[2] * axis[2] * (1.0 - cos_a);

    vec3 old_pos, temp, temp2, n1, n2;
    ContactListEntry *cl_entry;

    for(int p = last_particle; p >= first_particle; --p)
    {
        // rotate coordinates
        old_pos[0] = sim->pos_old[X_COORD(p)];
        old_pos[1] = sim->pos_old[Y_COORD(p)];
        old_pos[2] = sim->pos_old[Z_COORD(p)];

        sim->pos_old[X_COORD(p)] = a11 * old_pos[0] + a12 * old_pos[1] + a13 * old_pos[2];
        sim->pos_old[Y_COORD(p)] = a21 * old_pos[0] + a22 * old_pos[1] + a23 * old_pos[2];
        sim->pos_old[Z_COORD(p)] = a31 * old_pos[0] + a32 * old_pos[1] + a33 * old_pos[2];

        // rotate contact pointers
        cl_entry = sim->contact_list[p];

        while(cl_entry)
        {
            cl_entry->contact->getCurrentN1(&n1);
            cl_entry->contact->getCurrentN2(&n2);

            // rotate contact pointer of particle p
            temp[0] = old_pos[0] + n1[0];
            temp[1] = old_pos[1] + n1[1];
            temp[2] = old_pos[2] + n1[2];

            temp2[0] = a11 * temp[0] + a12 * temp[1] + a13 * temp[2];
            temp2[1] = a21 * temp[0] + a22 * temp[1] + a23 * temp[2];
            temp2[2] = a31 * temp[0] + a32 * temp[1] + a33 * temp[2];

            n1[0] = temp2[0] - sim->pos_old[X_COORD(p)];
            n1[1] = temp2[1] - sim->pos_old[Y_COORD(p)];
            n1[2] = temp2[2] - sim->pos_old[Z_COORD(p)];
            cl_entry->contact->updateN1Initial(n1);

            // rotate contact pointer of the other particle forming the contact
            temp[0] = sim->pos_old[X_COORD(cl_entry->id)] + n2[0];
            temp[1] = sim->pos_old[Y_COORD(cl_entry->id)] + n2[1];
            temp[2] = sim->pos_old[Z_COORD(cl_entry->id)] + n2[2];

            temp2[0] = a11 * temp[0] + a12 * temp[1] + a13 * temp[2];
            temp2[1] = a21 * temp[0] + a22 * temp[1] + a23 * temp[2];
            temp2[2] = a31 * temp[0] + a32 * temp[1] + a33 * temp[2];

            // calculate rotated pos of other particle (has not been done yet)
            temp[0] = a11 * sim->pos_old[X_COORD(cl_entry->id)] + a12 * sim->pos_old[Y_COORD(cl_entry->id)] + a13 * sim->pos_old[Z_COORD(cl_entry->id)];
            temp[1] = a21 * sim->pos_old[X_COORD(cl_entry->id)] + a22 * sim->pos_old[Y_COORD(cl_entry->id)] + a23 * sim->pos_old[Z_COORD(cl_entry->id)];
            temp[2] = a31 * sim->pos_old[X_COORD(cl_entry->id)] + a32 * sim->pos_old[Y_COORD(cl_entry->id)] + a33 * sim->pos_old[Z_COORD(cl_entry->id)];

            n2[0] = temp2[0] - temp[0];
            n2[1] = temp2[1] - temp[1];
            n2[2] = temp2[2] - temp[2];
            cl_entry->contact->updateN2Initial(n2);

            cl_entry = cl_entry->next;
        }
    }

    sim->grid.resetGrid();
    sim->grid.addParticles(sim->pos_old, sim->number_of_particles);
}

void SimLib::rotateParticles(Simulation *sim, double angle1, double angle2, double angle3, int first_particle, int last_particle)
{
    double cos1 = cos(angle1);
    double sin1 = sin(angle1);
    double cos2 = cos(angle2);
    double sin2 = sin(angle2);
    double cos3 = cos(angle3);
    double sin3 = sin(angle3);

    // rotation matrix
    double a11 = - sin1 * cos3 + cos1 * cos2 * cos3;
    double a12 = - cos1 * sin3 - sin1 * cos2 * cos3;
    double a13 = sin1 * cos2;
    double a21 = sin1 * cos3 + cos1 * cos2 * sin3;
    double a22 = cos1 * cos2 * cos3 - sin1 * sin3;
    double a23 = - cos1 * sin2;
    double a31 = sin2 * sin3;
    double a32 = sin2 * cos3;
    double a33 = cos2;

    vec3 old_pos;

    for(int p = first_particle; p <= last_particle; ++p)
    {
        old_pos[0] = sim->pos_old[X_COORD(p)];
        old_pos[1] = sim->pos_old[Y_COORD(p)];
        old_pos[2] = sim->pos_old[Z_COORD(p)];

        sim->pos_old[X_COORD(p)] = a11 * old_pos[0] + a12 * old_pos[1] + a13 * old_pos[2];
        sim->pos_old[Y_COORD(p)] = a21 * old_pos[0] + a22 * old_pos[1] + a23 * old_pos[2];
        sim->pos_old[Z_COORD(p)] = a31 * old_pos[0] + a32 * old_pos[1] + a33 * old_pos[2];
    }
}

void SimLib::centerCMS(Simulation *sim)
{
    vec3 center_of_mass;

    SimLib::getCenterOfMass(&center_of_mass, *sim, 0, sim->number_of_particles-1);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        sim->pos_old[X_COORD(p)] -= center_of_mass[0];
        sim->pos_old[Y_COORD(p)] -= center_of_mass[1];
        sim->pos_old[Z_COORD(p)] -= center_of_mass[2];
    }

    sim->grid.resetGrid();
    sim->grid.addParticles(sim->pos_old, sim->number_of_particles);
}

void SimLib::getSize(const Simulation &sim, double *gyration_radius, double *outer_radius)
{
    if(sim.number_of_particles == 0)
    {
        if(outer_radius)
            *outer_radius = 0.0;

        if(gyration_radius)
            *gyration_radius = 0.0;

        return;
    }
    else if(sim.number_of_particles == 1)
    {
        if(gyration_radius)
            *gyration_radius = 0.0;

        if(outer_radius)
            *outer_radius = particle_radius;

        return;
    }

    // determine CMS
    vec3 cms, delta_pos;
    SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);

    // calc gyration radius
    double result = 0.0;
    double dist;
    double max_dist = 0.0;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        delta_pos[0] = cms[0] - sim.pos_old[X_COORD(p)];
        delta_pos[1] = cms[1] - sim.pos_old[Y_COORD(p)];
        delta_pos[2] = cms[2] - sim.pos_old[Z_COORD(p)];
        dist = norm_squared(delta_pos);

        if(dist > max_dist)
            max_dist = dist;

        result += dist;
    }

    result /= (double)sim.number_of_particles;

    if(gyration_radius)
        *gyration_radius = sqrt(result);

    if(outer_radius)
        *outer_radius = sqrt(max_dist);
}




/*
void SimLib::accelerate_inward(const Simulation &sim, double speed)
{
    if(sim.number_of_particles == 0)
        return;

    // determine CMS
    vec3 cms, delta_pos;
    SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);

    // calc gyration radius
    double result = 0.0;
    double dist;
    double max_dist = 0.0;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        delta_pos[0] = cms[0] - sim.pos_old[X_COORD(p)];
        delta_pos[1] = cms[1] - sim.pos_old[Y_COORD(p)];
        delta_pos[2] = cms[2] - sim.pos_old[Z_COORD(p)];
        dist = norm_squared(delta_pos);

        if(dist > max_dist)
            max_dist = dist;

        result += dist;
    }




    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        delta_pos[0] = cms[0] - sim.pos_old[X_COORD(p)];
        delta_pos[1] = cms[1] - sim.pos_old[Y_COORD(p)];
        delta_pos[2] = cms[2] - sim.pos_old[Z_COORD(p)];
        dist = norm(delta_pos);

        double vel = speed/max_dist;

        sim.vel[X_COORD(p)] += vel * delta_pos[0];
        sim.vel[Y_COORD(p)] += vel * delta_pos[1];
        sim.vel[Z_COORD(p)] += vel * delta_pos[2];
    }





}
*/





void SimLib::getCrossSection(Simulation &sim, double cell_size, unsigned int rotations, double *cross_section, double *sigma_cross_section, double radius_cutoff_factor)
{
    if(rotations == 0)
        return;

    SimLib::centerCMS(&sim);

    // determine max size of agglomerate
    double outer_radius;
    double dist_squared;
    getSize(sim, NULL, &outer_radius);

    double max_dist = radius_cutoff_factor * radius_cutoff_factor * outer_radius * outer_radius;

    double x_min = -(outer_radius + 2.0 * particle_radius);
    double x_max = outer_radius + 2.0 *particle_radius;
    double z_min = -(outer_radius + 2.0 * particle_radius);
    double z_max = outer_radius + 2.0 * particle_radius;

    // allocate sufficient number of cells
    int x_cells = (int)( (x_max - x_min) / cell_size );
    int z_cells = (int)( (z_max - z_min) / cell_size );
    char *cells = new char[x_cells*z_cells];

    // range of cells that have to be checked when a particle is added (assuming that a particle covers more than 1 cell)
    int check_cells = (int)(particle_radius / cell_size)+1;

    double pos_x, pos_z;
    int x_id, z_id, cell_count;
    double my_cross_section = 0.0;
    double my_sigma_cross_section = 0.0;

    // get positions of particles
    double *positions = new double[3*sim.number_of_particles];
    memcpy(positions, sim.pos_old, 3 * sim.number_of_particles * sizeof(double));

    // rotation axis
    vec3 axis;
    axis[0] = 0.0;
    axis[1] = 1.0;
    axis[2] = 0.0;

    // calculate rotation matrix
    double cos_a = cos(M_PI / (double)rotations);
    double sin_a = sin(M_PI / (double)rotations);
    double a11 = cos_a + axis[0] * axis[0] * (1.0 - cos_a);
    double a12 = axis[0] * axis[1] * (1.0 - cos_a) - axis[2] * sin_a;
    double a13 = axis[0] * axis[2] * (1.0 - cos_a) + axis[1] * sin_a;
    double a21 = axis[1] * axis[0] * (1.0 - cos_a) + axis[2] * sin_a;
    double a22 = cos_a + axis[1] * axis[1] * (1.0 - cos_a);
    double a23 = axis[1] * axis[2] * (1.0 - cos_a) - axis[0] * sin_a;
    double a31 = axis[2] * axis[0] * (1.0 - cos_a) - axis[1] * sin_a;
    double a32 = axis[2] * axis[1] * (1.0 - cos_a) + axis[0] * sin_a;
    double a33 = cos_a + axis[2] * axis[2] * (1.0 - cos_a);

    for(unsigned int i = 0; i < rotations; ++i)
    {
        // init empty cells
        memset(cells, 0, x_cells*z_cells * sizeof(char));
        cell_count = 0;

        //////////////////////////////////////////////////////////////////////////////////////
        // apply rotation
        //////////////////////////////////////////////////////////////////////////////////////

        vec3 old_pos;

        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            // rotate coordinates
            old_pos[0] = positions[X_COORD(p)];
            old_pos[1] = positions[Y_COORD(p)];
            old_pos[2] = positions[Z_COORD(p)];

            positions[X_COORD(p)] = a11 * old_pos[0] + a12 * old_pos[1] + a13 * old_pos[2];
            positions[Y_COORD(p)] = a21 * old_pos[0] + a22 * old_pos[1] + a23 * old_pos[2];
            positions[Z_COORD(p)] = a31 * old_pos[0] + a32 * old_pos[1] + a33 * old_pos[2];

        }
/* Original by Alex
        // determine cross section
        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            // determine distance to CMS
            double dist = positions[X_COORD(p)]*positions[X_COORD(p)] + positions[Y_COORD(p)]*positions[Y_COORD(p)] + positions[Z_COORD(p)]*positions[Z_COORD(p)];

            // only take paticles into account that are within desired dist to cms
            if(dist < max_dist)
            {
                // determine id of the cell where the center of the particle is located
                x_id = (int)( (positions[X_COORD(p)] - x_min) / cell_size );
                z_id = (int)( (positions[Z_COORD(p)] - z_min) / cell_size );

                // check surrounding cells
                pos_z = z_min + (double)(z_id - check_cells) * cell_size + 0.5 * cell_size;

                for(int z = z_id - check_cells; z < z_id + check_cells; ++z)
                {
                    pos_x = x_min + (double)(x_id - check_cells) * cell_size + 0.5 * cell_size;

                    for(int x = x_id - check_cells; x < x_id + check_cells; ++x)
                    {
                        // cell is considered to be covered by the particle if its center lies within the particle radius
                        dist_squared = (pos_x - positions[X_COORD(p)])*(pos_x - positions[X_COORD(p)]) + (pos_z - positions[Z_COORD(p)])*(pos_z - positions[Z_COORD(p)]);

                        if(dist_squared < particle_radius*particle_radius)
                        {
                            // only count cells that have not been marked as occupied yet
                            if(cells[x + x_cells * z] == 0)
                            {
                                ++cell_count;
                                cells[x + x_cells * z] = 1;
                            }
                        }

                        pos_x += cell_size;
                    }

                    pos_z += cell_size;
                }
            }
        }
*/


        // determine cross section
        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            // determine id of the cell where the center of the particle is located
            x_id = (int)( (sim.pos_old[X_COORD(p)] - x_min) / cell_size );
            z_id = (int)( (sim.pos_old[Z_COORD(p)] - z_min) / cell_size );

            // check surrounding cells
            pos_z = z_min + (double)(z_id - check_cells) * cell_size + 0.5 * cell_size;

            for(int z = z_id - check_cells; z < z_id + check_cells; ++z)
            {
                pos_x = x_min + (double)(x_id - check_cells) * cell_size + 0.5 * cell_size;

                for(int x = x_id - check_cells; x < x_id + check_cells; ++x)
                {
                    // cell is considered to be covered by the particle if its center lies within the particle radius
                    dist_squared = (pos_x - sim.pos_old[X_COORD(p)])*(pos_x - sim.pos_old[X_COORD(p)]) + (pos_z - sim.pos_old[Z_COORD(p)])*(pos_z - sim.pos_old[Z_COORD(p)]);

                    if(dist_squared < (particle_radius)*(particle_radius)) // +1.27324*cell_size is correction factor from testing
                    {
                        // only count cells that have not been marked as occupied yet
                        if(cells[x + x_cells * z] == 0)
                        {
                            ++cell_count;
                            cells[x + x_cells * z] = 1;
                        }
                    }

                    pos_x += cell_size;
                }

                pos_z += cell_size;
            }
        }


        my_cross_section += (double)cell_count * cell_size*cell_size;
        my_sigma_cross_section += (double)cell_count * cell_size*cell_size * (double)cell_count * cell_size*cell_size;
    }

    delete [] cells;
    delete [] positions;

    if(cross_section)
        *cross_section = my_cross_section / (double)rotations;

    if(sigma_cross_section)
    {
        if(sim.number_of_particles > 1)
        {
            my_sigma_cross_section /= (double)rotations;
            *sigma_cross_section = sqrt(fabs(my_sigma_cross_section - my_cross_section * my_cross_section / (double)(rotations * rotations)));
        }
        else
            *sigma_cross_section = 0.0;
    }
}




void SimLib::getCrossSectionNoRotation(const Simulation &sim, const double cell_size, double &cross_section)
{

    vec3 cms;
    SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);

    // determine max size of agglomerate
    double outer_radius;
    double dist_squared;
    SimLib::getSize(sim, NULL, &outer_radius);


    const double x_min = cms[0] - (outer_radius + 2.0 * particle_radius);
    const double x_max = cms[0] + outer_radius + 2.0 *particle_radius;
    const double z_min = cms[2] - (outer_radius + 2.0 * particle_radius);
    const double z_max = cms[2] + outer_radius + 2.0 * particle_radius;


    // allocate sufficient number of cells
    const int x_cells = (int)( (x_max - x_min) / cell_size );
    const int z_cells = (int)( (z_max - z_min) / cell_size );
    char *cells = new char[x_cells*z_cells];

    // range of cells that have to be checked when a particle is added (assuming that a particle covers more than 1 cell)
    int check_cells = (int)(particle_radius / cell_size)+1;

    double pos_x, pos_z;
    int x_id, z_id, cell_count;

    // init empty cells
    memset(cells, 0, x_cells*z_cells * sizeof(char));
    cell_count = 0;


    // determine cross section
    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        // determine id of the cell where the center of the particle is located
        x_id = (int)( (sim.pos_old[X_COORD(p)] - x_min) / cell_size );
        z_id = (int)( (sim.pos_old[Z_COORD(p)] - z_min) / cell_size );

        // check surrounding cells
        pos_z = z_min + (double)(z_id - check_cells) * cell_size + 0.5 * cell_size;

        for(int z = z_id - check_cells; z < z_id + check_cells; ++z)
        {
            pos_x = x_min + (double)(x_id - check_cells) * cell_size + 0.5 * cell_size;

            for(int x = x_id - check_cells; x < x_id + check_cells; ++x)
            {
                // cell is considered to be covered by the particle if its center lies within the particle radius
                dist_squared = (pos_x - sim.pos_old[X_COORD(p)])*(pos_x - sim.pos_old[X_COORD(p)]) + (pos_z - sim.pos_old[Z_COORD(p)])*(pos_z - sim.pos_old[Z_COORD(p)]);

                if(dist_squared < (particle_radius+1.27324*cell_size)*(particle_radius)) // +1.27324*cell_size is correction factor from testing
                {
                    // only count cells that have not been marked as occupied yet
                    if(cells[x + x_cells * z] == 0)
                    {
                        ++cell_count;
                        cells[x + x_cells * z] = 1;
                    }
                }

                pos_x += cell_size;
            }

            pos_z += cell_size;
        }
    }

    cross_section = (double)cell_count * cell_size*cell_size;


    delete [] cells;

    return;

}








bool SimLib::inContactWithWall(Simulation *sim, int particle_id)
{
    ContactListEntry *cl_entry = sim->contact_list[particle_id];

    while(cl_entry)
    {
        if(cl_entry->id < 0)
            return true;
        else
            cl_entry = cl_entry->next;
    }

    return false;
}

void SimLib::printContactList(Simulation &sim, const char* filename)
{
    ContactListEntry *cl_entry;

    // build contact list
    /*std::vector< std::list<int> > contacts(sim->number_of_particles);

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        cl_entry = sim->contact_list[p];

        while(cl_entry)
        {
            contacts[p].push_back(cl_entry->id);
            contacts[cl_entry->id].push_back(p);

            cl_entry = cl_entry->next;
        }
    }*/

    FILE *file = fopen(filename, "w+");

    if(file)
    {
        /*for(int p = 0; p < sim->number_of_particles; ++p)
        {
            for(std::list<int>::iterator contact = contacts[p].begin(); contact != contacts[p].end(); ++contact)
                fprintf(file, "%i %i\n", p, *contact, );
        }*/

        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            cl_entry = sim.contact_list[p];

            while(cl_entry)
            {
                //fprintf(file, "%i %i %lf %lf %lf %lf\n", p, cl_entry->id, cl_entry->contact->rot1.e0, cl_entry->contact->rot1.e1, cl_entry->contact->rot1.e2, cl_entry->contact->rot1.e3);
                fprintf(file, "%i %i\n", p, cl_entry->id);
                cl_entry = cl_entry->next;
            }
        }

        fclose(file);
    }
}

void SimLib::printEnergy(Simulation* sim, const char *filename)
{

    vec3 temp_vec;
    vec3 temp_vec2;
    vec3 n_c;
    vec3 n1;
    vec3 n2;
    double dist;

    double E_kin = 0.0;
    double E_rot = 0.0;

    double V_normal = 0.0;
    double V_roll = 0.0;
    double V_slide = 0.0;
    double V_twist = 0.0;

    // calculate kinetic + rotational energy
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        E_kin += 0.5 * mass * (sim->vel[X_COORD(p)] * sim->vel[X_COORD(p)] + sim->vel[Y_COORD(p)] * sim->vel[Y_COORD(p)] + sim->vel[Z_COORD(p)] * sim->vel[Z_COORD(p)]);
        E_rot += 0.5 * moment_of_inertia * (sim->vel_angular[X_COORD(p)] * sim->vel_angular[X_COORD(p)] + sim->vel_angular[Y_COORD(p)] * sim->vel_angular[Y_COORD(p)] + sim->vel_angular[Z_COORD(p)] * sim->vel_angular[Z_COORD(p)]);
    }

    E_kin *= ENERGY_UNIT;
    E_rot *= ENERGY_UNIT;

    //
    if(sim->sim_info.sim_type == SIM_TYPE_DYNAMIC_COMPRESSION)
    {
        for(int w = 0; w < sim->number_of_walls; ++w)
            E_kin += 0.5 * sim->walls[w].mass * norm_squared(sim->walls[w].velocity) * ENERGY_UNIT;
    }

    ContactListEntry *cl_entry;

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        cl_entry = sim->contact_list[p];

        while(cl_entry)
        {
            if(cl_entry->id < 0)
            {
                // calculate current contact pointers
                cl_entry->contact->getCurrentN1(&n1);
                memcpy(n2, sim->walls[-cl_entry->id-1].normal, sizeof(vec3));

                temp_vec[0] = sim->pos_old[X_COORD(p)];
                temp_vec[1] = sim->pos_old[Y_COORD(p)];
                temp_vec[2] = sim->pos_old[Z_COORD(p)];
                //V_normal += sim->normal_interaction.getJKRWallPotentialEnergy(2.0 * particle_radius - sim->walls[-cl_entry->id-1].getDistanceTo(temp_vec));

#ifdef WALL_ROLLING
                temp_vec[0] = particle_radius * (n1[0] + n2[0]);
                temp_vec[1] = particle_radius * (n1[1] + n2[1]);
                temp_vec[2] = particle_radius * (n1[2] + n2[2]);
                V_roll += 0.5 * wall_k_r * norm_squared(temp_vec);

#endif

#ifdef WALL_SLIDING
                //temp_vec[0] = sim->pos_old[X_COORD(p)] + particle_radius * n1[0] - cl_entry->contact->n2_initial[0] - particle_radius * n2[0];
                //temp_vec[1] = sim->pos_old[Y_COORD(p)] + particle_radius * n1[1] - cl_entry->contact->n2_initial[1] - particle_radius * n2[1];
                //temp_vec[2] = sim->pos_old[Z_COORD(p)] + particle_radius * n1[2] - cl_entry->contact->n2_initial[2] - particle_radius * n2[2];



                vec3 wall_contact_pos;
                // n2_initial is the position where the contact was made minus the position of the wall at that time
                wall_contact_pos[0] = cl_entry->contact->n2_initial[0] + sim->walls[-cl_entry->id-1].pos[0];
                wall_contact_pos[1] = cl_entry->contact->n2_initial[1] + sim->walls[-cl_entry->id-1].pos[1];
                wall_contact_pos[2] = cl_entry->contact->n2_initial[2] + sim->walls[-cl_entry->id-1].pos[2];


                temp_vec[0] = sim->pos_old[X_COORD(p)] - wall_contact_pos[0] + particle_radius * n1[0];
                temp_vec[1] = sim->pos_old[Y_COORD(p)] - wall_contact_pos[1] + particle_radius * n1[1];
                temp_vec[2] = sim->pos_old[Z_COORD(p)] - wall_contact_pos[2] + particle_radius * n1[2];



                double temp = dot_product(temp_vec, n2);
                temp_vec2[0] = temp_vec[0] - temp * n2[0];
                temp_vec2[1] = temp_vec[1] - temp * n2[1];
                temp_vec2[2] = temp_vec[2] - temp * n2[2];

                V_slide += 0.5 * wall_k_s * norm_squared(temp_vec2);

#endif

#ifdef WALL_TWISTING
                V_twist += 0.5 * wall_k_t * cl_entry->contact->twisting_displacement * cl_entry->contact->twisting_displacement;
#endif
            }
            else
            {
                // calculate current contact pointers
                cl_entry->contact->getCurrentN1(&n1);
                cl_entry->contact->getCurrentN2(&n2);

                // determine distance between particles
                n_c[0] = sim->pos_old[X_COORD(p)] - sim->pos_old[X_COORD(cl_entry->id)];
                n_c[1] = sim->pos_old[Y_COORD(p)] - sim->pos_old[Y_COORD(cl_entry->id)];
                n_c[2] = sim->pos_old[Z_COORD(p)] - sim->pos_old[Z_COORD(cl_entry->id)];

                dist = norm(n_c);
                n_c[0] /= dist;
                n_c[1] /= dist;
                n_c[2] /= dist;

                V_normal += sim->normal_interaction.getJKRPotentialEnergy(2.0 * particle_radius - dist);
#ifdef ROLLING
                temp_vec[0] = reduced_radius * (n1[0] + n2[0]);
                temp_vec[1] = reduced_radius * (n1[1] + n2[1]);
                temp_vec[2] = reduced_radius * (n1[2] + n2[2]);
                V_roll += 0.5 * k_r * norm_squared(temp_vec);
#endif

#ifdef SLIDING

                temp_vec[0] = (n1[0] - n2[0] + 2.0 * n_c[0]);
                temp_vec[1] = (n1[1] - n2[1] + 2.0 * n_c[1]);
                temp_vec[2] = (n1[2] - n2[2] + 2.0 * n_c[2]);

                double temp = dot_product(temp_vec, n_c);

                temp_vec2[0] = temp_vec[0] - temp * n_c[0];
                temp_vec2[1] = temp_vec[1] - temp * n_c[1];
                temp_vec2[2] = temp_vec[2] - temp * n_c[2];
                V_slide += 0.5 * k_s * particle_radius * particle_radius * norm_squared(temp_vec2);

#endif

#ifdef TWISTING
                V_twist += 0.5 * k_t * cl_entry->contact->twisting_displacement * cl_entry->contact->twisting_displacement;
#endif
            }

            cl_entry = cl_entry->next;
        }
    }

    V_normal *= ENERGY_UNIT;
    V_roll *= ENERGY_UNIT;
    V_slide *= ENERGY_UNIT;
    V_twist *= ENERGY_UNIT;



    // sum up
    double V_tot = V_normal + V_roll + V_slide + V_twist;
    double E_diss = sim->dissipated_contact_energy + sim->dissipated_rolling_energy + sim->dissipated_sliding_energy + sim->dissipated_twisting_energy + sim->dissipated_wall_energy + sim->dissipated_damping_energy;
    double E_tot = E_kin + E_rot + V_tot + E_diss;



    if(!filename)
    {

        //printf("#V_tot, V_normal, V_roll, V_slide, V_twist\n");
        printf("%.2g %.2g %.2g %.2g %.2g %.2g\n",
                V_tot,
                V_normal,
                V_roll,
                V_slide,
                V_twist,
                sim->dissipated_damping_energy
                );

        return;
    }

    else
    {
        FILE *file = fopen(filename, "a+");

        if(file)
        {

            fprintf(file, "%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g\n",
                    sim->current_time,
                    E_tot,
                    E_kin,
                    E_rot,
                    V_tot,
                    V_normal,
                    V_roll,
                    V_slide,
                    V_twist,
                    E_diss,
                    sim->dissipated_contact_energy,
                    sim->dissipated_rolling_energy,
                    sim->dissipated_sliding_energy,
                    sim->dissipated_twisting_energy,
                    sim->dissipated_wall_energy,
                    sim->dissipated_damping_energy
                    );


            fclose(file);
            return;
        }
    }
}


/*
double SimLib::getTotalEnergy(Simulation* sim)
{
    vec3 temp_vec;
    vec3 temp_vec2;
    vec3 n_c;
    vec3 n1;
    vec3 n2;
    double dist;

    double E_kin = 0.0;
    double E_rot = 0.0;

    double V_normal = 0.0;
    double V_roll = 0.0;
    double V_slide = 0.0;
    double V_twist = 0.0;

    // calculate kinetic + rotational energy
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        E_kin += 0.5 * mass * (sim->vel[X_COORD(p)] * sim->vel[X_COORD(p)] + sim->vel[Y_COORD(p)] * sim->vel[Y_COORD(p)] + sim->vel[Z_COORD(p)] * sim->vel[Z_COORD(p)]);
        E_rot += 0.5 * moment_of_inertia * (sim->vel_angular[X_COORD(p)] * sim->vel_angular[X_COORD(p)] + sim->vel_angular[Y_COORD(p)] * sim->vel_angular[Y_COORD(p)] + sim->vel_angular[Z_COORD(p)] * sim->vel_angular[Z_COORD(p)]);
    }

    E_kin *= ENERGY_UNIT;
    E_rot *= ENERGY_UNIT;

    //
    if(sim->sim_info.sim_type == SIM_TYPE_DYNAMIC_COMPRESSION)
    {
        for(int w = 0; w < sim->number_of_walls; ++w)
            E_kin += 0.5 * sim->walls[w].mass * norm_squared(sim->walls[w].velocity) * ENERGY_UNIT;
    }

    ContactListEntry *cl_entry;

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        cl_entry = sim->contact_list[p];

        while(cl_entry)
        {
            if(cl_entry->id < 0)
            {
                // calculate current contact pointers
                cl_entry->contact->getCurrentN1(&n1);
                memcpy(n2, sim->walls[-cl_entry->id-1].normal, sizeof(vec3));

                temp_vec[0] = sim->pos_old[X_COORD(p)];
                temp_vec[1] = sim->pos_old[Y_COORD(p)];
                temp_vec[2] = sim->pos_old[Z_COORD(p)];
                V_normal += sim->normal_interaction.getJKRWallPotentialEnergy(2.0 * particle_radius - sim->walls[-cl_entry->id-1].getDistanceTo(temp_vec));

#ifdef WALL_ROLLING
                temp_vec[0] = particle_radius * (n1[0] + n2[0]);
                temp_vec[1] = particle_radius * (n1[1] + n2[1]);
                temp_vec[2] = particle_radius * (n1[2] + n2[2]);
                V_roll += 0.5 * k_r * norm_squared(temp_vec);
#endif

#ifdef WALL_SLIDING

                vec3 wall_contact_pos;
                // n2_initial is the position where the contact was made minus the position of the wall at that time
                wall_contact_pos[0] = cl_entry->contact->n2_initial[0] + sim->walls[-cl_entry->id-1].pos[0];
                wall_contact_pos[1] = cl_entry->contact->n2_initial[1] + sim->walls[-cl_entry->id-1].pos[1];
                wall_contact_pos[2] = cl_entry->contact->n2_initial[2] + sim->walls[-cl_entry->id-1].pos[2];


                temp_vec[0] = sim->pos_old[X_COORD(p)] - wall_contact_pos[0] + particle_radius * n1[0];
                temp_vec[1] = sim->pos_old[Y_COORD(p)] - wall_contact_pos[1] + particle_radius * n1[1];
                temp_vec[2] = sim->pos_old[Z_COORD(p)] - wall_contact_pos[2] + particle_radius * n1[2];

                double temp = dot_product(temp_vec, n2);
                temp_vec2[0] = temp_vec[0] - temp * n2[0];
                temp_vec2[1] = temp_vec[1] - temp * n2[1];
                temp_vec2[2] = temp_vec[2] - temp * n2[2];

                V_slide += 0.5 * k_s * norm_squared(temp_vec2);
#endif

#ifdef WALL_TWISTING
                V_twist += 0.5 * k_t * cl_entry->contact->twisting_displacement * cl_entry->contact->twisting_displacement;
#endif
            }
            else
            {
                // calculate current contact pointers
                cl_entry->contact->getCurrentN1(&n1);
                cl_entry->contact->getCurrentN2(&n2);

                // determine distance between particles
                n_c[0] = sim->pos_old[X_COORD(p)] - sim->pos_old[X_COORD(cl_entry->id)];
                n_c[1] = sim->pos_old[Y_COORD(p)] - sim->pos_old[Y_COORD(cl_entry->id)];
                n_c[2] = sim->pos_old[Z_COORD(p)] - sim->pos_old[Z_COORD(cl_entry->id)];

                dist = norm(n_c);
                n_c[0] /= dist;
                n_c[1] /= dist;
                n_c[2] /= dist;

                V_normal += sim->normal_interaction.getJKRPotentialEnergy(2.0 * particle_radius - dist);

#ifdef ROLLING
                temp_vec[0] = reduced_radius * (n1[0] + n2[0]);
                temp_vec[1] = reduced_radius * (n1[1] + n2[1]);
                temp_vec[2] = reduced_radius * (n1[2] + n2[2]);
                V_roll += 0.5 * k_r * norm_squared(temp_vec);
#endif

#ifdef SLIDING
                temp_vec[0] = particle_radius * (n1[0] - n2[0] + 2.0 * n_c[0]);
                temp_vec[1] = particle_radius * (n1[1] - n2[1] + 2.0 * n_c[1]);
                temp_vec[2] = particle_radius * (n1[2] - n2[2] + 2.0 * n_c[2]);

                double temp = dot_product(temp_vec, n_c);
                temp_vec2[0] = temp_vec[0] - temp * n_c[0];
                temp_vec2[1] = temp_vec[1] - temp * n_c[1];
                temp_vec2[2] = temp_vec[2] - temp * n_c[2];
                V_slide += 0.5 * k_s * norm_squared(temp_vec2);
#endif

#ifdef TWISTING
                V_twist += 0.5 * k_t * cl_entry->contact->twisting_displacement * cl_entry->contact->twisting_displacement;
#endif
            }

            cl_entry = cl_entry->next;
        }
    }

    V_normal *= ENERGY_UNIT;
    V_roll *= ENERGY_UNIT;
    V_slide *= ENERGY_UNIT;
    V_twist *= ENERGY_UNIT;

    // sum up
    double V_tot = V_normal + V_roll + V_slide + V_twist;
    double E_diss = sim->dissipated_contact_energy + sim->dissipated_rolling_energy + sim->dissipated_sliding_energy + sim->dissipated_twisting_energy + sim->dissipated_wall_energy + sim->dissipated_damping_energy;
    double E_tot = E_kin + E_rot + V_tot + E_diss;

    //printf("E_tot = %e  %e  %e  %e  %e\n", E_tot, E_kin, E_rot, V_tot, E_diss);

    if(E_kin + E_rot < E_tot*1.e-14)
    {
        //printf("Stopping simulation because all energy is dissipated!\n");
        sim->stop_simulation = true;
    }
    return E_tot;
}
*/



double SimLib::getTotalEnergy(Simulation* sim)
{

    vec3 temp_vec;
    vec3 temp_vec2;
    vec3 n_c;
    vec3 n1;
    vec3 n2;
    double dist;

    double E_kin = 0.0;
    double E_rot = 0.0;

    double V_normal = 0.0;
    double V_roll = 0.0;
    double V_slide = 0.0;
    double V_twist = 0.0;

    // calculate kinetic + rotational energy
    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        E_kin += 0.5 * mass * (sim->vel[X_COORD(p)] * sim->vel[X_COORD(p)] + sim->vel[Y_COORD(p)] * sim->vel[Y_COORD(p)] + sim->vel[Z_COORD(p)] * sim->vel[Z_COORD(p)]);
        E_rot += 0.5 * moment_of_inertia * (sim->vel_angular[X_COORD(p)] * sim->vel_angular[X_COORD(p)] + sim->vel_angular[Y_COORD(p)] * sim->vel_angular[Y_COORD(p)] + sim->vel_angular[Z_COORD(p)] * sim->vel_angular[Z_COORD(p)]);
    }

    E_kin *= ENERGY_UNIT;
    E_rot *= ENERGY_UNIT;

    //
    if(sim->sim_info.sim_type == SIM_TYPE_DYNAMIC_COMPRESSION)
    {
        for(int w = 0; w < sim->number_of_walls; ++w)
            E_kin += 0.5 * sim->walls[w].mass * norm_squared(sim->walls[w].velocity) * ENERGY_UNIT;
    }

    ContactListEntry *cl_entry;

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        cl_entry = sim->contact_list[p];

        while(cl_entry)
        {
            if(cl_entry->id < 0)
            {
                // calculate current contact pointers
                cl_entry->contact->getCurrentN1(&n1);
                memcpy(n2, sim->walls[-cl_entry->id-1].normal, sizeof(vec3));

                temp_vec[0] = sim->pos_old[X_COORD(p)];
                temp_vec[1] = sim->pos_old[Y_COORD(p)];
                temp_vec[2] = sim->pos_old[Z_COORD(p)];
                //V_normal += sim->normal_interaction.getJKRWallPotentialEnergy(2.0 * particle_radius - sim->walls[-cl_entry->id-1].getDistanceTo(temp_vec));

#ifdef WALL_ROLLING
                temp_vec[0] = particle_radius * (n1[0] + n2[0]);
                temp_vec[1] = particle_radius * (n1[1] + n2[1]);
                temp_vec[2] = particle_radius * (n1[2] + n2[2]);
                V_roll += 0.5 * wall_k_r * norm_squared(temp_vec);

#endif

#ifdef WALL_SLIDING
                //temp_vec[0] = sim->pos_old[X_COORD(p)] + particle_radius * n1[0] - cl_entry->contact->n2_initial[0] - particle_radius * n2[0];
                //temp_vec[1] = sim->pos_old[Y_COORD(p)] + particle_radius * n1[1] - cl_entry->contact->n2_initial[1] - particle_radius * n2[1];
                //temp_vec[2] = sim->pos_old[Z_COORD(p)] + particle_radius * n1[2] - cl_entry->contact->n2_initial[2] - particle_radius * n2[2];



                vec3 wall_contact_pos;
                // n2_initial is the position where the contact was made minus the position of the wall at that time
                wall_contact_pos[0] = cl_entry->contact->n2_initial[0] + sim->walls[-cl_entry->id-1].pos[0];
                wall_contact_pos[1] = cl_entry->contact->n2_initial[1] + sim->walls[-cl_entry->id-1].pos[1];
                wall_contact_pos[2] = cl_entry->contact->n2_initial[2] + sim->walls[-cl_entry->id-1].pos[2];


                temp_vec[0] = sim->pos_old[X_COORD(p)] - wall_contact_pos[0] + particle_radius * n1[0];
                temp_vec[1] = sim->pos_old[Y_COORD(p)] - wall_contact_pos[1] + particle_radius * n1[1];
                temp_vec[2] = sim->pos_old[Z_COORD(p)] - wall_contact_pos[2] + particle_radius * n1[2];



                double temp = dot_product(temp_vec, n2);
                temp_vec2[0] = temp_vec[0] - temp * n2[0];
                temp_vec2[1] = temp_vec[1] - temp * n2[1];
                temp_vec2[2] = temp_vec[2] - temp * n2[2];

                V_slide += 0.5 * wall_k_s * norm_squared(temp_vec2);

#endif

#ifdef WALL_TWISTING
                V_twist += 0.5 * wall_k_t * cl_entry->contact->twisting_displacement * cl_entry->contact->twisting_displacement;
#endif
            }
            else
            {
                // calculate current contact pointers
                cl_entry->contact->getCurrentN1(&n1);
                cl_entry->contact->getCurrentN2(&n2);

                // determine distance between particles
                n_c[0] = sim->pos_old[X_COORD(p)] - sim->pos_old[X_COORD(cl_entry->id)];
                n_c[1] = sim->pos_old[Y_COORD(p)] - sim->pos_old[Y_COORD(cl_entry->id)];
                n_c[2] = sim->pos_old[Z_COORD(p)] - sim->pos_old[Z_COORD(cl_entry->id)];

                dist = norm(n_c);
                n_c[0] /= dist;
                n_c[1] /= dist;
                n_c[2] /= dist;

                V_normal += sim->normal_interaction.getJKRPotentialEnergy(2.0 * particle_radius - dist);
#ifdef ROLLING
                temp_vec[0] = reduced_radius * (n1[0] + n2[0]);
                temp_vec[1] = reduced_radius * (n1[1] + n2[1]);
                temp_vec[2] = reduced_radius * (n1[2] + n2[2]);
                V_roll += 0.5 * k_r * norm_squared(temp_vec);
#endif

#ifdef SLIDING

                temp_vec[0] = (n1[0] - n2[0] + 2.0 * n_c[0]);
                temp_vec[1] = (n1[1] - n2[1] + 2.0 * n_c[1]);
                temp_vec[2] = (n1[2] - n2[2] + 2.0 * n_c[2]);

                double temp = dot_product(temp_vec, n_c);

                temp_vec2[0] = temp_vec[0] - temp * n_c[0];
                temp_vec2[1] = temp_vec[1] - temp * n_c[1];
                temp_vec2[2] = temp_vec[2] - temp * n_c[2];
                V_slide += 0.5 * k_s * particle_radius * particle_radius * norm_squared(temp_vec2);

#endif

#ifdef TWISTING
                V_twist += 0.5 * k_t * cl_entry->contact->twisting_displacement * cl_entry->contact->twisting_displacement;
#endif
            }

            cl_entry = cl_entry->next;
        }
    }

    V_normal *= ENERGY_UNIT;
    V_roll *= ENERGY_UNIT;
    V_slide *= ENERGY_UNIT;
    V_twist *= ENERGY_UNIT;



    // sum up
    double V_tot = V_normal + V_roll + V_slide + V_twist;
    double E_diss = sim->dissipated_contact_energy + sim->dissipated_rolling_energy + sim->dissipated_sliding_energy + sim->dissipated_twisting_energy + sim->dissipated_wall_energy + sim->dissipated_damping_energy;
    double E_tot = E_kin + E_rot + V_tot + E_diss;

    if(E_kin + E_rot < E_tot*1.e-14)
    {
        //printf("Stopping simulation because all energy is dissipated!\n");
        sim->stop_simulation = true;
    }
    return E_tot;
}







#if defined(TRACK_DISSIPATED_ENERGY_PER_PARTICLE) && defined(TRACK_CONTACTS_PER_PARTICLE)
void SimLib::saveParticleEnergyDissipation(Simulation &sim, char* filename)
{
    std::vector<double> dissipated_energy(12, 0.0);
    std::vector<int> num_particles(12, 0);

    // try to open specified file
    FILE *file = fopen(filename, "w+");

    if(file == NULL)
        return;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        dissipated_energy[ sim.initial_number_of_contacts_of_particle[p] ] += sim.dissipated_energy_of_particle[p];
        num_particles[ sim.initial_number_of_contacts_of_particle[p] ] += 1;
    }

    fprintf(file, "#average dissipated energy with respect to initial contacts\nnumber of initial contacts      avg dissipated energy       sigma dissipated energy\n");

    for(int i = 0; i < 12; ++i)
    {
        if(num_particles[i] > 0)
            dissipated_energy[i] /= (double)num_particles[i];

        fprintf(file, "%i %lf %lf\n", i, dissipated_energy[i], dissipated_energy[i]);
    }

    fclose(file);
}
#endif

void SimLib::printForces(Simulation &sim, const char *filename, bool new_forces)
{
    FILE *file = fopen(filename, "w+");
    //fprintf(file, "# Particles: %i\n\n", sim.number_of_particles);
    fprintf(file, "# Interaction:" );
#ifdef ROLLING
    fprintf(file, " Rolling");
#endif
#ifdef SLIDING
    fprintf(file, " Sliding");
#endif
#ifdef TWISTING
    fprintf(file, " Twisting");
#endif
    fprintf(file, "\n");

    double *torques;
    double *forces;

    if(new_forces)
    {
        forces = sim.force_new;
        torques = sim.torque_new;
    }
    else
    {
        forces = sim.force_old;
        torques = sim.torque_old;
    }

    if(file)
    {
        for(int p = 0; p < sim.number_of_particles; ++p)
            fprintf(file, "%i %g %g %g %g %g %g\n", p, forces[X_COORD(p)], forces[Y_COORD(p)], forces[Z_COORD(p)], torques[X_COORD(p)], torques[Y_COORD(p)], torques[Z_COORD(p)]);

        fclose(file);
    }
}

void SimLib::printContactPointers(Simulation &sim, const char *filename)
{
    FILE *file = fopen(filename, "w+");

    if(!file)
        return;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        ContactListEntry *entry = sim.contact_list[p];

        while(entry)
        {
            vec3 n1, n2;
            entry->contact->getCurrentN1(&n1);
            entry->contact->getCurrentN2(&n2);

            /*fprintf(file, "%i %i: n1 = ( %lf , %lf  , %lf ), n2 = ( %lf , %lf  , %lf )", p, entry->id, n1[0], n1[1], n1[2], n2[0], n2[1], n2[2]);

            fprintf(file, ", e1 = ( %lf , %lf  , %lf , %lf ), e2 = ( %lf , %lf  , %lf , %lf )\n",
                entry->contact->rot1.e0, entry->contact->rot1.e1, entry->contact->rot1.e2, entry->contact->rot1.e3,
                entry->contact->rot2.e0, entry->contact->rot2.e1, entry->contact->rot2.e2, entry->contact->rot2.e3);
            */

            /*vec3 n_c_old, n_c_new;

            n_c_new[0] = sim.pos_new[X_COORD(p)] - sim.pos_new[X_COORD(entry->id)];
            n_c_new[1] = sim.pos_new[Y_COORD(p)] - sim.pos_new[Y_COORD(entry->id)];
            n_c_new[2] = sim.pos_new[Z_COORD(p)] - sim.pos_new[Z_COORD(entry->id)];
            double d_new = norm(n_c_new);
            n_c_new[0] /= d_new;
            n_c_new[1] /= d_new;
            n_c_new[2] /= d_new;

            n_c_old[0] = sim.pos_old[X_COORD(p)] - sim.pos_old[X_COORD(entry->id)];
            n_c_old[1] = sim.pos_old[Y_COORD(p)] - sim.pos_old[Y_COORD(entry->id)];
            n_c_old[2] = sim.pos_old[Z_COORD(p)] - sim.pos_old[Z_COORD(entry->id)];
            double d_old = norm(n_c_old);
            n_c_old[0] /= d_old;
            n_c_old[1] /= d_old;
            n_c_old[2] /= d_old;

            fprintf(file, "    n_c_old = ( %lf , %lf , %lf ), d_old = %lf,  n_c_new = ( %lf , %lf , %lf ), d_new = %lf\n", n_c_old[0], n_c_old[1], n_c_old[2], d_old, n_c_new[0], n_c_new[1], n_c_new[2], d_new);
            */

            vec3 n_c;
            vec3 delta_n;

            delta_n[0] = n1[0] - n2[0];
            delta_n[1] = n1[1] - n2[1];
            delta_n[2] = n1[2] - n2[2];

            n_c[0] = sim.pos_new[X_COORD(p)] - sim.pos_new[X_COORD(entry->id)];
            n_c[1] = sim.pos_new[Y_COORD(p)] - sim.pos_new[Y_COORD(entry->id)];
            n_c[2] = sim.pos_new[Z_COORD(p)] - sim.pos_new[Z_COORD(entry->id)];
            double particle_distance = norm(n_c);

            n_c[0] /= particle_distance;
            n_c[1] /= particle_distance;
            n_c[2] /= particle_distance;

            double temp = dot_product(n_c, delta_n);
            //double a = delta_n[0];
            //double b = temp * n_c[0];
            //double c = delta_n[0] - (temp * n_c[0]);

            fprintf(file, "%i %i: n1 = ( %lg , %lg  , %lg ), n2 = ( %lg , %lg  , %lg )\n", p, entry->id, k_s * particle_radius * particle_radius * temp / particle_distance * (delta_n[0] - temp * n_c[0]), temp, n_c[0], delta_n[0], delta_n[1], k_s * particle_radius * particle_radius * temp / particle_distance);
            //fprintf(file, "%i %i: n1 = ( %lf , %lf  , %lf ), n2 = ( %lf , %lf  , %lf )\n", p, entry->id, n1[0], n2[0], delta_n[0],
            //                                                  a, b, c);

            entry = entry->next;
        }
    }
}

void SimLib::printContactNormals(Simulation &sim, const char *filename)
{
    FILE *file = fopen(filename, "w+");

    if(!file)
        return;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        ContactListEntry *entry = sim.contact_list[p];

        while(entry)
        {
            double timestep = 3e-10;

            vec3 n_c;
            n_c[0] = sim.pos_old[X_COORD(p)] - sim.pos_old[X_COORD(entry->id)];
            n_c[1] = sim.pos_old[Y_COORD(p)] - sim.pos_old[Y_COORD(entry->id)];
            n_c[2] = sim.pos_old[Z_COORD(p)] - sim.pos_old[Z_COORD(entry->id)];

            double dist = sqrt(n_c[0]*n_c[0] + n_c[1]*n_c[1] + n_c[2]*n_c[2]);
            n_c[0] /= dist;
            n_c[1] /= dist;
            n_c[2] /= dist;

            /*vec3 delta_v;
            delta_v[0] = sim.vel[X_COORD(p)] - 0.5 * mass_inv * timestep * (sim.force_new[X_COORD(p)] + sim.force_old[X_COORD(p)]) - sim.vel[X_COORD(entry->id)] + 0.5 * mass_inv * timestep * (sim.force_new[X_COORD(entry->id)] + sim.force_old[X_COORD(entry->id)]);
            delta_v[1] = sim.vel[Y_COORD(p)] - 0.5 * mass_inv * timestep * (sim.force_new[Y_COORD(p)] + sim.force_old[Y_COORD(p)]) - sim.vel[Y_COORD(entry->id)] + 0.5 * mass_inv * timestep * (sim.force_new[Y_COORD(entry->id)] + sim.force_old[Y_COORD(entry->id)]);
            delta_v[2] = sim.vel[Z_COORD(p)] - 0.5 * mass_inv * timestep * (sim.force_new[Z_COORD(p)] + sim.force_old[Z_COORD(p)]) - sim.vel[Z_COORD(entry->id)] + 0.5 * mass_inv * timestep * (sim.force_new[Z_COORD(entry->id)] + sim.force_old[Z_COORD(entry->id)]);

            // calculate normal component of velocity & damping force
            double v_rel = dot_product(delta_v, n_c);*/

            vec3 n1, n2, delta_n;
            entry->contact->getCurrentN1(&n1);
            entry->contact->getCurrentN2(&n2);
            delta_n[0] = n1[0] - n2[0];
            delta_n[1] = n1[1] - n2[1];
            delta_n[2] = n1[2] - n2[2];

            double temp = dot_product(delta_n, n_c);

            /*
            // force sliding
            vec3 force;
            force[0] = k_s * temp / dist * (delta_n[0] - temp * n_c[0]);
            force[1] = k_s * temp / dist * (delta_n[1] - temp * n_c[1]);
            force[2] = k_s * temp / dist * (delta_n[2] - temp * n_c[2]);
            */

            // torque (sliding)
            vec3 temp_vec;
            temp_vec[0] = delta_n[0] - temp * n_c[0];
            temp_vec[1] = delta_n[1] - temp * n_c[1];
            temp_vec[2] = delta_n[2] - temp * n_c[2];

            vec3 M_s;
            M_s[0] = n1[2] - n2[2];
            M_s[1] = temp * n_c[2];
            M_s[2] = temp_vec[2];

            //double contact_radius = sim.normal_interaction.getJKRContactRadius(2.0 * particle_radius - dist);
            //double force = sim.normal_interaction.getJKRForce(contact_radius);

            fprintf(file, "%i %i: n1 = ( %.18lg, %.18lg, %.18lg )  n2 = ( %.18lg, %.18lg, %.18lg ), n_c_new = ( %.20lf , %.20lf , %.20lf ), dist = %.18lg\n",
                        p, entry->id, M_s[0], M_s[1], M_s[2], n2[0], n2[1], n2[2], n_c[0], n_c[1], n_c[2], dist);

            entry = entry->next;
        }
    }

    fclose(file);
}

void SimLib::printSimData(Simulation &sim, const char *filename)
{

    FILE *file = fopen(filename, "w+");

    if(!file)
        return;

    fprintf(file, "Id   pos     vel     omega   force   torque\n");
    //fprintf(file, "Id   pos   omega    torque\n");

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        fprintf(file, "%i: ( %.8lf, %.8lf, %.8lf)   ( %.8lf, %.8lf, %.8lf)   ( %.8lf, %.8lf, %.8lf)   ( %.8lg, %.8lg, %.8lg)   ( %.8lg, %.8lg, %.8lg)\n",
                p,
                sim.pos_old[X_COORD(p)], sim.pos_old[Y_COORD(p)], sim.pos_old[Z_COORD(p)],
                sim.vel[X_COORD(p)], sim.vel[Y_COORD(p)], sim.vel[Z_COORD(p)],
                sim.vel_angular[X_COORD(p)], sim.vel_angular[Y_COORD(p)], sim.vel_angular[Z_COORD(p)],
                sim.force_old[X_COORD(p)], sim.force_old[Y_COORD(p)], sim.force_old[Z_COORD(p)],
                sim.torque_old[X_COORD(p)], sim.torque_old[Y_COORD(p)], sim.torque_old[Z_COORD(p)]);
    }





    fprintf(file, "Id   new pos new force   new torque\n");
    //fprintf(file, "Id   pos   omega    torque\n");

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        fprintf(file, "%i: ( %.8lf, %.8lf, %.8lf)   ( %.8lg, %.8lg, %.8lg)   ( %.8lg, %.8lg, %.8lg)\n",
                p,
                sim.pos_new[X_COORD(p)], sim.pos_new[Y_COORD(p)], sim.pos_new[Z_COORD(p)],
                sim.force_new[X_COORD(p)], sim.force_new[Y_COORD(p)], sim.force_new[Z_COORD(p)],
                sim.torque_new[X_COORD(p)], sim.torque_new[Y_COORD(p)], sim.torque_new[Z_COORD(p)]);
    }




    fprintf(file, "\n\nContacts: id1 id2  n1 n2  n1_init n2_init rot1 rot2\n");

    vec3 n1, n2;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        ContactListEntry *entry = sim.contact_list[p];

        while(entry)
        {

            vec3 n_c;
            if(entry->contact->id2 >= 0)
            {
                n_c[0] = sim.pos_old[X_COORD(p)] - sim.pos_old[X_COORD(entry->id)];
                n_c[1] = sim.pos_old[Y_COORD(p)] - sim.pos_old[Y_COORD(entry->id)];
                n_c[2] = sim.pos_old[Z_COORD(p)] - sim.pos_old[Z_COORD(entry->id)];
                entry->contact->getCurrentN1(&n1);
                entry->contact->getCurrentN2(&n2);
            }
            else
            {
                n_c[0] = sim.pos_old[X_COORD(p)] - sim.walls[WALL_ID(entry->id)].pos[0];
                n_c[1] = sim.pos_old[Y_COORD(p)] - sim.walls[WALL_ID(entry->id)].pos[1];
                n_c[2] = sim.pos_old[Z_COORD(p)] - sim.walls[WALL_ID(entry->id)].pos[2];
                entry->contact->getCurrentN1(&n1);
                n2[0] = 0.0;
                n2[1] = 0.0;
                n2[2] = 0.0;
            }


            double dist = sqrt(n_c[0]*n_c[0] + n_c[1]*n_c[1] + n_c[2]*n_c[2]);
            n_c[0] /= dist;
            n_c[1] /= dist;
            n_c[2] /= dist;

            fprintf(file,
                    "%i %i: twist: %.8lg,\
                    n_c_new = (%.8lg, %.8lg, %.8lg),\
                    n_c_old = (%.8lg, %.8lg, %.8lg),\
                    dist = %.8lg,\
                    n1 = (%.8lf, %.8lf, %.8lf)\
                    n2 = (%.8lf, %.8lf, %.8lf)\
                    n1_i = (%.8lf, %.8lf, %.8lf)\
                    n2_i = (%.8lf, %.8lf, %.8lf)\
                    rot1 = (%.8lg,%.8lg,%.8lg,%.8lg)\
                    rot2 = (%.8lg,%.8lg,%.8lg,%.8lg)",
                p, entry->id, entry->contact->twisting_displacement,
                n_c[0], n_c[1], n_c[2],
                entry->contact->old_contact_normal[0], entry->contact->old_contact_normal[1], entry->contact->old_contact_normal[2],
                dist,
                n1[0], n1[1], n1[2],
                n2[0], n2[1], n2[2],
                entry->contact->n1_initial[0], entry->contact->n1_initial[1], entry->contact->n1_initial[2],
                entry->contact->n2_initial[0], entry->contact->n2_initial[1], entry->contact->n2_initial[2],
                entry->contact->rot1.e0, entry->contact->rot1.e1, entry->contact->rot1.e2, entry->contact->rot1.e3,
                entry->contact->rot2.e0, entry->contact->rot2.e1, entry->contact->rot2.e2, entry->contact->rot2.e3);

            vec3 temp_vec;
            temp_vec[0] = n1[0] - n2[0];
            temp_vec[1] = n1[1] - n2[1];
            temp_vec[2] = n1[2] - n2[2];

            // sliding force
            double temp = temp_vec[0] * n_c[0] + temp_vec[1] * n_c[1] + temp_vec[2] * n_c[2];
            double interaction_strength = k_s * temp / dist;

            temp_vec[0] -= temp * n_c[0];
            temp_vec[1] -= temp * n_c[1];
            temp_vec[2] -= temp * n_c[2];

            vec3 force;
            force[0] = interaction_strength * temp_vec[0];
            force[1] = interaction_strength * temp_vec[1];
            force[2] = interaction_strength * temp_vec[2];

            // sliding torque
            vec3 torque;
            torque[0] = -k_s * (n1[1] * temp_vec[2] - n1[2] * temp_vec[1]);
            torque[1] = -k_s * (n1[2] * temp_vec[0] - n1[0] * temp_vec[2]);
            torque[2] = -k_s * (n1[0] * temp_vec[1] - n1[1] * temp_vec[0]);

            fprintf(file, " F_S = (%.8lg, %.8lg, %.8lg), M_S = (%.8lg, %.8lg, %.8lg)\n", force[0], force[1], force[2], torque[0], torque[1], torque[2]);
            entry = entry->next;
        }
    }

    double E_kin = 0.0;
    double E_rot = 0.0;

    // calculate kinetic + rotational energy
    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        E_kin += 0.5 * mass * (sim.vel[X_COORD(p)] * sim.vel[X_COORD(p)] + sim.vel[Y_COORD(p)] * sim.vel[Y_COORD(p)] + sim.vel[Z_COORD(p)] * sim.vel[Z_COORD(p)]);
        E_rot += 0.5 * moment_of_inertia * (sim.vel_angular[X_COORD(p)] * sim.vel_angular[X_COORD(p)] + sim.vel_angular[Y_COORD(p)] * sim.vel_angular[Y_COORD(p)] + sim.vel_angular[Z_COORD(p)] * sim.vel_angular[Z_COORD(p)]);
    }

    E_kin *= ENERGY_UNIT;
    E_rot *= ENERGY_UNIT;

    fprintf(file, "E_kin = %g, E_rot = %g\n", E_kin, E_rot);

    fclose(file);
}

/*void SimLib::printPositions(Simulation &sim, const char *filename, bool new_positions, bool orientation)
{
    FILE *file = fopen(filename, "w+");

    if(!file)
        return;

    fprintf(file, "# particle radius: %g cm\n", particle_radius);

    for(int p = 0; p < sim.number_of_particles; ++p)
        fprintf(file, "%g %g %g %g %g %g %g %g %g\n", sim.pos_old[X_COORD(p)], sim.pos_old[Y_COORD(p)], sim.pos_old[Z_COORD(p)], sim.vel[X_COORD(p)], sim.vel[Y_COORD(p)], sim.vel[Z_COORD(p)], sim.vel_angular[X_COORD(p)], sim.vel_angular[Y_COORD(p)], sim.vel_angular[Z_COORD(p)]);

    fclose(file);
}*/


void SimLib::printPositions(Simulation &sim, const char *filename, bool new_positions, bool orientation)
{
    FILE *file = fopen(filename, "w+");

    if(!file)
        return;

    fprintf(file, "# particle radius: %g cm\n", particle_radius);
    fprintf(file, "# walls: %i\n", sim.number_of_walls);

    for(int w = 0; w < sim.number_of_walls; ++w)
    {
        sim.walls[w].updateEdges();
        fprintf(file, "%g %g %g %g %g %g %g\n", sim.walls[w].alpha, sim.walls[w].edges[0], sim.walls[w].edges[1], sim.walls[w].edges[2], sim.walls[w].edges[18], sim.walls[w].edges[19], sim.walls[w].edges[20]);
    }

    bool *show_particle = new bool[sim.number_of_particles];


#ifdef TRACK_CONTACTS_PER_PARTICLE
    // get number of contacts per particle
    std::vector< int > contacts(sim.number_of_particles, 0);
    ContactListEntry *cl_entry = NULL;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        cl_entry = sim.contact_list[p];

        while(cl_entry)
        {
            // ignore contacts with walls
            if(cl_entry->id >= 0)
            {
                ++contacts[p];
                ++contacts[cl_entry->id];
            }

            cl_entry = cl_entry->next;
        }
    }

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        //if(sim.broken_contacts_of_particle[p] > 0 && contacts[p] < sim.initial_number_of_contacts_of_particle[p])
        if(sim.broken_contacts_of_particle[p] > 0 && sim.initial_number_of_contacts_of_particle[p] >= 6)
            show_particle[p] = true;
        else
            show_particle[p] = false;
    }
#else
    for(int p = 0; p < sim.number_of_particles; ++p)
        show_particle[p] = true;
#endif

#ifdef TRACK_PARTICLE_ORIENTATION
    if(orientation)
    {
        vec3 v;

        if(new_positions)
        {
            for(int p = 0; p < sim.number_of_particles; ++p)
            {
                if(show_particle[p])
                {
                    v[0] = 2.0 * (sim.orientation[4*p] * sim.orientation[4*p] + sim.orientation[4*p+1] * sim.orientation[4*p+1] - 0.5);
                    v[1] = 2.0 * (sim.orientation[4*p+1] * sim.orientation[4*p+2] + sim.orientation[4*p+3] * sim.orientation[4*p]);
                    v[2] = 2.0 * (sim.orientation[4*p+1] * sim.orientation[4*p+3] - sim.orientation[4*p+2] * sim.orientation[4*p]);

                    fprintf(file, "%g %g %g %g %g %g\n", sim.pos_new[X_COORD(p)],  sim.pos_new[Y_COORD(p)], sim.pos_new[Z_COORD(p)], v[0], v[1], v[2]);
                }
            }
        }
        else
        {
            for(int p = 0; p < sim.number_of_particles; ++p)
            {
                if(show_particle[p])
                {
                    v[0] = 2.0 * (sim.orientation[4*p] * sim.orientation[4*p] + sim.orientation[4*p+1] * sim.orientation[4*p+1] - 0.5);
                    v[1] = 2.0 * (sim.orientation[4*p+1] * sim.orientation[4*p+2] + sim.orientation[4*p+3] * sim.orientation[4*p]);
                    v[2] = 2.0 * (sim.orientation[4*p+1] * sim.orientation[4*p+3] - sim.orientation[4*p+2] * sim.orientation[4*p]);

                    fprintf(file, "%g %g %g %g %g %g\n", sim.pos_old[X_COORD(p)],  sim.pos_old[Y_COORD(p)], sim.pos_old[Z_COORD(p)], v[0], v[1], v[2]);
                }
            }
        }

        fclose(file);
        return;
    }
#endif

    if(new_positions)
    {
        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            if(show_particle[p])
                fprintf(file, "%g %g %g\n", sim.pos_new[X_COORD(p)],  sim.pos_new[Y_COORD(p)], sim.pos_new[Z_COORD(p)]);
        }
    }
    else
    {
        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            if(show_particle[p])
                fprintf(file, "%g %g %g\n", sim.pos_old[X_COORD(p)],  sim.pos_old[Y_COORD(p)], sim.pos_old[Z_COORD(p)]);
        }
    }

    fclose(file);

    delete [] show_particle;
}

void SimLib::getContactHistogram(Simulation &sim, int *particles, bool only_inner_particles)
{
    // get number of contacts per particle
    std::vector< int > contacts(sim.number_of_particles);
    ContactListEntry *cl_entry = NULL;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        cl_entry = sim.contact_list[p];

        while(cl_entry)
        {
            // ignore contacts with walls
            if(cl_entry->id >= 0)
            {
                ++contacts[p];
                ++contacts[cl_entry->id];
            }

            cl_entry = cl_entry->next;
        }
    }

    // determine how many particles have a certain number of contacts
    memset(particles, 0, 13 * sizeof(int));

    if(only_inner_particles)
    {
        vec3 lower_pos, upper_pos;
        sim.getEnclosingBox(&lower_pos, &upper_pos);

        lower_pos[0] += 2.0 * particle_radius;
        lower_pos[1] += 2.0 * particle_radius;
        lower_pos[2] += 2.0 * particle_radius;
        upper_pos[0] += 2.0 * particle_radius;
        upper_pos[1] += 2.0 * particle_radius;
        upper_pos[2] += 2.0 * particle_radius;

        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            if(sim.pos_old[X_COORD(p)] > lower_pos[0] && sim.pos_old[X_COORD(p)] < upper_pos[0]
            && sim.pos_old[Y_COORD(p)] > lower_pos[1] && sim.pos_old[Y_COORD(p)] < upper_pos[1]
            && sim.pos_old[Z_COORD(p)] > lower_pos[2] && sim.pos_old[Z_COORD(p)] < upper_pos[2])
                ++(particles[contacts[p]]);
        }
    }
    else
    {
        for(int p = 0; p < sim.number_of_particles; ++p)
            ++(particles[contacts[p]]);
    }
}

ErrorCode SimLib::printContactHistogram(Simulation &sim, const char *filename, bool only_inner_particles)
{
    FILE *file = fopen(filename, "w+");

    if(file)
    {
        int particles[13];

        getContactHistogram(sim, particles, only_inner_particles);
        fprintf(file, "# number of contacts         number of particles\n");

        for(int i = 0; i < 13; ++i)
            fprintf(file, "%i %i\n", i, particles[i]);

        fclose(file);
    }
    else
        return EC_FILE_NOT_FOUND;

    return EC_OK;
}

ErrorCode SimLib::printFragmentVelocities(Simulation &sim, const char *filename)
{
    std::vector<int> fragment_ids;                          // array storing the fragment id of every particle
    std::vector<int> size_of_fragment;                      // number of particles of the fragment
    std::vector< std::list<int> > particles_of_fragment;    // ids of the particles of a specific fragment

    SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, &particles_of_fragment);

    int fragments = (int)size_of_fragment.size();
    double *fragment_velocities = new double[3*fragments];
    memset(fragment_velocities, 0, 3 * fragments * sizeof(double));

    // determine largest and second largest fragment
    int largest_fragment = 0;
    int second_largest_fragment = 0;

    FILE *file = fopen(filename, "w+");
    fprintf(file, "# impact speed: %g\n", sim.sim_info.info_storage[2]);
    fprintf(file, "# number of monomers     velocity (x,y,z)    speed\n");

    for(int f = 0; f < fragments; ++f)
    {
        size_t fragment_size = particles_of_fragment[f].size();

        if(fragment_size > particles_of_fragment[largest_fragment].size())
        {
            second_largest_fragment = largest_fragment;
            largest_fragment = f;
        }
        else
        {
            if(fragment_size > particles_of_fragment[second_largest_fragment].size() || second_largest_fragment == largest_fragment)
                second_largest_fragment = f;
        }

        // calculate velocity of fragment
        for(std::list<int>::iterator p = particles_of_fragment[f].begin(); p != particles_of_fragment[f].end(); ++p)
        {
            fragment_velocities[3*f+0] += sim.vel[X_COORD(*p)];
            fragment_velocities[3*f+1] += sim.vel[Y_COORD(*p)];
            fragment_velocities[3*f+2] += sim.vel[Z_COORD(*p)];
        }

        fragment_velocities[3*f+0] /= (double)fragment_size;
        fragment_velocities[3*f+1] /= (double)fragment_size;
        fragment_velocities[3*f+2] /= (double)fragment_size;

        double speed = sqrt(fragment_velocities[3*f+0]*fragment_velocities[3*f+0] + fragment_velocities[3*f+1]*fragment_velocities[3*f+1] + fragment_velocities[3*f+2]*fragment_velocities[3*f+2]);

        fprintf(file, "%i %g %g %g %g\n", (int)fragment_size, fragment_velocities[3*f+0], fragment_velocities[3*f+1], fragment_velocities[3*f+2], speed);
    }

    if(largest_fragment != second_largest_fragment )
    {
        if( SimLib::getCollisionResult((int)particles_of_fragment[largest_fragment].size(), (int)particles_of_fragment[second_largest_fragment].size(), sim.number_of_particles) == COLLISION_RESULT_BOUNCING)
        {
            vec3 delta_v;
            delta_v[0] = fragment_velocities[3*largest_fragment+0] - fragment_velocities[3*second_largest_fragment+0];
            delta_v[1] = fragment_velocities[3*largest_fragment+1] - fragment_velocities[3*second_largest_fragment+1];
            delta_v[2] = fragment_velocities[3*largest_fragment+2] - fragment_velocities[3*second_largest_fragment+2];
            double speed = sqrt(delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2]);

            fprintf(file, "# largest / second largest fragment: %i / %i\n", largest_fragment, second_largest_fragment);
            fprintf(file, "# delta v = (%g, %g, %g) => ||delta v|| = %g\n", delta_v[0], delta_v[1], delta_v[2], speed);

            if(sim.sim_info.info_storage[2] > 0)
                fprintf(file, "# => coefficient of restitution = %g\n", speed / sim.sim_info.info_storage[2]);
        }
    }

    fclose(file);
    delete [] fragment_velocities;

    return EC_OK;
}

ErrorCode SimLib::printElasticCharge(Simulation &sim, FILE *file)
{
    int charged_contacts = 0;
    double elastic_charge = SimLib::getElasticCharge(sim, &charged_contacts);
    double kinetic_energy = SimLib::getKineticEnergy(sim);

    /*fprintf(file, "# number of contacts: %i\n", sim.number_of_contacts);
    fprintf(file, "# number of charged contacts: %i\n", charged_contacts);
    fprintf(file, "# elastic charge: %g (in F_c * delta_c)\n", elastic_charge *  ENEGRY_UNIT);*/

    if(file)
        fprintf(file, "%i %i %g %g\n", sim.number_of_contacts, charged_contacts, elastic_charge * ENERGY_UNIT, kinetic_energy * ENERGY_UNIT);

    return EC_OK;
}

double SimLib::getElasticCharge(Simulation &sim, int *charged_contacts)
{
    double elastic_charge = 0;
    ContactListEntry *cl_entry;
    vec3 n1;
    vec3 n2;
    vec3 delta_n;
    vec3 n_c;
    vec3 displacement;
    double rolling_displacement;
    double sliding_displacement;

    *charged_contacts = 0;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        cl_entry = sim.contact_list[p];

        while(cl_entry)
        {
            bool inelastic = false;

            // calculate current contact pointers
            cl_entry->contact->getCurrentN1(&n1);
            cl_entry->contact->getCurrentN2(&n2);
            delta_n[0] = n1[0] - n2[0];
            delta_n[1] = n1[1] - n2[1];
            delta_n[2] = n1[2] - n2[2];

            // determine distance between particles & contact normal
            n_c[0] = sim.pos_new[X_COORD(p)] - sim.pos_new[X_COORD(cl_entry->id)];
            n_c[1] = sim.pos_new[Y_COORD(p)] - sim.pos_new[Y_COORD(cl_entry->id)];
            n_c[2] = sim.pos_new[Z_COORD(p)] - sim.pos_new[Z_COORD(cl_entry->id)];
            double particle_distance = norm(n_c);

            n_c[0] /= particle_distance;
            n_c[1] /= particle_distance;
            n_c[2] /= particle_distance;

            // check inelastic rolling
            displacement[0] = reduced_radius * (n1[0] + n2[0]);
            displacement[1] = reduced_radius * (n1[1] + n2[1]);
            displacement[2] = reduced_radius * (n1[2] + n2[2]);
            rolling_displacement = norm_squared(displacement);

            if(rolling_displacement > 0.99 * crit_rolling_displacement_squared)
                inelastic = true;

            // check inelastic sliding
            double temp = dot_product(delta_n,  n_c);
            displacement[0] = particle_radius * (delta_n[0]  - temp * n_c[0]);
            displacement[1] = particle_radius * (delta_n[1]  - temp * n_c[1]);
            displacement[2] = particle_radius * (delta_n[2]  - temp * n_c[2]);
            sliding_displacement = norm_squared(displacement);

            // check if we are in the inelastic regime
            if(sliding_displacement > 0.99 * crit_sliding_displacement_squared)
                inelastic = true;

            if(!inelastic)
            {
                *charged_contacts += 1;
                elastic_charge += 0.5 * k_r * rolling_displacement;
                elastic_charge += 0.5 * k_s * sliding_displacement;
            }

            cl_entry = cl_entry->next;
        }
    }

    return elastic_charge;
}

double SimLib::getKineticEnergy(Simulation &sim)
{
    double E_kin = 0.0;
    double E_rot = 0.0;

    // calculate kinetic + rotational energy
    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        E_kin += 0.5 * mass * (sim.vel[X_COORD(p)] * sim.vel[X_COORD(p)] + sim.vel[Y_COORD(p)] * sim.vel[Y_COORD(p)] + sim.vel[Z_COORD(p)] * sim.vel[Z_COORD(p)]);
        E_rot += 0.5 * moment_of_inertia * (sim.vel_angular[X_COORD(p)] * sim.vel_angular[X_COORD(p)] + sim.vel_angular[Y_COORD(p)] * sim.vel_angular[Y_COORD(p)] + sim.vel_angular[Z_COORD(p)] * sim.vel_angular[Z_COORD(p)]);
    }

    return (E_kin + E_rot);
}

void SimLib::getOrthoVector(Simulation *sim, vec3 *ortho_vec, vec3 &vec)
{
    // determine two random angles
    double phi = sim->get_random_zero_twoPi();
    double cos_theta = sim->get_random_cos_theta();
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

    // calculate vector
    vec3 vec2;
    vec2[0] = sin_theta * cos(phi);
    vec2[1] = sin_theta * sin(phi);
    vec2[2] = cos_theta;

    (*ortho_vec)[0] = vec[1]*vec2[2] - vec[2]*vec2[1];
    (*ortho_vec)[1] = vec[2]*vec2[0] - vec[0]*vec2[2];
    (*ortho_vec)[2] = vec[0]*vec2[1] - vec[1]*vec2[0];

    normalize(ortho_vec);
}

double SimLib::getMaxSpeed(Simulation *sim)
{
    double max_speed = 0, speed;

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        speed = sim->vel[X_COORD(p)]*sim->vel[X_COORD(p)] + sim->vel[Y_COORD(p)]*sim->vel[Y_COORD(p)] + sim->vel[Z_COORD(p)]*sim->vel[Z_COORD(p)];

        if(speed > max_speed)
            max_speed = speed;
    }

    return sqrt(max_speed);
}



double SimLib::getFillingFactorOfBox(Simulation &sim, vec3 &lower_pos, vec3 &upper_pos)
{
    double V_particle = 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius;
    double box_volume = (upper_pos[0] - lower_pos[0]) * (upper_pos[1] - lower_pos[1]) * (upper_pos[2] - lower_pos[2]);

    if(box_volume < 0)
        return 0;

    double filling_factor = 0.0;
    double x_lower_dist;
    double x_upper_dist;
    double y_lower_dist;
    double y_upper_dist;
    double z_lower_dist;
    double z_upper_dist;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        x_lower_dist = sim.pos_old[X_COORD(p)] - lower_pos[0];
        x_upper_dist = upper_pos[0] - sim.pos_old[X_COORD(p)];
        y_lower_dist = sim.pos_old[Y_COORD(p)] - lower_pos[1];
        y_upper_dist = upper_pos[1] - sim.pos_old[Y_COORD(p)];
        z_lower_dist = sim.pos_old[Z_COORD(p)] - lower_pos[2];
        z_upper_dist = upper_pos[2] - sim.pos_old[Z_COORD(p)];

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check if particle is located within box
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(x_lower_dist > 0.0 && x_upper_dist > 0.0 && y_lower_dist > 0.0 && y_upper_dist > 0.0 && z_lower_dist > 0.0 && z_upper_dist > 0.0)
        {
            filling_factor += V_particle;

            // detect intersections with edges
            if(x_lower_dist < particle_radius)
            {
                double height = particle_radius - x_lower_dist;
                filling_factor -=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(x_upper_dist < particle_radius)
            {
                double height = particle_radius - x_upper_dist;
                filling_factor -=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(y_lower_dist < particle_radius)
            {
                double height = particle_radius - y_lower_dist;
                filling_factor -=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(y_upper_dist < particle_radius)
            {
                double height = particle_radius - y_upper_dist;
                filling_factor -=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(z_lower_dist < particle_radius)
            {
                double height = particle_radius - z_lower_dist;
                filling_factor -=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(z_upper_dist < particle_radius)
            {
                double height = particle_radius - z_upper_dist;
                filling_factor -=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check if particle is located within less than a particle radius outside of the box
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        else if(x_lower_dist > -particle_radius && x_upper_dist > -particle_radius && y_lower_dist > -particle_radius && y_upper_dist > -particle_radius && z_lower_dist > -particle_radius && z_upper_dist > -particle_radius)
        {
            // detect intersections with edges
            if(x_lower_dist < 0.0)
            {
                double height = particle_radius + x_lower_dist;
                filling_factor +=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(x_upper_dist < 0.0)
            {
                double height = particle_radius + x_upper_dist;
                filling_factor +=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(y_lower_dist < 0.0)
            {
                double height = particle_radius + y_lower_dist;
                filling_factor +=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(y_upper_dist < 0.0)
            {
                double height = particle_radius + y_upper_dist;
                filling_factor +=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(z_lower_dist < 0.0)
            {
                double height = particle_radius + z_lower_dist;
                filling_factor +=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }

            if(z_upper_dist < 0.0)
            {
                double height = particle_radius + z_upper_dist;
                filling_factor +=  M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
            }
        }
    }

    return (filling_factor / box_volume);
}


double SimLib::getFillingFactorOfSlice(Simulation &sim, double y_lower, double y_upper, vec3 &lower_pos, vec3 &upper_pos)
{
    double V_particle = 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius;
    double slice_volume = (upper_pos[0] - lower_pos[0]) * (upper_pos[2] - lower_pos[2]) * (y_upper - y_lower);

    if(slice_volume < 0.0)
        return 0;

    double filling_factor = 0.0;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        if(sim.pos_old[Y_COORD(p)] > y_lower - particle_radius)
        {
            if(sim.pos_old[Y_COORD(p)] < y_lower)
            {
                double height = particle_radius - y_lower + sim.pos_old[Y_COORD(p)];
                double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                filling_factor += V;
            }
            else if(sim.pos_old[Y_COORD(p)] < y_upper)
            {
                if(sim.pos_old[Y_COORD(p)] < y_lower + particle_radius)
                {
                    double height = particle_radius - sim.pos_old[Y_COORD(p)] + y_lower;
                    double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                    filling_factor += (V_particle - V);
                }
                else if(sim.pos_old[Y_COORD(p)] > y_upper - particle_radius)
                {
                    double height = particle_radius - y_upper + sim.pos_old[Y_COORD(p)];
                    double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                    filling_factor += (V_particle - V);
                }
                else
                    filling_factor += V_particle;
            }
            else if(sim.pos_old[Y_COORD(p)] < y_upper + particle_radius)
            {
                double height = particle_radius - sim.pos_old[Y_COORD(p)] + y_upper;
                double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                filling_factor += V;
            }
        }
    }

    filling_factor /= slice_volume;
    return filling_factor;
}



void SimLib::getFillingFactorOfSphere(Simulation &sim, double r_start, double r_end, double &filling_fact, double &coord_number)
{
    double V_particle = 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius;
    double layer_volume = 4.0 / 3.0 * M_PI * (r_end*r_end*r_end - r_start*r_start*r_start);

    if(layer_volume < 0.0)
    {
        filling_fact = 0.0;
        coord_number = 0.0;
        return;
    }

    int number_of_contacts = 0;
    int number_of_particles = 0;

    ContactListEntry *cl_entry;
    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        double r = sqrt(sim.pos_old[X_COORD(p)]*sim.pos_old[X_COORD(p)]
                + sim.pos_old[Y_COORD(p)]*sim.pos_old[Y_COORD(p)]
                + sim.pos_old[Z_COORD(p)]*sim.pos_old[Z_COORD(p)]);

        if(r > r_start && r < r_end)
        {

            ++number_of_particles;

            cl_entry = sim.contact_list[p];

            while(cl_entry)
            {
                if(cl_entry->id >= 0)
                    ++number_of_contacts;

                cl_entry = cl_entry->next;
            }


        }
    }

    if(number_of_particles > 0)
        coord_number = 2.0*(double)number_of_contacts/(double)number_of_particles;
    else
        coord_number = 0.0;


    double filling_factor = 0.0;
    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        double r = sqrt(sim.pos_old[X_COORD(p)]*sim.pos_old[X_COORD(p)]
                + sim.pos_old[Y_COORD(p)]*sim.pos_old[Y_COORD(p)]
                + sim.pos_old[Z_COORD(p)]*sim.pos_old[Z_COORD(p)]);

        if(r > r_start - particle_radius)
        {
            if(r < r_start)
            {
                double height = particle_radius - r_start + r;
                double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                filling_factor += V;
            }
            else if(r < r_end)
            {
                if(r < r_start + particle_radius)
                {
                    double height = particle_radius - r + r_start;
                    double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                    filling_factor += (V_particle - V);
                }
                else if(r > r_end - particle_radius)
                {
                    double height = particle_radius - r_end + r;
                    double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                    filling_factor += (V_particle - V);
                }
                else
                    filling_factor += V_particle;
            }
            else if(r < r_end + particle_radius)
            {
                double height = particle_radius - r + r_end;
                double V = M_PI / 3.0 * height*height * (3.0 * particle_radius - height);
                filling_factor += V;
            }
        }
    }




    filling_factor /= layer_volume;


    filling_fact = filling_factor;
}


void SimLib::getParticlesInSlice(Simulation &sim, vec3 &lower_pos, vec3 &upper_pos, int *particles, int *contacts)
{
    int num_particles = 0;
    int num_contacts = 0;

    double x_lower_dist, x_upper_dist, y_lower_dist, y_upper_dist, z_lower_dist, z_upper_dist;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        x_lower_dist = sim.pos_old[X_COORD(p)] - lower_pos[0];
        x_upper_dist = upper_pos[0] - sim.pos_old[X_COORD(p)];
        y_lower_dist = sim.pos_old[Y_COORD(p)] - lower_pos[1];
        y_upper_dist = upper_pos[1] - sim.pos_old[Y_COORD(p)];
        z_lower_dist = sim.pos_old[Z_COORD(p)] - lower_pos[2];
        z_upper_dist = upper_pos[2] - sim.pos_old[Z_COORD(p)];

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check if particle is located within box
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(x_lower_dist > 0 && x_upper_dist > 0 && y_lower_dist > 0 && y_upper_dist > 0 && z_lower_dist > 0 && z_upper_dist > 0)
        {
            ++num_particles;

            if(contacts)
            {
                ContactListEntry *cl_entry = sim.contact_list[p];

                while(cl_entry)
                {
                    ++num_contacts;
                    cl_entry = cl_entry->next;
                }
            }
        }
    }

    if(particles)
        *particles = num_particles;

    if(contacts)
        *contacts = num_contacts;
}

void SimLib::printFillingFactorProfile(Simulation &sim, const char *filename, double slice_height)
{
    vec3 lower_pos, upper_pos;
    sim.getEnclosingBox(&lower_pos, &upper_pos);

    FILE *file = fopen(filename, "w+");

    if(file)
    {
        fprintf(file, "# height  filling factor\n");
        fprintf(file, "# current time: %g\n", sim.current_time);
        double start_height = lower_pos[1] + particle_radius;
        double end_height = upper_pos[1] - particle_radius - slice_height;
        double y_pos = lower_pos[1] + particle_radius;

        while(y_pos < end_height)
        {
            fprintf(file, "%lf %lf\n", y_pos - start_height, SimLib::getFillingFactorOfSlice(sim, y_pos, y_pos + slice_height, lower_pos, upper_pos));
            y_pos += slice_height;
        }

        fclose(file);
    }
}


void SimLib::printFillingFactorProfileSphere(Simulation &sim, const char *filename, double layer_width)
{
    SimLib::centerCMS(&sim);

    double sphere_radius = 0.0;
    double sphere_gradius = 0.0;

    SimLib::getSize(sim, &sphere_gradius, &sphere_radius);

    double cross_section1;
    double sigma_cross_section;

    SimLib::getCrossSection(sim, 0.1*particle_radius, 12, &cross_section1, &sigma_cross_section);


    double radius1 = sqrt(cross_section1 / M_PI);


    FILE *file = fopen(filename, "w+");

    if(file)
    {
        fprintf(file, "# radius in µm  filling factor   coordination number\n");
        fprintf(file, "%e\n", radius1*1.e+4);

        sphere_radius = sphere_radius + particle_radius;
        double filling_factor = sim.number_of_particles*particle_radius*particle_radius*particle_radius/(sphere_radius*sphere_radius*sphere_radius);
        double coord_number = 2.0 * double(sim.getNumberOfContacts()) / double(sim.number_of_particles);
        fprintf(file, "%lf %lf  %lf\n", sphere_radius*1.e+4, filling_factor, coord_number);


        double start_radius = 0.0;
        double end_radius = particle_radius * 5.0;
        SimLib::getFillingFactorOfSphere(sim, start_radius, end_radius, filling_factor, coord_number);

        while(filling_factor != 0.0)
        {
            fprintf(file, "%lf %lf  %lf\n", end_radius*1.e+4, filling_factor, coord_number);
            start_radius = end_radius;
            end_radius += layer_width;
            SimLib::getFillingFactorOfSphere(sim, start_radius, end_radius, filling_factor, coord_number);
        }

        fprintf(file, "%lf %lf  %lf\n", end_radius*1.e+4, filling_factor, coord_number);

        fclose(file);
    }
}

ErrorCode SimLib::printFractalDimension(Simulation &sim, const char *filename, double radius_increment, double inner_cutoff, double outer_cutoff)
{
    vec3 cms;
    SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);

    double outer_radius;
    double gyration_radius;
    SimLib::getSize(sim, &gyration_radius, &outer_radius);

    // determine number of shells
    int shells = (int) ceil( outer_radius / radius_increment );
    std::vector<int> particles_in_shell(shells, 0);

    // sort particles in shells
    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        double dist = (cms[0] - sim.pos_old[X_COORD(p)]) * (cms[0] - sim.pos_old[X_COORD(p)]) + (cms[1] - sim.pos_old[Y_COORD(p)]) * (cms[1] - sim.pos_old[Y_COORD(p)]) + (cms[2] - sim.pos_old[Z_COORD(p)]) * (cms[2] - sim.pos_old[Z_COORD(p)]);
        dist = sqrt(dist);

        int shell = (int)floor(dist / radius_increment);
        ++particles_in_shell[shell];
    }

    // calculate & print results
    FILE* file = fopen(filename, "w+");
    if(!file)
        return EC_FILE_NOT_FOUND;

    fprintf(file, "# monomers: %i   outer_radius: %g   gyration_radius: %g\n", sim.number_of_particles, outer_radius, gyration_radius);
    fprintf(file, "# sphere radius    number of particles   particles/volume\n");

    for(int s = 0; s < shells; ++s)
    {
        double radius = (double)(s+1) * radius_increment;

        if(s > 0)
        {
            particles_in_shell[s] += particles_in_shell[s-1];
            int diff = particles_in_shell[s] - particles_in_shell[s-1];
            double old_radius = (double)s * radius_increment;

            if(radius > outer_radius * inner_cutoff && radius < outer_radius * outer_cutoff)
                fprintf(file, "%g %i %g\n", radius, particles_in_shell[s], (double)diff * particle_radius*particle_radius*particle_radius / (radius*radius*radius - old_radius*old_radius*old_radius));
        }
        else
        {
            if(radius > outer_radius * inner_cutoff && radius < outer_radius * outer_cutoff)
                fprintf(file, "%g %i %g\n", radius, particles_in_shell[s], (double)particles_in_shell[s] * particle_radius*particle_radius*particle_radius / (radius*radius*radius));
        }
    }

    fclose(file);

    return EC_OK;
}

void SimLib::checkEquiblibriumDistances(Simulation &sim, const char *filename)
{
    FILE *file = fopen(filename, "w+");

    if(file)
    {
        ContactListEntry *entry;
        double particle_distance;
        double equilibrium_distance = 2.0 * particle_radius - delta_0;
        vec3 n_c;

        double deviation;
        double max_deviation = 0.0;
        double avg_deviation = 0.0;

        for(int p = 0; p < sim.number_of_particles; ++p)
        {
            entry = sim.contact_list[p];

            while(entry)
            {
                if(entry->id >= 0)
                {
                    n_c[0] = sim.pos_old[X_COORD(p)] - sim.pos_old[X_COORD(entry->id)];
                    n_c[1] = sim.pos_old[Y_COORD(p)] - sim.pos_old[Y_COORD(entry->id)];
                    n_c[2] = sim.pos_old[Z_COORD(p)] - sim.pos_old[Z_COORD(entry->id)];
                    particle_distance = norm(n_c);

                    deviation = fabs(1.0 - particle_distance / equilibrium_distance);
                    fprintf(file, "%i <-> %i: %g\n", p , entry->id, deviation);

                    if(deviation > max_deviation)
                        max_deviation = deviation;

                    avg_deviation += deviation;
                }

                entry = entry->next;
            }
        }

        avg_deviation /= (double)sim.number_of_contacts;

        fprintf(file, "\navg_deviation: %g   max_deviation: %g\n", avg_deviation, max_deviation);

        fclose(file);
    }
}

CollisionResult SimLib::getCollisionResult(int size_largest_fragment, int size_second_largest_fragment, int number_of_particles, double sticking_threshhold)
{
    if( size_largest_fragment > (int)(sticking_threshhold * (double)number_of_particles) )
        return COLLISION_RESULT_STICKING;
    else if( size_second_largest_fragment > (int)(0.2 * (double)number_of_particles)
        && (size_largest_fragment + size_second_largest_fragment) > (int)(0.95 * (double)number_of_particles) )
        return COLLISION_RESULT_BOUNCING;
    else
        return COLLISION_RESULT_FRAGMENTATION;
}



ErrorCode SimLib::importParticlePositions(Simulation *sim, const char *filename)
{

    std::setlocale(LC_ALL, "C"); // for some reasong local was set to comma as digit separator, causing the files to be read incorrectly


    std::vector<double> vec_x;
    std::vector<double> vec_y;
    std::vector<double> vec_z;

    FILE *file = fopen(filename, "r");

    if(!file)
        return EC_FILE_NOT_FOUND;


    int number_of_walls;
    double radius;

    fscanf(file, "# particle radius: %lf cm\n", &radius);
    fscanf(file, "# walls: %i\n", &number_of_walls);

    char buf[1024];

    Wall* walls = new Wall[6];

    printf("num walls = %d  particle radius = %e\n", number_of_walls, radius);
    for(int w = 0; w < number_of_walls; ++w)
    {
        fgets(buf, sizeof(buf), file);
        printf("    %s", buf);
        char* p = buf;


        p = strtok (p," ");


        while (p != NULL)
        {
          int err = 0;
          double val;
          while(err == 0)
          {
            err = sscanf(p, "%lf", &val);
            printf("erg = |%s|    %f\n", p, val);
          }
          p = strtok (NULL, " ");
        }


   }
    printf("walls loaded\n");

    particle_radius = radius;


    while (fgets(buf, sizeof(buf), file) != NULL) {
        double x;
        double y;
        double z;
        sscanf(buf, "%lf %lf %lf\n", &x, &y, &z);
        vec_x.push_back(x);
        vec_y.push_back(y);
        vec_z.push_back(z);
    }


    sim->resizeArrays(vec_x.size(), number_of_walls);

    printf("number of particles loaded %d\n", sim->number_of_particles);

    for(int w = 0; w < number_of_walls; ++w)
    {
        sim->walls[w].alpha = walls[w].alpha;
        sim->walls[w].edges[0] = walls[w].edges[0];
        sim->walls[w].edges[1] = walls[w].edges[1];
        sim->walls[w].edges[2] = walls[w].edges[2];
        sim->walls[w].edges[18] = walls[w].edges[18];
        sim->walls[w].edges[19] = walls[w].edges[19];
        sim->walls[w].edges[20] = walls[w].edges[20];
    }


    delete[] walls;

    for(int p = 0; p < vec_x.size(); ++p)
    {

        sim->pos_old[X_COORD(p)] = vec_x[p];
        sim->pos_old[Y_COORD(p)] = vec_y[p];
        sim->pos_old[Z_COORD(p)] = vec_z[p];

    }

    // set default values
    sim->sim_info.info_storage[0] = 0.0;
    sim->sim_info.info_storage[1] = 0.0;
    sim->sim_info.info_storage[2] = 0.0;
    sim->sim_info.info_storage[3] = 0.0;
    sim->sim_info.info_storage[4] = 0.0;
    sim->sim_info.info_storage[5] = 0.0;
    sim->sim_info.info_storage[6] = 0.0;
    sim->sim_info.info_storage[7] = 0.0;

    fclose(file);
    return EC_OK;


}



/*
ErrorCode SimLib::importParticlePositions(Simulation *sim, const char *filename)
{
    FILE *file = fopen(filename, "r");

    if(!file)
        return EC_FILE_NOT_FOUND;

    char version[200];
    fscanf(file, "%s", version);

    if(strcmp(version, "IMPORT_TYPE_1") == 0)
    {
        int number_of_particles;
        double radius;

        fscanf(file, "%i %lf", &number_of_particles, &radius);

        if(fabs(radius - particle_radius) > 1e-10)
            return EC_DIFFERING_MATERIAL;

        sim->resizeArrays(number_of_particles, 0);

        // set default values
        sim->sim_info.info_storage[0] = 0.0;
        sim->sim_info.info_storage[1] = 0.0;
        sim->sim_info.info_storage[2] = 0.0;
        sim->sim_info.info_storage[3] = 0.0;
        sim->sim_info.info_storage[4] = 0.0;
        sim->sim_info.info_storage[5] = 0.0;
        sim->sim_info.info_storage[6] = 0.0;
        sim->sim_info.info_storage[7] = 0.0;


        for(int p = 0; p < number_of_particles; ++p)
            fscanf(file, "%lf %lf %lf\n", &(sim->pos_old[X_COORD(p)]), &(sim->pos_old[Y_COORD(p)]), &(sim->pos_old[Z_COORD(p)]));

        fclose(file);
        return EC_OK;
    }
    else if(strcmp(version, "IMPORT_TYPE_2") == 0)
    {
        int number_of_particles;
        double radius;

        fscanf(file, "%i %lf", &number_of_particles, &radius);

        if(fabs(radius - particle_radius) > 1e-10)
            return EC_DIFFERING_MATERIAL;

        sim->resizeArrays(number_of_particles, 0);

        // set default values
        sim->sim_info.info_storage[0] = 0.0;
        sim->sim_info.info_storage[1] = 0.0;
        sim->sim_info.info_storage[2] = 0.0;
        sim->sim_info.info_storage[3] = 0.0;
        sim->sim_info.info_storage[4] = 0.0;
        sim->sim_info.info_storage[5] = 0.0;
        sim->sim_info.info_storage[6] = 0.0;
        sim->sim_info.info_storage[7] = 0.0;


        double dummy;

        for(int p = 0; p < number_of_particles; ++p)
        {
            fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &(sim->pos_old[X_COORD(p)]), &(sim->pos_old[Y_COORD(p)]), &(sim->pos_old[Z_COORD(p)]), &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy);
            sim->pos_old[X_COORD(p)] *= 100.0;
            sim->pos_old[Y_COORD(p)] *= 100.0;
            sim->pos_old[Z_COORD(p)] *= 100.0;
        }

        fclose(file);
        return EC_OK;
    }
    else
        return EC_INVALID_FILE_VERSION;
}
*/

#ifdef ENABLE_FIXED_PARTICLES
ErrorCode SimLib::fixateParticles(Simulation *sim, int first_id, int last_id)
{
    if(sim->number_of_particles <= last_id)
        return EC_INVALID_PARAMETER;

    if(!sim->fixed_particles)
        sim->fixed_particles = new bool[sim->number_of_particles];

    memset(sim->fixed_particles, 0, sim->number_of_particles * sizeof(bool));

    for(int p = 0; p < sim->number_of_particles; ++p)
    {
        if(p >= first_id && p <= last_id)
            sim->fixed_particles[p] = true;
    }

    return EC_OK;
}
#endif

ErrorCode SimLib::resetContacts(Simulation *sim)
{


    sim->deleteContacts();

    memcpy(sim->pos_new, sim->pos_old, sim->number_of_particles * 3 * sizeof(double));
    sim->updateSticking();
    //sim->resetContacts();
    return EC_OK;
}

