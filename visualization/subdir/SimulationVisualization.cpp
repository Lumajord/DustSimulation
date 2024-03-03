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
 *   Free Software Foundation, Inc.											*
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.				*
 ***************************************************************************/

#include "SimulationVisualization.h"
#include "SimulationLib.h"
#include <math.h>
#include <string.h>
#include <numeric>	// std::accumulate

extern double particle_radius;
extern double reduced_radius;
extern double mass;
extern double moment_of_inertia;
extern double k_r;
extern double k_s;
extern double k_t;
extern double ENERGY_UNIT;
extern double delta_0;
extern double F_c;
extern double v_c;
extern double surface_energy;
extern double crit_rolling_displacement;

extern int vis_neighbour_search_dist;
extern int vis_neighbour_max_particles;
extern int vis_neighbour_min_offset;
extern int vis_neighbour_max_offset;
extern double vis_density_max_filling_factor;
extern bool vis_density_consider_walls;
extern double vis_dislocation_min_value;
extern double vis_dislocation_max_value;

SimulationVis::SimulationVis(void) : Simulation()
{
	max_neighbours = 0;
	max_diss_energy = 0;
	max_speed = 0;
	avg_speed[0] = avg_speed[1] = avg_speed[2] = 0;

	fragment_ids = NULL;
	fragment_id_offset = 0;
	
	force_avg_counter = 0;
	force_avg_size = 0;
	forces_storage = NULL;

	initial_positions = NULL;
}

SimulationVis::~SimulationVis(void)
{
	delete [] fragment_ids;
	fragment_ids = NULL;

	delete [] initial_positions;
	initial_positions = NULL;

	if(forces_storage)
	{
		delete [] forces_storage;
		forces_storage = NULL;
		force_avg_size = 0;
	}

}

void SimulationVis::update()
{
	Simulation::update();

	if(box && forces_storage)
	{
		if(sim_info.sim_type == SIM_TYPE_SHOCKWAVE || sim_info.sim_type == SIM_TYPE_OPEN_BOX)
			forces_storage[force_avg_counter] = walls[box->bottom_wall_id].getWallForce();
		else if(sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
			forces_storage[force_avg_counter] = walls[box->top_wall_id].getWallForce();
		else
			forces_storage[force_avg_counter] = walls[box->top_wall_id].getWallForce();

		++force_avg_counter;

		if(force_avg_counter >= force_avg_size)
			force_avg_counter = 0;
	}
}

void SimulationVis::predictor(double timestep)
{
	Simulation::predictor(timestep);

	// update wall edges (for visualization only)
	for(int w = 0; w < walls.size(); ++w)
		walls[w].updateEdges();
}

ErrorCode SimulationVis::initSimData(SimType type, double collision_speed, double stop_filling_factor, double stop_dissipation_factor, double initial_kinetic_energy, double shockwave_penetration_depth)
{
    ErrorCode status = Simulation::initSimData(type, collision_speed, stop_filling_factor, stop_dissipation_factor, initial_kinetic_energy, shockwave_penetration_depth);

	if(status != EC_OK)
		return status;
	else
	{
		memcpy(initial_positions, pos_old, 3*number_of_particles*sizeof(double));

		fragment_id_offset = SimLib::detectFragments(*this, fragment_ids);

		prepareColorScaling();

		return EC_OK;
	}
}

ErrorCode SimulationVis::resizeArrays(int new_number_of_particles, int new_number_of_walls)
{
	Simulation::resizeArrays(new_number_of_particles, new_number_of_walls);

	delete [] initial_positions;
	initial_positions = new double[3*new_number_of_particles];
	memcpy(initial_positions, pos_old, 3 * new_number_of_particles * sizeof(double));

	delete [] fragment_ids;
	fragment_ids = new int[new_number_of_particles];

	return EC_OK;
}

ErrorCode SimulationVis::addParticlesFromSim(Simulation &sim)
{
	int new_number_of_particles = number_of_particles + sim.number_of_particles;

	delete [] initial_positions;
	initial_positions = new double[3*new_number_of_particles];

	delete [] fragment_ids;
	fragment_ids = new int[new_number_of_particles];

	return Simulation::addParticlesFromSim(sim);
}

ErrorCode SimulationVis::addParticles(int particles)
{
	int new_number_of_particles = number_of_particles + particles;

	delete [] initial_positions;
	initial_positions = new double[3*new_number_of_particles];

	delete [] fragment_ids;
	fragment_ids = new int[new_number_of_particles];

	return Simulation::addParticles(particles);
}

void SimulationVis::prepareColorScaling()
{
	int particles;
	double speed;

	max_neighbours = 1;
	max_speed = 0;
	avg_speed[0] = avg_speed[1] = avg_speed[2] = 0;
	max_diss_energy = 0;

	// determine max number of neighbours and average particle speed
	for(int p = 0; p < number_of_particles; ++p)
	{
		particles = grid.getNeighbourCount(p, vis_neighbour_search_dist);

		if(particles > max_neighbours)
			max_neighbours = particles;

		avg_speed[0] += vel[X_COORD(p)];
		avg_speed[1] += vel[Y_COORD(p)];
		avg_speed[2] += vel[Z_COORD(p)];
	}

	// update max/avg velocity
	if(number_of_particles > 0)
	{
		avg_speed[0] /= (double)number_of_particles;
		avg_speed[1] /= (double)number_of_particles;
		avg_speed[2] /= (double)number_of_particles;

		for(int p = 0; p < number_of_particles; ++p)
		{
			// determine speed of particle
			vec3 velocity;
			velocity[0]  = vel[X_COORD(p)] - avg_speed[0];
			velocity[1]  = vel[Y_COORD(p)] - avg_speed[1];
			velocity[2]  = vel[Z_COORD(p)] - avg_speed[2];

			speed = norm(velocity);
			speed -= 0.1 * v_c;

			if(speed > max_speed)
				max_speed = speed;
		}

		max_speed += 0.2 * v_c;
	}
	else
		max_speed = 0.2 * v_c;
}

float SimulationVis::getDepthColor(int particle_id)
{
	// determine number of neighbouring particles and update max. number
	int particles = grid.getNeighbourCount(particle_id, vis_neighbour_search_dist);

	if(particles > max_neighbours_new)
		max_neighbours_new = particles;

	return calcDepthColor(particles);
}

float SimulationVis::calcDepthColor(int number_of_neighbours)
{
	// apply scaling
	number_of_neighbours -= vis_neighbour_min_offset;

	if(number_of_neighbours < 0 || max_neighbours <= vis_neighbour_min_offset)
		return 1.0f;
	else if(number_of_neighbours >= (max_neighbours - vis_neighbour_min_offset + vis_neighbour_max_offset))
		return 0.2f;
	else
		return (1.0f - 0.8f * (float)(number_of_neighbours) / (float)(max_neighbours - vis_neighbour_min_offset + vis_neighbour_max_offset) );	
}

void SimulationVis::storeParticlePos(double *dest)
{
	size_t size = number_of_particles*3*sizeof(double);
	memcpy(dest, pos_old, size);
}

void SimulationVis::storeWalls(double *new_pos, double *normals, float *alpha_values)
{
	for(int w = 0; w < number_of_walls; ++w)
	{
		alpha_values[w] = walls[w].alpha;
		memcpy(&(new_pos[24*w]), &(walls[w].edges), 24 * sizeof(double));
		normals[3*w] = walls[w].normal[0];
		normals[3*w+1] = walls[w].normal[1];
		normals[3*w+2] = walls[w].normal[2];
	}
}

void SimulationVis::storeDissipatedContactEnergy(float *dest, bool depth) 
{
	double E_roll = 6.0 * M_PI * M_PI * surface_energy * reduced_radius * crit_rolling_displacement * ENERGY_UNIT;
    double scale = (dissipated_contact_energy + dissipated_rolling_energy + dissipated_twisting_energy) / ( (double)number_of_particles * E_roll);
	float color = 1.0f;
	double max_diss_energy_new = 0;

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
	double dissipation;
	float red, green;
#endif
	
	if(depth)
		max_neighbours_new = 1;

	if(scale < 0.1)
		scale = 0.1;
	else if(scale > 1.0)
		scale = 1.0;

	scale = 1.0 / (scale * E_roll);

	for(int p = 0; p < number_of_particles; ++p)
	{
		if(depth)
			color = getDepthColor(p);

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
		if(dissipated_energy_of_particle[p] > max_diss_energy_new)
			max_diss_energy_new = dissipated_energy_of_particle[p];

		if(max_diss_energy > 0)
			dissipation = 0.1 + 0.9 * dissipated_energy_of_particle[p] / max_diss_energy;
		else
			dissipation = 0;

		if(dissipation < 0.5)
		{
			red = 0.2f + 1.6f * dissipation;
			green = 1.0f;
		}
		else
		{
			red = 1.0f;
			green = 1.0f - 1.6f * (dissipation - 0.5);
		}

		dest[X_COORD_(p)] = red * color;
		dest[Y_COORD_(p)] = green * color;
		dest[Z_COORD_(p)] = 0.3f * color;
#else
		dest[X_COORD_(p)] = color;
		dest[Y_COORD_(p)] = color;
		dest[Z_COORD_(p)] = color;
#endif
	}

	max_diss_energy = max_diss_energy_new + 0.5 * E_roll;

	if(depth)
		max_neighbours = max_neighbours_new;
}

void SimulationVis::storeColorParticleContacts(float *dest, bool depth) 
{
	// determine number of neighbours of every particle
	float *contacts_coefficient = new float[number_of_particles];

	for(int p = 0; p < number_of_particles; ++p)
	{
		ContactListEntry *cl_entry = contact_list[p];

		while(cl_entry)
		{
			contacts_coefficient[p] += 0.1f;

			if(cl_entry->id >= 0)
				contacts_coefficient[cl_entry->id] += 0.1f;

			cl_entry = cl_entry->next;
		}
	}

	// determine color
	float color = 1.0f;

	if(depth)
		max_neighbours_new = 1;

	for(int p = 0; p < number_of_particles; ++p)
	{
		getColorCoordinationNumber(&dest[X_COORD_(p)], contacts_coefficient[p]); 

		if(depth)
		{
			color = getDepthColor(p);
			dest[X_COORD_(p)] *= color;
			dest[Y_COORD_(p)] *= color;
			dest[Z_COORD_(p)] *= color;
		}
	}

	delete [] contacts_coefficient;

	if(depth)
		max_neighbours = max_neighbours_new;
}

void SimulationVis::storeDensity(float *dest, bool depth) 
{
	float density;
	float color = 1.0f;
	int particles;

	// volume that is searched for particles
	double V_search = grid.x_width * grid.y_width *  grid.z_width * pow( (double) (2*vis_neighbour_search_dist+1), 3);
	double V_particle = 4.0/3.0 * M_PI * particle_radius*particle_radius*particle_radius;

	double volume_factor = V_particle / V_search;

	max_neighbours_new = 1;

	// maximum distance of another particle from a given cell to be taken into account
	double max_search_dist = (0.5 + (double)vis_neighbour_search_dist) * grid.x_width;
	std::vector<double> wall_dist(walls.size());
	int close_walls;
	vec3 cell_center;

	bool check_for_walls = (walls.size() > 0) && vis_density_consider_walls;

	for(int p = 0; p < number_of_particles; ++p)
	{
		particles = grid.getNeighbourCount(p, vis_neighbour_search_dist);

		if(depth)
			color = calcDepthColor(particles);

		if(particles > max_neighbours_new)
			max_neighbours_new = particles;

		double filling_factor = (double)particles * volume_factor;

		// if particles are close to walls the volume where neighbours may be found is reduced
		// without correction, particles close to walls will be assigned a lower filling factor
		if(check_for_walls)
		{
			// check if particle is close to a wall
			close_walls = 0;
			grid.getCellCenter(&cell_center, pos_old[X_COORD(p)], pos_old[Y_COORD(p)], pos_old[Z_COORD(p)]);

			for(int w = 0; w < walls.size(); ++w)
			{
				double dist = walls[w].getDistanceTo(cell_center);

				if(dist < max_search_dist)
				{
					wall_dist[close_walls] = dist;
					++close_walls;
				}
			}

			// correct filling factor
			for(int w = 0; w < close_walls; ++w)
				filling_factor *= (1.75 - 0.75 * wall_dist[w] / max_search_dist);
		}

		if(filling_factor < vis_density_max_filling_factor)
			density = (float) (filling_factor / vis_density_max_filling_factor);
		else
			density = 1.0f;

		// determine color
		getColorDensity(&dest[X_COORD_(p)], density);

		if(depth)
		{
			dest[X_COORD_(p)] *= color;
			dest[Y_COORD_(p)] *= color;
			dest[Z_COORD_(p)] *= color;
		}
	}

	max_neighbours = max_neighbours_new;
}

void SimulationVis::storeParticleSpeed(float *dest, bool depth)
{
	if(number_of_particles < 1)
		return;

	float color = 1.0f;
	double max_speed_new = 0;
	vec3 avg_speed_new;
	avg_speed_new[0] = avg_speed_new[1] = avg_speed_new[2] = 0;
	vec3 velocity;
	double speed;
	bool relative_speed = true;

	if(depth)
		max_neighbours_new = 1;

	if(sim_info.sim_type == SIM_TYPE_COMPRESSION_WITH_SIDE_WALLS || sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS || sim_info.sim_type == SIM_TYPE_WALL_COLLISION)
		relative_speed = false;

	for(int p = 0; p < number_of_particles; ++p)
	{
		if(depth)
			color = getDepthColor(p);

		// update velocity of cms
		avg_speed_new[0] += vel[X_COORD(p)];
		avg_speed_new[1] += vel[Y_COORD(p)];
		avg_speed_new[2] += vel[Z_COORD(p)];

		// determine speed of particle
		if(relative_speed)
		{
			velocity[0] = vel[X_COORD(p)] - avg_speed[0];
			velocity[1] = vel[Y_COORD(p)] - avg_speed[1];
			velocity[2] = vel[Z_COORD(p)] - avg_speed[2];

			speed = norm(velocity);
		}
		else
			speed = sqrt(vel[X_COORD(p)]*vel[X_COORD(p)] + vel[Y_COORD(p)]*vel[Y_COORD(p)] + vel[Z_COORD(p)]*vel[Z_COORD(p)]);

		// filter normal oscillations
		speed -= 0.1 * v_c;

		if(speed > max_speed_new)
			max_speed_new = speed;

		// scale to 0-1
		speed /= (max_speed+0.01);

		// particles below average speed will have default color
		if(speed > 0)
		{
			if(speed < 0.5)
			{
				dest[X_COORD_(p)] = (0.2f + 1.6f * (float)speed) * color;
				dest[Y_COORD_(p)] = color;
			}
			if(speed < 0.5)
			{
				dest[X_COORD_(p)] = (0.2f + 1.6f * (float)speed) * color;
				dest[Y_COORD_(p)] = color;
			}
			else if(speed < 1.0)
			{
				dest[X_COORD_(p)] = color;
				dest[Y_COORD_(p)] = (1.0f - 1.6f * ((float)speed - 0.5f)) * color;
			}
			else
			{
				dest[X_COORD_(p)] = color;
				dest[Y_COORD_(p)] = 0.2f * color;
			}

			dest[Z_COORD_(p)] = 0.3f * color;
		}
		else
		{
			dest[X_COORD_(p)] = color;
			dest[Y_COORD_(p)] = color;
			dest[Z_COORD_(p)] = color;
		}
	}

	if(depth)
		max_neighbours = max_neighbours_new;

	avg_speed_new[0] /= (double)number_of_particles;
	avg_speed_new[1] /= (double)number_of_particles;
	avg_speed_new[2] /= (double)number_of_particles;

	max_speed = max_speed_new + 0.2 * v_c;
	memcpy(avg_speed, avg_speed_new, sizeof(vec3));
}

void SimulationVis::storeColorWallContacts(float *dest, bool depth)
{
	ContactListEntry *cl_entry;
	float color = 1.0f;

	if(depth)
		max_neighbours_new = 1;

	for(int p = 0; p < number_of_particles; ++p)
	{
		if(depth)
			color = getDepthColor(p);

		dest[X_COORD_(p)] = color;
		dest[Y_COORD_(p)] = color;
		dest[Z_COORD_(p)] = color;

		cl_entry = contact_list[p];

		while(cl_entry)
		{
			if(cl_entry->id < 0)
			{
				if(sim_info.sim_type != SIM_TYPE_PULL_STRENGTH_TEST)
				{
					dest[Y_COORD_(p)] = 0.3f * color;
					dest[Z_COORD_(p)] = 0.3f * color;
					break;
				}
				else if(box && (WALL_ID(cl_entry->id) == box->top_wall_id || WALL_ID(cl_entry->id) == box->bottom_wall_id))
				{
					dest[Y_COORD_(p)] = 0.3f * color;
					dest[Z_COORD_(p)] = 0.3f * color;
					break;
				}
			}

			cl_entry = cl_entry->next;
		}

#ifdef ENABLE_WALL_GLUE
		if(sim_info.sim_type == SIM_TYPE_PULL_STRENGTH_TEST || sim_info.sim_type == SIM_TYPE_SHEAR_STRENGTH_TEST)
		{
			double wall_dist = 2.0 * particle_radius * wall_glue_distance;

			double temp = walls[box->top_wall_id].getDistanceTo(pos_new[X_COORD(p)], pos_new[Y_COORD(p)], pos_new[Z_COORD(p)]);
			if(temp < wall_dist)
				wall_dist = temp;

			temp = walls[box->bottom_wall_id].getDistanceTo(pos_new[X_COORD(p)], pos_new[Y_COORD(p)], pos_new[Z_COORD(p)]);
			if(temp < wall_dist)
				wall_dist = temp;

			if(wall_dist < particle_radius * wall_glue_distance)
			{
				double k = (particle_radius * wall_glue_distance - wall_dist) / (particle_radius * wall_glue_distance);

				dest[Y_COORD_(p)] = (float)(1.0 - 0.7 * k) * color;
				dest[Z_COORD_(p)] = (float)(1.0 - 0.7 * k) * color;
			}
		}
#endif
	}

	if(depth)
		max_neighbours = max_neighbours_new;
}

void SimulationVis::storeWallAngle(float *dest, bool depth)
{
	if(walls.size() == 0)
	{
		storeDepth(dest, depth);
		return;
	}

	float color = 1.0f;
	float angle;

	if(depth)
		max_neighbours_new = 1;

	for(int p = 0; p < number_of_particles; ++p)
	{
		if(depth)
			color = getDepthColor(p);

		// calculate angle between velocity and wall
		angle = (float) acos( ( walls[0].normal[0] * vel[X_COORD(p)] + walls[0].normal[1] * vel[Y_COORD(p)] + walls[0].normal[2] * vel[Z_COORD(p)]) * sqrt(vel[X_COORD(p)]*vel[X_COORD(p)] + vel[Y_COORD(p)]*vel[Y_COORD(p)] + vel[Z_COORD(p)]*vel[Z_COORD(p)]) );
		angle /= (float)M_PI;

		if(angle < 0.25f)
		{
			dest[X_COORD_(p)] = (0.2f + 3.2f * angle) * color;
			dest[Y_COORD_(p)] = color;
			dest[Z_COORD_(p)] = 0.3f * color;
		}
		else if(angle < 0.5f)
		{
			dest[X_COORD_(p)] = color;
			dest[Y_COORD_(p)] = (1.0f - 3.2f * (angle - 0.25f)) * color;
			dest[Z_COORD_(p)] = 0.3f * color;
		}
		else // use standard color if particle is moving towards the wall
		{
			dest[X_COORD_(p)] = color;
			dest[Y_COORD_(p)] = color;
			dest[Z_COORD_(p)] = color;
		}
	}

	if(depth)
		max_neighbours = max_neighbours_new;
}

void SimulationVis::storeFragments(float *dest, bool depth)
{
	float color = 1.0f;
	int color_id = 0;

	if(depth)
		max_neighbours_new = 1;
	
	for(int p = 0; p < number_of_particles; ++p)
	{
		if(depth)
			color = getDepthColor(p);

		if(fragment_ids != NULL)
			color_id = (fragment_ids[p] - fragment_id_offset)%FRAGMENT_COLORS;

		dest[X_COORD_(p)] = color * FRAGMENT_COLOR_R[color_id];
		dest[Y_COORD_(p)] = color * FRAGMENT_COLOR_G[color_id];
		dest[Z_COORD_(p)] = color * FRAGMENT_COLOR_B[color_id];
	}

	if(depth)
		max_neighbours = max_neighbours_new;
}

void SimulationVis::storeDislocation(float *dest, bool depth)
{
	if(depth)
		max_neighbours_new = 1;

	double dist;
	double min_dist = (vis_dislocation_min_value * delta_0) * (vis_dislocation_min_value * delta_0);
	double max_dist = (vis_dislocation_max_value * particle_radius) * (vis_dislocation_max_value * particle_radius);

	for(int p = 0; p < number_of_particles; ++p)
	{
		dist = (initial_positions[X_COORD(p)] - pos_old[X_COORD(p)]) * (initial_positions[X_COORD(p)] - pos_old[X_COORD(p)])
			+ (initial_positions[Y_COORD(p)] - pos_old[Y_COORD(p)]) * (initial_positions[Y_COORD(p)] - pos_old[Y_COORD(p)]) 
			+ (initial_positions[Z_COORD(p)] - pos_old[Z_COORD(p)]) * (initial_positions[Z_COORD(p)] - pos_old[Z_COORD(p)]);

		dist -= min_dist;
		dist /= (max_dist - min_dist);
		getColorDislocation(&dest[X_COORD(p)], dist);

		if(depth)
		{
			float color = getDepthColor(p);
			dest[X_COORD_(p)] *= color;
			dest[Y_COORD_(p)] *= color;
			dest[Z_COORD_(p)] *= color;
		}
	}

	if(depth)
		max_neighbours = max_neighbours_new;
}

void SimulationVis::store50(float *dest, bool depth)
{
    if(depth)
        max_neighbours_new = 1;

    for(int p = 0; p < number_of_particles; ++p)
    {
        getColor50(&dest[X_COORD(p)], p, number_of_particles);

        if(depth)
        {
            float color = getDepthColor(p);
            dest[X_COORD_(p)] *= color;
            dest[Y_COORD_(p)] *= color;
            dest[Z_COORD_(p)] *= color;
        }
    }

    if(depth)
        max_neighbours = max_neighbours_new;
}

void SimulationVis::storeDepth(float *dest, bool depth)
{
	float color = 1.0f;

	if(depth)
		max_neighbours_new = 1;

	for(int p = 0; p < number_of_particles; ++p)
	{
		if(depth)
			color = getDepthColor(p);

		dest[X_COORD_(p)] = color;
		dest[Y_COORD_(p)] = color;
		dest[Z_COORD_(p)] = color;
	}

	if(depth)
		max_neighbours = max_neighbours_new;
}

////////////////////////////////////////////////////////////////////////////////////////
// color functions 
////////////////////////////////////////////////////////////////////////////////////////

void SimulationVis::getColorDensity(float *color, float density)
{
	if(density < 0.5f)
	{
		color[0] = (0.2f + 1.6f * density);
		color[1] = 1.0f;
		color[2] = 0.3f;
	}
	else
	{
		color[0] = 1.0f;
		color[1] = (1.8f - 1.6f * density);
		color[2] = 0.3f;
	}
}

void SimulationVis::getColorCoordinationNumber(float *color, float k)
{
	if(k == 0)
	{
		color[0] = 1.0f;
		color[1] = 1.0f;
		color[2] = 1.0f;
	}
	else if(k < 0.5f)
	{
		color[0] = (0.2f + 1.6f * k);
		color[1] = 1.0f;
		color[2] = 0.3f;
	}
	else if(k < 1.0f)
	{
		color[0] = 1.0f;
		color[1] = (1.8f - 1.6f * k);
		color[2] = 0.3f;
	}
	else
	{
		color[0] = 1.0f;
		color[1] = 0.2f;
		color[2] = 0.3f;
	}
}

void SimulationVis::getColorDislocation(float *color, float k)
{
	if(k <= 0)
	{
		color[0] = 1.0f;
		color[1] = 1.0f;
		color[2] = 1.0f;
	}
	else if(k < 0.5f)
	{
		color[0] = (0.2f + 1.6f * k);
		color[1] = 1.0f;
		color[2] = 0.3f;
	}
	else if(k < 1.0f)
	{
		color[0] = 1.0f;
		color[1] = (1.8f - 1.6f * k);
		color[2] = 0.3f;
	}
	else
	{
		color[0] = 1.0f;
		color[1] = 0.2f;
		color[2] = 0.3f;
	}
}


void SimulationVis::getColor50(float *color, int i, int N)
{
    if(i < N/2)
    {
        color[0] = 1.0f;
        color[1] = 0.3f;
        color[2] = 0.3f;
    }
    else
    {
        color[0] = 0.3f;
        color[1] = 0.3f;
        color[2] = 1.0f;
    }
}


void SimulationVis::initPressureAveraging(int averaging_steps)
{
	if(forces_storage)
		delete [] forces_storage;

	forces_storage = new double[averaging_steps];
	memset(forces_storage, 0, averaging_steps * sizeof(double));
	force_avg_size = averaging_steps;
	force_avg_counter = 0;
}

double SimulationVis::getAvgWallPressure()
{
	if(box && force_avg_size > 0)
	{
		double avg_force = 0;
		for(int i = 0; i < force_avg_size; ++i)
			avg_force += forces_storage[i];

		avg_force /= (double)force_avg_size;

		// 0.1 to convert from CGS to SI
		double pressure = 0.1 * avg_force;
			
		if(sim_info.sim_type == SIM_TYPE_SHOCKWAVE || sim_info.sim_type == SIM_TYPE_OPEN_BOX)
			pressure /= box->base;
		else if(sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
        {
            double cross_section = 0.0;
             SimLib::getCrossSectionNoRotation(*this, 0.2*particle_radius, cross_section);

             if(cross_section != 0.0)
                pressure /= cross_section;
             else
                pressure = 0.0;
        }
		else
			pressure /= box->base;
			
		return pressure;
	}

	return -1.0;
}

double SimulationVis::getAvgWallForce()
{
	if(box && force_avg_size > 0)
	{
		double avg_force = 0;
		for(int i = 0; i < force_avg_size; ++i)
			avg_force += forces_storage[i];

		avg_force /= (double)force_avg_size;

		// convert from CGS to SI
		avg_force *= 1e-5;

		return avg_force;
	}

	return -1.0;
}
