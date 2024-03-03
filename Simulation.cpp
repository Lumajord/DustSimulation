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

#include "Simulation.h"
#include <cfloat>
#include <cstring>
#include <ctime>
#include <list>
#include "SimulationLib.h"

#include <locale.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(USE_JKR) && defined(USE_DMT)
#error JKR and DMT interaction model enabled - select only one!
#endif

#if defined(USE_CONSTANT_DAMPING) && defined(USE_VISCOELASTIC_DAMPING)
#error Both constant and viscoelastic damping enabled - select only one!
#endif

// global variables that do not change during the sim (material parameters, etc.)
bool material_initialized = false; // prevent multiple instances from reinitializing the material constants
char material_name[200];
double particle_radius;		// in cm
double particle_radius_inv;	// in cm^-1
double reduced_radius;		// in cm

// material parameters
double surface_energy;		// surface energy / area (in mJ / m^2)
double nu;					// poisson number
double young_mod;			// Young's modulus (in Ba = g / (cm s^2) )
double young_mod_reduced;
double shear_mod;			// shear modulus (in Ba = g / (cm s^2) )
double shear_mod_reduced;
double bulk_mod;			// bulk modulus (in Ba = g / (cm s^2) )
double density;				// in g/cm^3
double mass;				// in g
double mass_inv;			// in g^-1
double moment_of_inertia;	// in g cm^2
double moment_of_inertia_inv;
double osc_damping_factor;	// additional damping for normal oscillations
double T_vis;
double viscoelastic_damping_constant;  // viscoelastic damping constant as proposed by Sebastian Krijt
double viscous_constant = 8e-10;

double equilibrium_contact_radius;	// radius of contact area when no normal is acting on the particles
double equilibrium_distance;		// distance when no normal is acting on the particles
double delta_0;				// in cm
double delta_c;				// in cm
double F_c;					// in 10^-5 N = 10^-2 mJ/m
double t_c;					// characteristic time (in s)
double T_normal;			// period of normal oscillations (in s)
double bulk_sound_speed;	// sound speed of the bulk material (in cm/s)
double v_c;					// characteristic stiking velocity (in cm/s)
double ENERGY_UNIT;

double k_s;	// constant for sliding potential
double k_r;	// constant for rolling potential
double k_t;	// constant for twisting potential

double contact_breaking_dist;			// in cm
double contact_breaking_dist_squared;	// in cm^2
double contact_making_dist;				// in cm
double contact_making_dist_squared;		// in cm^2

double crit_rolling_displacement;			// in cm
double crit_rolling_displacement_squared;	// in cm^2
double crit_sliding_displacement;
double crit_sliding_displacement_squared;
double crit_twisting_displacement;
double crit_twisting_displacement_squared;

// for wall interaction
double wall_reduced_radius;
double wall_equilibrium_contact_radius;
double wall_delta_c;
double wall_delta_0;
double wall_contact_breaking_dist;
double wall_F_c;
double wall_k_s;
double wall_k_r;
double wall_k_t;
double crit_wall_sliding_displacement;
double crit_wall_sliding_displacement_squared;
double wall_acceleration;
double wall_inertia = 1.0;

double azimuthal_acceleration = 0.0;

// the following modifiers can used to tweak various aspects of the D&T particle-particle interaction
double rolling_modifier = 1.0;
double sliding_modifier = 1.0;
double twisting_modifier = 1.0;
double crit_sliding_displacement_modifier = 1.0;
double crit_wall_sliding_displacement_modifier = 1.0;

// direction and strength of gravity
#ifdef ENABLE_GRAVITY
    double gravity_strength = 9.81e2; // in cm/s^2
    vec3 gravity_direction = {0.0, -1.0, 0.0};
    bool gravity_enabled = false;
    double gravity_modifier = 1.0;
#endif

double damping_factor = 0.7;
bool damping_enabled = false;


// returns a random number in [0.0, 2*pi)
double Simulation::get_random_zero_twoPi()
{
    return dist_zero_twoPi(rand_generator);
}
// returns a random number in [0.0, 1.0)
double Simulation::get_random_zero_one()
{
    return dist_zero_one(rand_generator);
}
// returns a random number in [-1.0, 1.0]
double Simulation::get_random_cos_theta()
{
    return dist_cos_theta(rand_generator);
}
// returns a random number in [0.0, 1.0]
double Simulation::get_random_zero_one_incl()
{
    return dist_zero_one_incl(rand_generator);
}

Simulation::Simulation(void)
{
    initialized = false;

    rand_generator.seed(time(NULL));
    dist_cos_theta = std::uniform_real_distribution<double>(-1.0, std::nextafter(1.0, DBL_MAX));	// [-1, 1]
    dist_zero_twoPi = std::uniform_real_distribution<double>(0.0, 2.0*M_PI);						// [0, 2*pi)
    dist_zero_one = std::uniform_real_distribution<double>(0.0, 1.0);								// [0, 1)
    dist_zero_one_incl = std::uniform_real_distribution<double>(0.0, std::nextafter(1.0, DBL_MAX));	// [0, 1]

    timestep = 0;
    current_time = 0;
    end_time = 0;
    min_end_time = 0;

    number_of_particles = 0;
    number_of_contacts = 0;
    number_of_walls = 0;
    number_of_wall_contacts = 0;

    initialized = false;


#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
    dissipated_energy_of_particle = NULL;
#endif

#ifdef TRACK_CONTACTS_PER_PARTICLE
    initial_number_of_contacts_of_particle = NULL;
    broken_contacts_of_particle = NULL;
    created_contacts_of_particle = NULL;
#endif

    pos_old = NULL;
    pos_new = NULL;
    vel = NULL;
    vel_angular = NULL;
    force_new = NULL;
    force_old = NULL;
    torque_new = NULL;
    torque_old = NULL;

#ifdef TRACK_PARTICLE_ORIENTATION
    orientation = NULL;
#endif

#ifdef ENABLE_FIXED_PARTICLES
    fixed_particles = NULL;
#endif

    contact_list = NULL;
    energies_filename = NULL;
    positions_path = NULL;

    // init default material if something went wrong
    if(!material_initialized)
    {
        setMaterialConstants(6e-5, 2.65, 25.0, 0.17, 5.4e11, 2e-7, 1e-6, 1e-11);
        strcpy(material_name, "Silicate");

#ifdef INTERPOLATE_NORMAL_FORCES
        normal_interaction.init(NORMAL_FORCES_INTERPOLATION_POINTS);
#else
        normal_interaction.init();
#endif

        material_initialized = true;
    }

    box = NULL;

    sim_info.info_storage[0] = 0.0;
    sim_info.info_storage[1] = 0.0;
    sim_info.info_storage[2] = 0.0;
    sim_info.info_storage[3] = 0.0;
    sim_info.info_storage[4] = 0.0;
    sim_info.info_storage[5] = 0.0;
    sim_info.info_storage[6] = 0.0;
    sim_info.info_storage[7] = 0.0;

    sim_info.sim_type = SIM_TYPE_GENERAL;

    check_potential_variation_interval = 0;
    check_potential_variation_counter = 0;
    potential_variation_stop_threshold = 0;
    last_dissipated_energy = 0;
}

Simulation::~Simulation(void)
{
    cleanUp();
}

void Simulation::cleanUp()
{

    if(box)
    {
        delete box;
        box = NULL;
    }

    deleteContacts();


    delete [] contact_list;
    contact_list = NULL;

    number_of_particles = 0;

    // clean up arrays
    delete [] pos_old;
    delete [] pos_new;
    delete [] vel;
    delete [] vel_angular;
    delete [] force_new;
    delete [] force_old;
    delete [] torque_new;
    delete [] torque_old;

    pos_old = NULL;
    pos_new = NULL;
    vel = NULL;
    vel_angular = NULL;
    force_new = NULL;
    force_old = NULL;
    torque_new = NULL;
    torque_old = NULL;


#ifdef TRACK_PARTICLE_ORIENTATION
    delete [] orientation;
    orientation = NULL;
#endif


#ifdef ENABLE_FIXED_PARTICLES
    if(fixed_particles)
    {
        delete [] fixed_particles;
        fixed_particles = NULL;
    }
#endif


    if(energies_filename)
    {
        delete [] energies_filename;
        energies_filename = NULL;
    }


    if(positions_path)
    {
        delete [] positions_path;
        positions_path = NULL;
    }


#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
    if(dissipated_energy_of_particle)
    {
        delete [] dissipated_energy_of_particle;
        dissipated_energy_of_particle = NULL;
    }
#endif


#ifdef TRACK_CONTACTS_PER_PARTICLE
    if(initial_number_of_contacts_of_particle)
    {
        delete [] initial_number_of_contacts_of_particle;
        initial_number_of_contacts_of_particle = NULL;
    }

    if(broken_contacts_of_particle)
    {
        delete [] broken_contacts_of_particle;
        broken_contacts_of_particle = NULL;
    }

    if(created_contacts_of_particle)
    {
        delete [] created_contacts_of_particle;
        created_contacts_of_particle = NULL;
    }
#endif

    walls.clear();
    number_of_walls = 0;

    current_time = 0;
    end_time = 0;
    stop_simulation = true;
}

void Simulation::deleteContacts()
{
    if(contact_list)
    {
        ContactListEntry *cl_entry, *last_entry;

        for(int p = 0; p < number_of_particles; ++p)
        {
            cl_entry = contact_list[p];

            while(cl_entry)
            {
                last_entry = cl_entry;
                cl_entry = cl_entry->next;

                delete last_entry->contact;
                delete last_entry;
            }

            contact_list[p] = NULL;
        }
    }

    number_of_contacts = 0;
    number_of_wall_contacts = 0;
}

void Simulation::resetContacts()
{
    if(contact_list)
    {
        ContactListEntry *cl_entry;

        for(int p = 0; p < number_of_particles; ++p)
        {
            cl_entry = contact_list[p];

            while(cl_entry)
            {


                int id1 = cl_entry->contact->id1;
                int id2 = cl_entry->contact->id2;
                vec3 n2_initial;
                n2_initial[0] = pos_new[X_COORD(id1)] - pos_new[X_COORD(id2)];
                n2_initial[1] = pos_new[Y_COORD(id1)] - pos_new[Y_COORD(id2)];
                n2_initial[2] = pos_new[Z_COORD(id1)] - pos_new[Z_COORD(id2)];

                normalize(&n2_initial);


                cl_entry->contact->n1_initial[0] = -n2_initial[0];
                cl_entry->contact->n1_initial[1] = -n2_initial[1];
                cl_entry->contact->n1_initial[2] = -n2_initial[2];

                cl_entry->contact->n2_initial[0] = n2_initial[0];
                cl_entry->contact->n2_initial[1] = n2_initial[1];
                cl_entry->contact->n2_initial[2] = n2_initial[2];

                cl_entry->contact->old_contact_normal[0] = n2_initial[0];
                cl_entry->contact->old_contact_normal[1] = n2_initial[1];
                cl_entry->contact->old_contact_normal[2] = n2_initial[2];

                cl_entry->contact->compression_length = 0;
                cl_entry->contact->twisting_displacement = 0;

                cl_entry = cl_entry->next;

            }
        }
    }

}


ErrorCode Simulation::resizeArrays(int new_number_of_particles, int new_number_of_walls)
{
    // clean up first
    cleanUp();

    number_of_particles = new_number_of_particles;
    number_of_walls = new_number_of_walls;

    // allocate memory
    pos_new = new double[3 * number_of_particles];
    pos_old = new double[3 * number_of_particles];
    vel = new double[3 * number_of_particles];
    vel_angular = new double[3 * number_of_particles];
    force_new = new double[3 * number_of_particles];
    force_old = new double[3 * number_of_particles];
    torque_new = new double[3 * number_of_particles];
    torque_old = new double[3 * number_of_particles];

#ifdef TRACK_PARTICLE_ORIENTATION
    orientation = new double[4 * number_of_particles];
    memset(orientation, 0, 4 * sizeof(double) * number_of_particles);
#endif

    // set default values
    memset(pos_new, 0, 3 * sizeof(double) * number_of_particles);
    memset(pos_old, 0, 3 * sizeof(double) * number_of_particles);
    memset(vel, 0, 3 * sizeof(double) * number_of_particles);
    memset(vel_angular, 0, 3 * sizeof(double) * number_of_particles);
    memset(force_new, 0, 3 * sizeof(double) * number_of_particles);
    memset(force_old, 0, 3 * sizeof(double) * number_of_particles);
    memset(torque_new, 0, 3 * sizeof(double) * number_of_particles);
    memset(torque_old, 0, 3 * sizeof(double) * number_of_particles);

    // init new contact list
    contact_list = new ContactListEntry*[number_of_particles];
    memset(contact_list, 0, number_of_particles * sizeof(ContactListEntry*));

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
    dissipated_energy_of_particle = new double[number_of_particles];
    memset(dissipated_energy_of_particle, 0.0, sizeof(double) * number_of_particles);
#endif

#ifdef TRACK_CONTACTS_PER_PARTICLE
    initial_number_of_contacts_of_particle = new int[number_of_particles];
    memset(initial_number_of_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    broken_contacts_of_particle = new int[number_of_particles];
    memset(broken_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    created_contacts_of_particle = new int[number_of_particles];
    memset(created_contacts_of_particle, 0, sizeof(int) * number_of_particles);
#endif

    grid.init(GRID_CELLS, 2.1 * particle_radius, number_of_particles);

#ifdef TRACK_PARTICLE_ORIENTATION
    initParticleOrientation();
#endif

    walls.resize(number_of_walls);
    return EC_OK;
}

ErrorCode Simulation::addParticlesFromSim(Simulation &sim)
{
    if(number_of_walls > 0 || sim.number_of_walls > 0)
        return EC_INVALID_BOUNDARIES;

    int new_number_of_particles = number_of_particles + sim.number_of_particles;

    // create new arrays
    double *n_pos_old = new double[3 * new_number_of_particles];
    double *n_pos_new = new double[3 * new_number_of_particles];
    double *n_vel = new double[3 * new_number_of_particles];
    double *n_vel_angular = new double[3 * new_number_of_particles];
    double *n_force_old = new double[3 * new_number_of_particles];
    double *n_torque_old = new double[3 * new_number_of_particles];
    double *n_force_new = new double[3 * new_number_of_particles];
    double *n_torque_new = new double[3 * new_number_of_particles];

#ifdef TRACK_PARTICLE_ORIENTATION
    double *n_orientation = new double[4 * new_number_of_particles];
#endif

#ifdef ENABLE_FIXED_PARTICLES
    bool *n_fixed_particles = NULL;

    if(fixed_particles || sim.fixed_particles)
    {
        n_fixed_particles = new bool[new_number_of_particles];
        memset(n_fixed_particles, 0, new_number_of_particles * sizeof(bool));

        if(fixed_particles)
            memcpy(n_fixed_particles, fixed_particles, number_of_particles * sizeof(bool));

        if(sim.fixed_particles)
            memcpy(n_fixed_particles+number_of_particles, sim.fixed_particles, sim.number_of_particles * sizeof(bool));
    }
#endif

    // copy data of all existing particles
    if(number_of_particles > 0)
    {
        memcpy(n_pos_old, pos_old, number_of_particles * sizeof(vec3) );
        memcpy(n_pos_new, pos_new, number_of_particles * sizeof(vec3) );
        memcpy(n_vel, vel, number_of_particles * sizeof(vec3) );
        memcpy(n_vel_angular, vel_angular, number_of_particles * sizeof(vec3) );
        memcpy(n_force_old, force_old, number_of_particles * sizeof(vec3) );
        memcpy(n_force_new, force_new, number_of_particles * sizeof(vec3) );
        memcpy(n_torque_old, torque_old, number_of_particles * sizeof(vec3) );
        memcpy(n_torque_new, torque_new, number_of_particles * sizeof(vec3) );

#ifdef TRACK_PARTICLE_ORIENTATION
        memcpy(n_orientation, orientation, 4 * number_of_particles * sizeof(double) );
#endif
    }

    // append new particles
    memcpy(n_pos_old+3*number_of_particles, sim.pos_old, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_pos_new+3*number_of_particles, sim.pos_new, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_vel+3*number_of_particles, sim.vel, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_vel_angular+3*number_of_particles, sim.vel_angular, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_force_old+3*number_of_particles, sim.force_old, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_force_new+3*number_of_particles, sim.force_new, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_torque_old+3*number_of_particles, sim.torque_old, sim.number_of_particles * sizeof(vec3) );
    memcpy(n_torque_new+3*number_of_particles, sim.torque_new, sim.number_of_particles * sizeof(vec3) );

#ifdef TRACK_PARTICLE_ORIENTATION
    memcpy(n_orientation+4*number_of_particles, sim.orientation, sim.number_of_particles * 4 * sizeof(double) );
#endif

    // new contact list
    ContactListEntry **n_contact_list = new ContactListEntry*[new_number_of_particles];

    // copy existing contacts
    for(int p = 0; p < number_of_particles; ++p)
        n_contact_list[p] = contact_list[p];

    // prevent contact list entries from being deleted when calling cleanUp() later
    delete [] contact_list;
    contact_list = NULL;

    // add new contacts (copy contact list and adjust particle ids)
    ContactListEntry *new_cl_entry, *cl_entry;

    for(int p = 0; p < sim.number_of_particles; ++p)
    {
        // start with empty contact list of current particle
        n_contact_list[number_of_particles + p] = NULL;

        new_cl_entry = NULL;
        cl_entry = sim.contact_list[p];

        while(cl_entry)
        {
            Contact *new_contact = new Contact();
            *new_contact = *(cl_entry->contact);

            new_contact->id1 += number_of_particles;
            new_contact->id2 += number_of_particles;

            if(n_contact_list[number_of_particles + p] == NULL) // particle has no other contacts yet
            {
                // create contact list entry
                n_contact_list[number_of_particles + p] = new ContactListEntry;
                n_contact_list[number_of_particles + p]->next = NULL;
                n_contact_list[number_of_particles + p]->id = cl_entry->id + number_of_particles;
                n_contact_list[number_of_particles + p]->contact = new_contact;
                new_cl_entry = n_contact_list[number_of_particles + p];
            }
            else
            {
                new_cl_entry->next = new ContactListEntry;
                new_cl_entry->next->next = NULL;
                new_cl_entry->next->id = cl_entry->id + number_of_particles;
                new_cl_entry->next->contact = new_contact;
                new_cl_entry = new_cl_entry->next;
            }

            cl_entry = cl_entry->next;
        }
    }

    // clean up
    cleanUp();

    // replace
    pos_old = n_pos_old;
    pos_new = n_pos_new;
    vel = n_vel;
    vel_angular = n_vel_angular;
    force_old = n_force_old;
    force_new = n_force_new;
    torque_old = n_torque_old;
    torque_new = n_torque_new;
    contact_list = n_contact_list;

#ifdef TRACK_PARTICLE_ORIENTATION
    orientation = n_orientation;
#endif

#ifdef ENABLE_FIXED_PARTICLES
    fixed_particles = n_fixed_particles;
#endif

    // set new number of particles and init grid
    number_of_particles = new_number_of_particles;
    number_of_contacts = getNumberOfContacts();
    number_of_wall_contacts = getNumberOfWallContacts();

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
    dissipated_energy_of_particle = new double[number_of_particles];
    memset(dissipated_energy_of_particle, 0.0, sizeof(double) * number_of_particles);
#endif

#ifdef TRACK_CONTACTS_PER_PARTICLE
    initial_number_of_contacts_of_particle = new int[number_of_particles];
    memset(initial_number_of_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    broken_contacts_of_particle = new int[number_of_particles];
    memset(broken_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    created_contacts_of_particle = new int[number_of_particles];
    memset(created_contacts_of_particle, 0, sizeof(int) * number_of_particles);
#endif

    grid.init(256, 2.1 * particle_radius, number_of_particles);
    grid.addParticles(pos_old, number_of_particles);

    walls.resize(0);

    return EC_OK;
}

ErrorCode Simulation::addParticles(int particles)
{
    int new_number_of_particles = number_of_particles + particles;

    // create new arrays
    double *n_pos_old = new double[3 * new_number_of_particles];
    double *n_pos_new = new double[3 * new_number_of_particles];
    double *n_vel = new double[3 * new_number_of_particles];
    double *n_vel_angular = new double[3 * new_number_of_particles];
    double *n_force_old = new double[3 * new_number_of_particles];
    double *n_torque_old = new double[3 * new_number_of_particles];
    double *n_force_new = new double[3 * new_number_of_particles];
    double *n_torque_new = new double[3 * new_number_of_particles];

    // set default values
    memset(n_vel, 0, 3 * sizeof(double) * new_number_of_particles);
    memset(n_vel_angular, 0, 3 * sizeof(double) * new_number_of_particles);
    memset(n_force_new, 0, 3 * sizeof(double) * new_number_of_particles);
    memset(n_force_old, 0, 3 * sizeof(double) * new_number_of_particles);
    memset(n_torque_new, 0, 3 * sizeof(double) * new_number_of_particles);
    memset(n_torque_old, 0, 3 * sizeof(double) * new_number_of_particles);

#ifdef TRACK_PARTICLE_ORIENTATION
    double *n_orientation = new double[4 * new_number_of_particles];
#endif

#ifdef ENABLE_FIXED_PARTICLES
    bool *n_fixed_particles = NULL;

    if(fixed_particles)
    {
        n_fixed_particles = new bool[new_number_of_particles];
        memset(n_fixed_particles, 0, new_number_of_particles * sizeof(bool));

        memcpy(n_fixed_particles, fixed_particles, number_of_particles * sizeof(bool));
    }
#endif

    // copy data of all existing particles
    memcpy(n_pos_old, pos_old, number_of_particles * sizeof(vec3) );
    memcpy(n_pos_new, pos_new, number_of_particles * sizeof(vec3) );
    memcpy(n_vel, vel, number_of_particles * sizeof(vec3) );
    memcpy(n_vel_angular, vel_angular, number_of_particles * sizeof(vec3) );
    memcpy(n_force_old, force_old, number_of_particles * sizeof(vec3) );
    memcpy(n_force_new, force_new, number_of_particles * sizeof(vec3) );
    memcpy(n_torque_old, torque_old, number_of_particles * sizeof(vec3) );
    memcpy(n_torque_new, torque_new, number_of_particles * sizeof(vec3) );

    // new contact list
    ContactListEntry **n_contact_list = new ContactListEntry*[new_number_of_particles];

    // copy existing contacts
    for(int p = 0; p < number_of_particles; ++p)
        n_contact_list[p] = contact_list[p];

    // prevent contact list entries from being deleted when calling cleanUp() later
    delete [] contact_list;
    contact_list = NULL;

    // start with empty contact list of current particle
    for(int p = 0; p < particles; ++p)
        n_contact_list[number_of_particles + p] = NULL;

    // save walls & box info
    std::vector<Wall> n_walls;
    n_walls = walls;

    WallBox *n_box = box;
    box = NULL;

    // clean up
    cleanUp();

    // replace
    pos_old = n_pos_old;
    pos_new = n_pos_new;
    vel = n_vel;
    vel_angular = n_vel_angular;
    force_old = n_force_old;
    force_new = n_force_new;
    torque_old = n_torque_old;
    torque_new = n_torque_new;
    contact_list = n_contact_list;

#ifdef TRACK_PARTICLE_ORIENTATION
    orientation = n_orientation;
#endif

#ifdef ENABLE_FIXED_PARTICLES
    fixed_particles = n_fixed_particles;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
    dissipated_energy_of_particle = new double[new_number_of_particles];
    memset(dissipated_energy_of_particle, 0.0, sizeof(double) * new_number_of_particles);
#endif

#ifdef TRACK_CONTACTS_PER_PARTICLE
    initial_number_of_contacts_of_particle = new int[new_number_of_particles];
    memset(initial_number_of_contacts_of_particle, 0, sizeof(int) * new_number_of_particles);

    broken_contacts_of_particle = new int[new_number_of_particles];
    memset(broken_contacts_of_particle, 0, sizeof(int) * new_number_of_particles);

    created_contacts_of_particle = new int[new_number_of_particles];
    memset(created_contacts_of_particle, 0, sizeof(int) * new_number_of_particles);
#endif

    grid.init(256, 2.1 * particle_radius, new_number_of_particles);
    vec3 pos;
    for(int i = 0; i < new_number_of_particles - particles; ++i)
    {
        pos[0] = pos_old[X_COORD(i)];
        pos[1] = pos_old[Y_COORD(i)];
        pos[2] = pos_old[Z_COORD(i)];
        grid.addParticle(pos, i);
    }

    walls = n_walls;
    number_of_walls = walls.size();

    box = n_box;

    // set new number of particles
    number_of_particles = new_number_of_particles;
    number_of_contacts = getNumberOfContacts();
    number_of_wall_contacts = getNumberOfWallContacts();

#ifdef TRACK_PARTICLE_ORIENTATION
    initParticleOrientation();
#endif

    return EC_OK;
}

void Simulation::removeParticles(bool *remove_list, int removed_particles)
{
    int new_number_of_particles = number_of_particles - removed_particles;

    int *new_id = new int[number_of_particles];

    // create new arrays
    double *n_pos_old = new double[3 * new_number_of_particles ];
    double *n_pos_new = new double[3 * new_number_of_particles ];
    double *n_vel = new double[3 * new_number_of_particles ];
    double *n_vel_angular = new double[3 * new_number_of_particles ];
    double *n_force_old = new double[3 * new_number_of_particles ];
    double *n_torque_old = new double[3 * new_number_of_particles ];
    double *n_force_new = new double[3 * new_number_of_particles ];
    double *n_torque_new = new double[3 * new_number_of_particles ];

#ifdef TRACK_PARTICLE_ORIENTATION
    double *n_orientation = new double[4 * new_number_of_particles ];
#endif

#ifdef ENABLE_FIXED_PARTICLES
    bool *n_fixed_particles = NULL;

    if(fixed_particles)
        n_fixed_particles = new bool[new_number_of_particles];
#endif

    // copy data of all particles that are kept to the new arrays
    int next_id = 0;

    for(int p = 0; p < number_of_particles; ++p)
    {
        // copy particle
        if(!remove_list[p])
        {
            n_pos_old[X_COORD(next_id)] = pos_old[X_COORD(p)];
            n_pos_old[Y_COORD(next_id)] = pos_old[Y_COORD(p)];
            n_pos_old[Z_COORD(next_id)] = pos_old[Z_COORD(p)];

            n_pos_new[X_COORD(next_id)] = pos_new[X_COORD(p)];
            n_pos_new[Y_COORD(next_id)] = pos_new[Y_COORD(p)];
            n_pos_new[Z_COORD(next_id)] = pos_new[Z_COORD(p)];

            n_vel[X_COORD(next_id)] = vel[X_COORD(p)];
            n_vel[Y_COORD(next_id)] = vel[Y_COORD(p)];
            n_vel[Z_COORD(next_id)] = vel[Z_COORD(p)];

            n_vel_angular[X_COORD(next_id)] = vel_angular[X_COORD(p)];
            n_vel_angular[Y_COORD(next_id)] = vel_angular[Y_COORD(p)];
            n_vel_angular[Z_COORD(next_id)] = vel_angular[Z_COORD(p)];

            memcpy(n_force_old+3*next_id, force_old+3*p, sizeof(vec3));
            memcpy(n_torque_old+3*next_id, torque_old+3*p, sizeof(vec3));

#ifdef TRACK_PARTICLE_ORIENTATION
            n_orientation[4*next_id] = orientation[4*p];
            n_orientation[4*next_id+1] = orientation[4*p+1];
            n_orientation[4*next_id+2] = orientation[4*p+2];
            n_orientation[4*next_id+3] = orientation[4*p+3];
#endif

#ifdef ENABLE_FIXED_PARTICLES
            if(fixed_particles)
                n_fixed_particles[next_id] = fixed_particles[p];
#endif

            new_id[p] = next_id;
            ++next_id;
        }
        else
            new_id[p] = -1;
    }

    // copy contact list
    ContactListEntry **n_contact_list = new ContactListEntry*[new_number_of_particles];
    ContactListEntry *cl_entry = NULL, *last_cl_entry = NULL;

    for(int p = 0; p < number_of_particles; ++p)
    {
        // particle is going to be deleted -> delete contacts
        if(new_id[p] == -1)
        {
            cl_entry = contact_list[p];

            while(cl_entry)
            {
                last_cl_entry = cl_entry;
                cl_entry = cl_entry->next;

                delete last_cl_entry->contact;
                delete last_cl_entry;
            }
        }
        // particle is kept -> update id of its contacts with new values
        else
        {
            n_contact_list[new_id[p]] = contact_list[p];
            cl_entry = contact_list[p];

            while(cl_entry)
            {
                if(new_id[cl_entry->id] < 0)	// remove contact as other particle will be deleted
                {
                    // remove contact at the beginning of the contact list
                    if(cl_entry == n_contact_list[new_id[p]])
                    {
                        n_contact_list[new_id[p]] = cl_entry->next;
                        delete cl_entry->contact;
                        delete cl_entry;
                        cl_entry = n_contact_list[new_id[p]];
                    }
                    else // remove contact in the middle of the contact list
                    {
                        last_cl_entry->next = cl_entry->next;
                        delete cl_entry->contact;
                        delete cl_entry;
                        cl_entry = last_cl_entry->next;
                    }
                }
                else // keep contact
                {
                    cl_entry->id = new_id[cl_entry->id];
                    cl_entry->contact->id1 = new_id[cl_entry->contact->id1];
                    cl_entry->contact->id2 = new_id[cl_entry->contact->id2];

                    last_cl_entry = cl_entry;
                    cl_entry = cl_entry->next;
                }
            }
        }
    }

    delete [] contact_list;
    contact_list = NULL;

    // clean up existing
    cleanUp();

    // replace
    pos_old = n_pos_old;
    pos_new = n_pos_new;
    vel = n_vel;
    vel_angular = n_vel_angular;
    force_old = n_force_old;
    force_new = n_force_new;
    torque_old = n_torque_old;
    torque_new = n_torque_new;

#ifdef TRACK_PARTICLE_ORIENTATION
    orientation = n_orientation;
#endif

    contact_list = n_contact_list;

    // set new number of particles and init grid
    number_of_particles = new_number_of_particles;
    number_of_contacts = getNumberOfContacts();
    number_of_wall_contacts = getNumberOfWallContacts();

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
    dissipated_energy_of_particle = new double[number_of_particles];
    memset(dissipated_energy_of_particle, 0.0, sizeof(double) * number_of_particles);
#endif

#ifdef TRACK_CONTACTS_PER_PARTICLE
    initial_number_of_contacts_of_particle = new int[number_of_particles];
    memset(initial_number_of_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    broken_contacts_of_particle = new int[number_of_particles];
    memset(broken_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    created_contacts_of_particle = new int[number_of_particles];
    memset(created_contacts_of_particle, 0, sizeof(int) * number_of_particles);
#endif

    grid.init(256, 2.1 * particle_radius, number_of_particles);
    grid.addParticles(pos_old, number_of_particles);

    initSimData(SIM_TYPE_GENERAL);

    sim_info.info_storage[5] = 2.0 * (double)number_of_contacts / (double)number_of_particles;

    delete [] new_id;
}

void Simulation::removeWalls()
{
    if(walls.size() == 0)
        return;

    // delete walls
    walls.clear();
    number_of_walls = 0;

    // remove any wall contacts from contact list
    ContactListEntry *cl_entry = NULL, *last_cl_entry = NULL;

    for(int p = 0; p < number_of_particles; ++p)
    {
        cl_entry = contact_list[p];

        while(cl_entry)
        {
            if(cl_entry->id < 0)
            {
                // remove contact at the beginning of the contact list
                if(cl_entry == contact_list[p])
                {
                    contact_list[p] = cl_entry->next;
                    delete cl_entry->contact;
                    delete cl_entry;
                    cl_entry = contact_list[p];
                }
                else // remove contact in the middle of the contact list
                {
                    last_cl_entry->next = cl_entry->next;
                    delete cl_entry->contact;
                    delete cl_entry;
                    cl_entry = last_cl_entry->next;
                }
            }
            else
            {
                last_cl_entry = cl_entry;
                cl_entry = cl_entry->next;
            }
        }
    }

    // clean up box
    if(box)
    {
        delete box;
        box = NULL;
    }

    initSimData(SIM_TYPE_GENERAL);

    sim_info.info_storage[5] = 2.0 * (double)getNumberOfContacts() / (double)number_of_particles, 0.0;
}

ErrorCode Simulation::loadFromFile(const char* filename)
{

    // set default values
    sim_info.info_storage[0] = 0.0;
    sim_info.info_storage[1] = 0.0;
    sim_info.info_storage[2] = 0.0;
    sim_info.info_storage[3] = 0.0;
    sim_info.info_storage[4] = 0.0;
    sim_info.info_storage[5] = 0.0;
    sim_info.info_storage[6] = 0.0;
    sim_info.info_storage[7] = 0.0;

    // try to open specified file
    FILE *file = fopen(filename, "rb");

    if(!file)
        return EC_FILE_NOT_FOUND;

    ////////////////////////////////////////////////////////////////////////////////////////
    // determine file version
    ////////////////////////////////////////////////////////////////////////////////////////

    int file_version = 0;
    char version_buffer[23];
    fread(version_buffer, sizeof(char), 22, file);
    version_buffer[22] = 0;



    if(strcmp(version_buffer, "DATA_FILE_VERSION_1_10") == 0)
        file_version = 10;

    if(file_version == 0)
    {
        printf("%s - ", version_buffer);
        return EC_INVALID_FILE_VERSION;
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    // load description
    ////////////////////////////////////////////////////////////////////////////////////////

    int num_particles, num_walls, num_contacts, num_wall_contacts;
    unsigned int broken, created;

    // stores box info (will otherwise be deleted when resizing arrays)
    WallBox *temp_box = NULL;

    if(file_version >= 9)
    {
        DataFileOffsets offsets;
        SimInfo info;

        fread(&offsets, sizeof(DataFileOffsets), 1, file);
        fread(&num_particles, sizeof(int), 1, file);
        fread(&num_walls, sizeof(int), 1, file);
        fread(&num_contacts, sizeof(int), 1, file);
        fread(&broken, sizeof(unsigned int), 1, file);
        fread(&created, sizeof(unsigned int), 1, file);
        fread(&num_wall_contacts, sizeof(int), 1, file);
        fread(&info, sizeof(SimInfo), 1, file);

        if(strcmp(sim_info.material_identifier, info.material_identifier) == 0)
            memcpy(&sim_info, &info, sizeof(SimInfo));
        else
        {
            printf("Sim: %s   File: %s\n", sim_info.material_identifier, info.material_identifier);
            return EC_DIFFERING_MATERIAL;
        }


        // load box info
        if(sim_info.sim_type >= SIM_TYPE_COMPRESSION_NO_SIDE_WALLS && sim_info.sim_type <= SIM_TYPE_SHEAR_STRENGTH_TEST)
        {
            temp_box = new WallBox;
            fread(temp_box, sizeof(WallBox), 1, file);
        }
    }



    /////////////////////////////////////////////////////////////////////////////////////////////
    // reset values & resize arrays
    /////////////////////////////////////////////////////////////////////////////////////////////

    resizeArrays(num_particles, num_walls);
    memset(force_old, 0, number_of_particles * sizeof(vec3));
    memset(force_new, 0, number_of_particles * sizeof(vec3));
    memset(torque_old, 0, number_of_particles * sizeof(vec3));
    memset(torque_new, 0, number_of_particles * sizeof(vec3));

    // set box
    box = temp_box;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // read particle data
    /////////////////////////////////////////////////////////////////////////////////////////////

    if(file_version >= 9)
    {
        fread(pos_old, sizeof(double), 3*number_of_particles, file);
        fread(vel, sizeof(double), 3*number_of_particles, file);
        fread(vel_angular, sizeof(double), 3*number_of_particles, file);

        fread(force_old, sizeof(double), 3*number_of_particles, file);
        fread(torque_old, sizeof(double), 3*number_of_particles, file);
    }



    /////////////////////////////////////////////////////////////////////////////////////////////
    // read walls
    /////////////////////////////////////////////////////////////////////////////////////////////

    for(int w = 0; w < number_of_walls; ++w)
    {
        fread(&(walls[w]), sizeof(Wall), 1, file);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    // read contacts - use contact as buffer
    /////////////////////////////////////////////////////////////////////////////////////////////

    bool contact_list_corrupted = false;
    ContactListEntry *cl_entry;

    number_of_contacts = num_contacts;
    number_of_wall_contacts = num_wall_contacts;

    for(int c = 0; c < number_of_contacts + number_of_wall_contacts; ++c)
    {
        Contact *new_contact = new Contact;

        if(file_version == 10)
        {
            fread(&(new_contact->id1), sizeof(int), 1, file);
            fread(&(new_contact->id2), sizeof(int), 1, file);
            fread(&(new_contact->n1_initial), sizeof(vec3), 1, file);
            fread(&(new_contact->n2_initial), sizeof(vec3), 1, file);
            fread(&(new_contact->old_contact_normal), sizeof(vec3), 1, file);
            fread(&(new_contact->rot1), sizeof(RotParam), 1, file);
            fread(&(new_contact->rot2), sizeof(RotParam), 1, file);
            fread(&(new_contact->twisting_displacement), sizeof(double), 1, file);
            fread(&(new_contact->compression_length), sizeof(double), 1, file);
        }


        // FIX: filter out corrupted contacts (occured very rarely)
        if(new_contact->id2 < number_of_particles && new_contact->id2 > - (number_of_walls+1) )
        {

            /*  Lucas: this part alters the data for no apparent reason
            printf("%d  %d  %d\n", number_of_particles, new_contact->id2, (number_of_walls+1));

            if(new_contact->id2 >= 0)
            {
                vec3 n_c;
                n_c[0] = pos_old[X_COORD(new_contact->id1)] - pos_old[X_COORD(new_contact->id2)];
                n_c[1] = pos_old[Y_COORD(new_contact->id1)] - pos_old[Y_COORD(new_contact->id2)];
                n_c[2] = pos_old[Z_COORD(new_contact->id1)] - pos_old[Z_COORD(new_contact->id2)];
                normalize(&n_c);

                new_contact->old_contact_normal[0] = n_c[0];
                new_contact->old_contact_normal[1] = n_c[1];
                new_contact->old_contact_normal[2] = n_c[2];
            }
            else
            {
                if(file_version <= 9)
                {
                    int wall_id = WALL_ID(new_contact->id2);
                    new_contact->n2_initial[0] -= walls[wall_id].pos[0];
                    new_contact->n2_initial[1] -= walls[wall_id].pos[1];
                    new_contact->n2_initial[2] -= walls[wall_id].pos[2];
                }
            }
            */

            if(contact_list[new_contact->id1] == NULL) // particle has no other contacts yet
            {
                // create contact list entry
                contact_list[new_contact->id1] = new ContactListEntry;
                contact_list[new_contact->id1]->next = NULL;
                contact_list[new_contact->id1]->id = new_contact->id2;
                contact_list[new_contact->id1]->contact = new_contact;
            }
            else  // particle has other contacts -> append new contact at the end
            {
                // check if particles are not already in contact
                cl_entry = contact_list[new_contact->id1];

                while(cl_entry)
                {
                    if(cl_entry->next)
                        cl_entry = cl_entry->next;
                    else
                    {
                        // create contact list entry
                        cl_entry->next = new ContactListEntry;
                        cl_entry->next->next = NULL;
                        cl_entry->next->id = new_contact->id2;
                        cl_entry->next->contact = new_contact;

                        // contact has been appended -> abort loop
                        cl_entry = NULL;
                    }
                }
            }
        }
        else
        {
            contact_list_corrupted = true;
        }

    }







    fclose(file);

    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up data for sim
    /////////////////////////////////////////////////////////////////////////////////////////////

    memcpy(pos_new, pos_old, number_of_particles * sizeof(vec3));


    grid.init(256, 2.1 * particle_radius, number_of_particles);
    grid.addParticles(pos_old, number_of_particles);

    initSimData(sim_info.sim_type, sim_info.info_storage[2], sim_info.info_storage[0], sim_info.info_storage[4], sim_info.info_storage[1], sim_info.info_storage[3]);

    broken_contacts = broken;
    created_contacts = created;


    SimLib::printEnergy(this, NULL);

    if(contact_list_corrupted)
        return EC_CONTACT_LIST_CORRUPTED;
    else
        return EC_OK;
}

ErrorCode Simulation::saveToFile(const char* filename, bool save_contacts)
{

    if(number_of_particles < 1)
    {
        printf("Warning trying to save empty Simulation!\n");
        return EC_NO_PARTICLES;
    }

    // try to open specified file
    FILE *file = fopen(filename, "wb+");

    if(!file)
        return EC_FILE_NOT_FOUND;


    for(int i = 0; i < this->number_of_walls; ++i)
    {
        this->walls[i].updateEdges();
    }

    DataFileOffsets offsets;
    offsets.pos_offset = 0;
    offsets.vel_offset = 0;
    offsets.vel_angular_offset = 0;
    offsets.walls_offset = 0;
    offsets.contacts_offset = 0;

    int num_of_contacts = getNumberOfContacts();
    int num_of_wall_contacts = getNumberOfWallContacts();

    if(!save_contacts)
        num_of_contacts = 0;
    int num_of_walls = (unsigned int)walls.size();

    // write file version, simulation type, etc.
    fwrite(DATA_FILE_VERSION, sizeof(char), 22, file);
    fwrite(&offsets, sizeof(DataFileOffsets), 1, file);
    fwrite(&number_of_particles, sizeof(int), 1, file);
    fwrite(&num_of_walls, sizeof(int), 1, file);
    fwrite(&num_of_contacts, sizeof(int), 1, file);
    fwrite(&broken_contacts, sizeof(unsigned int), 1, file);
    fwrite(&created_contacts, sizeof(unsigned int), 1, file);
    fwrite(&num_of_wall_contacts, sizeof(int), 1, file);
    fwrite(&sim_info, sizeof(SimInfo), 1, file);

    // store box info
    if(sim_info.sim_type >= SIM_TYPE_COMPRESSION_NO_SIDE_WALLS && sim_info.sim_type <= SIM_TYPE_SHEAR_STRENGTH_TEST)
    {
        if(box)
            fwrite(box, sizeof(WallBox), 1, file);
        else
        {
            fclose(file);
            return EC_NO_BOX;
        }
    }




    fwrite(pos_old, sizeof(double), 3*number_of_particles, file);
    fwrite(vel, sizeof(double), 3*number_of_particles, file);
    fwrite(vel_angular, sizeof(double), 3*number_of_particles, file);

    fwrite(force_old, sizeof(double), 3*number_of_particles, file);
    fwrite(torque_old, sizeof(double), 3*number_of_particles, file);

    for(int w = 0; w < number_of_walls; ++w)
    {
        fwrite(&(walls[w]),  sizeof(Wall), 1, file);
    }


    if(save_contacts)
    {
        ContactListEntry *cl_entry = NULL;

        for(int p = 0; p < number_of_particles; ++p)
        {
            cl_entry = contact_list[p];

            while(cl_entry)
            {
                fwrite(&(cl_entry->contact->id1), sizeof(int), 1, file);
                fwrite(&(cl_entry->contact->id2), sizeof(int), 1, file);
                fwrite(&(cl_entry->contact->n1_initial), sizeof(vec3), 1, file);
                fwrite(&(cl_entry->contact->n2_initial), sizeof(vec3), 1, file);
                fwrite(&(cl_entry->contact->old_contact_normal), sizeof(vec3), 1, file);
                fwrite(&(cl_entry->contact->rot1), sizeof(RotParam), 1, file);
                fwrite(&(cl_entry->contact->rot2), sizeof(RotParam), 1, file);
                fwrite(&(cl_entry->contact->twisting_displacement), sizeof(double), 1, file);
                fwrite(&(cl_entry->contact->compression_length), sizeof(double), 1, file);

                //fwrite(cl_entry->contact, sizeof(Contact), 1, file);
                cl_entry = cl_entry->next;
            }
        }
    }

    fclose(file);

    return EC_OK;
}

#ifdef TRACK_CONTACTS_PER_PARTICLE
void Simulation::initInitialContacts()
{
    memset(initial_number_of_contacts_of_particle, 0, sizeof(int) * number_of_particles);
    memset(broken_contacts_of_particle, 0, sizeof(int) * number_of_particles);
    memset(created_contacts_of_particle, 0, sizeof(int) * number_of_particles);

    ContactListEntry *cl_entry;

    for(int p = 0; p < number_of_particles; ++p)
    {
        cl_entry = contact_list[p];

        while(cl_entry)
        {
            // ignore contacts with walls
            if(cl_entry->id >= 0)
            {
                ++initial_number_of_contacts_of_particle[cl_entry->id];
                ++initial_number_of_contacts_of_particle[p];
            }

            cl_entry = cl_entry->next;
        }
    }
}
#endif

ErrorCode Simulation::initSimData(SimType type, double collision_speed, double stop_filling_factor, double stop_dissipation_factor, double initial_kinetic_energy, double shockwave_penetration_depth)
{
    if(type >= SIM_TYPE_COMPRESSION_NO_SIDE_WALLS && type <= SIM_TYPE_SHOCKWAVE)
    {
        if(!box)
            return EC_NO_BOX;
    }

    sim_info.sim_type = type;
    sim_info.info_storage[0] = stop_filling_factor;
    sim_info.info_storage[1] = initial_kinetic_energy;
    sim_info.info_storage[2] = collision_speed;
    sim_info.info_storage[3] = shockwave_penetration_depth;
    sim_info.info_storage[4] = stop_dissipation_factor;


    dissipated_contact_energy = 0;
    dissipated_rolling_energy = 0;
    dissipated_sliding_energy = 0;
    dissipated_twisting_energy = 0;
    dissipated_wall_energy = 0;
    dissipated_damping_energy = 0;

    current_time = 0;
    end_time = 0;

    broken_contacts = 0;
    created_contacts = 0;

    broken_wall_contacts = 0;
    created_wall_contacts = 0;

#ifdef TRACK_CONTACTS_PER_PARTICLE
    initInitialContacts();
#endif

    grid.resetGrid();
    grid.addParticles(pos_old, number_of_particles);

    number_of_contacts = getNumberOfContacts();
    number_of_wall_contacts = getNumberOfWallContacts();

    initialized = true;
    stop_simulation = false;

    return EC_OK;
}

ErrorCode Simulation::loadMaterial(const char *filename)
{

    printf("Loading material %s\n", filename);
    cleanUp();

    setlocale(LC_ALL, "C"); // for some reasong local was set to comma as digit separator, causing the files to be read incorrectly

    FILE *file = fopen(filename, "r");

    // init with default values (dust) if opening of the file failed
    if(!file)
    {
        setMaterialConstants(6e-5, 2.65, 25.0, 0.17, 5.4e11, 4e-8, 1e-6, 1e-11);
        printf("EC_FILE_NOT_FOUND\n");
        return EC_FILE_NOT_FOUND;
    }

    // check if correct file version
    char buffer[200];
    fscanf(file, "%200s", buffer);

    if( strcmp(buffer, "MATERIAL_FILE_VERSION_1_0") && strcmp(buffer, "MATERIAL_FILE_VERSION_1_1") && strcmp(buffer, "MATERIAL_FILE_VERSION_1_2") && strcmp(buffer, "MATERIAL_FILE_VERSION_1_3") && strcmp(buffer, "MATERIAL_FILE_VERSION_1_4") )
    {
        setMaterialConstants(6e-5, 2.65, 25.0, 0.17, 5.4e11, 4e-8, 1e-6, 1e-11);
        printf("EC_INVALID_MATERIAL_FILE_VERSION\n");
        return EC_INVALID_MATERIAL_FILE_VERSION;
    }


#ifdef USE_JKR
    T_vis = 3.159667e-11; // Brisset r = 7.6 um
#endif // USE_JKR
#ifdef USE_DMT
    T_vis = 4.495807e-11; // Brisset + DMT r = 7.6 um
#endif // USE_DMT
    T_vis = 1.25e-11; // Alex paper
    //T_vis = 3.159667e-11; // Brisset r = 7.6 um   fittet to v_s = 105.869
    //T_vis = 7.354933e-12; // Brisset r = 2.5 um   fittet to v_s = 190
    //T_vis = 4.495807e-11; // Brisset + DMT r = 7.6 um "
    //T_vis = 1.076728e-11; // Brisset + DMT r = 2.5 um "
    rolling_modifier = 1.0;
    sliding_modifier = 1.0;
    twisting_modifier = 1.0;
    crit_sliding_displacement_modifier = 1.0;
    crit_wall_sliding_displacement_modifier = 1.0;

    // load description
    fscanf(file, "%s", material_name);

    // load values
    if(strcmp(buffer, "MATERIAL_FILE_VERSION_1_0") == 0)
        fscanf(file, "%lf %lf %lf %lf %lf %lf", &particle_radius, &density, &surface_energy, &nu, &young_mod, &crit_rolling_displacement);
    else if(strcmp(buffer, "MATERIAL_FILE_VERSION_1_1") == 0)
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &particle_radius, &density, &surface_energy, &nu, &young_mod, &crit_rolling_displacement, &rolling_modifier, &sliding_modifier, &twisting_modifier);
    else if(strcmp(buffer, "MATERIAL_FILE_VERSION_1_2") == 0)
    {
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &particle_radius,
               &density,
               &surface_energy,
               &nu,
               &young_mod,
               &crit_rolling_displacement,
               &osc_damping_factor,
               &rolling_modifier,
               &sliding_modifier,
               &twisting_modifier);

    }
    else if(strcmp(buffer, "MATERIAL_FILE_VERSION_1_3") == 0)
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &particle_radius, &density, &surface_energy, &nu, &young_mod, &crit_rolling_displacement, &osc_damping_factor, &rolling_modifier, &sliding_modifier, &twisting_modifier, &crit_sliding_displacement_modifier, &crit_wall_sliding_displacement_modifier);
    else if(strcmp(buffer, "MATERIAL_FILE_VERSION_1_4") == 0)
    {
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &particle_radius,
               &density,
               &surface_energy,
               &nu,
               &young_mod,
               &crit_rolling_displacement,
               &osc_damping_factor,
               &T_vis,
               &rolling_modifier,
               &sliding_modifier,
               &twisting_modifier);
    }
    setMaterialConstants(
                particle_radius,
                density,
                surface_energy,
                nu,
                young_mod,
                crit_rolling_displacement,
                osc_damping_factor,
                T_vis,
                rolling_modifier,
                sliding_modifier,
                twisting_modifier,
                crit_sliding_displacement_modifier,
                crit_wall_sliding_displacement_modifier
                );

    fclose(file);
    return EC_OK;
}

void Simulation::setMaterialConstants(
        const double particle_radius,
        const double density,
        const double surface_energy,
        const double nu,
        const double young_mod,
        const double crit_rolling_displacement,
        const double osc_damping_factor,
        const double T_vis_,
        const double rolling_modifier,
        const double sliding_modifier,
        const double twisting_modifier,
        const double crit_sliding_displacement_modifier,
        const double crit_wall_sliding_displacement_modifier
        )
{

    // set values
    ::particle_radius = particle_radius;
    ::density = density;
    ::surface_energy = surface_energy;
    ::nu = nu;
    ::young_mod = young_mod;
    ::crit_sliding_displacement_modifier = crit_sliding_displacement_modifier;
    ::crit_wall_sliding_displacement_modifier = crit_wall_sliding_displacement_modifier;
    ::osc_damping_factor = osc_damping_factor;
    ::T_vis = T_vis_;

    // calculate the other constants
    young_mod_reduced = 0.5 * young_mod / (1.0 - nu * nu);
    shear_mod = 0.5 * young_mod / (1.0 + nu);
    shear_mod_reduced = 0.5 * shear_mod / (2.0 - nu);
    bulk_mod = young_mod * shear_mod / (3.0 * (3.0 * shear_mod - young_mod) );

    particle_radius_inv = 1.0 / particle_radius;
    reduced_radius = 0.5 * particle_radius;
    equilibrium_contact_radius = pow(9.0 * M_PI * surface_energy * reduced_radius * reduced_radius / young_mod_reduced, 1.0/3.0);
    delta_0 = equilibrium_contact_radius * equilibrium_contact_radius / (3.0 * reduced_radius);
    equilibrium_distance = 2.0 * particle_radius - delta_0;
    delta_c = 0.5 * equilibrium_contact_radius * equilibrium_contact_radius / ( reduced_radius * pow(6.0, 1.0 / 3.0));
    F_c = 3.0 * M_PI * surface_energy * reduced_radius;
    mass = density * 4.0 / 3.0 * M_PI * particle_radius * particle_radius * particle_radius;
    mass_inv = 1.0 / mass;
    moment_of_inertia = 2.0 / 5.0 * mass * particle_radius * particle_radius;
    moment_of_inertia_inv = 1.0 /moment_of_inertia;

#ifdef USE_KRIJT_CRITICAL_ROLLING_DISPLACEMENT
    ::crit_rolling_displacement = equilibrium_contact_radius / 12.0 * crit_rolling_displacement;
#else
    ::crit_rolling_displacement = crit_rolling_displacement;
#endif



    viscoelastic_damping_constant =  2.0 * T_vis / (nu * nu) * young_mod_reduced;


#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
    k_r = rolling_modifier * 4.0 * F_c / reduced_radius;
    k_s = sliding_modifier * 8.0 * shear_mod_reduced; // * particle_radius * particle_radius;
    k_t = twisting_modifier * 16.0 / 3.0 * shear_mod_reduced;
#else
    k_r = rolling_modifier * 4.0 * F_c / reduced_radius;
    k_s = sliding_modifier * 8.0  * shear_mod_reduced * equilibrium_contact_radius; // * particle_radius * particle_radius;
    // bug from Alex: shear_mod_reduced = G*  != G = 0.5 * shear_mod
    //k_t = twisting_modifier * 16.0 / 3.0 * shear_mod_reduced * equilibrium_contact_radius * equilibrium_contact_radius * equilibrium_contact_radius;
    k_t = twisting_modifier * 16.0 / 3.0 * shear_mod*0.5 * equilibrium_contact_radius * equilibrium_contact_radius * equilibrium_contact_radius;
#endif

    printf("constants rolling sliding twisting  %f  %f  %f\n", k_r, k_s, k_t);

    t_c = sqrt(mass * delta_c / F_c);	// see Wada et al. 2007
    T_normal = 2.0 * M_PI * sqrt(5.0/6.0 * mass / (young_mod_reduced * equilibrium_contact_radius) );	// see eq. 34/35 D&T 2004 "Resistance to rolling in the adhesive contact of two elastic spheres"
    v_c = 1.07 * pow(surface_energy/reduced_radius, 5.0/6.0) * pow(density, -0.5) * pow(young_mod_reduced, - 1.0/3.0); // see eq. 3 Dominik&Paszun 2008 "Determination of material propertiesof dust cakes"
    bulk_sound_speed = sqrt( (bulk_mod + 4.0/3.0 * shear_mod) / density);
    ENERGY_UNIT = 1.0 / (F_c * delta_c);



#ifdef USE_JKR
    contact_breaking_dist = 2.0 * particle_radius + delta_c;
    contact_making_dist = 2.0 * particle_radius;
    //contact_making_dist = (2.0 * particle_radius + 0.465345 * delta_c); // dissipate less energy on contact creation
    //contact_making_dist = (2.0 * particle_radius + delta_c); // dissipate no energy on contact creation
    contact_breaking_dist_squared = contact_breaking_dist * contact_breaking_dist;
    contact_making_dist_squared = contact_making_dist * contact_making_dist;
#endif


#ifdef USE_DMT
    contact_breaking_dist = 2.0 * particle_radius;
    contact_making_dist = 2.0 * particle_radius;
    contact_breaking_dist_squared = contact_breaking_dist * contact_breaking_dist;
    contact_making_dist_squared = contact_making_dist * contact_making_dist;
#endif



    // for wall
#ifdef USE_DEFAULT_WALL_NORMAL_INTERACTION
    wall_reduced_radius = reduced_radius;
#else
    wall_reduced_radius = 2.0 * reduced_radius;
#endif

    wall_equilibrium_contact_radius = pow(9.0 * M_PI * surface_energy * wall_reduced_radius * wall_reduced_radius / young_mod_reduced, 1.0/3.0);
    wall_delta_c = 0.5 * wall_equilibrium_contact_radius * wall_equilibrium_contact_radius / ( wall_reduced_radius * pow(6.0, 1.0 / 3.0));
    wall_delta_0 = wall_equilibrium_contact_radius * wall_equilibrium_contact_radius / (3.0 * wall_reduced_radius);
#ifdef USE_JKR
    wall_contact_breaking_dist = 2.0 * particle_radius + wall_delta_c; // wall_contact_breaking_dist has a factor of 2.0 because of the algorithm design by Alex
#endif
#ifdef USE_DMT
    wall_contact_breaking_dist = 2.0 * particle_radius;
#endif

    wall_F_c = 3.0 * M_PI * surface_energy * wall_reduced_radius;
    wall_k_s = sliding_modifier * 8.0 * wall_equilibrium_contact_radius * shear_mod_reduced;
    wall_k_r = rolling_modifier * 4.0 * wall_F_c / wall_reduced_radius;       // Alex orginal code : wall_k_r = rolling_modifier * 4.0 * wall_F_c / reduced_radius;
                                                                         // was reduced_radius instead of wall_reduced_radius a bug?? see A. Seizinger et al. 2012  Compression behavior of porous dust agglomerates
                                                                         // "...where kr,wall is equivalent to the rolling constant kr given in Eq. (17) taking the different reduced radius of the particle-wall interaction into account."

    wall_k_t = twisting_modifier * 16.0 / 3.0 * shear_mod_reduced * wall_equilibrium_contact_radius * wall_equilibrium_contact_radius * wall_equilibrium_contact_radius;
    crit_wall_sliding_displacement = crit_wall_sliding_displacement_modifier * (2.0 - nu)/(16.0 * M_PI) * wall_equilibrium_contact_radius;
    crit_wall_sliding_displacement_squared = crit_wall_sliding_displacement*crit_wall_sliding_displacement;

    wall_acceleration = 2.0 * particle_radius / pow(50.0*T_normal, 2);

    crit_rolling_displacement_squared = crit_rolling_displacement * crit_rolling_displacement;
    crit_sliding_displacement = crit_sliding_displacement_modifier * (2.0 - nu)/(16.0 * M_PI) * equilibrium_contact_radius;
    crit_sliding_displacement_squared = crit_sliding_displacement * crit_sliding_displacement;
    crit_twisting_displacement = 1.0 / (16.0 * M_PI);
    crit_twisting_displacement_squared = crit_twisting_displacement * crit_twisting_displacement;


#ifdef INTERPOLATE_NORMAL_FORCES
    normal_interaction.init(NORMAL_FORCES_INTERPOLATION_POINTS);
#else
    normal_interaction.init();
#endif

    // set material identifier
    if(strlen(material_name) >= 3)
    {
        sim_info.material_identifier[0] = material_name[0];
        sim_info.material_identifier[1] = material_name[1];
        sim_info.material_identifier[2] = material_name[2];
    }
    else
    {
        sim_info.material_identifier[0] = 'U';
        sim_info.material_identifier[1] = 'N';
        sim_info.material_identifier[2] = 'K';
    }

    char buffer[10];
    sprintf(buffer, "%1.2g", particle_radius);
    sim_info.material_identifier[3] = buffer[0];
    sim_info.material_identifier[4] = buffer[3];

    sprintf(buffer, "%1.2g", density);
    sim_info.material_identifier[5] = buffer[0];
    sim_info.material_identifier[6] = buffer[2];

    sprintf(buffer, "%1.2g", surface_energy);
    sim_info.material_identifier[7] = buffer[0];
    sim_info.material_identifier[8] = buffer[1];

    sprintf(buffer, "%1.2g", young_mod);
    sim_info.material_identifier[9] = buffer[0];
    sim_info.material_identifier[10] = buffer[2];

    sim_info.material_identifier[11] = 0;
}

void Simulation::setPotentialVariationStopCondition(double min_time, double potential_variation_stop_threshold, int check_potential_variation_interval)
{
    min_end_time = current_time + min_time;
    check_potential_variation_counter = 0;
    this->check_potential_variation_interval = check_potential_variation_interval;
    this->potential_variation_stop_threshold = potential_variation_stop_threshold;
    last_dissipated_energy = dissipated_contact_energy + dissipated_damping_energy + dissipated_rolling_energy + dissipated_sliding_energy + dissipated_twisting_energy;
}

ErrorCode Simulation::startSimulation(double duration, double timestep, int print_energies_interval, int print_positions_interval, const char *energy_filename, const char *print_positions_path, bool use_gravity)
{
    if(initialized)
    {
        check_potential_variation_counter = 0;

        // check timestep
#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
        double sliding_osc_period = 2.0 * M_PI * sqrt(mass / (k_s * equilibrium_contact_radius) );
#else
        double sliding_osc_period = 2.0 * M_PI * sqrt(mass / k_s);
#endif

        printf("P_sliding / dt = %.12g\n\n", sliding_osc_period/timestep);


        if(timestep > 0.05 * sliding_osc_period)
            return EC_TIMESTEP_TOO_LARGE;



        this->timestep = timestep;
        end_time = current_time + duration;
        stop_simulation = false;
        print_energies_counter = 0;
        print_positions_counter = 0;
        print_positions_file_counter = 0;
        this->print_energies_interval = print_energies_interval;
        this->print_positions_interval = print_positions_interval;

        if(energies_filename)
        {
            delete [] energies_filename;
            energies_filename = NULL;
        }

        if(positions_path)
        {
            delete [] positions_path;
            positions_path = NULL;
        }

        if(print_energies_interval > 0)
        {
            size_t length = strlen(energy_filename);
            energies_filename = new char[length+1];

            strcpy(energies_filename, energy_filename);

            FILE *file = fopen(energy_filename, "w+");

            if(file)
            {
                fprintf(file, "#time, E_tot, E_kin, E_rot, V_tot, V_normal, V_roll, V_slide, V_twist, E_diss, diss_contact, diss_rolling, diss_sliding, diss_twisting, diss_wall, diss_damping)\n");
                fclose(file);
            }
        }


        if(print_positions_path)
        {
            size_t length = strlen(print_positions_path);

            if(length > 0)
            {
                positions_path = new char[length+1];
                strcpy(positions_path, print_positions_path);
            }
        }


#ifdef ENABLE_GRAVITY
        gravity_enabled = use_gravity;
#endif

        return EC_OK;
    }
    else
        return EC_SIM_NOT_INITIALIZED;
}

void Simulation::update()
{

    // predictor step: calculate new positions/contact pointers
    predictor(timestep);

    updateContacts(timestep);

    // check for changes of contacts due to updated positions
    updateSticking();

    // evaluation step: calculate forces/torques based on predicted positions/contact pointers
    updateParticleInteraction(timestep);

    if(azimuthal_acceleration > 0)
        induceRotation(azimuthal_acceleration);

    // calculate velocity based on forces/torques from evaluation step
    corrector(timestep);

    if(damping_enabled)
        dampVelocities(damping_factor);

#ifdef TRACK_PARTICLE_ORIENTATION
    updateParticleOrientation(timestep);
#endif


    // change pointers - new value are now old values and the former old values will be overwritten with new values in the next integration step
    switchOldNewValues();


    // log data if necessary
    ++print_energies_counter;
    ++print_positions_counter;

    if(print_energies_interval &&  print_energies_counter >= print_energies_interval)
    {
        print_energies_counter = 0;

        if(energies_filename)
            SimLib::printEnergy(this, energies_filename);

        /*{
            double mean, sigma;
            calculateRotation(&mean, &sigma);

            fprintf(file_energies, "%g %g %g\n", current_time, mean, sigma);
        }*/
    }

    if(print_positions_interval && print_positions_counter >= print_positions_interval)
    {
        print_positions_counter = 0;
        ++print_positions_file_counter;

        char filename[200];
        char buf[200];

        // copy path
        strcpy(filename, positions_path);

        // determine filename
        sprintf(buf, "positions_%i.dat", print_positions_file_counter);

        // append filename to path
        strcat(filename, buf);

        SimLib::printPositions(*this, filename);
    }

    current_time += timestep;

    // check time
    if(current_time > end_time)
        stop_simulation = true;

    // check energy dissipation
    if(check_potential_variation_interval > 0)
    {
        ++check_potential_variation_counter;

        if(check_potential_variation_counter > check_potential_variation_interval)
        {
            check_potential_variation_counter = 0;

            // determine the ammount of dissipated energy since last check
            double dissipated_energy = dissipated_contact_energy + dissipated_damping_energy + dissipated_rolling_energy + dissipated_sliding_energy + dissipated_twisting_energy;

            if(dissipated_energy - last_dissipated_energy < potential_variation_stop_threshold && current_time > min_end_time)
                stop_simulation = true;
            else
                last_dissipated_energy = dissipated_energy;
        }
    }

    // check sim type specific behaviour
    if(sim_info.sim_type == SIM_TYPE_COMPRESSION_RELAXATION)
    {
        if(sim_info.info_storage[0] > 0)
        {
            // stop wall if info_storage[0] has been reached
            if(getBoxFillingFactor() > sim_info.info_storage[0])
            {
                sim_info.info_storage[0] = 0;
                walls[box->top_wall_id].velocity[0] = 0;
                walls[box->top_wall_id].velocity[1] = 0;
                walls[box->top_wall_id].velocity[2] = 0;
                sim_info.info_storage[1] = SimLib::getKineticEnergy(*this);
            }
        }
        else	// wall is already stopped -> check dissipation of kinetic energy
        {
            if(SimLib::getKineticEnergy(*this) / sim_info.info_storage[1] < sim_info.info_storage[4])
                stop_simulation = true;
        }
    }
    else if(sim_info.sim_type == SIM_TYPE_SHOCKWAVE)
    {
        if(sim_info.info_storage[3] > 0)
            //&& (box->lower_pos_y + box->height - walls[box->top_wall_id].pos[1]) > penetration_depth * particle_radius)
        {
            double bottom_pressure = 0.1 * walls[box->bottom_wall_id].getWallForce() / box->base;

            if(bottom_pressure > sim_info.info_storage[3])
                stop_simulation = true;
        }
    }
    else if(sim_info.sim_type == SIM_TYPE_COMPACTION || sim_info.sim_type == SIM_TYPE_COMPRESSION_WITH_SIDE_WALLS || sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
    {
        // stop sim if walls are too close
        if(box->height < 4.0 * particle_radius)
        {
            printf("Wall too low\n");
            stop_simulation = true;
        }

        /*
        if(sim_info.info_storage[0] > 0 && getBoxFillingFactor() > sim_info.info_storage[0])
        {
            printf("Stop filling factor reached %e > %e\n", getBoxFillingFactor(), sim_info.info_storage[0]);
            stop_simulation = true;
        }
        */

        if(sim_info.info_storage[7] > 0 && box->height < sim_info.info_storage[7])
        {
            printf("Stop height reached %e < %e\n", box->height, sim_info.info_storage[7]);
            stop_simulation = true;
        }

    }

}

void Simulation::predictor(double timestep)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // move particles
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_FIXED_PARTICLES
    if(fixed_particles)
    {
        // fixed particles are present
        for(int p = 0; p < number_of_particles; ++p)
        {
            if(fixed_particles[p])
            {
                pos_new[X_COORD(p)] = pos_old[X_COORD(p)] + timestep * vel[X_COORD(p)];
                pos_new[Y_COORD(p)] = pos_old[Y_COORD(p)] + timestep * vel[Y_COORD(p)];
                pos_new[Z_COORD(p)] = pos_old[Z_COORD(p)] + timestep * vel[Z_COORD(p)];
            }
            else
            {
                pos_new[X_COORD(p)] = pos_old[X_COORD(p)] + timestep * vel[X_COORD(p)] + 0.5 * mass_inv * timestep * timestep * force_old[X_COORD(p)];
                pos_new[Y_COORD(p)] = pos_old[Y_COORD(p)] + timestep * vel[Y_COORD(p)] + 0.5 * mass_inv * timestep * timestep * force_old[Y_COORD(p)];
                pos_new[Z_COORD(p)] = pos_old[Z_COORD(p)] + timestep * vel[Z_COORD(p)] + 0.5 * mass_inv * timestep * timestep * force_old[Z_COORD(p)];
            }
        }
    }
    else
    {
#endif
        for(int p = 0; p < number_of_particles; ++p)
        {
            pos_new[X_COORD(p)] = pos_old[X_COORD(p)] + timestep * vel[X_COORD(p)] + 0.5 * mass_inv * timestep * timestep * force_old[X_COORD(p)];
            pos_new[Y_COORD(p)] = pos_old[Y_COORD(p)] + timestep * vel[Y_COORD(p)] + 0.5 * mass_inv * timestep * timestep * force_old[Y_COORD(p)];
            pos_new[Z_COORD(p)] = pos_old[Z_COORD(p)] + timestep * vel[Z_COORD(p)] + 0.5 * mass_inv * timestep * timestep * force_old[Z_COORD(p)];
        }
#ifdef ENABLE_FIXED_PARTICLES
    }
#endif

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // move walls
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(unsigned int w = 0; w < walls.size(); ++w)
    {
        walls[w].pos[0] += timestep * walls[w].velocity[0];
        walls[w].pos[1] += timestep * walls[w].velocity[1];
        walls[w].pos[2] += timestep * walls[w].velocity[2];
    }

    updateBox();
}

void Simulation::updateContacts(double timestep)
{
    ContactListEntry *next_cl_entry;
    vec3 omega1, omega2, omega1_dot, omega2_dot;

    for(int p = 0; p < number_of_particles; ++p)
    {
        next_cl_entry = contact_list[p];

        while(next_cl_entry)
        {
            // contact between wall and particle
            if(next_cl_entry->id < 0)
            {
                omega1[0] = vel_angular[X_COORD(p)];
                omega1[1] = vel_angular[Y_COORD(p)];
                omega1[2] = vel_angular[Z_COORD(p)];

                omega1_dot[0] = moment_of_inertia_inv * torque_old[X_COORD(p)];
                omega1_dot[1] = moment_of_inertia_inv * torque_old[Y_COORD(p)];
                omega1_dot[2] = moment_of_inertia_inv * torque_old[Z_COORD(p)];

                next_cl_entry->contact->updateRotation(omega1, omega1_dot, timestep);
            }
            else
            {
                omega1_dot[0] = moment_of_inertia_inv * torque_old[X_COORD(p)];
                omega1_dot[1] = moment_of_inertia_inv * torque_old[Y_COORD(p)];
                omega1_dot[2] = moment_of_inertia_inv * torque_old[Z_COORD(p)];
                omega2_dot[0] = moment_of_inertia_inv * torque_old[X_COORD(next_cl_entry->id)];
                omega2_dot[1] = moment_of_inertia_inv * torque_old[Y_COORD(next_cl_entry->id)];
                omega2_dot[2] = moment_of_inertia_inv * torque_old[Z_COORD(next_cl_entry->id)];

                omega1[0] = vel_angular[X_COORD(p)];
                omega1[1] = vel_angular[Y_COORD(p)];
                omega1[2] = vel_angular[Z_COORD(p)];
                omega2[0] = vel_angular[X_COORD(next_cl_entry->id)];
                omega2[1] = vel_angular[Y_COORD(next_cl_entry->id)];
                omega2[2] = vel_angular[Z_COORD(next_cl_entry->id)];

                next_cl_entry->contact->updateRotation(omega1, omega2, omega1_dot, omega2_dot, timestep);
            }

            next_cl_entry = next_cl_entry->next;
        }
    }
}

void Simulation::corrector(double timestep)
{

#ifdef ENABLE_FIXED_PARTICLES
    if(fixed_particles)
    {
        // fixed particles are present
        for(int p = 0; p < number_of_particles; ++p)
        {
            if(!fixed_particles[p])
            {
                vel[X_COORD(p)] += 0.5 * mass_inv * timestep * (force_new[X_COORD(p)] + force_old[X_COORD(p)]);
                vel[Y_COORD(p)] += 0.5 * mass_inv * timestep * (force_new[Y_COORD(p)] + force_old[Y_COORD(p)]);
                vel[Z_COORD(p)] += 0.5 * mass_inv * timestep * (force_new[Z_COORD(p)] + force_old[Z_COORD(p)]);

                vel_angular[X_COORD(p)] += 0.5 * moment_of_inertia_inv * timestep * (torque_new[X_COORD(p)] + torque_old[X_COORD(p)]);
                vel_angular[Y_COORD(p)] += 0.5 * moment_of_inertia_inv * timestep * (torque_new[Y_COORD(p)] + torque_old[Y_COORD(p)]);
                vel_angular[Z_COORD(p)] += 0.5 * moment_of_inertia_inv * timestep * (torque_new[Z_COORD(p)] + torque_old[Z_COORD(p)]);
            }
            else
            {
                torque_new[X_COORD(p)] = 0;
                torque_new[Y_COORD(p)] = 0;
                torque_new[Z_COORD(p)] = 0;
            }
        }
    }
    else
    {
#endif
        for(int p = 0; p < number_of_particles; ++p)
        {
            vel[X_COORD(p)] += 0.5 * mass_inv * timestep * (force_new[X_COORD(p)] + force_old[X_COORD(p)]);
            vel[Y_COORD(p)] += 0.5 * mass_inv * timestep * (force_new[Y_COORD(p)] + force_old[Y_COORD(p)]);
            vel[Z_COORD(p)] += 0.5 * mass_inv * timestep * (force_new[Z_COORD(p)] + force_old[Z_COORD(p)]);

            vel_angular[X_COORD(p)] += 0.5 * moment_of_inertia_inv * timestep * (torque_new[X_COORD(p)] + torque_old[X_COORD(p)]);
            vel_angular[Y_COORD(p)] += 0.5 * moment_of_inertia_inv * timestep * (torque_new[Y_COORD(p)] + torque_old[Y_COORD(p)]);
            vel_angular[Z_COORD(p)] += 0.5 * moment_of_inertia_inv * timestep * (torque_new[Z_COORD(p)] + torque_old[Z_COORD(p)]);

        }
#ifdef ENABLE_FIXED_PARTICLES
    }
#endif

    for(unsigned int w = 0; w < walls.size(); ++w)
    {
        // solve equation of motion
        if(walls[w].mass > 0)
        {
            // total force on wall
            double F_tot = walls[w].normal[0] * walls[w].total_force[0] + walls[w].normal[1] * walls[w].total_force[1] + walls[w].normal[2] * walls[w].total_force[2];

            // net acceleration on wall
            double k = timestep * (1.0 + F_tot / walls[w].mass) * wall_inertia * wall_acceleration;

            walls[w].velocity[0] += k * walls[w].normal[0];
            walls[w].velocity[1] += k * walls[w].normal[1];
            walls[w].velocity[2] += k * walls[w].normal[2];
        }
    }
}

void Simulation::dampVelocities(double damping_factor)
{
    for(int p = 0; p < number_of_particles; ++p)
    {
        vel[X_COORD(p)] *= damping_factor;
        vel[Y_COORD(p)] *= damping_factor;
        vel[Z_COORD(p)] *= damping_factor;

        vel_angular[X_COORD(p)] *= damping_factor;
        vel_angular[Y_COORD(p)] *= damping_factor;
        vel_angular[Z_COORD(p)] *= damping_factor;
    }
}

void Simulation::induceRotation(double azimuthal_acceleration)
{
    vec3 e_r = {0.0, 0.0, 0.0};
    vec3 e_phi = {0.0, 1.0, 0.0};
    vec3 rot_axis = {0.0, 1.0, 0.0};

    for(int p = 0; p < number_of_particles; ++p)
    {
        double dist = sqrt(pos_old[X_COORD(p)]*pos_old[X_COORD(p)] + pos_old[Z_COORD(p)]*pos_old[Z_COORD(p)]);

        e_r[0] = pos_old[X_COORD(p)] / dist;
        e_r[2] = pos_old[Z_COORD(p)] / dist;

        e_phi[0] = rot_axis[1] * e_r[2] - rot_axis[2] * e_r[1];
        e_phi[1] = rot_axis[2] * e_r[0] - rot_axis[0] * e_r[2];
        e_phi[2] = rot_axis[0] * e_r[1] - rot_axis[1] * e_r[0];

        double strength = azimuthal_acceleration * mass * dist;

        force_new[X_COORD(p)] += strength * e_phi[0];
        force_new[Y_COORD(p)] += strength * e_phi[1];
        force_new[Z_COORD(p)] += strength * e_phi[2];
    }
}

void Simulation::updateSticking()
{
    double dist;
    vec3 delta_pos;
    std::list<int> particles;

    ContactListEntry *next_cl_entry, *last_cl_entry;

    GridEntry *next_entry;

    // update grid first
    grid.resetGrid();
    grid.addParticles(pos_new, number_of_particles);

    for(int p = 0; p < number_of_particles; ++p)
    {
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check if existing contacts have broken
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        next_cl_entry = contact_list[p];
        last_cl_entry = contact_list[p];

        while(next_cl_entry)
        {
            if(next_cl_entry->id < 0) // contact with wall
            {
                dist = walls[WALL_ID(next_cl_entry->id)].getDistanceTo(pos_new[X_COORD(p)], pos_new[Y_COORD(p)], pos_new[Z_COORD(p)]);

                // break contacts when particles leave the box while sliding/rolling along the side walls
                if(number_of_walls == 5 && WALL_ID(next_cl_entry->id) != box->bottom_wall_id)
                {
                    if(pos_new[Y_COORD(p)] > box->lower_pos_y + box->height)
                        dist = 1.001 * wall_contact_breaking_dist;
                }

                if(dist > wall_contact_breaking_dist)
                {
#ifdef TRACK_DISSIPATED_ENERGY
                    updateDissipatedBreakingEnergy(next_cl_entry->contact, dist);
#endif
                    ++broken_wall_contacts;
                    --number_of_wall_contacts;

                    // remove contact at the beginning of the contact list
                    if(next_cl_entry == contact_list[p])
                    {
                        contact_list[p] = next_cl_entry->next;
                        delete next_cl_entry->contact;
                        delete next_cl_entry;
                        next_cl_entry = contact_list[p];
                    }
                    else // remove contact in the middle of the contact list
                    {
                        last_cl_entry->next = next_cl_entry->next;
                        delete next_cl_entry->contact;
                        delete next_cl_entry;
                        next_cl_entry = last_cl_entry->next;
                    }
                }
                else
                {
                    last_cl_entry = next_cl_entry;
                    next_cl_entry = next_cl_entry->next;
                }
            }
            else	// contact with other particle
            {
                dist = (pos_new[X_COORD(p)] - pos_new[X_COORD(next_cl_entry->id)]) * (pos_new[X_COORD(p)] - pos_new[X_COORD(next_cl_entry->id)])
                     + (pos_new[Y_COORD(p)] - pos_new[Y_COORD(next_cl_entry->id)]) * (pos_new[Y_COORD(p)] - pos_new[Y_COORD(next_cl_entry->id)])
                     + (pos_new[Z_COORD(p)] - pos_new[Z_COORD(next_cl_entry->id)]) * (pos_new[Z_COORD(p)] - pos_new[Z_COORD(next_cl_entry->id)]);

                if( dist > contact_breaking_dist_squared )
                {

#ifdef TRACK_CONTACTS_PER_PARTICLE
                    ++broken_contacts_of_particle[p];
                    ++broken_contacts_of_particle[next_cl_entry->id];
#endif

#ifdef TRACK_DISSIPATED_ENERGY
                    dist = 1.0 / sqrt(dist);
                    updateDissipatedBreakingEnergy(next_cl_entry->contact, dist);
#endif
                    ++broken_contacts;
                    --number_of_contacts;

                    // remove contact at the beginning of the contact list
                    if(next_cl_entry == contact_list[p])
                    {
                        contact_list[p] = next_cl_entry->next;
                        delete next_cl_entry->contact;
                        delete next_cl_entry;
                        next_cl_entry = contact_list[p];
                    }
                    else // remove contact in the middle of the contact list
                    {
                        last_cl_entry->next = next_cl_entry->next;
                        delete next_cl_entry->contact;
                        delete next_cl_entry;
                        next_cl_entry = last_cl_entry->next;
                    }
                }
                else
                {
                    last_cl_entry = next_cl_entry;
                    next_cl_entry = next_cl_entry->next;
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check if new contacts with other particles occurred
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ContactListEntry *cl_entry;

        for(int x = -1; x <= 1; ++x)
        {
            for(int y = -1; y <= 1; ++y)
            {
                for(int z = -1; z <= 1; ++z)
                {
                    int cell_id = grid.cell_of_particle[p] - x - y * grid.x_cells - z * grid.x_cells * grid.y_cells;

                    next_entry = grid.particles_in_cell[cell_id];

                    while(next_entry)
                    {
                        if(next_entry->particle < p)
                        {
                            delta_pos[0] = pos_new[X_COORD(p)] - pos_new[X_COORD(next_entry->particle)];
                            delta_pos[1] = pos_new[Y_COORD(p)] - pos_new[Y_COORD(next_entry->particle)];
                            delta_pos[2] = pos_new[Z_COORD(p)] - pos_new[Z_COORD(next_entry->particle)];
                            dist = norm_squared(delta_pos);

                            if(dist < contact_making_dist_squared)
                            {
                                // check if particles are not already in contact
                                cl_entry = getContactInsertPos(p, next_entry->particle);

                                if(cl_entry)
                                {

                                    // create new contact
                                    double i_dist = 1.0 / sqrt(dist);
                                    delta_pos[0] *= i_dist;
                                    delta_pos[1] *= i_dist;
                                    delta_pos[2] *= i_dist;

                                    Contact *contact = new Contact(p, next_entry->particle, delta_pos);

                                    // add to contact list
                                    cl_entry->contact = contact;

                                    ++created_contacts;
                                    ++number_of_contacts;

#ifdef TRACK_CONTACTS_PER_PARTICLE
                                    ++created_contacts_of_particle[p];
                                    ++created_contacts_of_particle[next_entry->particle];
#endif

#ifdef TRACK_DISSIPATED_ENERGY
                                    dissipated_contact_energy -= normal_interaction.getJKRPotentialEnergy(0.0) * ENERGY_UNIT;
#endif
                                }
                            }
                        }

                        next_entry = next_entry->next;
                    }
                }
            }
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // check if new contacts with walls occurred
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////

        for(int w = 0; w < number_of_walls; ++w)
        {
            dist = walls[w].getDistanceTo(pos_new[X_COORD(p)], pos_new[Y_COORD(p)], pos_new[Z_COORD(p)]);


            // filter out particles that are above the side walls if box is open
            if(number_of_walls == 5 && w != box->bottom_wall_id)
            {
                if(pos_new[Y_COORD(p)] > box->lower_pos_y + box->height)
                    dist = 2.0 * contact_making_dist;
            }

            if(dist < contact_making_dist)
            {

                // check if particle is not already in contact with wall
                cl_entry = getContactInsertPos(p, -w-1);

                if(cl_entry)
                {


                    // delta_pos stores the location, where the contact has been established
                    // this is needed as an anchor to keep track of moving walls
                    delta_pos[0] = pos_new[X_COORD(p)] - dist * walls[w].normal[0] - walls[w].pos[0];
                    delta_pos[1] = pos_new[Y_COORD(p)] - dist * walls[w].normal[1] - walls[w].pos[1];
                    delta_pos[2] = pos_new[Z_COORD(p)] - dist * walls[w].normal[2] - walls[w].pos[2];

                    // create new contact
                    Contact *contact = new Contact(p, w, walls[w].normal, delta_pos);

                #ifdef SURFACE_ASPERITIES
                    surfaceAsperitiesDampingWall(p, w);
                #endif

                    // append contact
                    cl_entry->contact = contact;

                    ++created_wall_contacts;
                    ++number_of_wall_contacts;

#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_wall_energy -= normal_interaction.getJKRWallPotentialEnergy(0.0) * ENERGY_UNIT;
#endif
                }
            }
        }
    }
}

ContactListEntry* Simulation::getContactInsertPos(int id1, int id2)
{
    // append at the beginning
    if(contact_list[id1] == NULL)
    {
        contact_list[id1] = new ContactListEntry;
        contact_list[id1]->next = NULL;
        contact_list[id1]->id = id2;
        return contact_list[id1];
    }

    ContactListEntry *cl_entry = contact_list[id1];

    // check if already in contact
    while(true)
    {
        // already in contact -> do not append
        if(cl_entry->id == id2)
            return NULL;
        else if(cl_entry->next == NULL)
            break;
        else
            cl_entry = cl_entry->next;
    }

    // append at the end
    cl_entry->next = new ContactListEntry;
    cl_entry->next->next = NULL;
    cl_entry->next->id = id2;
    return cl_entry->next;
}

#ifdef SURFACE_ASPERITIES
/* this function seems out of date
void Simulation::surfaceAsperitiesDamping(int id1, int id2)
{
    vec3 v_normal, n_c;

    // calculate relative velocity
    v_normal[0] = vel[4*id1] - vel[4*id2];
    v_normal[1] = vel[4*id1+1] - vel[4*id2+1];
    v_normal[2] = vel[4*id1+2] - vel[4*id2+2];

    // calculate normal
    n_c[0] = pos_new[4*id1] - pos_new[4*id2];
    n_c[1] = pos_new[4*id1+1] - pos_new[4*id2+1];
    n_c[2] = pos_new[4*id1+2] - pos_new[4*id2+2];
    normalize(&n_c);

    // calculate normal component of velocity
    double temp = dot_product(v_normal, n_c);
    v_normal[0] = temp * n_c[0];
    v_normal[1] = temp * n_c[1];
    v_normal[2] = temp * n_c[2];

    // impact energy
    double E_impact = 0.5 * mass * norm_squared(v_normal);

    // calculate maximum contact radius
    double F_impact = pow(15.625 * E_impact*E_impact*E_impact * reduced_radius * young_mod_reduced*young_mod_reduced, 0.2);//pow(2.5 * E, 0.6) * pow(reduced_radius, 0.2) * pow(young_mod_reduced, 0.4);
    temp = 6.0 * M_PI * surface_energy * reduced_radius;

    double a_max = pow( ( F_impact + temp + sqrt(temp*temp + 2.0 * temp * F_impact) ) * 0.75 * reduced_radius / young_mod_reduced, 1.0/3.0);

    // determine total volume of surface bumps in contact area
    double bump_radius = 1e-7; //0.002 * particle_radius;
    double V_tot = 4.0 / 3.0 * M_PI * bump_radius * a_max*a_max * surface_bumpiness;

    // calculate dissipated energy (*10 to convert yield_strength from Pa to Barye)
    double E_damp = 10.0 * yield_strength * V_tot;
    double E_kin = 0.5 * mass * ( vel[4*id1]*vel[4*id1] + vel[4*id1+1]*vel[4*id1+1] + vel[4*id1+2]*vel[4*id1+2] + vel[4*id2]*vel[4*id2] + vel[4*id2+1]*vel[4*id2+1] + vel[4*id2*vel[4*id2);

    if(E_kin <= 0)
        return;

#ifdef TRACK_DISSIPATED_ENERGY
    dissipated_damping_energy += (E_damp * ENERGY_UNIT);
#endif

    double k = sqrt(1.0 - E_damp / E_kin );

    // calculate damped velocity vectors
    vel[4*id1] *= k;
    vel[4*id1+1] *= k;
    vel[4*id1+2] *= k;

    vel[4*id2] *= k;
    vel[4*id2+1] *= k;
    vel[4*id2+2] *= k;
}
*/

void Simulation::surfaceAsperitiesDamping(int id1, int id2)
{
    vec3 v_normal, n_c;

    // calculate relative velocity
    v_normal[0] = vel[X_COORD(id1)] - vel[X_COORD(id2)];
    v_normal[1] = vel[Y_COORD(id1)] - vel[Y_COORD(id2)];
    v_normal[2] = vel[Z_COORD(id1)] - vel[Z_COORD(id2)];

    // calculate normal
    n_c[0] = pos_new[X_COORD(id1)] - pos_new[X_COORD(id2)];
    n_c[1] = pos_new[Y_COORD(id1)] - pos_new[Y_COORD(id2)];
    n_c[2] = pos_new[Z_COORD(id1)] - pos_new[Z_COORD(id2)];
    normalize(&n_c);

    // calculate normal component of velocity
    double temp = dot_product(v_normal, n_c);
    v_normal[0] = temp * n_c[0];
    v_normal[1] = temp * n_c[1];
    v_normal[2] = temp * n_c[2];

    // impact energy
    double E_impact = 0.5 * mass * norm_squared(v_normal);

    // calculate maximum contact radius
    double F_impact = pow(15.625 * E_impact*E_impact*E_impact * reduced_radius * young_mod_reduced*young_mod_reduced, 0.2);//pow(2.5 * E, 0.6) * pow(reduced_radius, 0.2) * pow(young_mod_reduced, 0.4);
    temp = 6.0 * M_PI * surface_energy * reduced_radius;

    double a_max = pow( ( F_impact + temp + sqrt(temp*temp + 2.0 * temp * F_impact) ) * 0.75 * reduced_radius / young_mod_reduced, 1.0/3.0);

    // determine total volume of surface bumps in contact area
    double bump_radius = 1e-7; //0.002 * particle_radius;
    double surface_bumpiness = 1.0;
    double V_tot = 4.0 / 3.0 * M_PI * bump_radius * a_max*a_max * surface_bumpiness;

    // calculate dissipated energy (*10 to convert yield_strength from Pa to Barye)
    double yield_strength = 104.e6;
    double E_damp = 10.0 * yield_strength * V_tot;
    double E_kin = 0.5 * mass * (
                vel[X_COORD(id1)]*vel[X_COORD(id1)] +
                vel[Y_COORD(id1)]*vel[Y_COORD(id1)] +
                vel[Z_COORD(id1)]*vel[Z_COORD(id1)] +

                vel[X_COORD(id2)]*vel[X_COORD(id2)] +
                vel[Y_COORD(id2)]*vel[Y_COORD(id2)] +
                vel[Z_COORD(id2)]*vel[Z_COORD(id2)]);

    if(E_kin <= 0)
        return;

#ifdef TRACK_DISSIPATED_ENERGY
    dissipated_damping_energy += (E_damp * ENERGY_UNIT);
#endif

    double k = sqrt(1.0 - E_damp / E_kin );

    // calculate damped velocity vectors
    vel[X_COORD(id1)] *= k;
    vel[Y_COORD(id1)] *= k;
    vel[Z_COORD(id1)] *= k;

    vel[X_COORD(id2)] *= k;
    vel[Y_COORD(id2)] *= k;
    vel[Z_COORD(id2)] *= k;

}

void Simulation::surfaceAsperitiesDampingWall(int pid, int wid)
{
    vec3 v_normal;

    // calculate relative velocity
    v_normal[0] = vel[X_COORD(pid)] - walls[wid].velocity[0];
    v_normal[1] = vel[Y_COORD(pid)] - walls[wid].velocity[1];
    v_normal[2] = vel[Z_COORD(pid)] - walls[wid].velocity[2];



    // calculate normal component of velocity
    double v_n = dot_product(v_normal, walls[wid].normal);


    // impact energy
    double E_impact = 0.5 * mass * norm_squared(v_normal);

    // calculate maximum contact radius
    double F_impact = pow(15.625 * E_impact*E_impact*E_impact * wall_reduced_radius * young_mod_reduced*young_mod_reduced, 0.2);//pow(2.5 * E, 0.6) * pow(reduced_radius, 0.2) * pow(young_mod_reduced, 0.4);
    double temp = 6.0 * M_PI * surface_energy * wall_reduced_radius;

    double a_max = pow( ( F_impact + temp + sqrt(temp*temp + 2.0 * temp * F_impact) ) * 0.75 * wall_reduced_radius / young_mod_reduced, 1.0/3.0);

    // determine total volume of surface bumps in contact area
    double bump_radius = 1e-7; //0.002 * particle_radius;
    double surface_bumpiness = 10.e9;
    double area = M_PI*a_max*a_max;


    double V_tot = 4.0 / 3.0 * M_PI * bump_radius * a_max*a_max;
    double V_tot2 = 2.0/3.0 * M_PI * bump_radius * 0.02 * a_max*a_max;
    printf("V_tot = %e  %e\n", V_tot, V_tot2);

    // calculate dissipated energy (*10 to convert yield_strength from Pa to Barye)
    double yield_strength = 104.e6;
    double E_damp = 10.0 * yield_strength * V_tot;
    double E_kin = 0.5 * mass * v_n * v_n;


    double k = 0.0;

    if(E_kin <= E_damp)
        k = 0.0;
    else
    {

#ifdef TRACK_DISSIPATED_ENERGY
     dissipated_damping_energy += (E_damp * ENERGY_UNIT);
#endif


     k = sqrt(1.0 - E_damp / E_kin );
    }
    printf("total area = %e k = %e\n", area*surface_bumpiness, k);


    // calculate damped velocity vectors
    vel[X_COORD(pid)] *= k;
    vel[Y_COORD(pid)] *= k;
    vel[Z_COORD(pid)] *= k;

}

#endif

void Simulation::updateDissipatedBreakingEnergy(Contact *contact, double particle_dist)
{
    vec3 n1, n2, n_c, vec1, vec2;

    // check if contact with wall
    if(contact->id2 < 0)
    {
        // calculate current contact pointers
        contact->getCurrentN1(&n1);
        memcpy(n2, walls[-contact->id2-1].normal,sizeof(vec3));

        dissipated_wall_energy += normal_interaction.getJKRWallPotentialEnergy(- wall_delta_c) * ENERGY_UNIT;

#ifdef WALL_ROLLING
        vec1[0] = n1[0] + n2[0];
        vec1[1] = n1[1] + n2[1];
        vec1[2] = n1[2] + n2[2];
        dissipated_wall_energy += 0.5 * wall_k_r * particle_radius*particle_radius * norm_squared(vec1) * ENERGY_UNIT;
#endif

#ifdef WALL_SLIDING
        vec1[0] = pos_new[X_COORD(contact->id1)] - (contact->n2_initial[0] + walls[-contact->id2-1].pos[0]) + particle_radius * n1[0];
        vec1[1] = pos_new[Y_COORD(contact->id1)] - (contact->n2_initial[1] + walls[-contact->id2-1].pos[1]) + particle_radius * n1[1];
        vec1[2] = pos_new[Z_COORD(contact->id1)] - (contact->n2_initial[2] + walls[-contact->id2-1].pos[2]) + particle_radius * n1[2];

        double temp = dot_product(vec1, n2);
        vec2[0] = vec1[0] - temp * n2[0];
        vec2[1] = vec1[1] - temp * n2[1];
        vec2[2] = vec1[2] - temp * n2[2];

        dissipated_wall_energy += 0.5 * wall_k_s * norm_squared(vec2) * ENERGY_UNIT;
#endif

#ifdef WALL_TWISTING
        dissipated_wall_energy += 0.5 * wall_k_t * contact->twisting_displacement * contact->twisting_displacement * ENERGY_UNIT;
#endif
    }
    else // contact with other particle
    {
        contact->getCurrentN1(&n1);
        contact->getCurrentN2(&n2);

        n_c[0] = particle_dist * (pos_new[X_COORD(contact->id1)] - pos_new[X_COORD(contact->id2)]);
        n_c[1] = particle_dist * (pos_new[Y_COORD(contact->id1)] - pos_new[Y_COORD(contact->id2)]);
        n_c[2] = particle_dist * (pos_new[Z_COORD(contact->id1)] - pos_new[Z_COORD(contact->id2)]);

        dissipated_contact_energy += normal_interaction.getJKRPotentialEnergy(- delta_c) * ENERGY_UNIT;

#ifdef ROLLING
        vec1[0] = n1[0] + n2[0];
        vec1[1] = n1[1] + n2[1];
        vec1[2] = n1[2] + n2[2];
        dissipated_contact_energy += 0.5 * k_r * reduced_radius * reduced_radius * norm_squared(vec1) * ENERGY_UNIT;
#endif

#ifdef SLIDING
        vec3 delta_n;
        /* Lucas: probably error by Alex, correction below
        delta_n[0] = (n1[0] - n2[0]);
        delta_n[1] = (n1[1] - n2[1]);
        delta_n[2] = (n1[2] - n2[2]);
        */

        delta_n[0] = (n1[0] - n2[0] + 2.0*n_c[0]); // == \xi_0, helper vector for sliding displacement
        delta_n[1] = (n1[1] - n2[1] + 2.0*n_c[1]);
        delta_n[2] = (n1[2] - n2[2] + 2.0*n_c[2]);


        double temp = dot_product(delta_n, n_c);
        vec2[0] = delta_n[0] - temp * n_c[0]; // == sliding displacement \xi
        vec2[1] = delta_n[1] - temp * n_c[1];
        vec2[2] = delta_n[2] - temp * n_c[2];

        dissipated_contact_energy += 0.5 * k_s * particle_radius * particle_radius * norm_squared(vec2) * ENERGY_UNIT;
#endif

#ifdef TWISTING
        dissipated_contact_energy += 0.5 * k_t * contact->twisting_displacement * contact->twisting_displacement * ENERGY_UNIT;
#endif
    }
}

void Simulation::updateParticleInteraction(double timestep)
{
    // init forces & torques with 0
    memset(force_new, 0, 3*number_of_particles * sizeof(double));
    memset(torque_new, 0, 3*number_of_particles * sizeof(double));

    for(int w = 0; w < number_of_walls; ++w)
    {
        memset(walls[w].total_force, 0, sizeof(vec3));
        memset(walls[w].box_volume_force, 0, sizeof(vec3));
        memset(walls[w].total_torque, 0, sizeof(vec3));
    }

#ifdef ENABLE_GRAVITY
    if(gravity_enabled)
    {
        for(int p = 0; p < number_of_particles; ++p)
        {
            force_new[X_COORD(p)] += mass * gravity_modifier * gravity_strength * gravity_direction[0];
            force_new[Y_COORD(p)] += mass * gravity_modifier * gravity_strength * gravity_direction[1];
            force_new[Z_COORD(p)] += mass * gravity_modifier * gravity_strength * gravity_direction[2];
        }

        for(int w = 0; w < number_of_walls; ++w)
        {
            walls[w].total_force[0] = walls[w].mass * gravity_modifier * gravity_strength * gravity_direction[0];
            walls[w].total_force[1] = walls[w].mass * gravity_modifier * gravity_strength * gravity_direction[1];
            walls[w].total_force[2] = walls[w].mass * gravity_modifier * gravity_strength * gravity_direction[2];
        }
    }
#endif

    double particle_distance;
    double displacement_norm;
    double compression_length;
    double force;
    double alpha1, alpha2, inv_norm, temp;
    vec3 n1, n2, delta_n, n_c;
    vec3 temp_vec, displacement, displacement_correction;
    ContactListEntry *next_cl_entry;


#ifdef TRACK_FORCES
    char filename[1024];
    sprintf(filename, "../surface_energy_%d_timestep_%d.dat", int(surface_energy), int(timestep/1.e-14));
    FILE* track_forces = fopen(filename, "a");

    if(!track_forces)
        return;

#endif // TRACK_FORCES

    //#pragma omp for private(next_cl_entry, n1, n2, n_c, delta_n, displacement, displacement_correction, temp_vec, alpha1, alpha2, inv_norm, temp, particle_distance, displacement_norm)
    for(int i = 0; i < number_of_particles; ++i)
    {

        next_cl_entry = contact_list[i];

        while(next_cl_entry)
        {
            if(next_cl_entry->id >= 0)
            {
                // calculate current contact pointers
                next_cl_entry->contact->getCurrentN1(&n1);
                next_cl_entry->contact->getCurrentN2(&n2);
                delta_n[0] = n1[0] - n2[0];
                delta_n[1] = n1[1] - n2[1];
                delta_n[2] = n1[2] - n2[2];

                // determine distance between particles & contact normal
                n_c[0] = pos_new[X_COORD(i)] - pos_new[X_COORD(next_cl_entry->id)];
                n_c[1] = pos_new[Y_COORD(i)] - pos_new[Y_COORD(next_cl_entry->id)];
                n_c[2] = pos_new[Z_COORD(i)] - pos_new[Z_COORD(next_cl_entry->id)];
                particle_distance = norm(n_c);

                double tmp = 1.0 / particle_distance;

                n_c[0] *= tmp;
                n_c[1] *= tmp;
                n_c[2] *= tmp;

                // determine compression length
                compression_length = 2.0 * particle_radius - particle_distance;

                double contact_radius = normal_interaction.getJKRContactRadius(compression_length);


#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
                double k_r2 = k_r * pow(contact_radius / equilibrium_contact_radius, 1.5);
                double k_s2 = k_s * contact_radius;
                double k_t2 = k_t * contact_radius * contact_radius * contact_radius;
                //double k_r2 = k_r * pow(equilibrium_contact_radius / equilibrium_contact_radius, 1.5);
                //double k_s2 = k_s * equilibrium_contact_radius;
                //double k_t2 = k_t * equilibrium_contact_radius * equilibrium_contact_radius * equilibrium_contact_radius;
#endif

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // force in normal direction (compression/sticking)
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // get force for current compression length according to JKR/DMT model
#ifdef USE_DMT
                force = normal_interaction.getDMTForce(compression_length);
#endif

#ifdef USE_JKR
  #ifdef INTERPOLATE_NORMAL_FORCES
                force = normal_interaction.getJKRApproxForce(compression_length);
  #else
                force = normal_interaction.getJKRForce(contact_radius);

                /*
                double P_mean = force / (M_PI * contact_radius * contact_radius);
                double P_max = 6.1e9 * 10; // *10 to convert from SI to cgs

                if(P_mean > P_max*0.1)
                {
                    printf("%d  %e  \n", i, double(P_mean/P_max));
                }
                */

  #endif
#endif


#ifdef TRACK_FORCES
                if(i == 0 || i == 1)
                {
                    fprintf(track_forces, "%e	%e %e    ", current_time, compression_length/delta_c, force/F_c);
                }
#endif // TRACK_FORCES

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // additional gluing of particles hat are close to walls
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef ENABLE_WALL_GLUE
                if(sim_info.sim_type == SIM_TYPE_PULL_STRENGTH_TEST || sim_info.sim_type == SIM_TYPE_SHEAR_STRENGTH_TEST)
                {
                    double wall_dist = 2.0 * particle_radius * wall_glue_distance;

                    double temp = walls[box->top_wall_id].getDistanceTo(pos_new[X_COORD(i)], pos_new[Y_COORD(i)], pos_new[Z_COORD(i)]);
                    if(temp < wall_dist)
                        wall_dist = temp;

                    temp = walls[box->bottom_wall_id].getDistanceTo(pos_new[X_COORD(i)], pos_new[Y_COORD(i)], pos_new[Z_COORD(i)]);
                    if(temp < wall_dist)
                        wall_dist = temp;

                    if(wall_dist < particle_radius * wall_glue_distance)
                    {
                        double k = wall_glue_strength * (particle_radius * wall_glue_distance - wall_dist) / (particle_radius * wall_glue_distance);
                        force *= (1.0 + k);
                    }
                }
#endif

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // damping of normal oscillations
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_CONSTANT_DAMPING
                double v_rel = (compression_length - next_cl_entry->contact->compression_length) /  timestep;
                force += osc_damping_factor * v_rel;

#ifdef TRACK_DISSIPATED_ENERGY
                dissipated_damping_energy += osc_damping_factor * v_rel*v_rel * timestep * ENERGY_UNIT;
#endif
#endif

#ifdef USE_VISCOELASTIC_DAMPING
                double v_rel = (compression_length - next_cl_entry->contact->compression_length) /  timestep;
                force += viscoelastic_damping_constant * contact_radius * v_rel;

#ifdef TRACK_DISSIPATED_ENERGY
                dissipated_damping_energy += viscoelastic_damping_constant * contact_radius * v_rel*v_rel * timestep * ENERGY_UNIT;
#endif

#ifdef TRACK_FORCES
                if(i == 0 || i == 1)
                {
                    fprintf(track_forces, "%e %e    %e  ", force/F_c, v_rel, contact_radius);
                }
#endif // TRACK_FORCES

#endif


                force_new[X_COORD(i)] += force * n_c[0];
                force_new[Y_COORD(i)] += force * n_c[1];
                force_new[Z_COORD(i)] += force * n_c[2];


                force_new[X_COORD(next_cl_entry->id)] -=  force * n_c[0];
                force_new[Y_COORD(next_cl_entry->id)] -=  force * n_c[1];
                force_new[Z_COORD(next_cl_entry->id)] -=  force * n_c[2];

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // apply special treatment for inelastic rolling, sliding & twisting
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if defined(INELASTIC_ROLLING) || defined(INELASTIC_SLIDING)
                bool update_contact_pointers = false;
#endif

#ifdef INELASTIC_SLIDING
                temp = dot_product(delta_n,  n_c);

                displacement[0] = particle_radius * (delta_n[0]  - temp * n_c[0]);
                displacement[1] = particle_radius * (delta_n[1]  - temp * n_c[1]);
                displacement[2] = particle_radius * (delta_n[2]  - temp * n_c[2]);

                displacement_norm = norm_squared(displacement);

                // check if we are in the inelastic regime
                if(displacement_norm > crit_sliding_displacement_squared)
                {
                    displacement_norm = sqrt(displacement_norm);

                    // determine correction of contact pointers
                    displacement_correction[0] = (1.0 - crit_sliding_displacement / displacement_norm) * displacement[0];
                    displacement_correction[1] = (1.0 - crit_sliding_displacement / displacement_norm) * displacement[1];
                    displacement_correction[2] = (1.0 - crit_sliding_displacement / displacement_norm) * displacement[2];

                    // calculate correction factor (see Wada et al. 2007 appendix for details)
                    inv_norm = 1.0 / norm_squared(displacement_correction);

                    temp = dot_product(n1, displacement_correction);
                    alpha1 = 1.0 / (1.0 - temp * temp * inv_norm);

                    temp = dot_product(n2, displacement_correction);
                    alpha2 = 1.0 / (1.0 - temp * temp * inv_norm);

                    // correct contact pointers
                    n1[0] -= 0.5 * particle_radius_inv * alpha1 * displacement_correction[0];
                    n1[1] -= 0.5 * particle_radius_inv * alpha1 * displacement_correction[1];
                    n1[2] -= 0.5 * particle_radius_inv * alpha1 * displacement_correction[2];

                    n2[0] += 0.5 * particle_radius_inv * alpha2 * displacement_correction[0];
                    n2[1] += 0.5 * particle_radius_inv * alpha2 * displacement_correction[1];
                    n2[2] += 0.5 * particle_radius_inv * alpha2 * displacement_correction[2];

                    normalize(&n1);
                    normalize(&n2);

                    update_contact_pointers = true;

#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_sliding_energy += k_s * crit_sliding_displacement * (displacement_norm - crit_sliding_displacement) * ENERGY_UNIT;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
                    dissipated_energy_of_particle[i] += k_s * crit_sliding_displacement * (displacement_norm - crit_sliding_displacement) * ENERGY_UNIT;
                    dissipated_energy_of_particle[next_cl_entry->id] += k_s * crit_sliding_displacement * (displacement_norm - crit_sliding_displacement) * ENERGY_UNIT;
#endif
                }
#endif	// INELASTIC_SLIDING

#ifdef INELASTIC_ROLLING
                displacement[0] = reduced_radius * (n1[0] + n2[0]);
                displacement[1] = reduced_radius * (n1[1] + n2[1]);
                displacement[2] = reduced_radius * (n1[2] + n2[2]);
                displacement_norm = norm_squared(displacement);

                // special treatment of contact pointers if we are in the inelastic regime
                if(displacement_norm > crit_rolling_displacement_squared)
                {
                    displacement_norm = sqrt(displacement_norm);

                    // determine correction of contact pointers
                    displacement_correction[0] = (1.0 - crit_rolling_displacement / displacement_norm) * displacement[0];
                    displacement_correction[1] = (1.0 - crit_rolling_displacement / displacement_norm) * displacement[1];
                    displacement_correction[2] = (1.0 - crit_rolling_displacement / displacement_norm) * displacement[2];

                    // calculate correction factor (see Wada et al. 2007 appendix for details)
                    inv_norm = 1.0 / norm_squared(displacement_correction);

                    temp = dot_product(n1, displacement_correction);
                    alpha1 = 1.0 / (1.0 - temp * temp * inv_norm);

                    temp = dot_product(n2, displacement_correction);
                    alpha2 = 1.0 / (1.0 - temp * temp * inv_norm);

                    // correct contact pointers
                    n1[0] -= particle_radius_inv * alpha1 * displacement_correction[0];
                    n1[1] -= particle_radius_inv * alpha1 * displacement_correction[1];
                    n1[2] -= particle_radius_inv * alpha1 * displacement_correction[2];
                    n2[0] -= particle_radius_inv * alpha2 * displacement_correction[0];
                    n2[1] -= particle_radius_inv * alpha2 * displacement_correction[1];
                    n2[2] -= particle_radius_inv * alpha2 * displacement_correction[2];

                    normalize(&n1);
                    normalize(&n2);

                    update_contact_pointers = true;

#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_rolling_energy += k_r * crit_rolling_displacement * ( displacement_norm - crit_rolling_displacement) * ENERGY_UNIT;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
                    dissipated_energy_of_particle[i] += k_r * crit_rolling_displacement * (displacement_norm - crit_rolling_displacement) * ENERGY_UNIT;
                    dissipated_energy_of_particle[next_cl_entry->id] += k_r * crit_rolling_displacement * (displacement_norm - crit_rolling_displacement) * ENERGY_UNIT;
#endif
                }
#endif	// INELASTIC_ROLLING

#if defined(INELASTIC_ROLLING) || defined(INELASTIC_SLIDING)
                if(update_contact_pointers)
                {
                    next_cl_entry->contact->updateN1Initial(n1);
                    next_cl_entry->contact->updateN2Initial(n2);

                    delta_n[0] = n1[0] - n2[0];
                    delta_n[1] = n1[1] - n2[1];
                    delta_n[2] = n1[2] - n2[2];
                }
#endif

                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // calculate forces/torques according to the new positions/rotations
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef SLIDING
                temp = dot_product(delta_n, n_c);
                temp_vec[0] = delta_n[0] - temp * n_c[0];
                temp_vec[1] = delta_n[1] - temp * n_c[1];
                temp_vec[2] = delta_n[2] - temp * n_c[2];

                // force in tangential direction (sliding)
                double interaction_strength = K_S * particle_radius * particle_radius / particle_distance * temp;


#ifdef TRACK_FORCES
                if(i == 0 || i == 1)
                {
                    fprintf(track_forces, "%e", interaction_strength * norm(temp_vec));
                }
#endif // TRACK_FORCES



                force_new[X_COORD(i)] += interaction_strength * temp_vec[0];
                force_new[Y_COORD(i)] += interaction_strength * temp_vec[1];
                force_new[Z_COORD(i)] += interaction_strength * temp_vec[2];

                force_new[X_COORD(next_cl_entry->id)] -= interaction_strength * temp_vec[0];
                force_new[Y_COORD(next_cl_entry->id)] -= interaction_strength * temp_vec[1];
                force_new[Z_COORD(next_cl_entry->id)] -= interaction_strength * temp_vec[2];

                // torque (sliding)
                interaction_strength = K_S * particle_radius * particle_radius;

                torque_new[X_COORD(i)] -= interaction_strength * (n1[1] * temp_vec[2] - n1[2] * temp_vec[1]);
                torque_new[Y_COORD(i)] -= interaction_strength * (n1[2] * temp_vec[0] - n1[0] * temp_vec[2]);
                torque_new[Z_COORD(i)] -= interaction_strength * (n1[0] * temp_vec[1] - n1[1] * temp_vec[0]);

                torque_new[X_COORD(next_cl_entry->id)] += interaction_strength * (n2[1] * temp_vec[2] - n2[2] * temp_vec[1]);
                torque_new[Y_COORD(next_cl_entry->id)] += interaction_strength * (n2[2] * temp_vec[0] - n2[0] * temp_vec[2]);
                torque_new[Z_COORD(next_cl_entry->id)] += interaction_strength * (n2[0] * temp_vec[1] - n2[1] * temp_vec[0]);

#endif

#ifdef ROLLING
                temp_vec[0] = K_R * reduced_radius * reduced_radius * (n1[1] * n2[2] - n1[2] * n2[1]);
                temp_vec[1] = K_R * reduced_radius * reduced_radius * (n1[2] * n2[0] - n1[0] * n2[2]);
                temp_vec[2] = K_R * reduced_radius * reduced_radius * (n1[0] * n2[1] - n1[1] * n2[0]);

                torque_new[X_COORD(i)] -= temp_vec[0];
                torque_new[Y_COORD(i)] -= temp_vec[1];
                torque_new[Z_COORD(i)] -= temp_vec[2];

                torque_new[X_COORD(next_cl_entry->id)] += temp_vec[0];
                torque_new[Y_COORD(next_cl_entry->id)] += temp_vec[1];
                torque_new[Z_COORD(next_cl_entry->id)] += temp_vec[2];

#endif

#ifdef TWISTING
                // update twisting displacement - use second order integration: Phi^n+1 = Phi^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
                vec3 delta_omega_old, delta_omega_new;

                delta_omega_old[0] = vel_angular[X_COORD(i)] - vel_angular[X_COORD(next_cl_entry->id)];
                delta_omega_old[1] = vel_angular[Y_COORD(i)] - vel_angular[Y_COORD(next_cl_entry->id)];
                delta_omega_old[2] = vel_angular[Z_COORD(i)] - vel_angular[Z_COORD(next_cl_entry->id)];

                delta_omega_new[0] = delta_omega_old[0] + timestep * moment_of_inertia_inv * (torque_old[X_COORD(i)] - torque_old[X_COORD(next_cl_entry->id)]);
                delta_omega_new[1] = delta_omega_old[1] + timestep * moment_of_inertia_inv * (torque_old[Y_COORD(i)] - torque_old[Y_COORD(next_cl_entry->id)]);
                delta_omega_new[2] = delta_omega_old[2] + timestep * moment_of_inertia_inv * (torque_old[Z_COORD(i)] - torque_old[Z_COORD(next_cl_entry->id)]);


                next_cl_entry->contact->twisting_displacement += 0.5 * timestep * (dot_product(delta_omega_old, next_cl_entry->contact->old_contact_normal) + dot_product(delta_omega_new, n_c));



#ifdef INELASTIC_TWISTING
                if(next_cl_entry->contact->twisting_displacement > crit_twisting_displacement)
                {
#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_twisting_energy += k_t * crit_twisting_displacement * (next_cl_entry->contact->twisting_displacement - crit_twisting_displacement) * ENERGY_UNIT;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
                    dissipated_energy_of_particle[i] += k_t * crit_twisting_displacement * (next_cl_entry->contact->twisting_displacement - crit_twisting_displacement) * ENERGY_UNIT;
                    dissipated_energy_of_particle[next_cl_entry->id] += k_t * crit_twisting_displacement * (next_cl_entry->contact->twisting_displacement - crit_twisting_displacement) * ENERGY_UNIT;
#endif

                    next_cl_entry->contact->twisting_displacement = crit_twisting_displacement;
                }
#endif // INELASTIC_TWISTING

                torque_new[X_COORD(i)] -= K_T * next_cl_entry->contact->twisting_displacement * n_c[0];
                torque_new[Y_COORD(i)] -= K_T * next_cl_entry->contact->twisting_displacement * n_c[1];
                torque_new[Z_COORD(i)] -= K_T * next_cl_entry->contact->twisting_displacement * n_c[2];
                torque_new[X_COORD(next_cl_entry->id)] += K_T * next_cl_entry->contact->twisting_displacement * n_c[0];
                torque_new[Y_COORD(next_cl_entry->id)] += K_T * next_cl_entry->contact->twisting_displacement * n_c[1];
                torque_new[Z_COORD(next_cl_entry->id)] += K_T * next_cl_entry->contact->twisting_displacement * n_c[2];

#endif // TWISTING


                // store current contact normal for next update step
                next_cl_entry->contact->old_contact_normal[0] = n_c[0];
                next_cl_entry->contact->old_contact_normal[1] = n_c[1];
                next_cl_entry->contact->old_contact_normal[2] = n_c[2];
                next_cl_entry->contact->compression_length = compression_length;

            }
            else
            {
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                // calculate forces/torques caused by collisions with walls
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                int wall_id = WALL_ID(next_cl_entry->id);

                next_cl_entry->contact->getCurrentN1(&n1);
                memcpy(n2, walls[wall_id].normal, sizeof(vec3));


                vec3 wall_contact_pos;
                // n2_initial is the position where the contact was made minus the position of the wall at that time
                wall_contact_pos[0] = next_cl_entry->contact->n2_initial[0] + walls[wall_id].pos[0];
                wall_contact_pos[1] = next_cl_entry->contact->n2_initial[1] + walls[wall_id].pos[1];
                wall_contact_pos[2] = next_cl_entry->contact->n2_initial[2] + walls[wall_id].pos[2];


                n_c[0] = pos_new[X_COORD(i)] - wall_contact_pos[0];
                n_c[1] = pos_new[Y_COORD(i)] - wall_contact_pos[1];
                n_c[2] = pos_new[Z_COORD(i)] - wall_contact_pos[2];
                normalize(&n_c);

                // force in normal direction
                compression_length = 2.0 * particle_radius - walls[wall_id].getDistanceTo(pos_new[X_COORD(i)], pos_new[Y_COORD(i)], pos_new[Z_COORD(i)]);


#ifdef USE_DMT
                double contact_radius = normal_interaction.getDMTWallContactRadius(compression_length);
                force = walls[wall_id].compression_modifier * normal_interaction.getDMTWallForce(compression_length);
#endif

#ifdef USE_JKR
                double contact_radius = normal_interaction.getJKRWallContactRadius(compression_length);
                force = walls[wall_id].compression_modifier * normal_interaction.getJKRWallForce(contact_radius);
#endif

#ifdef TRACK_FORCES
                if(i == 0)
                {
                    fprintf(track_forces, "%e	%e %e    ", current_time, compression_length, force);
                }
#endif // TRACK_FORCES


                // track force on wall before oscillation damping is applied
                walls[wall_id].total_force[0] -= force * walls[wall_id].normal[0];
                walls[wall_id].total_force[1] -= force * walls[wall_id].normal[1];
                walls[wall_id].total_force[2] -= force * walls[wall_id].normal[2];


                if(sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
                {
                    if(particleInBox(i))
                    {
                        walls[wall_id].box_volume_force[0] -= force * walls[wall_id].normal[0];
                        walls[wall_id].box_volume_force[1] -= force * walls[wall_id].normal[1];
                        walls[wall_id].box_volume_force[2] -= force * walls[wall_id].normal[2];
                    }
                }


#ifdef USE_VISCOELASTIC_DAMPING
    // calculate relative velocity
    temp_vec[0] = vel[X_COORD(i)] - walls[wall_id].velocity[0];
    temp_vec[1] = vel[Y_COORD(i)] - walls[wall_id].velocity[1];
    temp_vec[2] = vel[Z_COORD(i)] - walls[wall_id].velocity[2];

    // calculate normal component of velocity & damping force
    double v_rel = dot_product(temp_vec, walls[wall_id].normal);


    force -= viscoelastic_damping_constant * contact_radius * v_rel;


#ifdef TRACK_DISSIPATED_ENERGY
                double v_particle = vel[X_COORD(i)] * walls[wall_id].normal[0] + vel[Y_COORD(i)] * walls[wall_id].normal[1] + vel[Z_COORD(i)] * walls[wall_id].normal[2];
                dissipated_damping_energy += viscoelastic_damping_constant * contact_radius * v_rel * v_particle * timestep * ENERGY_UNIT;
#endif  // TRACK_DISSIPATED_ENERGY
#ifdef TRACK_FORCES
                if(i == 0)
                {
                    fprintf(track_forces, "%e %e    %e  ", v_rel, contact_radius, viscoelastic_damping_constant * contact_radius * v_rel);
                }
#endif // TRACK_FORCES


#endif // USE_VISCOELASTIC_DAMPING



#ifdef DAMP_WALL_OSCILLATIONS
                // calculate relative velocity
                temp_vec[0] = vel[X_COORD(i)] - walls[wall_id].velocity[0];
                temp_vec[1] = vel[Y_COORD(i)] - walls[wall_id].velocity[1];
                temp_vec[2] = vel[Z_COORD(i)] - walls[wall_id].velocity[2];

                // calculate normal component of velocity & damping force
                double v_rel = dot_product(temp_vec, walls[wall_id].normal);

                force -= osc_damping_factor * v_rel;

#ifdef TRACK_DISSIPATED_ENERGY
                double v_particle = vel[X_COORD(i)] * walls[wall_id].normal[0] + vel[Y_COORD(i)] * walls[wall_id].normal[1] + vel[Z_COORD(i)] * walls[wall_id].normal[2];
                dissipated_damping_energy += osc_damping_factor * v_rel * v_particle * timestep * ENERGY_UNIT;
#endif  // TRACK_DISSIPATED_ENERGY
#endif	// DAMP_WALL_OSCILLATIONS

                force_new[X_COORD(i)] += force * walls[wall_id].normal[0];
                force_new[Y_COORD(i)] += force * walls[wall_id].normal[1];
                force_new[Z_COORD(i)] += force * walls[wall_id].normal[2];


#ifdef INELASTIC_WALL_ROLLING
                displacement[0] = particle_radius * (n1[0] + n2[0]);
                displacement[1] = particle_radius * (n1[1] + n2[1]);
                displacement[2] = particle_radius * (n1[2] + n2[2]);
                displacement_norm = norm_squared(displacement);

                // special treatment of contact pointers if we are in the inelastic regime
                if(displacement_norm > crit_rolling_displacement_squared)
                {
                    displacement_norm = sqrt(displacement_norm);

                    // determine correction of contact pointers
                    displacement_correction[0] = (1.0 - crit_rolling_displacement / displacement_norm) * displacement[0];
                    displacement_correction[1] = (1.0 - crit_rolling_displacement / displacement_norm) * displacement[1];
                    displacement_correction[2] = (1.0 - crit_rolling_displacement / displacement_norm) * displacement[2];

                    // calculate correction factor (see Wada et al. 2007 appendix for details)
                    temp = dot_product(n1, displacement_correction);
                    alpha1 = 1.0 / (1.0 - temp * temp / norm_squared(displacement_correction));

                    // store n1 for sliding correction
                    memcpy(temp_vec, n1, sizeof(vec3));

                    // correct contact pointers
                    n1[0] -= particle_radius_inv * alpha1 * displacement_correction[0];
                    n1[1] -= particle_radius_inv * alpha1 * displacement_correction[1];
                    n1[2] -= particle_radius_inv * alpha1 * displacement_correction[2];
                    normalize(&n1);
                    next_cl_entry->contact->updateN1Initial(n1);

                    // additional correction required to keep sliding displacement invariant under change of contact pointer
                    delta_n[0] = n1[0] - temp_vec[0];
                    delta_n[1] = n1[1] - temp_vec[1];
                    delta_n[2] = n1[2] - temp_vec[2];

                    temp = dot_product(delta_n, n2);
                    temp_vec[0] = particle_radius * (delta_n[0] + temp * n2[0]);
                    temp_vec[1] = particle_radius * (delta_n[1] + temp * n2[1]);
                    temp_vec[2] = particle_radius * (delta_n[2] + temp * n2[2]);

                    next_cl_entry->contact->n2_initial[0] += temp_vec[0];
                    next_cl_entry->contact->n2_initial[1] += temp_vec[1];
                    next_cl_entry->contact->n2_initial[2] += temp_vec[2];

#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_wall_energy += walls[wall_id].rolling_modifier * wall_k_r * crit_rolling_displacement * (displacement_norm - crit_rolling_displacement) * ENERGY_UNIT;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
                    dissipated_energy_of_particle[i] += walls[wall_id].rolling_modifier * wall_k_r * crit_rolling_displacement * (displacement_norm - crit_rolling_displacement) * ENERGY_UNIT;
#endif
                }
#endif	// INELASTIC_WALL_ROLLING

#ifdef WALL_TWISTING
                // update twisting displacement - use second order integration: Phi^n+1 = Phi^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
                // omega1_dot, omega2_dot are used to store delta_omega^n and delta_omega^n+1 respectively
                vec3 n_c_old;
                n_c_old[0] = pos_old[X_COORD(i)] - wall_contact_pos[0];
                n_c_old[1] = pos_old[Y_COORD(i)] - wall_contact_pos[1];
                n_c_old[2] = pos_old[Z_COORD(i)] - wall_contact_pos[2];
                normalize(&n_c_old);

                vec3 n_c_new;
                n_c_new[0] = pos_new[X_COORD(i)] - wall_contact_pos[0];
                n_c_new[1] = pos_new[Y_COORD(i)] - wall_contact_pos[1];
                n_c_new[2] = pos_new[Z_COORD(i)] - wall_contact_pos[2];
                normalize(&n_c_new);

                vec3 omega2_dot;
                omega2_dot[0] = vel_angular[X_COORD(i)] + timestep * moment_of_inertia_inv * torque_old[X_COORD(i)];
                omega2_dot[1] = vel_angular[Y_COORD(i)] + timestep * moment_of_inertia_inv * torque_old[Y_COORD(i)];
                omega2_dot[2] = vel_angular[Z_COORD(i)] + timestep * moment_of_inertia_inv * torque_old[Z_COORD(i)];

                next_cl_entry->contact->twisting_displacement += 0.5 * timestep *
                    (
                    vel_angular[X_COORD(i)] * n_c_old[0] +
                    vel_angular[Y_COORD(i)] * n_c_old[1] +
                    vel_angular[Z_COORD(i)] * n_c_old[2] +
                    omega2_dot[0] * n_c_new[0] +
                    omega2_dot[1] * n_c_new[1] +
                    omega2_dot[2] * n_c_new[2]
                    );


#ifdef INELASTIC_WALL_TWISTING
                if(next_cl_entry->contact->twisting_displacement > crit_twisting_displacement)
                {
#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_wall_energy += walls[wall_id].twisting_modifier * wall_k_t * crit_twisting_displacement * (next_cl_entry->contact->twisting_displacement - crit_twisting_displacement) * ENERGY_UNIT;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
                    dissipated_energy_of_particle[i] += walls[wall_id].twisting_modifier * wall_k_t * crit_twisting_displacement * (next_cl_entry->contact->twisting_displacement - crit_twisting_displacement) * ENERGY_UNIT;
#endif
                    next_cl_entry->contact->twisting_displacement = crit_twisting_displacement;
                }
#endif // INELASTIC_WALL_TWISTING
#endif // WALL_TWISTING

#ifdef WALL_ROLLING
                // torque rolling
                temp_vec[0] = walls[wall_id].rolling_modifier * wall_k_r * particle_radius*particle_radius * (n1[1] * n2[2] - n1[2] * n2[1]);
                temp_vec[1] = walls[wall_id].rolling_modifier * wall_k_r * particle_radius*particle_radius * (n1[2] * n2[0] - n1[0] * n2[2]);
                temp_vec[2] = walls[wall_id].rolling_modifier * wall_k_r * particle_radius*particle_radius * (n1[0] * n2[1] - n1[1] * n2[0]);

                torque_new[X_COORD(i)] -= temp_vec[0];
                torque_new[Y_COORD(i)] -= temp_vec[1];
                torque_new[Z_COORD(i)] -= temp_vec[2];


                //walls[WALL_ID(next_cl_entry->id)].total_torque[0] += temp_vec[0];
                //walls[WALL_ID(next_cl_entry->id)].total_torque[1] += temp_vec[1];
                //walls[WALL_ID(next_cl_entry->id)].total_torque[2] += temp_vec[2];
#endif

#ifdef WALL_SLIDING
                // calculate displacement
                temp_vec[0] = pos_new[X_COORD(i)] - wall_contact_pos[0] + particle_radius * n1[0];
                temp_vec[1] = pos_new[Y_COORD(i)] - wall_contact_pos[1] + particle_radius * n1[1];
                temp_vec[2] = pos_new[Z_COORD(i)] - wall_contact_pos[2] + particle_radius * n1[2];

                temp = dot_product(temp_vec, n2);
                displacement[0] = temp_vec[0] - temp * n2[0];
                displacement[1] = temp_vec[1] - temp * n2[1];
                displacement[2] = temp_vec[2] - temp * n2[2];


#ifdef INELASTIC_WALL_SLIDING
                displacement_norm = norm_squared(displacement);

                if(displacement_norm > crit_wall_sliding_displacement_squared)
                {
                    displacement_norm = sqrt(displacement_norm);

                    // correct the position where the contact has been established
                    next_cl_entry->contact->n2_initial[0] += (1.0 - crit_wall_sliding_displacement / displacement_norm) * displacement[0];
                    next_cl_entry->contact->n2_initial[1] += (1.0 - crit_wall_sliding_displacement / displacement_norm) * displacement[1];
                    next_cl_entry->contact->n2_initial[2] += (1.0 - crit_wall_sliding_displacement / displacement_norm) * displacement[2];

                    // correct current displacement displacement
                    displacement[0] *= (crit_wall_sliding_displacement / displacement_norm);
                    displacement[1] *= (crit_wall_sliding_displacement / displacement_norm);
                    displacement[2] *= (crit_wall_sliding_displacement / displacement_norm);

#ifdef TRACK_DISSIPATED_ENERGY
                    dissipated_wall_energy += walls[wall_id].sliding_modifier * wall_k_s * crit_sliding_displacement * (displacement_norm - crit_sliding_displacement) * ENERGY_UNIT;
#endif

#ifdef TRACK_DISSIPATED_ENERGY_PER_PARTICLE
                    dissipated_energy_of_particle[i] += walls[wall_id].sliding_modifier * wall_k_s * crit_sliding_displacement * (displacement_norm - crit_sliding_displacement) * ENERGY_UNIT;
#endif
                }
#endif // INELASTIC_WALL_SLIDING

                // force
                force_new[X_COORD(i)] -= (walls[wall_id].sliding_modifier * wall_k_s) * displacement[0];
                force_new[Y_COORD(i)] -= (walls[wall_id].sliding_modifier * wall_k_s) * displacement[1];
                force_new[Z_COORD(i)] -= (walls[wall_id].sliding_modifier * wall_k_s) * displacement[2];


#ifdef TRACK_FORCES
                if(i == 0)
                {
                    fprintf(track_forces, "%e", (walls[wall_id].sliding_modifier * wall_k_s) * norm(displacement));
                }
#endif // TRACK_FORCES


                walls[wall_id].total_force[0] += (walls[wall_id].sliding_modifier * wall_k_s) * displacement[0];
                walls[wall_id].total_force[1] += (walls[wall_id].sliding_modifier * wall_k_s) * displacement[1];
                walls[wall_id].total_force[2] += (walls[wall_id].sliding_modifier * wall_k_s) * displacement[2];



                if(sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
                {
                    if(particleInBox(i))
                    {
                        walls[wall_id].box_volume_force[0] += (walls[wall_id].sliding_modifier * wall_k_s) * displacement[0];
                        walls[wall_id].box_volume_force[1] += (walls[wall_id].sliding_modifier * wall_k_s) * displacement[1];
                        walls[wall_id].box_volume_force[2] += (walls[wall_id].sliding_modifier * wall_k_s) * displacement[2];
                    }
                }

                // torque
                torque_new[X_COORD(i)] -= walls[wall_id].sliding_modifier * wall_k_s * particle_radius * (n1[1] * displacement[2] - n1[2] * displacement[1]);
                torque_new[Y_COORD(i)] -= walls[wall_id].sliding_modifier * wall_k_s * particle_radius * (n1[2] * displacement[0] - n1[0] * displacement[2]);
                torque_new[Z_COORD(i)] -= walls[wall_id].sliding_modifier * wall_k_s * particle_radius * (n1[0] * displacement[1] - n1[1] * displacement[0]);

                //walls[-next_cl_entry->id-1].total_torque.x += k_s * (n2.y * displacement.z - n2.z * displacement.y);
                //walls[-next_cl_entry->id-1].total_torque.y += k_s * (n2.z * displacement.x - n2.x * displacement.z);
                //walls[-next_cl_entry->id-1].total_torque.z += k_s * (n2.x * displacement.y - n2.y * displacement.x);
#endif // WALL_SLIDING

#ifdef WALL_TWISTING
                torque_new[X_COORD(i)] -= walls[wall_id].twisting_modifier * wall_k_t * next_cl_entry->contact->twisting_displacement * n_c[0];
                torque_new[Y_COORD(i)] -= walls[wall_id].twisting_modifier * wall_k_t * next_cl_entry->contact->twisting_displacement * n_c[1];
                torque_new[Z_COORD(i)] -= walls[wall_id].twisting_modifier * wall_k_t * next_cl_entry->contact->twisting_displacement * n_c[2];
                //walls[-next_cl_entry->id-1].total_torque[0] += k_t * next_cl_entry->contact->twisting_displacement * n_c[0];
                //walls[-next_cl_entry->id-1].total_torque[1] += k_t * next_cl_entry->contact->twisting_displacement * n_c[1];
                //walls[-next_cl_entry->id-1].total_torque[2] += k_t * next_cl_entry->contact->twisting_displacement * n_c[2];

#endif
            }

#ifdef TRACK_FORCES
                if(i == 0 || i == 1)
                {
                    fprintf(track_forces, "\n");
                }

#endif // TRACK_FORCES

            next_cl_entry = next_cl_entry->next;
        }

    }

#ifdef TRACK_FORCES
                    fclose(track_forces);
#endif // TRACK_FORCES


}

void Simulation::switchOldNewValues()
{
    double *temp;


    temp = pos_new;
    pos_new = pos_old;
    pos_old = temp;


    temp = force_new;
    force_new =	force_old;
    force_old = temp;


    temp = torque_new;
    torque_new = torque_old;
    torque_old = temp;
}

double Simulation::getCrossSection()
{
    double cell_size = 2.1 * particle_radius;

    vec3 cms;
    SimLib::getCenterOfMass(&cms, *this, 0, number_of_particles-1);

    // the center of mass will be located in the middle of the area covered by the grid
    cms[0] -= (double)(GRID_CELLS/2) * cell_size;
    cms[2] -= (double)(GRID_CELLS/2) * cell_size;

    memset(cross_section_cells, 0, sizeof(int) * GRID_CELLS * GRID_CELLS);

    int x_id, z_id, cell_count = 0;

    for(int p = 0; p < number_of_particles; ++p)
    {
        x_id = (int)( (pos_old[X_COORD(p)] - cms[0])/cell_size );
        z_id = (int)( (pos_old[Z_COORD(p)] - cms[2])/cell_size );

        if(x_id < 0)
            x_id = 0;
        else if(x_id >= GRID_CELLS)
            x_id = GRID_CELLS-1;

        if(z_id < 0)
            z_id = 0;
        else if(z_id >= GRID_CELLS)
            z_id = GRID_CELLS-1;

        // only count cells that have not been marked as occupied yet
        if(cross_section_cells[x_id + GRID_CELLS * z_id] == 0)
        {
            ++cell_count;
            cross_section_cells[x_id + GRID_CELLS * z_id] = 1;
        }
    }

    return (double)cell_count * cell_size*cell_size;
}

double Simulation::getCrossSectionFillingFactor()
{
    double cross_section = 0.0;
    SimLib::getCrossSectionNoRotation(*this, 0.1*particle_radius, cross_section);

    if(cross_section > 0 && box)
        return (number_of_particles * 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius) / (cross_section * box->height);
    else
        return 0;
}

double Simulation::getBoxFillingFactor()
{
    if(box)
    {
        if(number_of_walls == 2 )
        {
            // determine number of particles within the volume
            /*
            int particles = 0;

            for(int p = 0; p < number_of_particles; ++p)
            {
                if(particleInBox(p))
                    ++particles;
            }

            if(box->base > 0 && box->height > 0)
                return (particles * 4.0 / 3.0 * M_PI * particle_radius*particle_radius*particle_radius) / (box->base * box->height);
            */
            return getCrossSectionFillingFactor();
        }
        else
        {



            // determine box
            vec3 lower_pos, upper_pos;

            // determine filling factor of a box inside our simulation box to reduce boundary effects
            lower_pos[0] = walls[LEFT_WALL_ID].pos[0] + 3.0 * particle_radius;
            lower_pos[1] = walls[box->bottom_wall_id].pos[1] + 3.0 * particle_radius;
            lower_pos[2] = walls[BACK_WALL_ID].pos[2] + 3.0 * particle_radius;
            upper_pos[0] = walls[RIGHT_WALL_ID].pos[0] - 3.0 * particle_radius;
            upper_pos[1] = walls[box->top_wall_id].pos[1] - 3.0 * particle_radius;
            upper_pos[2] = walls[FRONT_WALL_ID].pos[2] - 3.0 * particle_radius;


            return SimLib::getFillingFactorOfBox(*this, lower_pos, upper_pos);


        }
    }
    else
    {
        vec3 lower_pos, upper_pos;
        getEnclosingBox(&lower_pos, &upper_pos);

        // determine filling factor of a box inside our enclosing box to reduce boundary effects
        lower_pos[0] += 2.0 * particle_radius;
        lower_pos[1] += 2.0 * particle_radius;
        lower_pos[2] += 2.0 * particle_radius;
        upper_pos[0] -= 2.0 * particle_radius;
        upper_pos[1] -= 2.0 * particle_radius;
        upper_pos[2] -= 2.0 * particle_radius;

        return SimLib::getFillingFactorOfBox(*this, lower_pos, upper_pos);
    }

    return 0;
}



ErrorCode Simulation::getEnclosingBox(vec3 *min_pos, vec3 *max_pos)
{
    if(number_of_particles <= 0)
        return EC_NO_PARTICLES;

    (*min_pos)[0] = pos_old[X_COORD(0)];
    (*min_pos)[1] = pos_old[Y_COORD(0)];
    (*min_pos)[2] = pos_old[Z_COORD(0)];

    (*max_pos)[0] = pos_old[X_COORD(0)];
    (*max_pos)[1] = pos_old[Y_COORD(0)];
    (*max_pos)[2] = pos_old[Z_COORD(0)];

    for(int p = 1; p < number_of_particles; ++p)
    {
        if(pos_old[X_COORD(p)] < (*min_pos)[0])
            (*min_pos)[0] = pos_old[X_COORD(p)];

        if(pos_old[Y_COORD(p)] < (*min_pos)[1])
            (*min_pos)[1] = pos_old[Y_COORD(p)];

        if(pos_old[Z_COORD(p)] < (*min_pos)[2])
            (*min_pos)[2] = pos_old[Z_COORD(p)];

        if(pos_old[X_COORD(p)] > (*max_pos)[0])
            (*max_pos)[0] = pos_old[X_COORD(p)];

        if(pos_old[Y_COORD(p)] > (*max_pos)[1])
            (*max_pos)[1] = pos_old[Y_COORD(p)];

        if(pos_old[Z_COORD(p)] > (*max_pos)[2])
            (*max_pos)[2] = pos_old[Z_COORD(p)];
    }

    (*min_pos)[0] -= particle_radius;
    (*min_pos)[1] -= particle_radius;
    (*min_pos)[2] -= particle_radius;

    (*max_pos)[0] += particle_radius;
    (*max_pos)[1] += particle_radius;
    (*max_pos)[2] += particle_radius;

    return EC_OK;
}

void Simulation::initBox(double height, double x_size, double y_size, double lower_pos_x, double lower_pos_y, double lower_pos_z, double bottom_wall_id, double top_wall_id)
{
    // clean up old box
    if(box)
    {
        delete box;
        box = NULL;
    }

    box = new WallBox();

    box->height = height;
    box->x_size = x_size;
    box->y_size = y_size;
    box->base = x_size*y_size;

    box->lower_pos_x = lower_pos_x;
    box->lower_pos_y = lower_pos_y;
    box->lower_pos_z = lower_pos_z;

    box->bottom_wall_id = bottom_wall_id;
    box->top_wall_id = top_wall_id;
}

void Simulation::updateBox()
{
    if(box)
    {
        // check height
        if(box->top_wall_id != -1)
            box->height = walls[box->top_wall_id].pos[1] - walls[box->bottom_wall_id].pos[1] - 2.0 * particle_radius;

        // check size
        if(number_of_walls == 6)
        {
            box->lower_pos_x = walls[LEFT_WALL_ID].pos[0] + particle_radius;
            box->lower_pos_y = walls[box->bottom_wall_id].pos[1] + particle_radius;
            box->lower_pos_z = walls[BACK_WALL_ID].pos[2] + particle_radius;

            box->x_size = walls[RIGHT_WALL_ID].pos[0] - walls[LEFT_WALL_ID].pos[0] - 2.0 * particle_radius;
            box->y_size = walls[FRONT_WALL_ID].pos[2] - walls[BACK_WALL_ID].pos[2] - 2.0 * particle_radius;

            box->base = box->x_size*box->y_size;
        }
    }
}

bool Simulation::particleInBox(int id)
{
    return (pos_old[X_COORD(id)] > box->lower_pos_x && pos_old[X_COORD(id)] < box->lower_pos_x + box->x_size && pos_old[Z_COORD(id)] > box->lower_pos_z && pos_old[Z_COORD(id)] < box->lower_pos_z + box->y_size);
}

void Simulation::calculateRotation(double *mean_omega, double *sigma_omega)
{
    *mean_omega = 0;
    *sigma_omega = 0;

    for(int p = 0; p < number_of_particles; ++p)
    {
        double r_squared = pos_new[X_COORD(p)]*pos_new[X_COORD(p)] + pos_new[Z_COORD(p)]*pos_new[Z_COORD(p)];

        double omega = fabs(- vel[X_COORD(p)] * pos_new[Z_COORD(p)] + vel[Z_COORD(p)] * pos_new[X_COORD(p)]) / r_squared;

        *mean_omega += omega;
        *sigma_omega += omega*omega;
    }

    *mean_omega /= (double)number_of_particles;
    *sigma_omega /= (double)number_of_particles;

    *sigma_omega = sqrt(*sigma_omega  - *mean_omega * *mean_omega);
}

void Simulation::getErrorMessage(ErrorCode error_code, char *message)
{
    if(error_code == EC_OK)
        sprintf(message, "OK");
    else if(error_code == EC_INTERNAL_ERROR)
        sprintf(message, "Internal error");
    else if(error_code == EC_FILE_NOT_FOUND)
        sprintf(message, "Could not access file");
    else if(error_code == EC_INVALID_FILE_VERSION)
        sprintf(message, "Invalid file version");
    else if(error_code == EC_INVALID_MATERIAL_FILE_VERSION)
        sprintf(message, "Invalid material file version");
    else if(error_code == EC_NO_BOX)
        sprintf(message, "Box not initialized");
    else if(error_code == EC_INVALID_BOX_TYPE)
        sprintf(message, "Invalid box type");
    else if(error_code == EC_SIM_NOT_INITIALIZED)
        sprintf(message, "Simulation not initialized");
    else if(error_code == EC_NO_PARTICLES)
        sprintf(message, "Simulation has zero particles");
    else if(error_code == EC_DIFFERING_MATERIAL)
        sprintf(message, "Different material of data file");
    else if(error_code == EC_INVALID_PARAMETER)
        sprintf(message, "Invalid parameter");
    else if(error_code == EC_NO_CUDA_DEVICE)
        sprintf(message, "No suitable CUDA devide found");
    else if(error_code == EC_CUDA_DEVICE_UNAVAILABLE)
        sprintf(message, "CUDA device could not be accessed");
    else if(error_code == EC_INSUFFICIENT_GPU_MEMORY)
        sprintf(message, "Insufficient GPU memory");
    else if(error_code == EC_CUDA_NOT_INITIALIZED)
        sprintf(message, "Cuda not initiliazed");
    else if(error_code == EC_CUDA_INVALID_BOUNDARIES)
        sprintf(message, "Number of walls not supported on GPU");
    else if(error_code == EC_INVALID_BOUNDARIES)
        sprintf(message, "Invalid number of walls");
    else if(error_code == EC_TIMESTEP_TOO_LARGE)
        sprintf(message, "Timestep may be too large");
    else if(error_code == EC_NO_HIT_AND_STICK_POS_FOUND)
        sprintf(message, "No hit and stick configuration found");
    else if(error_code == EC_FIXED_PARTICLES_NOT_INITIALIZED)
        sprintf(message, "Fixed particles are not initialized");
    else if(error_code == EC_CONTACT_LIST_CORRUPTED)
        sprintf(message, "Contact list corrupted");
    else
        sprintf(message, "Unknown error code");
}

#ifdef TRACK_PARTICLE_ORIENTATION

void Simulation::initParticleOrientation()
{
    double phi;
    double cos_theta;
    double sin_theta;
    double rot_angle;
    vec3 axis;

    for(int p = 0; p < number_of_particles; ++p)
    {
        // determine random axis/rotation angle
        phi = get_random_zero_twoPi();
        cos_theta = get_random_cos_theta();
        sin_theta = sqrt(1.0 - cos_theta*cos_theta);
        rot_angle = get_random_zero_twoPi();

        axis[0] = sin_theta * cos(phi);
        axis[1] = sin_theta * sin(phi);
        axis[2] = cos_theta;

        orientation[4*p] = cos(0.5 * rot_angle);
        orientation[4*p+1] = axis[0] * sin(0.5 * rot_angle);
        orientation[4*p+2] = axis[1] * sin(0.5 * rot_angle);
        orientation[4*p+3] = axis[2] * sin(0.5 * rot_angle);
    }
}

void Simulation::updateParticleOrientation(double timestep)
{
    double e_dot[4];
    double e_ddot[4];
    vec3 omega_dot;
    double omega_squared, norm_inv;

    for(int p = 0; p < number_of_particles; ++p)
    {
        omega_squared = vel_angular[X_COORD(p)]*vel_angular[X_COORD(p)] + vel_angular[Y_COORD(p)]*vel_angular[Y_COORD(p)] + vel_angular[Z_COORD(p)]*vel_angular[Z_COORD(p)];
        omega_dot[0] = moment_of_inertia_inv * torque_old[X_COORD(p)];
        omega_dot[1] = moment_of_inertia_inv * torque_old[Y_COORD(p)];
        omega_dot[2] = moment_of_inertia_inv * torque_old[Z_COORD(p)];

        e_dot[0] = - 0.5 * (orientation[4*p+1] * vel_angular[X_COORD(p)] + orientation[4*p+2] * vel_angular[Y_COORD(p)] + orientation[4*p+3] * vel_angular[Z_COORD(p)]);
        e_dot[1] = 0.5 * (orientation[4*p] * vel_angular[X_COORD(p)] - orientation[4*p+2] * vel_angular[Z_COORD(p)] + orientation[4*p+3] * vel_angular[Y_COORD(p)]);
        e_dot[2] = 0.5 * (orientation[4*p] * vel_angular[Y_COORD(p)] - orientation[4*p+3] * vel_angular[X_COORD(p)] + orientation[4*p+1] * vel_angular[Z_COORD(p)]);
        e_dot[3] = 0.5 * (orientation[4*p] * vel_angular[Z_COORD(p)] - orientation[4*p+1] * vel_angular[Y_COORD(p)] + orientation[4*p+2] * vel_angular[X_COORD(p)]);

        double temp = - 0.25 * (orientation[4*p+1] * vel_angular[X_COORD(p)] + orientation[4*p+2] * vel_angular[Y_COORD(p)] + orientation[4*p+3] * vel_angular[Z_COORD(p)]);

        e_ddot[0] = - 0.25 * (orientation[4*p] * omega_squared + 2.0 * (orientation[4*p+1] * omega_dot[0] + orientation[4*p+2] * omega_dot[1] + orientation[4*p+3] * omega_dot[2]));
        e_ddot[1] = temp * vel_angular[X_COORD(p)] + 0.5 * (orientation[4*p] * omega_dot[0] - orientation[4*p+2] * omega_dot[2] + orientation[4*p+3] * omega_dot[1]);
        e_ddot[2] = temp * vel_angular[Y_COORD(p)] + 0.5 * (orientation[4*p] * omega_dot[1] - orientation[4*p+3] * omega_dot[0] + orientation[4*p+1] * omega_dot[2]);
        e_ddot[3] = temp * vel_angular[Z_COORD(p)] + 0.5 * (orientation[4*p] * omega_dot[2] - orientation[4*p+1] * omega_dot[1] + orientation[4*p+2] * omega_dot[0]);

        orientation[4*p] += timestep * e_dot[0] + 0.5 * timestep * timestep * e_ddot[0];
        orientation[4*p+1] += timestep * e_dot[1] + 0.5 * timestep * timestep * e_ddot[1];
        orientation[4*p+2] += timestep * e_dot[2] + 0.5 * timestep * timestep * e_ddot[2];
        orientation[4*p+3] += timestep * e_dot[3] + 0.5 * timestep * timestep * e_ddot[3];

        // make sure that e0^2 + .. + e3^2 = 1
        norm_inv = 1.0 / sqrt(orientation[4*p]*orientation[4*p] + orientation[4*p+1]*orientation[4*p+1] + orientation[4*p+2]*orientation[4*p+2] + orientation[4*p+3]*orientation[4*p+3]);

        orientation[4*p] *= norm_inv;
        orientation[4*p+1] *= norm_inv;
        orientation[4*p+2] *= norm_inv;
        orientation[4*p+3] *= norm_inv;
    }
}
#endif

#ifdef ENABLE_FIXED_PARTICLES
void Simulation::initFixedPaticles(bool fixed)
{
    if(fixed_particles)
        delete [] fixed_particles;

    fixed_particles = new bool[number_of_particles];

    for(int p = 0; p < number_of_particles; ++p)
        fixed_particles[p] = fixed;
}

ErrorCode Simulation::fixateParticle(int particle_id, bool fixed)
{
    if(fixed_particles)
    {
        fixed_particles[particle_id] = fixed;
        return EC_OK;
    }
    else
        return EC_FIXED_PARTICLES_NOT_INITIALIZED;
}
#endif

int Simulation::getNumberOfContacts()
{
    int number_of_contacts = 0;
    ContactListEntry *cl_entry;

    for(int p = 0; p < number_of_particles; ++p)
    {
        cl_entry = contact_list[p];

        while(cl_entry)
        {
            if(cl_entry->id >= 0)
                ++number_of_contacts;

            cl_entry = cl_entry->next;
        }
    }

    return number_of_contacts;
}
int Simulation::getNumberOfWallContacts()
{
    int number_of_wall_contacts = 0;
    ContactListEntry *cl_entry;

    for(int p = 0; p < number_of_particles; ++p)
    {
        cl_entry = contact_list[p];

        while(cl_entry)
        {
            if(cl_entry->id < 0)
                ++number_of_wall_contacts;

            cl_entry = cl_entry->next;
        }
    }

    return number_of_wall_contacts;
}

void Simulation::debugDump(const char *filename)
{
    FILE *file = fopen(filename, "w+");

    if(!file)
        return;

    fprintf(file, "# number of particles: %i   number of contacts: %i\n", number_of_particles, number_of_contacts);
    fprintf(file, "# Enabled interaction: Normal interaction");

#ifdef ROLLING
    fprintf(file, "  Rolling");
#endif
#ifdef SLIDING
    fprintf(file, "  Sliding");
#endif
#ifdef TWISTING
    fprintf(file, "  Twisting");
#endif
#ifdef  INELASTIC_ROLLING
    fprintf(file, "  Inelastic Rolling");
#endif
#ifdef  INELASTIC_SLIDING
    fprintf(file, "  Inelastic Sliding");
#endif
#ifdef  INELASTIC_TWISTING
    fprintf(file, "  Inelastic Twisting");
#endif

    fprintf(file, "\n# Rolling / sliding / twisting constant: %g / %g / %g\n", k_r, k_s, k_t);

    fprintf(file, "# Timestep: %g   Elapsed time: %g\n", timestep, current_time);

    fprintf(file, "\n# New positions / velocities / angular velocities\n");
    for(int p = 0; p < number_of_particles; ++p)
    {
        fprintf(file, "%.20g %.20g %.20g  /  %.20g %.20g %.20g  /  %.20g %.20g %.20g\n", pos_new[X_COORD(p)], pos_new[Y_COORD(p)], pos_new[Z_COORD(p)], vel[X_COORD(p)], vel[Y_COORD(p)], vel[Z_COORD(p)], vel_angular[X_COORD(p)], vel_angular[Y_COORD(p)], vel_angular[Z_COORD(p)]);
        printf("%.20g %.20g %.20g\n", pos_new[X_COORD(p)], pos_new[Y_COORD(p)], pos_new[Z_COORD(p)]);
    }
    fprintf(file, "\n# New forces / torques\n");
    for(int p = 0; p < number_of_particles; ++p)
        fprintf(file, "%.20g %.20g %.20g  /  %.20g %.20g %.20g\n", force_new[X_COORD(p)], force_new[Y_COORD(p)], force_new[Z_COORD(p)], torque_new[X_COORD(p)], torque_new[Y_COORD(p)], torque_new[Z_COORD(p)]);

    fprintf(file, "\n# Contacts: id1-id2: e1_0 e1_1 e1_2 e1_3 / e2_0 e2_1 e2_2 e2_3\n");

    ContactListEntry *cl_entry;

    for(int p = 0; p < number_of_particles; ++p)
    {
        cl_entry = contact_list[p];

        while(cl_entry)
        {
            Contact *contact = cl_entry->contact;
            fprintf(file, "%i-%i: %.20g %.20g %.20g %.20g  / %.20g %.20g %.20g %.20g\n", contact->id1, contact->id2, contact->rot1.e0, contact->rot1.e1, contact->rot1.e2, contact->rot1.e3, contact->rot2.e0, contact->rot2.e1, contact->rot2.e2, contact->rot2.e3);

            vec3 n1, n2;
            contact->getCurrentN1(&n1);
            contact->getCurrentN2(&n2);

            fprintf(file, "   (%.20g, %.20g, %.20g)  (%.20g, %.20g, %.20g)\n", n1[0], n1[1], n1[2], n2[0], n2[1], n2[2]);

            cl_entry = cl_entry->next;
        }
    }

    fclose(file);
}
