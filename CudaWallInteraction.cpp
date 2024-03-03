#include "CudaWallInteraction.h"
#include "CudaDefines.cuh"
#include "CudaWallKernelsWrapper.cuh"

#include "Wall.h"
#include <vector>

#include <cstdio>

#if defined _WIN32 || _WIN64
#define _USE_MATH_DEFINES
#include <math.h>
#else
#include <cmath>
#endif



// global variables that do not change during the sim (material parameters, etc.)
extern double particle_radius;		// in cm
extern double particle_radius_inv;	// in cm^-1


// material parameters
extern double surface_energy;		// surface energy / area (in mJ / m^2)

extern double moment_of_inertia;	// in g cm^2
extern double moment_of_inertia_inv;
extern double osc_damping_factor;	// additional damping for normal oscillations
extern double T_vis;
extern double viscoelastic_damping_constant;  // viscoelastic damping constant as proposed by Sebastian Krijt

extern double crit_rolling_displacement;			// in cm
extern double crit_rolling_displacement_squared;	// in cm^2


// for wall interaction
extern double wall_reduced_radius;
extern double wall_equilibrium_contact_radius;
extern double wall_delta_c;
extern double wall_contact_breaking_dist;
extern double wall_F_c;
extern double wall_k_s;
extern double wall_k_r;
extern double wall_k_t;
extern double crit_wall_sliding_displacement;
extern double crit_wall_sliding_displacement_squared;


particleWallInteraction::particleWallInteraction()
{

    m_wall_interaction_initialized = false;
    m_wall_interaction_allocated = false;
    m_wall_gpu_constants_initialized = false;
    m_wall_box_loaded = false;


}


particleWallInteraction::~particleWallInteraction()
{

    if(m_wall_interaction_allocated)
    {
        deleteParticleWallInteraction();
    }

}



void particleWallInteraction::updateWallContactNormals(
        const double* RESTRICT positions
        )

{

    if(m_wall_interaction_initialized)
    {

        updateWallContactNormalsWrapper(
                    m_gpu_wall_particle_ids,
                    m_gpu_dist_to_wall,
                    m_gpu_wall_contact_vector_new,
                    positions,
                    m_gpu_wall_initial_contact,
                    m_gpu_wall_positions,
                    m_number_of_particles
                    );

    }
    else
    {
        if(!m_wall_gpu_constants_initialized)
        {
            printf("Error trying to start updateWallContacts without wall material constants initialized\n");
        }
        if(!m_wall_interaction_allocated)
        {
            printf("Error trying to start updateWallContacts without having arrays allocated\n");
        }
    }

}


void particleWallInteraction::updateWallSticking(
        const double* RESTRICT positions,
        const double* RESTRICT gpu_angular_velocities,
        const double* RESTRICT gpu_torques,
        const double timestep
        )
{

    if(m_wall_interaction_initialized)
    {

        m_cpu_wall_positions[GPU_WALL_TOP].y += timestep * m_wallConstants.top_wall_speed;

        copyToGPU((void*)m_gpu_wall_positions, (void*)m_cpu_wall_positions, 6 * sizeof(double3));


        updateWallStickingWrapper(
            m_gpu_wall_particle_ids,
            positions,
            m_gpu_dist_to_wall,
            m_gpu_wall_contact_vector_new,
            m_gpu_wall_contact_vector_old,
            m_gpu_wall_contact_n_initial,
            m_gpu_wall_contact_rot,
            m_gpu_wall_contact_n,
            m_gpu_wall_initial_contact,
            m_gpu_wall_twisting_displacement,
            m_gpu_wall_positions,
            m_number_of_particles,

            gpu_angular_velocities,
            gpu_torques,
            timestep
                );
    }
    else
    {
        if(!m_wall_gpu_constants_initialized)
        {
            printf("Error trying to start updateWallSticking without wall material constants initialized\n");
        }
        if(!m_wall_interaction_allocated)
        {
            printf("Error trying to start updateWallSticking without having arrays allocated\n");
        }
    }

}




void particleWallInteraction::updateWallStickingNoSw(
        const double* RESTRICT positions,
        const double* RESTRICT gpu_angular_velocities,
        const double* RESTRICT gpu_torques,
        const double timestep
        )
{

    if(m_wall_interaction_initialized)
    {

        m_cpu_wall_positions[GPU_WALL_TOP].y += timestep * m_wallConstants.top_wall_speed;

        copyToGPU((void*)m_gpu_wall_positions, (void*)m_cpu_wall_positions, 2 * sizeof(double3));

        updateWallStickingWrapperNoSw(
            m_gpu_wall_particle_ids,
            positions,
            m_gpu_dist_to_wall,
            m_gpu_wall_contact_vector_new,
            m_gpu_wall_contact_vector_old,
            m_gpu_wall_contact_n_initial,
            m_gpu_wall_contact_rot,
            m_gpu_wall_contact_n,
            m_gpu_wall_initial_contact,
            m_gpu_wall_twisting_displacement,
            m_gpu_wall_positions,
            m_number_of_particles,

            gpu_angular_velocities,
            gpu_torques,
            timestep
                );
    }
    else
    {
        if(!m_wall_gpu_constants_initialized)
        {
            printf("Error trying to start updateWallStickingNoSw without wall material constants initialized\n");
        }
        if(!m_wall_interaction_allocated)
        {
            printf("Error trying to start updateWallStickingNoSw without having arrays allocated\n");
        }
    }

}





void particleWallInteraction::updateWallContacts(
        double * RESTRICT gpu_forces_new,
        double * RESTRICT gpu_torques_new,
        const double* RESTRICT gpu_velocities,
        const double* RESTRICT gpu_angular_velocities,
        const double* RESTRICT gpu_torques,
        const double timestep
    )

{

    if(m_wall_interaction_initialized)
    {

        updateWallContactsWrapper(
            m_gpu_wall_contact_rot,
            m_gpu_wall_contact_n_initial,
            m_gpu_wall_contact_n,
            gpu_forces_new,
            gpu_torques_new,
            m_gpu_top_wall_force, // force the particle exert on the top wall
            m_gpu_bot_wall_force, // force the particle exert on the bot wall
            m_gpu_dist_to_wall,
            m_gpu_wall_contact_vector_new,
            m_gpu_wall_contact_vector_old,
            m_gpu_wall_initial_contact, // = n2_initial in alex code
            gpu_velocities,
            gpu_angular_velocities,
            gpu_torques,
            m_gpu_wall_particle_ids,
            m_gpu_wall_twisting_displacement,
            timestep,
            m_number_of_particles
            );


        double3 *temp2;
        temp2 = m_gpu_wall_contact_vector_new;
        m_gpu_wall_contact_vector_new = m_gpu_wall_contact_vector_old;
        m_gpu_wall_contact_vector_old = temp2;

    }
    else
    {
        if(!m_wall_gpu_constants_initialized)
        {
            printf("Error trying to start updateWallContacts without wall material constants initialized\n");
        }
        if(!m_wall_interaction_allocated)
        {
            printf("Error trying to start updateWallContacts without having arrays allocated\n");
        }
    }
}


void particleWallInteraction::initParticleWallInteraction(const int number_of_particles)
{


    if(m_wall_interaction_allocated)
    {
        deleteParticleWallInteraction();
    }

    m_number_of_particles = number_of_particles;

    // assuming the box is bigger than a single particle, the particle can be in contact with 3 walls at most if it's directly in a corner
    int max_number_of_wall_contacts = 3*m_number_of_particles;

    // allocate gpu arrays

    allocateArray((void**) &m_gpu_wall_particle_ids, m_number_of_particles * sizeof(int));			// stores the ids of the other particles a certain particle is in contact with
    allocateArray((void**) &m_gpu_top_wall_force, m_number_of_particles * sizeof(double));
    allocateArray((void**) &m_gpu_bot_wall_force, m_number_of_particles * sizeof(double));



    allocateArray((void**) &m_gpu_dist_to_wall, max_number_of_wall_contacts * sizeof(double));

    allocateArray((void**) &m_gpu_wall_contact_vector_new, max_number_of_wall_contacts * sizeof(double3));		// stores contact normal
    allocateArray((void**) &m_gpu_wall_contact_vector_old, max_number_of_wall_contacts * sizeof(double3));
    allocateArray((void**) &m_gpu_wall_contact_n_initial, max_number_of_wall_contacts * sizeof(double3));		// initial contact pointer
    allocateArray((void**) &m_gpu_wall_contact_rot, max_number_of_wall_contacts * sizeof(double4));				// rotation state
    allocateArray((void**) &m_gpu_wall_contact_n, max_number_of_wall_contacts * sizeof(double3));				// current n1
    allocateArray((void**) &m_gpu_wall_initial_contact, max_number_of_wall_contacts * sizeof(double3));

    allocateArray((void**) &m_gpu_wall_twisting_displacement, max_number_of_wall_contacts * sizeof(double));

    allocateArray((void**) &m_gpu_wall_positions, 6 * sizeof(double3));
    m_cpu_wall_positions = new double3[6];


    memsetArray(m_gpu_wall_particle_ids, 0, m_number_of_particles*sizeof(int));
    memsetArray(m_gpu_top_wall_force, 0, m_number_of_particles * sizeof(double));
    memsetArray(m_gpu_bot_wall_force, 0, m_number_of_particles * sizeof(double));


    m_wallConstants;

    m_wall_interaction_allocated = true;
    m_wall_interaction_initialized = m_wall_interaction_allocated && m_wall_gpu_constants_initialized & m_wall_box_loaded;

}



void particleWallInteraction::deleteParticleWallInteraction()
{

    if(m_wall_interaction_allocated)
    {



        freeArray(m_gpu_wall_particle_ids);			// stores the ids of the other particles a certain particle is in contact with

        freeArray(m_gpu_top_wall_force);
        freeArray(m_gpu_bot_wall_force);

        freeArray(m_gpu_dist_to_wall);
        freeArray(m_gpu_wall_contact_vector_new);
        freeArray(m_gpu_wall_contact_vector_old);
        freeArray(m_gpu_wall_contact_n_initial);
        freeArray(m_gpu_wall_contact_rot);
        freeArray(m_gpu_wall_contact_n);
        freeArray(m_gpu_wall_initial_contact);

        freeArray(m_gpu_wall_twisting_displacement);

        freeArray(m_gpu_wall_positions);
        delete[]  m_cpu_wall_positions;


        m_wall_interaction_allocated = false;
        m_wall_interaction_initialized = false;
        m_wall_box_loaded = false;
    }
    else
    {
        printf("Warning: trying to delete particle wall interaction without having it initialized first\n");
    }

}




static double getWallContactRadiusBuild(const double compression_length, double reduced_particle_radius, double c1_contact_radius, double c2_contact_radius, double particle_radius)
{
    // contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
    // use equilibrium contact radius as starting value
    const double k = reduced_particle_radius * compression_length / 3.0;
    double x_pow3;
    double x_new;
    double x_old = c1_contact_radius;

    // use Newton-Raphson method to find root
    for (int i = 0; i < 200; ++i)
    {
        x_pow3 = x_old * x_old * x_old;
        x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - c2_contact_radius);

        if (fabs(x_new - x_old) / particle_radius < 1.e-14)
            return x_new * x_new;

        x_old = x_new;
    }

    return x_new * x_new;
}


void particleWallInteraction::initMaterialConstants(
        const NormalInteraction &normal_interaction
        )
{


    double wall_delta_0 = wall_equilibrium_contact_radius * wall_equilibrium_contact_radius / (3.0 * wall_reduced_radius);

    double get_contact_radius_c1 = getWallContactRadiusBuild(
                -wall_delta_c,
                wall_reduced_radius,
                normal_interaction.c1_contact_radius_wall,
                normal_interaction.c2_contact_radius_wall,
                particle_radius);

    double get_contact_radius_c2 = (wall_equilibrium_contact_radius - get_contact_radius_c1) / sqrt(wall_delta_c + wall_delta_0);

    m_wallConstants.c1_contact_radius = normal_interaction.c1_contact_radius_wall;
    m_wallConstants.c2_contact_radius = normal_interaction.c2_contact_radius_wall;

    m_wallConstants.c1_normal_force = normal_interaction.k1_wall;
    m_wallConstants.c2_normal_force = -normal_interaction.k2_wall;

    m_wallConstants.eq_contact_radius = wall_equilibrium_contact_radius;
    m_wallConstants.delta_c = wall_delta_c;

    m_wallConstants.contact_breaking_dist = wall_contact_breaking_dist;

    m_wallConstants.particle_radius = particle_radius;
    m_wallConstants.reduced_particle_radius = wall_reduced_radius;

    m_wallConstants.particle_radius_inv = particle_radius_inv;
    m_wallConstants.moment_of_inertia_inv = moment_of_inertia_inv;


    m_wallConstants.k_r = wall_k_r * particle_radius * particle_radius;
    m_wallConstants.k_s = wall_k_s;
    m_wallConstants.k_t = wall_k_t;

    m_wallConstants.crit_rolling_displacement = crit_rolling_displacement;
    m_wallConstants.crit_rolling_displacement_squared = crit_rolling_displacement*crit_rolling_displacement;

    m_wallConstants.crit_twisting_displacement = 1.0 / (16.0 * M_PI);

    m_wallConstants.crit_sliding_displacement = crit_wall_sliding_displacement;
    m_wallConstants.crit_sliding_displacement_squared = crit_wall_sliding_displacement_squared;

    m_wallConstants.get_contact_radius_c1 = get_contact_radius_c1;
    m_wallConstants.get_contact_radius_c2 = get_contact_radius_c2;

#ifdef GPU_WALL_USE_VISCOELASTIC_DAMPING
    m_wallConstants.osc_damping_factor = viscoelastic_damping_constant;
#endif // GPU_WALL_USE_VISCOELASTIC_DAMPING

#ifdef GPU_WALL_USE_CONSTANT_DAMPING
    m_wallConstants.osc_damping_factor = osc_damping_factor;
#endif // GPU_WALL_USE_CONSTANT_DAMPING







    m_wallConstants.top_wall_speed = 0.0;

    m_wallConstants.wall_normals[GPU_WALL_BOTTOM] = make_double3(0.0);
    m_wallConstants.wall_normals[GPU_WALL_LEFT  ] = make_double3(0.0);
    m_wallConstants.wall_normals[GPU_WALL_RIGHT ] = make_double3(0.0);
    m_wallConstants.wall_normals[GPU_WALL_FRONT ] = make_double3(0.0);
    m_wallConstants.wall_normals[GPU_WALL_BACK  ] = make_double3(0.0);
    m_wallConstants.wall_normals[GPU_WALL_TOP   ] = make_double3(0.0);


    m_wallConstants.wall_compression_modifier[GPU_WALL_BOTTOM] = 0.0;
    m_wallConstants.wall_compression_modifier[GPU_WALL_LEFT  ] = 0.0;
    m_wallConstants.wall_compression_modifier[GPU_WALL_RIGHT ] = 0.0;
    m_wallConstants.wall_compression_modifier[GPU_WALL_FRONT ] = 0.0;
    m_wallConstants.wall_compression_modifier[GPU_WALL_BACK  ] = 0.0;
    m_wallConstants.wall_compression_modifier[GPU_WALL_TOP   ] = 0.0;


    m_wallConstants.wall_sliding_modifier[GPU_WALL_BOTTOM]  = 0.0;
    m_wallConstants.wall_sliding_modifier[GPU_WALL_LEFT  ]  = 0.0;
    m_wallConstants.wall_sliding_modifier[GPU_WALL_RIGHT ]  = 0.0;
    m_wallConstants.wall_sliding_modifier[GPU_WALL_FRONT ]  = 0.0;
    m_wallConstants.wall_sliding_modifier[GPU_WALL_BACK  ]  = 0.0;
    m_wallConstants.wall_sliding_modifier[GPU_WALL_TOP   ]  = 0.0;


    m_wallConstants.rolling_modifier[GPU_WALL_BOTTOM] = 0.0;
    m_wallConstants.rolling_modifier[GPU_WALL_LEFT  ] = 0.0;
    m_wallConstants.rolling_modifier[GPU_WALL_RIGHT ] = 0.0;
    m_wallConstants.rolling_modifier[GPU_WALL_FRONT ] = 0.0;
    m_wallConstants.rolling_modifier[GPU_WALL_BACK  ] = 0.0;
    m_wallConstants.rolling_modifier[GPU_WALL_TOP   ] = 0.0;


    m_wallConstants.twisting_modifier[GPU_WALL_BOTTOM] = 0.0;
    m_wallConstants.twisting_modifier[GPU_WALL_LEFT  ] = 0.0;
    m_wallConstants.twisting_modifier[GPU_WALL_RIGHT ] = 0.0;
    m_wallConstants.twisting_modifier[GPU_WALL_FRONT ] = 0.0;
    m_wallConstants.twisting_modifier[GPU_WALL_BACK  ] = 0.0;
    m_wallConstants.twisting_modifier[GPU_WALL_TOP   ] = 0.0;


    // load wall positions from simulation walls vector to array in particleWallInteraction
    m_cpu_wall_positions[GPU_WALL_BOTTOM] = make_double3(0.0);
    m_cpu_wall_positions[GPU_WALL_LEFT  ] = make_double3(0.0);
    m_cpu_wall_positions[GPU_WALL_RIGHT ] = make_double3(0.0);
    m_cpu_wall_positions[GPU_WALL_FRONT ] = make_double3(0.0);
    m_cpu_wall_positions[GPU_WALL_BACK  ] = make_double3(0.0);
    m_cpu_wall_positions[GPU_WALL_TOP   ] = make_double3(0.0);







    init_GPUWallConstansWrapper(m_wallConstants);

    m_wall_gpu_constants_initialized = true;
    m_wall_box_loaded = false;
    m_wall_interaction_initialized = m_wall_interaction_allocated && m_wall_gpu_constants_initialized   && m_wall_box_loaded;

}




double particleWallInteraction::getGPUTopWallForce(CubPlan &cubPlan)
{
    double F = 0.0;

    F = cubPlan.reduce(m_gpu_top_wall_force, m_number_of_particles);

    memsetArray(m_gpu_top_wall_force, 0, m_number_of_particles * sizeof(double));


    return F;
}

double particleWallInteraction::getGPUBotWallForce(CubPlan &cubPlan)
{
    double F = 0.0;

    F = cubPlan.reduce(m_gpu_bot_wall_force, m_number_of_particles);

    memsetArray(m_gpu_bot_wall_force, 0, m_number_of_particles * sizeof(double));

    return F;
}



int wall_id_GPU_to_CPU(int gpu_wall_id, int num_walls)
{

    int cpu_wall_id = -1;

    if(num_walls == 6)
    {
        switch(gpu_wall_id)
        {
            case GPU_WALL_BOTTOM:   cpu_wall_id = 0; break;
            case GPU_WALL_LEFT:     cpu_wall_id = 1; break;
            case GPU_WALL_RIGHT:    cpu_wall_id = 2; break;
            case GPU_WALL_FRONT:    cpu_wall_id = 3; break;
            case GPU_WALL_BACK:     cpu_wall_id = 4; break;
            case GPU_WALL_TOP:      cpu_wall_id = 5; break;
        }

    }
    else if(num_walls == 2)
    {
        switch(gpu_wall_id)
        {
            case GPU_WALL_BOTTOM:   cpu_wall_id = 0; break;
            case GPU_WALL_TOP:      cpu_wall_id = 1; break;
        }
    }



    return cpu_wall_id;

}


int wall_id_CPU_to_GPU(int cpu_wall_id, int num_walls)
{

    int gpu_wall_id = -1;

    if(num_walls == 6)
    {
        switch(cpu_wall_id)
        {
            case 0: gpu_wall_id = GPU_WALL_BOTTOM;  break;
            case 1: gpu_wall_id = GPU_WALL_LEFT;    break;
            case 2: gpu_wall_id = GPU_WALL_RIGHT;   break;
            case 3: gpu_wall_id = GPU_WALL_FRONT;   break;
            case 4: gpu_wall_id = GPU_WALL_BACK;    break;
            case 5: gpu_wall_id = GPU_WALL_TOP;     break;
        }

    }
    else if(num_walls == 2)
    {
        switch(cpu_wall_id)
        {
            case 0: gpu_wall_id = GPU_WALL_BOTTOM;  break;
            case 1: gpu_wall_id = GPU_WALL_TOP;     break;
        }
    }

    return gpu_wall_id;

}



int cpu_find_free_wall_contact_id(int bitmap)
{
    int wall_id_id = 0;

    // choose free spot in the bitmap, empty byte means free spot
    // 65280 = 255<<8
    if ((bitmap & 65280) == 0) // try if second byte is empty
    {
        wall_id_id = 1;
    }
    // 16711680 = 255<<16
    else if ((bitmap & 16711680) == 0) // else try if third byte is empty
    {
        wall_id_id = 2;
    }
    // 4278190080 = 255<<24
    else if ((bitmap & 4278190080) == 0) // else try if fourth byte is empty
    {
        wall_id_id = 3;
    }

    return wall_id_id;
}






void particleWallInteraction::load_box_to_gpu(const Simulation &sim)
{



    // Load box to GPU
    if(sim.walls.size() != 6 && sim.walls.size() != 2)
    {
        printf("Error: not six or two walls!\n");
        m_wall_interaction_initialized = false;
        return;
    }


    m_wallConstants.top_wall_speed = sim.walls[sim.box->top_wall_id].velocity[1];


    if(sim.walls.size() == 6)
    {

        m_wallConstants.wall_normals[GPU_WALL_BOTTOM] = make_double3(sim.walls[0].normal[0], sim.walls[0].normal[1], sim.walls[0].normal[2]);
        m_wallConstants.wall_normals[GPU_WALL_LEFT  ] = make_double3(sim.walls[1].normal[0], sim.walls[1].normal[1], sim.walls[1].normal[2]);
        m_wallConstants.wall_normals[GPU_WALL_RIGHT ] = make_double3(sim.walls[2].normal[0], sim.walls[2].normal[1], sim.walls[2].normal[2]);
        m_wallConstants.wall_normals[GPU_WALL_FRONT ] = make_double3(sim.walls[3].normal[0], sim.walls[3].normal[1], sim.walls[3].normal[2]);
        m_wallConstants.wall_normals[GPU_WALL_BACK  ] = make_double3(sim.walls[4].normal[0], sim.walls[4].normal[1], sim.walls[4].normal[2]);
        m_wallConstants.wall_normals[GPU_WALL_TOP   ] = make_double3(sim.walls[5].normal[0], sim.walls[5].normal[1], sim.walls[5].normal[2]);


        m_wallConstants.wall_compression_modifier[GPU_WALL_BOTTOM] = sim.walls[0].compression_modifier;
        m_wallConstants.wall_compression_modifier[GPU_WALL_LEFT  ] = sim.walls[1].compression_modifier;
        m_wallConstants.wall_compression_modifier[GPU_WALL_RIGHT ] = sim.walls[2].compression_modifier;
        m_wallConstants.wall_compression_modifier[GPU_WALL_FRONT ] = sim.walls[3].compression_modifier;
        m_wallConstants.wall_compression_modifier[GPU_WALL_BACK  ] = sim.walls[4].compression_modifier;
        m_wallConstants.wall_compression_modifier[GPU_WALL_TOP   ] = sim.walls[5].compression_modifier;


        m_wallConstants.wall_sliding_modifier[GPU_WALL_BOTTOM] = sim.walls[0].sliding_modifier;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_LEFT  ] = sim.walls[1].sliding_modifier;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_RIGHT ] = sim.walls[2].sliding_modifier;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_FRONT ] = sim.walls[3].sliding_modifier;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_BACK  ] = sim.walls[4].sliding_modifier;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_TOP   ] = sim.walls[5].sliding_modifier;


        m_wallConstants.rolling_modifier[GPU_WALL_BOTTOM] = sim.walls[0].rolling_modifier;
        m_wallConstants.rolling_modifier[GPU_WALL_LEFT  ] = sim.walls[1].rolling_modifier;
        m_wallConstants.rolling_modifier[GPU_WALL_RIGHT ] = sim.walls[2].rolling_modifier;
        m_wallConstants.rolling_modifier[GPU_WALL_FRONT ] = sim.walls[3].rolling_modifier;
        m_wallConstants.rolling_modifier[GPU_WALL_BACK  ] = sim.walls[4].rolling_modifier;
        m_wallConstants.rolling_modifier[GPU_WALL_TOP   ] = sim.walls[5].rolling_modifier;


        m_wallConstants.twisting_modifier[GPU_WALL_BOTTOM] = sim.walls[0].twisting_modifier;
        m_wallConstants.twisting_modifier[GPU_WALL_LEFT  ] = sim.walls[1].twisting_modifier;
        m_wallConstants.twisting_modifier[GPU_WALL_RIGHT ] = sim.walls[2].twisting_modifier;
        m_wallConstants.twisting_modifier[GPU_WALL_FRONT ] = sim.walls[3].twisting_modifier;
        m_wallConstants.twisting_modifier[GPU_WALL_BACK  ] = sim.walls[4].twisting_modifier;
        m_wallConstants.twisting_modifier[GPU_WALL_TOP   ] = sim.walls[5].twisting_modifier;


        // load wall positions from simulation walls vector to array in particleWallInteraction
        m_cpu_wall_positions[GPU_WALL_BOTTOM] = make_double3(sim.walls[0].pos[0], sim.walls[0].pos[1], sim.walls[0].pos[2]);
        m_cpu_wall_positions[GPU_WALL_LEFT  ] = make_double3(sim.walls[1].pos[0], sim.walls[1].pos[1], sim.walls[1].pos[2]);
        m_cpu_wall_positions[GPU_WALL_RIGHT ] = make_double3(sim.walls[2].pos[0], sim.walls[2].pos[1], sim.walls[2].pos[2]);
        m_cpu_wall_positions[GPU_WALL_FRONT ] = make_double3(sim.walls[3].pos[0], sim.walls[3].pos[1], sim.walls[3].pos[2]);
        m_cpu_wall_positions[GPU_WALL_BACK  ] = make_double3(sim.walls[4].pos[0], sim.walls[4].pos[1], sim.walls[4].pos[2]);
        m_cpu_wall_positions[GPU_WALL_TOP   ] = make_double3(sim.walls[5].pos[0], sim.walls[5].pos[1], sim.walls[5].pos[2]);
    }

    else if(sim.walls.size() == 2)
    {

        m_wallConstants.wall_normals[GPU_WALL_BOTTOM] = make_double3(sim.walls[0].normal[0], sim.walls[0].normal[1], sim.walls[0].normal[2]);
        m_wallConstants.wall_normals[GPU_WALL_LEFT  ] = make_double3(0.0, 0.0, 0.0);
        m_wallConstants.wall_normals[GPU_WALL_RIGHT ] = make_double3(0.0, 0.0, 0.0);
        m_wallConstants.wall_normals[GPU_WALL_FRONT ] = make_double3(0.0, 0.0, 0.0);
        m_wallConstants.wall_normals[GPU_WALL_BACK  ] = make_double3(0.0, 0.0, 0.0);
        m_wallConstants.wall_normals[GPU_WALL_TOP   ] = make_double3(sim.walls[1].normal[0], sim.walls[1].normal[1], sim.walls[1].normal[2]);


        m_wallConstants.wall_compression_modifier[GPU_WALL_BOTTOM]  = sim.walls[0].compression_modifier;
        m_wallConstants.wall_compression_modifier[GPU_WALL_LEFT  ]  = 0.0;
        m_wallConstants.wall_compression_modifier[GPU_WALL_RIGHT ]  = 0.0;
        m_wallConstants.wall_compression_modifier[GPU_WALL_FRONT ]  = 0.0;
        m_wallConstants.wall_compression_modifier[GPU_WALL_BACK  ]  = 0.0;
        m_wallConstants.wall_compression_modifier[GPU_WALL_TOP   ]  = sim.walls[1].compression_modifier;


        m_wallConstants.wall_sliding_modifier[GPU_WALL_BOTTOM] = sim.walls[0].sliding_modifier;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_LEFT  ] = 0.0;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_RIGHT ] = 0.0;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_FRONT ] = 0.0;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_BACK  ] = 0.0;
        m_wallConstants.wall_sliding_modifier[GPU_WALL_TOP   ] = sim.walls[1].sliding_modifier;


        m_wallConstants.rolling_modifier[GPU_WALL_BOTTOM] = sim.walls[0].rolling_modifier;
        m_wallConstants.rolling_modifier[GPU_WALL_LEFT  ] = 0.0;
        m_wallConstants.rolling_modifier[GPU_WALL_RIGHT ] = 0.0;
        m_wallConstants.rolling_modifier[GPU_WALL_FRONT ] = 0.0;
        m_wallConstants.rolling_modifier[GPU_WALL_BACK  ] = 0.0;
        m_wallConstants.rolling_modifier[GPU_WALL_TOP   ] = sim.walls[1].rolling_modifier;


        m_wallConstants.twisting_modifier[GPU_WALL_BOTTOM] = sim.walls[0].twisting_modifier;
        m_wallConstants.twisting_modifier[GPU_WALL_LEFT  ] = 0.0;
        m_wallConstants.twisting_modifier[GPU_WALL_RIGHT ] = 0.0;
        m_wallConstants.twisting_modifier[GPU_WALL_FRONT ] = 0.0;
        m_wallConstants.twisting_modifier[GPU_WALL_BACK  ] = 0.0;
        m_wallConstants.twisting_modifier[GPU_WALL_TOP   ] = sim.walls[1].twisting_modifier;


        // load wall positions from simulation walls vector to array in particleWallInteraction
        m_cpu_wall_positions[GPU_WALL_BOTTOM] = make_double3(sim.walls[0].pos[0], sim.walls[0].pos[1], sim.walls[0].pos[2]);
        m_cpu_wall_positions[GPU_WALL_LEFT] = make_double3(0.0, 0.0, 0.0);
        m_cpu_wall_positions[GPU_WALL_RIGHT] = make_double3(0.0, 0.0, 0.0);
        m_cpu_wall_positions[GPU_WALL_FRONT] = make_double3(0.0, 0.0, 0.0);
        m_cpu_wall_positions[GPU_WALL_BACK] =make_double3(0.0, 0.0, 0.0);
        m_cpu_wall_positions[GPU_WALL_TOP] = make_double3(sim.walls[1].pos[0], sim.walls[1].pos[1], sim.walls[1].pos[2]);

    }







    // Load existing wall contacts to GPU
    int max_number_of_wall_contacts = 3*m_number_of_particles;

    int* cpu_wall_particle_ids = new int[m_number_of_particles];
    double * cpu_dist_to_wall = new double[max_number_of_wall_contacts];
    double3* cpu_wall_contact_vector_new = new double3[max_number_of_wall_contacts];
    double3* cpu_wall_contact_vector_old = new double3[max_number_of_wall_contacts];
    double3* cpu_wall_contact_n_initial  = new double3[max_number_of_wall_contacts];
    double4* cpu_wall_contact_rot        = new double4[max_number_of_wall_contacts];
    double3* cpu_wall_contact_n          = new double3[max_number_of_wall_contacts];
    double3* cpu_wall_initial_contact    = new double3[max_number_of_wall_contacts];
    double * cpu_wall_twisting_displacement=new double[max_number_of_wall_contacts];


    memset(cpu_wall_particle_ids,           0, m_number_of_particles    *  sizeof(int    ));
    memset(cpu_dist_to_wall,                0, max_number_of_wall_contacts*sizeof(double ));
    memset(cpu_wall_contact_vector_new,     0, max_number_of_wall_contacts*sizeof(double3));
    memset(cpu_wall_contact_vector_old,     0, max_number_of_wall_contacts*sizeof(double3));
    memset(cpu_wall_contact_n_initial,      0, max_number_of_wall_contacts*sizeof(double3));
    memset(cpu_wall_contact_rot,            0, max_number_of_wall_contacts*sizeof(double4));
    memset(cpu_wall_contact_n,              0, max_number_of_wall_contacts*sizeof(double3));
    memset(cpu_wall_initial_contact,        0, max_number_of_wall_contacts*sizeof(double3));
    memset(cpu_wall_twisting_displacement,  0, max_number_of_wall_contacts*sizeof(double ));




    for(int p = 0; p < m_number_of_particles; ++p)
    {

        const ContactListEntry *cl_entry = sim.contact_list[p];

        while(cl_entry)
        {

            if(cl_entry->id < 0) // contact with wall
            {

                int bitmap = cpu_wall_particle_ids[p];

                int wall_id = wall_id_CPU_to_GPU(WALL_ID(cl_entry->id), sim.number_of_walls);

                int wall_id_id = cpu_find_free_wall_contact_id(bitmap);

                if(wall_id_id == 0)
                {
                    printf("Error finding free wall id in load_contacts_to_gpu \n");
                }


                bitmap |= 1 << (wall_id);
                bitmap |= 1 << (wall_id + 8 * wall_id_id);

                int contact_id = p + (wall_id_id-1)*m_number_of_particles;

                double  dist_to_wall = sim.walls[WALL_ID(cl_entry->id)].getDistanceTo(sim.pos_new[X_COORD(p)], sim.pos_new[Y_COORD(p)], sim.pos_new[Z_COORD(p)]);


                double3 wall_initial_contact =
                        make_double3
                        (
                            cl_entry->contact->n2_initial[0],
                            cl_entry->contact->n2_initial[1],
                            cl_entry->contact->n2_initial[2]
                            );


                double3 wall_pos =
                        make_double3
                        (
                            sim.walls[WALL_ID(cl_entry->id)].pos[0],
                            sim.walls[WALL_ID(cl_entry->id)].pos[1],
                            sim.walls[WALL_ID(cl_entry->id)].pos[2]
                        );


                double3 pos_old =
                        make_double3
                        (
                        sim.pos_old[X_COORD(p)],
                        sim.pos_old[Y_COORD(p)],
                        sim.pos_old[Z_COORD(p)]
                            );


                double3 pos_new =
                        make_double3
                        (
                        sim.pos_new[X_COORD(p)],
                        sim.pos_new[Y_COORD(p)],
                        sim.pos_new[Z_COORD(p)]
                            );

                double3 wall_contact_pos = wall_pos + wall_initial_contact;

                double3 wall_contact_vector_new = pos_new - wall_contact_pos;
                double3 wall_contact_vector_old = pos_old - wall_contact_pos;


                double3 wall_contact_n_initial =
                        make_double3
                        (
                            cl_entry->contact->n1_initial[0],
                            cl_entry->contact->n1_initial[1],
                            cl_entry->contact->n1_initial[2]
                        );


                double4 wall_contact_rot =
                        make_double4
                        (
                            cl_entry->contact->rot1.e1,
                            cl_entry->contact->rot1.e2,
                            cl_entry->contact->rot1.e3,
                            cl_entry->contact->rot1.e0
                        );


                vec3 n1;
                cl_entry->contact->getCurrentN1(&n1);
                double3 wall_contact_n = make_double3(n1[0], n1[1], n1[2]);


                double wall_twisting_displacement = cl_entry->contact->twisting_displacement;

                cpu_wall_particle_ids[p] = bitmap;
                cpu_dist_to_wall[contact_id] = dist_to_wall;
                cpu_wall_contact_vector_new[contact_id] = wall_contact_vector_new;
                cpu_wall_contact_vector_old[contact_id] = wall_contact_vector_old;
                cpu_wall_contact_n_initial[contact_id] = wall_contact_n_initial;
                cpu_wall_contact_rot[contact_id] = wall_contact_rot;
                cpu_wall_contact_n[contact_id] = wall_contact_n;
                cpu_wall_initial_contact[contact_id] = wall_initial_contact;
                cpu_wall_twisting_displacement[contact_id] = wall_twisting_displacement;

            }

            cl_entry = cl_entry->next;
        }

    }

    copyToGPU((void*)m_gpu_wall_particle_ids,  (void*)cpu_wall_particle_ids, m_number_of_particles * sizeof(int));
    copyToGPU((void*)m_gpu_dist_to_wall, (void*)cpu_dist_to_wall, max_number_of_wall_contacts * sizeof(double));
    copyToGPU((void*)m_gpu_wall_contact_vector_new, (void*)cpu_wall_contact_vector_new, max_number_of_wall_contacts * sizeof(double3));
    copyToGPU((void*)m_gpu_wall_contact_vector_old, (void*)cpu_wall_contact_vector_old, max_number_of_wall_contacts * sizeof(double3));
    copyToGPU((void*)m_gpu_wall_contact_n_initial, (void*)cpu_wall_contact_n_initial, max_number_of_wall_contacts * sizeof(double3));
    copyToGPU((void*)m_gpu_wall_contact_rot, (void*)cpu_wall_contact_rot, max_number_of_wall_contacts * sizeof(double4));
    copyToGPU((void*)m_gpu_wall_contact_n, (void*)cpu_wall_contact_n, max_number_of_wall_contacts * sizeof(double3));
    copyToGPU((void*)m_gpu_wall_initial_contact, (void*)cpu_wall_initial_contact, max_number_of_wall_contacts * sizeof(double3));
    copyToGPU((void*)m_gpu_wall_twisting_displacement, (void*)cpu_wall_twisting_displacement, max_number_of_wall_contacts * sizeof(double));



    delete[] cpu_wall_particle_ids;
    delete[] cpu_dist_to_wall;
    delete[] cpu_wall_contact_vector_new;
    delete[] cpu_wall_contact_vector_old;
    delete[] cpu_wall_contact_n_initial;
    delete[] cpu_wall_contact_rot;
    delete[] cpu_wall_contact_n;
    delete[] cpu_wall_initial_contact;
    delete[] cpu_wall_twisting_displacement;


    init_GPUWallConstansWrapper(m_wallConstants);

    m_wall_box_loaded = true;
    m_wall_interaction_initialized = m_wall_interaction_allocated && m_wall_gpu_constants_initialized   && m_wall_box_loaded;




    return;


}




void particleWallInteraction::get_box_pos_from_gpu(Simulation &sim)
{
    // Load box to CPU
    sim.walls[sim.box->top_wall_id].pos[0] = m_cpu_wall_positions[GPU_WALL_TOP].x;
    sim.walls[sim.box->top_wall_id].pos[1] = m_cpu_wall_positions[GPU_WALL_TOP].y;
    sim.walls[sim.box->top_wall_id].pos[2] = m_cpu_wall_positions[GPU_WALL_TOP].z;

    sim.walls[sim.box->bottom_wall_id].pos[0] = m_cpu_wall_positions[GPU_WALL_BOTTOM].x;
    sim.walls[sim.box->bottom_wall_id].pos[1] = m_cpu_wall_positions[GPU_WALL_BOTTOM].y;
    sim.walls[sim.box->bottom_wall_id].pos[2] = m_cpu_wall_positions[GPU_WALL_BOTTOM].z;

    if(sim.walls.size() == 6)
    {
        sim.walls[1].pos[0] = m_cpu_wall_positions[GPU_WALL_LEFT ].x;
        sim.walls[1].pos[1] = m_cpu_wall_positions[GPU_WALL_LEFT ].y;
        sim.walls[1].pos[2] = m_cpu_wall_positions[GPU_WALL_LEFT ].z;

        sim.walls[2].pos[0] = m_cpu_wall_positions[GPU_WALL_RIGHT].x;
        sim.walls[2].pos[1] = m_cpu_wall_positions[GPU_WALL_RIGHT].y;
        sim.walls[2].pos[2] = m_cpu_wall_positions[GPU_WALL_RIGHT].z;

        sim.walls[3].pos[0] = m_cpu_wall_positions[GPU_WALL_FRONT].x;
        sim.walls[3].pos[1] = m_cpu_wall_positions[GPU_WALL_FRONT].y;
        sim.walls[3].pos[2] = m_cpu_wall_positions[GPU_WALL_FRONT].z;

        sim.walls[4].pos[0] = m_cpu_wall_positions[GPU_WALL_BACK ].x;
        sim.walls[4].pos[1] = m_cpu_wall_positions[GPU_WALL_BACK ].y;
        sim.walls[4].pos[2] = m_cpu_wall_positions[GPU_WALL_BACK ].z;
    }


    sim.updateBox();

    return;
}





void particleWallInteraction::load_box_from_gpu(Simulation &sim)
{



    // Load box to CPU
    sim.walls[sim.box->top_wall_id].pos[0] = m_cpu_wall_positions[GPU_WALL_TOP].x;
    sim.walls[sim.box->top_wall_id].pos[1] = m_cpu_wall_positions[GPU_WALL_TOP].y;
    sim.walls[sim.box->top_wall_id].pos[2] = m_cpu_wall_positions[GPU_WALL_TOP].z;

    sim.walls[sim.box->bottom_wall_id].pos[0] = m_cpu_wall_positions[GPU_WALL_BOTTOM].x;
    sim.walls[sim.box->bottom_wall_id].pos[1] = m_cpu_wall_positions[GPU_WALL_BOTTOM].y;
    sim.walls[sim.box->bottom_wall_id].pos[2] = m_cpu_wall_positions[GPU_WALL_BOTTOM].z;

    if(sim.walls.size() == 6)
    {
        sim.walls[1].pos[0] = m_cpu_wall_positions[GPU_WALL_LEFT ].x;
        sim.walls[1].pos[1] = m_cpu_wall_positions[GPU_WALL_LEFT ].y;
        sim.walls[1].pos[2] = m_cpu_wall_positions[GPU_WALL_LEFT ].z;

        sim.walls[2].pos[0] = m_cpu_wall_positions[GPU_WALL_RIGHT].x;
        sim.walls[2].pos[1] = m_cpu_wall_positions[GPU_WALL_RIGHT].y;
        sim.walls[2].pos[2] = m_cpu_wall_positions[GPU_WALL_RIGHT].z;

        sim.walls[3].pos[0] = m_cpu_wall_positions[GPU_WALL_FRONT].x;
        sim.walls[3].pos[1] = m_cpu_wall_positions[GPU_WALL_FRONT].y;
        sim.walls[3].pos[2] = m_cpu_wall_positions[GPU_WALL_FRONT].z;

        sim.walls[4].pos[0] = m_cpu_wall_positions[GPU_WALL_BACK ].x;
        sim.walls[4].pos[1] = m_cpu_wall_positions[GPU_WALL_BACK ].y;
        sim.walls[4].pos[2] = m_cpu_wall_positions[GPU_WALL_BACK ].z;
    }


    sim.updateBox();

    // Load existing wall contacts to CPU
    int max_number_of_wall_contacts = 3*m_number_of_particles;

    int* cpu_wall_particle_ids = new int[m_number_of_particles];
    double3* cpu_wall_contact_n_initial  = new double3[max_number_of_wall_contacts];
    double4* cpu_wall_contact_rot        = new double4[max_number_of_wall_contacts];
    double3* cpu_wall_contact_n          = new double3[max_number_of_wall_contacts];
    double3* cpu_wall_initial_contact    = new double3[max_number_of_wall_contacts];
    double * cpu_wall_twisting_displacement=new double[max_number_of_wall_contacts];


    copyFromGPU((void*)cpu_wall_particle_ids, (void*)m_gpu_wall_particle_ids,  m_number_of_particles * sizeof(int));
    copyFromGPU((void*)cpu_wall_contact_n_initial, (void*)m_gpu_wall_contact_n_initial, max_number_of_wall_contacts * sizeof(double3));
    copyFromGPU((void*)cpu_wall_contact_rot, (void*)m_gpu_wall_contact_rot, max_number_of_wall_contacts * sizeof(double4));
    copyFromGPU((void*)cpu_wall_contact_n, (void*)m_gpu_wall_contact_n, max_number_of_wall_contacts * sizeof(double3));
    copyFromGPU((void*)cpu_wall_initial_contact, (void*)m_gpu_wall_initial_contact, max_number_of_wall_contacts * sizeof(double3));
    copyFromGPU((void*)cpu_wall_twisting_displacement, (void*)m_gpu_wall_twisting_displacement, max_number_of_wall_contacts * sizeof(double));


    for(int p = 0; p < m_number_of_particles; ++p)
    {
        const int contact_bitmap = cpu_wall_particle_ids[p];

        if((contact_bitmap & 255) != 0) // if particle is in contact with any wall at all
        {
            for(int wall_id_id = 1; wall_id_id < 4; ++wall_id_id)
            {

                // get index of the wall
                int wall_id = contact_bitmap >> (8 * wall_id_id); // ingore the bytes before  wall_id_id

                wall_id = wall_id & 255; //ingore the bytes after wall_id_id

                if (wall_id == 0) // no contact with any wall stored at this position
                    continue;

                int targetlevel = 0;
                while (wall_id >>= 1) ++targetlevel;

                wall_id = targetlevel;


                int contact_id = p + m_number_of_particles * (wall_id_id-1);


                int cpu_wall_id = wall_id_GPU_to_CPU(wall_id, sim.number_of_walls);


                Contact *new_contact = new Contact;


                new_contact->id1 = p;
                new_contact->id2 = -cpu_wall_id-1;

                new_contact->n1_initial[0] = cpu_wall_contact_n_initial[contact_id].x;
                new_contact->n1_initial[1] = cpu_wall_contact_n_initial[contact_id].y;
                new_contact->n1_initial[2] = cpu_wall_contact_n_initial[contact_id].z;

                new_contact->n2_initial[0] = cpu_wall_initial_contact[contact_id].x;
                new_contact->n2_initial[1] = cpu_wall_initial_contact[contact_id].y;
                new_contact->n2_initial[2] = cpu_wall_initial_contact[contact_id].z;

                new_contact->old_contact_normal[0] = 0.0; // unsued for wall contacts
                new_contact->old_contact_normal[1] = 0.0;
                new_contact->old_contact_normal[2] = 0.0;

                new_contact->rot1.e0 = cpu_wall_contact_rot[contact_id].w;
                new_contact->rot1.e1 = cpu_wall_contact_rot[contact_id].x;
                new_contact->rot1.e2 = cpu_wall_contact_rot[contact_id].y;
                new_contact->rot1.e3 = cpu_wall_contact_rot[contact_id].z;

                new_contact->rot2.e0 = 1.0; // unused for wall contacts
                new_contact->rot2.e1 = 0.0;
                new_contact->rot2.e2 = 0.0;
                new_contact->rot2.e3 = 0.0;

                new_contact->twisting_displacement = cpu_wall_twisting_displacement[contact_id];


                // check if particle is not already in contact with wall
                ContactListEntry *cl_entry = sim.getContactInsertPos(p, -cpu_wall_id-1);

                if(cl_entry)
                {

                    // append contact
                    cl_entry->contact = new_contact;

                }
            }
        }
    }


    delete[] cpu_wall_particle_ids;
    delete[] cpu_wall_contact_n_initial;
    delete[] cpu_wall_contact_rot;
    delete[] cpu_wall_contact_n;
    delete[] cpu_wall_initial_contact;
    delete[] cpu_wall_twisting_displacement;



}






void particleWallInteraction::storeToFile(FILE* file)
{


    int max_number_of_wall_contacts = 3*m_number_of_particles;


    fwrite(&m_number_of_particles, sizeof(int), 1, file);

    fwrite(&m_wall_interaction_allocated, sizeof(bool), 1, file);
    fwrite(&m_wall_gpu_constants_initialized, sizeof(bool), 1, file);
    fwrite(&m_wall_interaction_initialized, sizeof(bool), 1, file);
    fwrite(&m_wall_box_loaded, sizeof(bool), 1, file);


    fwrite(&m_wallConstants, sizeof(GPUWallConstants), 1, file);





    int* int_buffer = new int[m_number_of_particles];
    copyFromGPU(int_buffer, m_gpu_wall_particle_ids, m_number_of_particles*sizeof(int));
    fwrite(int_buffer, sizeof(int), m_number_of_particles, file);

    delete[] int_buffer;


    double* double_buffer = new double[m_number_of_particles];

    copyFromGPU(double_buffer, m_gpu_top_wall_force, m_number_of_particles*sizeof(double));
    fwrite(double_buffer, sizeof(double), m_number_of_particles, file);

    copyFromGPU(double_buffer, m_gpu_bot_wall_force, m_number_of_particles*sizeof(double));
    fwrite(double_buffer, sizeof(double), m_number_of_particles, file);


    double* double_contact_buffer = new double[max_number_of_wall_contacts];

    copyFromGPU(double_contact_buffer, m_gpu_dist_to_wall, max_number_of_wall_contacts*sizeof(double));
    fwrite(double_contact_buffer, sizeof(double), max_number_of_wall_contacts, file);

    copyFromGPU(double_contact_buffer, m_gpu_wall_twisting_displacement, max_number_of_wall_contacts*sizeof(double));
    fwrite(double_contact_buffer, sizeof(double), max_number_of_wall_contacts, file);

    delete[] double_contact_buffer;


    double3* double3_buffer = new double3[max_number_of_wall_contacts];

    copyFromGPU(double3_buffer, m_gpu_wall_contact_vector_new, max_number_of_wall_contacts*sizeof(double3));
    fwrite(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);

    copyFromGPU(double3_buffer, m_gpu_wall_contact_vector_old, max_number_of_wall_contacts*sizeof(double3));
    fwrite(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);

    copyFromGPU(double3_buffer, m_gpu_wall_contact_n_initial, max_number_of_wall_contacts*sizeof(double3));
    fwrite(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);

    copyFromGPU(double3_buffer, m_gpu_wall_contact_n, max_number_of_wall_contacts*sizeof(double3));
    fwrite(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);

    copyFromGPU(double3_buffer, m_gpu_wall_initial_contact, max_number_of_wall_contacts*sizeof(double3));
    fwrite(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);


    delete[] double3_buffer;

    double4* double4_buffer = new double4[max_number_of_wall_contacts];

    copyFromGPU(double4_buffer, m_gpu_wall_contact_rot, max_number_of_wall_contacts*sizeof(double4));
    fwrite(double4_buffer, sizeof(double4), max_number_of_wall_contacts, file);

    fwrite(m_cpu_wall_positions, sizeof(double3), 6, file);

    delete[] double4_buffer;

    return;

}




void particleWallInteraction::loadFromFile(FILE* file)
{


    fread(&m_number_of_particles, sizeof(int), 1, file);

    this->initParticleWallInteraction(m_number_of_particles);

    fread(&m_wall_interaction_allocated, sizeof(bool), 1, file);
    fread(&m_wall_gpu_constants_initialized, sizeof(bool), 1, file);
    fread(&m_wall_interaction_initialized, sizeof(bool), 1, file);
    fread(&m_wall_box_loaded, sizeof(bool), 1, file);

    int max_number_of_wall_contacts = 3*m_number_of_particles;

    fread(&m_wallConstants, sizeof(GPUWallConstants), 1, file);

    init_GPUWallConstansWrapper(m_wallConstants);


    int* int_buffer = new int[m_number_of_particles];
    fread(int_buffer, sizeof(int), m_number_of_particles, file);
    copyToGPU(m_gpu_wall_particle_ids, int_buffer, m_number_of_particles*sizeof(int));
    delete[] int_buffer;


    double* double_buffer = new double[m_number_of_particles];
    fread(double_buffer, sizeof(double), m_number_of_particles, file);
    copyToGPU(m_gpu_top_wall_force, double_buffer, m_number_of_particles*sizeof(double));

    fread(double_buffer, sizeof(double), m_number_of_particles, file);
    copyToGPU(m_gpu_bot_wall_force, double_buffer, m_number_of_particles*sizeof(double));



    double* double_contact_buffer = new double[max_number_of_wall_contacts];
    fread(double_contact_buffer, sizeof(double), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_dist_to_wall, double_contact_buffer, max_number_of_wall_contacts*sizeof(double));

    fread(double_contact_buffer, sizeof(double), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_twisting_displacement, double_contact_buffer, max_number_of_wall_contacts*sizeof(double));

    delete[] double_contact_buffer;


    double3* double3_buffer = new double3[max_number_of_wall_contacts];

    fread(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_contact_vector_new, double3_buffer, max_number_of_wall_contacts*sizeof(double3));

    fread(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_contact_vector_old, double3_buffer, max_number_of_wall_contacts*sizeof(double3));

    fread(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_contact_n_initial, double3_buffer, max_number_of_wall_contacts*sizeof(double3));

    fread(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_contact_n, double3_buffer, max_number_of_wall_contacts*sizeof(double3));

    fread(double3_buffer, sizeof(double3), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_initial_contact, double3_buffer, max_number_of_wall_contacts*sizeof(double3));



    delete[] double3_buffer;

    double4* double4_buffer = new double4[max_number_of_wall_contacts];

    fread(double4_buffer, sizeof(double4), max_number_of_wall_contacts, file);
    copyToGPU(m_gpu_wall_contact_rot, double4_buffer, max_number_of_wall_contacts*sizeof(double4));
    delete[] double4_buffer;

    fread(m_cpu_wall_positions, sizeof(double3), 6, file);

    return;

}
