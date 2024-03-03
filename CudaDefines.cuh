
#ifndef CUDA_DEFINES_H
#define CUDA_DEFINES_H


//#define CPU_EQUIVALENT // try to reproduce the result of the CPU version  use together with nvcc compile flag : "-fmad=false"

// particle interaction
#define GPU_ROLLING
#define GPU_SLIDING
#define GPU_TWISTING

#define GPU_INELASTIC_ROLLING
#define GPU_INELASTIC_SLIDING
#define GPU_INELASTIC_TWISTING


// particle - wall interaction


#define GPU_WALL_ROLLING
#define GPU_WALL_SLIDING
#define GPU_WALL_TWISTING

#define GPU_WALL_INELASTIC_ROLLING
#define GPU_WALL_INELASTIC_SLIDING
#define GPU_WALL_INELASTIC_TWISTING

// enable wall damping method: simple (weak) damping of normal oscillations or viscoelastic damping proposed by Sebastian Krijt
//#define GPU_WALL_USE_CONSTANT_DAMPING
#define GPU_WALL_USE_VISCOELASTIC_DAMPING



// enable damping method: simple (weak) damping of normal oscillations or viscoelastic damping proposed by Sebastian Krijt
//#define GPU_USE_CONSTANT_DAMPING
#define GPU_USE_VISCOELASTIC_DAMPING

// enables tracking of dissipated energy
// Warning: has bug and can cause segmentation fault
// #define GPU_TRACK_DISSIPATED_ENERGY


#define BLOCK_SIZE 128

// addressing scheme for vectors
#define X_COORD_GPU(index, number_of_particles) (index)
#define Y_COORD_GPU(index, number_of_particles) (number_of_particles + index)
#define Z_COORD_GPU(index, number_of_particles) (2*number_of_particles + index)

// contact list
#define MAX_CONTACTS 12
#define MAX_WALL_CONTACTS 3

// adressing scheme for local contact list
#define CONTACT_ID(contact_list_pos, particle_index, number_of_particles) ( (contact_list_pos) * (number_of_particles) + (particle_index) )

struct GPUConstants
{
    double c1_contact_radius;
    double c2_contact_radius;
    double c1_normal_force;
    double c2_normal_force;
    double eq_contact_radius;
    double delta_c;
    double delta_0;
    double contact_breaking_dist;
    double contact_breaking_dist_squared;

    double particle_radius;
    double mass;
    double moment_of_inertia;

    double particle_radius_inv;
    double mass_inv;
    double moment_of_inertia_inv;

    double F_c;

    double k_r;                             // actually k_r * reduced_radius^2
    double k_s;                             // actually k_s * particle_radius^2
    double k_t;

    double crit_rolling_displacement;
    double crit_rolling_displacement_squared;
    double crit_sliding_displacement;
    double crit_sliding_displacement_squared;
    double crit_twisting_displacement;

    double osc_damping_factor;                  // actually used for krijt's viscoelastic damping constant

    // grid specifications
    unsigned int gpu_num_grid_cells;
    double gpu_grid_cell_width_inv;
    double gpu_grid_shift;

    double get_contact_radius_c1;
    double get_contact_radius_c2;

    double contact_making_dist_squared;
    double k_hertz;

#ifdef CPU_EQUIVALENT
    double reduced_radius;
#endif
};


struct GPUWallConstants
{
        double c1_contact_radius;
        double c2_contact_radius;
        double c1_normal_force;
        double c2_normal_force;
        double eq_contact_radius;
        double delta_c;
        double contact_breaking_dist;

        double particle_radius;
        double reduced_particle_radius;

        double particle_radius_inv;
        double moment_of_inertia_inv;


        double k_r;
        double k_s;
        double k_t;

        double crit_rolling_displacement;
        double crit_rolling_displacement_squared;
        double crit_sliding_displacement;
        double crit_sliding_displacement_squared;
        double crit_twisting_displacement;

        double osc_damping_factor;              // actually used for krijt's viscoelastic damping constant

        double top_wall_speed;
        double3 wall_normals[6];                // normals of the walls
        double wall_compression_modifier[6];    // modifier the the force in normal direction the walls exert on the particles
        double wall_sliding_modifier[6];
        double rolling_modifier[6];
        double twisting_modifier[6];

        double get_contact_radius_c1;
        double get_contact_radius_c2;

};


// Grid
#define NO_PARTICLE 0xffffffff


// Other stuff
#define RESTRICT __restrict


// right handed coordinate system:
//  the positive x and y axes point right and up,
//	and the negative z axis points forward
#define GPU_WALL_BOTTOM	0 // -y
#define GPU_WALL_TOP	1 // +y
#define GPU_WALL_LEFT	2 // -x
#define GPU_WALL_RIGHT	3 // +x
#define GPU_WALL_FRONT	4 // +z
#define GPU_WALL_BACK	5 // -z




#endif // CUDA_DEFINES_H
