#ifndef __CUDA_WALL_KERNELS_H__
#define __CUDA_WALL_KERNELS_H__

#include "CudaDefines.cuh"
#include "my_device_helper_math.h"

#include "device_launch_parameters.h"
#include "cuda_runtime.h"

#include <cuda.h>
#include <stdio.h>


__device__ __constant__ GPUWallConstants gpu_wall_constants;




void init_GPUWallConstansKernel(const GPUWallConstants &cpu_constants)
{
    checkCudaErrors(cudaMemcpyToSymbol(gpu_wall_constants, &cpu_constants, sizeof(GPUWallConstants)));
}




////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contact Management
////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * updateWallContactNormalsKernel
 * goes through all contacts and calculates the contact normal
 * and checks wether a contact has been broken
 * if a contact has been broken it will be deletet
 */
__global__
void updateWallContactNormalsKernel(
    int* RESTRICT contact_wall_particle_ids,
    double * RESTRICT dist_to_wall,
    double3* RESTRICT wall_contact_vector_new,
    const double* RESTRICT positions,
    const double3* RESTRICT wall_initial_contact,
    const double3* RESTRICT wall_positions,
    const int number_of_particles
        )
{
    const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index >= MAX_WALL_CONTACTS * number_of_particles)
        return;

    // determine ids
    const int particle_id = index % number_of_particles;
    const int wall_id_id = index / number_of_particles + 1; // index of the wall_id in contact_bitmap

    int contact_bitmap = contact_wall_particle_ids[particle_id];

    // TODO: check again if this works
    int wall_id = contact_bitmap>>(8*wall_id_id); // remove bytes before the wall_id

    wall_id = wall_id & 255; // only consider the first byte and set the rest to 0

    if (wall_id == 0) // no contact with any wall stored at this position
        return;

    // wall_id = find which bit is set to 1
    wall_id = __ffs(wall_id)-1; // __ffs(int C) finds the position least significant bit in C   with __ffs(0) = 0;

    //////////////////////////////////////////////////////////////////////////////
    // update contact normals
    //////////////////////////////////////////////////////////////////////////////
    double3 pos = make_double3(positions, particle_id, number_of_particles);

    double dist = fabs(dot(pos - wall_positions[wall_id], gpu_wall_constants.wall_normals[wall_id]));

    if (dist > gpu_wall_constants.contact_breaking_dist) // if distance between particle is too big, break the contact and return
    {

        contact_bitmap &= ~(1 << (wall_id));                // set bit at wall_id in first byte to 0
        contact_bitmap &= ~(1 << (wall_id + 8*wall_id_id)); // set bit at wall_id in the wall_id_id byte to 0   ( wall_id_id byte is the 2nd 3rd or 4th byte )

        contact_wall_particle_ids[particle_id] = contact_bitmap;

        return;
    }
    
    double3 wall_contact_pos = wall_initial_contact[index] + wall_positions[wall_id];

    double3 deltaPos = pos - wall_contact_pos;


    dist_to_wall[index] = dist;
    wall_contact_vector_new[index] = deltaPos;

    return;
}


__device__
int find_free_wall_contact_id(int bitmap)
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




/*
 * updateWallStickingKernel
 *
 * takes the current position of the walls and of the particles
 * and tests each particle for a collision with each wall
 * if a collision occurs a contact will be initialized and stored
 *
 */
template<int number_of_walls> __global__
void updateWallStickingKernel(
    int* RESTRICT contact_wall_particle_ids,
    const double* RESTRICT positions,
    double * RESTRICT dist_to_wall,
    double3* RESTRICT wall_contact_vector_new,
    double3* RESTRICT wall_contact_vector_old,
    double3* RESTRICT wall_contact_n_intitial,
    double4* RESTRICT wall_contact_rot,
    double3* RESTRICT wall_contact_n,
    double3* RESTRICT wall_initial_contact,
    double * RESTRICT wall_twisting_displacement,
    const double3* RESTRICT wall_pos,
    const int number_of_particles
        )
{
    const uint index = blockIdx.x * blockDim.x + threadIdx.x;

    if(index >= number_of_particles)
        return;          // handle case when no. of particles not multiple of block size

    // load position
    const double3 pos = make_double3(positions, index, number_of_particles);

    int bitmap = contact_wall_particle_ids[index];
    bool bitmap_changed = false;



    for(int wall = 0; wall < number_of_walls; ++wall)
    {
        if(!((bitmap >> wall) & 1)) // if particle is not in contact with wall
        {
            double dist = fabs(dot(pos - wall_pos[wall], gpu_wall_constants.wall_normals[wall]));

            if(dist < 2.0 * gpu_wall_constants.particle_radius)
            {
                // choose free id spot in the bitmap
                int wall_id_id = find_free_wall_contact_id(bitmap);

                if(wall_id_id == 0)
                {
                    printf("Error finding free wall id: Pid = %d Wall = %d Bitmap = %d	dist = %e   pos = %e    %e  %e  wall pos = %e   %e  %e\n\n",
                           index,
                           wall,
                           bitmap,
                            dist / gpu_wall_constants.particle_radius,
                           pos.x,
                           pos.y,
                           pos.z,
                           wall_pos[wall].x,
                           wall_pos[wall].y,
                           wall_pos[wall].z);
                    return;
                }


                // set the particle-wall contact in bitmap
                bitmap |= 1 << (wall);
                bitmap |= 1 << (wall + 8*wall_id_id);
                bitmap_changed = true;

                // initialize contact
                int contact_id = index + (wall_id_id-1)*number_of_particles;


                wall_contact_rot[contact_id] = make_double4(0.0, 0.0, 0.0, 1.0);

                // n_initial is the contact normal from the particle's view
                wall_contact_n_intitial[contact_id] = - gpu_wall_constants.wall_normals[wall];

                // TODO hier alles negieren wall_contact_vector_new und wall_contact_vector_old
                // is contact normals_new and _old needed or is it always the same as the wall normal
                dist_to_wall[contact_id] = dist;
                wall_contact_vector_new[contact_id] = dist * gpu_wall_constants.wall_normals[wall];
                wall_contact_vector_old[contact_id] = gpu_wall_constants.particle_radius * 2.0 * gpu_wall_constants.wall_normals[wall];
                wall_contact_n[contact_id] = make_double3(0.0, 0.0, 0.0); // current normal
                wall_twisting_displacement[contact_id] = 0.0;


                // initial_contact is the particles position minus the position of the edge when the contact is made	see Simulation::updateSticking()
                // this is needed as an anchor to keep track of moving walls
                wall_initial_contact[contact_id] = pos - dist*gpu_wall_constants.wall_normals[wall] - wall_pos[wall];

            }

        }


        if(bitmap_changed)
        {
            contact_wall_particle_ids[index] = bitmap;
        }
    }
}


template __global__
void updateWallStickingKernel<6>(
    int* RESTRICT contact_wall_particle_ids,
    const double* RESTRICT positions,
    double * RESTRICT dist_to_wall,
    double3* RESTRICT wall_contact_vector_new,
    double3* RESTRICT wall_contact_vector_old,
    double3* RESTRICT wall_contact_n_intitial,
    double4* RESTRICT wall_contact_rot,
    double3* RESTRICT wall_contact_n,
    double3* RESTRICT wall_initial_contact,
    double * RESTRICT wall_twisting_displacement,
    const double3* RESTRICT wall_pos,
    const int number_of_particles
        );


template __global__
void updateWallStickingKernel<2>(
    int* RESTRICT contact_wall_particle_ids,
    const double* RESTRICT positions,
    double * RESTRICT dist_to_wall,
    double3* RESTRICT wall_contact_vector_new,
    double3* RESTRICT wall_contact_vector_old,
    double3* RESTRICT wall_contact_n_intitial,
    double4* RESTRICT wall_contact_rot,
    double3* RESTRICT wall_contact_n,
    double3* RESTRICT wall_initial_contact,
    double * RESTRICT wall_twisting_displacement,
    const double3* RESTRICT wall_pos,
    const int number_of_particles
        );












/*
 * Approximate the contact radius with a root function to speed up
 * the convergence of the newton raphson method
 */
__device__
double getApproxWallContactRadius(const double compression)
{
    double a = gpu_wall_constants.get_contact_radius_c1 + gpu_wall_constants.get_contact_radius_c2*sqrt(compression);
    return sqrt(a);
}



__device__
double getWallContactRadius(const double compression_length)
{
    // contact radius can be obtained by finding the root of a fourth order polynomial where x^2 = contact_radius
    // use equilibrium contact radius as starting value
    const double k = gpu_wall_constants.reduced_particle_radius * compression_length / 3.0;
    double x_pow3;
    double x_new;

    double compression = gpu_wall_constants.delta_c + compression_length;
    double x_old;

    if(compression  >= 0.0)
    {
#ifdef CPU_EQUIVALENT
        x_old = gpu_wall_constants.c1_contact_radius;
#else
        x_old = getApproxWallContactRadius(compression);
#endif
    }
    else
    {
        printf("Error getWallContactRadius: too long bond not correctly broken\n");
        return 0.0;
    }

    //double x_old = gpu_wall_constants.c1_contact_radius;

    // use Newton-Raphson method to find root
    for(int i = 0; i < 20; ++i)

    {
        x_pow3 = x_old * x_old * x_old;
        x_new = 0.75 * (x_pow3 * x_old + k) / (x_pow3 - gpu_wall_constants.c2_contact_radius);

        if (fabs(x_new - x_old) / x_new < 1.e-14)
            break;

        x_old = x_new;
    }

    return x_new * x_new;
}

__device__
double getWallNormalForce(const double compression_length, double &_contact_radius)
{

    const double contact_radius = getWallContactRadius(compression_length);

    double x = sqrt(contact_radius*contact_radius*contact_radius);

#ifdef GPU_USE_VISCOELASTIC_DAMPING
    _contact_radius = contact_radius;
#endif // GPU_USE_VISCOELASTIC_DAMPING

    return x * (gpu_wall_constants.c1_normal_force * x + gpu_wall_constants.c2_normal_force);

}





__global__
void updateRotation(
	double4* RESTRICT       rot,
	const double* RESTRICT  angular_velocities,
	const double* RESTRICT  torques,
	const int* RESTRICT     contact_wall_particle_ids,
	const double			timestep,
	const unsigned int		number_of_particles
)
{
	const unsigned int particle_id = blockIdx.x *blockDim.x + threadIdx.x;

	if (particle_id >= number_of_particles) // per particle and iterate over walls
		return;

	const int contact_bitmap = contact_wall_particle_ids[particle_id];


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// load torque and angular velocity of the particle
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	const double3 omega1 = make_double3(angular_velocities, particle_id, number_of_particles);

	const double3 omega1_dot = gpu_wall_constants.moment_of_inertia_inv * make_double3(torques, particle_id, number_of_particles);

	for (int index = particle_id;
		index < number_of_particles * MAX_WALL_CONTACTS;
		index += number_of_particles)
	{

		const int wall_id_id = index / number_of_particles + 1; // we start at wall_id_id = 1 which represents the second byte ; the first byte stores all the wall the particle is in contact with


        // get index of the wall
        int wall_id = contact_bitmap >> (8 * wall_id_id); // ingore the bytes before  wall_id_id

        wall_id = wall_id & 255; //ingore the bytes after wall_id_id

		if (wall_id == 0) // no contact with any wall stored at this position
			continue;


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// update rotation parameters
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		double4 e_dot;
		double4 e_ddot;
		double4 rot_param1 = rot[index];


		e_dot.w = -0.5 * dot(rot_param1, omega1);
        e_dot.x = 0.5 * (diff_of_products_simple(rot_param1.w, omega1.x, rot_param1.y, omega1.z) + rot_param1.z * omega1.y);
        e_dot.y = 0.5 * (diff_of_products_simple(rot_param1.w, omega1.y, rot_param1.z, omega1.x) + rot_param1.x * omega1.z);
        e_dot.z = 0.5 * (diff_of_products_simple(rot_param1.w, omega1.z, rot_param1.x, omega1.y) + rot_param1.y * omega1.x);

		double temp = 0.5 * e_dot.w;

		e_ddot.w = -0.25 * (rot_param1.w * dot(omega1, omega1) + 2.0 * dot(rot_param1, omega1_dot));
        e_ddot.x = temp * omega1.x + 0.5 * (diff_of_products_simple(rot_param1.w, omega1_dot.x, rot_param1.y, omega1_dot.z) + rot_param1.z * omega1_dot.y);
        e_ddot.y = temp * omega1.y + 0.5 * (diff_of_products_simple(rot_param1.w, omega1_dot.y, rot_param1.z, omega1_dot.x) + rot_param1.x * omega1_dot.z);
        e_ddot.z = temp * omega1.z + 0.5 * (diff_of_products_simple(rot_param1.w, omega1_dot.z, rot_param1.x, omega1_dot.y) + rot_param1.y * omega1_dot.x);

		//rot_param1 += timestep * e_dot + 0.5 * timestep * timestep * e_ddot; // einmal timestep ausklammern
		rot_param1 += timestep * (e_dot + 0.5 * timestep * e_ddot);

		// normalize
		temp = rsqrt(rot_param1.w*rot_param1.w + rot_param1.x*rot_param1.x + rot_param1.y*rot_param1.y + rot_param1.z*rot_param1.z);
		rot_param1 *= temp;


		rot[index] = rot_param1;
	
	}

}




__global__
void updateWallContactsKernel(
    double4* RESTRICT       rot,
    double3* RESTRICT       n1_initial,
    double3* RESTRICT       n_1,
    double * RESTRICT       forces_new,
    double * RESTRICT       torques_new,
    double * RESTRICT       wall_top_force, // force the particle exert on the top wall
    double * RESTRICT       wall_bot_force, // force the particle exert on the bottom wall
    const double * RESTRICT dist_to_wall,
    const double3* RESTRICT contact_vectors_new,
    const double3* RESTRICT contact_vectors_old,
          double3* RESTRICT initial_contact, // = n2_initial in alex code
    const double* RESTRICT  velocities,
    const double* RESTRICT  angular_velocities,
    const double* RESTRICT  torques,
    const int* RESTRICT     contact_wall_particle_ids,
    double* RESTRICT        twisting_displacement_arr,
    const double            timestep,
    const uint              number_of_particles
    )
{
    const unsigned int particle_id = blockIdx.x *blockDim.x + threadIdx.x;

    if (particle_id >= number_of_particles) // per particle and iterate over walls
        return;


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // load torque and angular velocity of the particle
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const double3 omega1 = make_double3(angular_velocities, particle_id, number_of_particles);

    const double3 omega1_dot = gpu_wall_constants.moment_of_inertia_inv * make_double3(torques, particle_id, number_of_particles);

    ///////////////////////////////////////////////////////////////////////////////
    // load the bitmap of the particle containing all the walls the particle is in contact with
    ///////////////////////////////////////////////////////////////////////////////

    const int contact_bitmap = contact_wall_particle_ids[particle_id];


    // variables for results
    double3 force = make_double3(0.0);
    double3 torque = make_double3(0.0);


    bool update_force_and_torque = false;

    //////////////////////////////////////////////////////////////////////////////
    // start the loop over all wall contacts of the particle
    //////////////////////////////////////////////////////////////////////////////

    for(int index = particle_id;
            index < number_of_particles * MAX_WALL_CONTACTS;
            index += number_of_particles)
    {

        const int wall_id_id = index/number_of_particles + 1; // we start at wall_id_id = 1 which represents the second byte ; the first byte stores all the wall the particle is in contact with


        // get index of the wall
        int wall_id = contact_bitmap>>(8*wall_id_id); // ignore the bytes before  wall_id_id

        wall_id = wall_id & 255; //ignore the bytes after wall_id_id

        if (wall_id == 0) // no contact with any wall stored at this position
            continue;


        update_force_and_torque = true;

        wall_id = __ffs(wall_id)-1; // wall_id is coded as the position of the first bit set to 1

        const double3 n2 = gpu_wall_constants.wall_normals[wall_id]; // normal of the wall the particle is in contact with
                                                                     // called n2 to be consinstent with Alex code in Simulation::updateParticleInteraction(double timestep)

		/*
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // update rotation parameters
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        double4 e_dot;
        double4 e_ddot;
        
        e_dot.w = -0.5 * dot(rot_param1, omega1);
        e_dot.x = 0.5 * (rot_param1.w * omega1.x - rot_param1.y * omega1.z + rot_param1.z * omega1.y);
        e_dot.y = 0.5 * (rot_param1.w * omega1.y - rot_param1.z * omega1.x + rot_param1.x * omega1.z);
        e_dot.z = 0.5 * (rot_param1.w * omega1.z - rot_param1.x * omega1.y + rot_param1.y * omega1.x);

        double temp = 0.5 * e_dot.w;

        e_ddot.w = -0.25 * (rot_param1.w * dot(omega1, omega1) + 2.0 * dot(rot_param1, omega1_dot));
        e_ddot.x = temp * omega1.x + 0.5 * (rot_param1.w * omega1_dot.x - rot_param1.y * omega1_dot.z + rot_param1.z * omega1_dot.y);
        e_ddot.y = temp * omega1.y + 0.5 * (rot_param1.w * omega1_dot.y - rot_param1.z * omega1_dot.x + rot_param1.x * omega1_dot.z);
        e_ddot.z = temp * omega1.z + 0.5 * (rot_param1.w * omega1_dot.z - rot_param1.x * omega1_dot.y + rot_param1.y * omega1_dot.x);

        //rot_param1 += timestep * e_dot + 0.5 * timestep * timestep * e_ddot; // einmal timestep ausklammern
        rot_param1 += timestep * (e_dot + 0.5 * timestep * e_ddot);

        // normalize
        temp = rsqrt(rot_param1.w*rot_param1.w + rot_param1.x*rot_param1.x + rot_param1.y*rot_param1.y + rot_param1.z*rot_param1.z);
        rot_param1 *= temp;

        rot[index] = rot_param1;
		*/

        //////////////////////////////////////////////////////////////////////////////
        // determine current contact pointer
        //////////////////////////////////////////////////////////////////////////////




#if defined(GPU_WALL_INELASTIC_ROLLING) || defined(GPU_WALL_SLIDING)

        double4 rot_param1 = rot[index];

        double3 n_initial = n1_initial[index]; // == Contact::getCurrentN1(...)
        double3 n1;
        n1.x = 2.0 * ((rot_param1.w * rot_param1.w + rot_param1.x * rot_param1.x - 0.5) * n_initial.x + diff_of_products_simple(rot_param1.x, rot_param1.y, rot_param1.z, rot_param1.w) * n_initial.y + (rot_param1.x * rot_param1.z + rot_param1.y * rot_param1.w) * n_initial.z);
        n1.y = 2.0 * ((rot_param1.x * rot_param1.y + rot_param1.z * rot_param1.w) * n_initial.x + (rot_param1.w * rot_param1.w + rot_param1.y * rot_param1.y - 0.5) * n_initial.y + diff_of_products_simple(rot_param1.y, rot_param1.z, rot_param1.x, rot_param1.w) * n_initial.z);
        n1.z = 2.0 * (diff_of_products_simple(rot_param1.x, rot_param1.z, rot_param1.y, rot_param1.w) * n_initial.x + (rot_param1.y * rot_param1.z + rot_param1.x * rot_param1.w) * n_initial.y + (rot_param1.w * rot_param1.w + rot_param1.z * rot_param1.z - 0.5) * n_initial.z);
#endif

		



        //////////////////////////////////////////////////////////////////////////////
        // update twisting displacement
        //////////////////////////////////////////////////////////////////////////////

        // use second order integration: Phi^n+1 = Phi^n + 0.5 (<delta_omega^n, n_c^n> + <delta_omega^n+1, n_c^n+1>)
        const double3 contact_vector = contact_vectors_new[index];       // n_c = pos_new - wall_contact_pos  | .w = distance between wall and center of particle
                                                                        // wall_contact_pos = wall_contact_pos_initial + wall_pos
        const double tmp = 1.0/sqrt(dot(make_double3(contact_vector), make_double3(contact_vector)));

        const double dist = dist_to_wall[index];
        const double3 n_c = contact_vector * tmp;


#if defined(GPU_WALL_TWISTING) || defined(GPU_WALL_USE_VISCOELASTIC_DAMPING)
        const double3 contact_vector_old = contact_vectors_old[index];	// n_c_old = pos_old - wall_contact_pos  | .w = old distance between wall and center of particle

        const double tmp2 = 1.0/sqrt(dot(make_double3(contact_vector_old), make_double3(contact_vector_old)));

        const double3 n_c_old = contact_vector_old * tmp2;
#endif


#ifdef GPU_WALL_TWISTING
        const double3 omega1_new = omega1 + timestep * omega1_dot;

        double twisting_displacement = twisting_displacement_arr[index] + 0.5 * timestep *
            (
                fma(omega1.x, n_c_old.x,
                fma(omega1.y, n_c_old.y,
                fma(omega1.z, n_c_old.z,
                fma(omega1_new.x, n_c.x,
                fma(omega1_new.y, n_c.y,
                omega1_new.z * n_c.z)))))
                );


    #ifdef GPU_WALL_INELASTIC_TWISTING
        if (twisting_displacement > gpu_wall_constants.crit_twisting_displacement)
            twisting_displacement = gpu_wall_constants.crit_twisting_displacement;
    #endif // GPU_WALL_INELASTIC_TWISTING

        torque -= gpu_wall_constants.twisting_modifier[wall_id] * gpu_wall_constants.k_t * twisting_displacement * make_double3(n_c);

        twisting_displacement_arr[index] = twisting_displacement;


#endif // GPU_WALL_TWISTING



        /////////////////////////////////////////////////////////////////////////////
        // force in normal direction
        /////////////////////////////////////////////////////////////////////////////

        const double compression_length = 2.0 * gpu_wall_constants.particle_radius - dist;

        double contact_radius;

        double force_normal = gpu_wall_constants.wall_compression_modifier[wall_id] * getWallNormalForce(compression_length, contact_radius);


        // track forces on wall before damping
        // we only want to track the forces on the top and bottom wall
        if(wall_id == GPU_WALL_TOP)
        {
            wall_top_force[particle_id] += force_normal;
        }

        if(wall_id == GPU_WALL_BOTTOM)
        {
            wall_bot_force[particle_id] -= force_normal;
        }

#ifdef GPU_WALL_USE_VISCOELASTIC_DAMPING

        // viscoelastic model from krijt 2013

        double v_rel = dot(make_double3(velocities, particle_id, number_of_particles), gpu_wall_constants.wall_normals[wall_id]);
        if(wall_id == GPU_WALL_TOP)
            v_rel += gpu_wall_constants.top_wall_speed;

        double damping_strength = gpu_wall_constants.osc_damping_factor * contact_radius;
        force_normal -= damping_strength * v_rel;
#endif // GPU_WALL_USE_VISCOELASTIC_DAMPING

#ifdef GPU_WALL_USE_CONSTANT_DAMPING
        double v_rel = dot(make_double3(velocities, particle_id, number_of_particles), gpu_wall_constants.wall_normals[wall_id]);
        if(wall_id == GPU_WALL_TOP)
            v_rel += gpu_wall_constants.top_wall_speed;

        force_normal -= gpu_wall_constants.osc_damping_factor * v_rel;

#endif // GPU_WALL_USE_CONSTANT_DAMPING


        force += force_normal * n2; // reminder: n2 is the normal of the wall


        //////////////////////////////////////////////////////////////////////////////
        // apply inelastic corrections
        //////////////////////////////////////////////////////////////////////////////


#ifdef GPU_WALL_INELASTIC_ROLLING
		{
			double3 displacement;
			displacement = gpu_wall_constants.particle_radius * (n1 + n2);


			double norm = dot(displacement, displacement);

			// special treatment of contact pointers if we are in the inelastic regime
			if (norm > gpu_wall_constants.crit_rolling_displacement_squared)
			{
				norm = sqrt(norm);

				// determine correction of contact pointers
				displacement *= (1.0 - gpu_wall_constants.crit_rolling_displacement / norm);


				norm = dot(displacement, displacement);

                if(norm > 0.0) // norm == 0 -> displacement == 0 -> division by 0 = nan -> crash
                {
                    // correct contact pointers
                    // calculate correction factor alpha (see Wada et al. 2007 appendix for details)
                    double temp = dot(n1, displacement);
                    double alpha = norm / (norm - temp * temp);

                    // store n1 for sliding correction
                    double3 temp_vec = n1;

                    n1 -= gpu_wall_constants.particle_radius_inv * alpha * displacement;

                    // normalize contact pointers
                    norm = rsqrt(n1.x*n1.x + n1.y*n1.y + n1.z*n1.z);
                    n1 *= norm;

                    // new contact pointer, write to global memory at end of function
                    n_initial.x = 2.0 * ((0.5 - rot_param1.y * rot_param1.y - rot_param1.z * rot_param1.z) * n1.x + (rot_param1.x * rot_param1.y + rot_param1.z * rot_param1.w) * n1.y + diff_of_products_simple(rot_param1.x, rot_param1.z, rot_param1.y, rot_param1.w) * n1.z);
                    n_initial.y = 2.0 * (diff_of_products_simple(rot_param1.x, rot_param1.y, rot_param1.z, rot_param1.w) * n1.x + (0.5 - rot_param1.x * rot_param1.x - rot_param1.z * rot_param1.z) * n1.y + (rot_param1.y * rot_param1.z + rot_param1.x * rot_param1.w) * n1.z);
                    n_initial.z = 2.0 * ((rot_param1.x * rot_param1.z + rot_param1.y * rot_param1.w) * n1.x + diff_of_products_simple(rot_param1.y, rot_param1.z, rot_param1.x, rot_param1.w) * n1.y + (0.5 - rot_param1.x * rot_param1.x - rot_param1.y * rot_param1.y) * n1.z);


                    double3 delta_n = n1 - temp_vec;
                    temp = dot(delta_n, n2);

                    temp_vec = gpu_wall_constants.particle_radius * (delta_n + temp * n2);

                    initial_contact[index] += temp_vec;
                    n1_initial[index] = n_initial;
                    n_1[index] = n1;
                }
			}
		}
#endif	// GPU_WALL_INELASTIC_ROLLING




#ifdef GPU_WALL_ROLLING
        // Note: k_r is actually k_r * particle_radius^2
        double3 temp_vec = gpu_wall_constants.rolling_modifier[wall_id] * gpu_wall_constants.k_r * cross(n1, n2);

        torque -= temp_vec;


#endif // GPU_WALL_ROLLING





		
#ifdef GPU_WALL_SLIDING // see void Simulation::updateParticleInteraction
		{
			double3 displacement;

            // contact_vector = pos_new - wall_contact_pos
            displacement = fma(n1, gpu_wall_constants.particle_radius, contact_vector);
                            //contact_vector + gpu_wall_constants.particle_radius * n1;

			double temp = dot(displacement, n2);

            //displacement -= temp * n2;
            displacement = fma(n2, -temp, displacement);

#ifdef GPU_WALL_INELASTIC_SLIDING

			double norm = dot(displacement, displacement);

			// check if we are in the inelastic regime
			if (norm > gpu_wall_constants.crit_sliding_displacement_squared)
			{
				norm = sqrt(norm);

				// determine correction of contact pointers
				initial_contact[index] += (1.0 - gpu_wall_constants.crit_sliding_displacement / norm) * displacement;

				displacement = gpu_wall_constants.crit_sliding_displacement / norm * displacement;

			}

#endif	// GPU_WALL_INELASTIC_SLIDING

			/////////////////////////////////////////////////////////////////////////////////////
			// calculate new forces and torques from sliding interaction
			/////////////////////////////////////////////////////////////////////////////////////

					// force
			double3 sliding_force = gpu_wall_constants.wall_sliding_modifier[wall_id] * gpu_wall_constants.k_s * displacement;

            force -= sliding_force;

            // TODO: sliding_force.y always = 0.0 ? remove then
            // we only want to track the forces on the top and bottom wall
            if (wall_id == GPU_WALL_TOP)
            {
                wall_top_force[particle_id] += sliding_force.y;

            }
            if(wall_id == GPU_WALL_BOTTOM)
            {
                wall_bot_force[particle_id] += sliding_force.y;
            }




			// torque
			torque -= gpu_wall_constants.wall_sliding_modifier[wall_id] * gpu_wall_constants.k_s * gpu_wall_constants.particle_radius * cross(n1, displacement);

		}
#endif	// GPU_WALL_SLIDING

    }

    if(update_force_and_torque)
    {

        add_double3(forces_new, particle_id, number_of_particles, force);
        add_double3(torques_new, particle_id, number_of_particles, torque);

    }

    return;
}













#endif // __CUDA_WALL_KERNELS_H__
