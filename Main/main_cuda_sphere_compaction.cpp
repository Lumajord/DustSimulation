#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "SimulationCuda.h"
#include "SimulationLib.h"

extern double k_s;
extern double mass;

int main(int argc, char **argv)
{
	/////////////////////////////////////////////////////////////////////////////////////////
	// get params
	/////////////////////////////////////////////////////////////////////////////////////////

    double timestep;
    double destination_filling_factor;
    int material_file_index = -1;

    double agg_mass = 0.0;
    double agg_size = 0.0;
    int seed_agg;
    int GPU_id = -1;

    if(argc == 7)
	{
        destination_filling_factor = atof(argv[1]);
        agg_mass = atof(argv[2]);
        agg_size = atof(argv[3]);
        material_file_index = 4;
        seed_agg = atoi(argv[5]);
        GPU_id = atoi(argv[6]);
	}
	else
	{
        printf("Wrong number of arguments!  %d || 7 Use:\n\
               -destination_filling_factor\
               -agg_mass\
               -agg_size\
               -material_filename\
               -seed_agg\
               - GPU_ID\n", argc);
		return EXIT_SUCCESS;
	}


    if(GPU_id >= 0)
        cudaSetDevice(GPU_id);

    /////////////////////////////////////////////////////////////////////////////////////
    // Files
    ////////////////////////////////////////////////////////////////////////////////////


    char* shpere1_file = new char[1024];
    char* shpere2_file = new char[1024];

#ifdef INIT
    sprintf(shpere1_file, "../data/init/RBDsphere_");
    sprintf(shpere2_file, "../data/init/CRBDsphere_");
#else
    sprintf(shpere2_file, "../data/init/RCRBDsphere_");
    sprintf(shpere1_file, "../data/init/CRBDsphere_");
#endif

    char buf[256];


    sprintf(buf, "_size_%d_%d_%d.dat", int(agg_mass*100), int(agg_size*1.e4), seed_agg);


    size_t len = strlen(argv[material_file_index]);
    strncat(shpere1_file, argv[material_file_index]+3, len-7);
    strncat(shpere2_file, argv[material_file_index]+3, len-7);

    strcat(shpere1_file, buf);
    strcat(shpere2_file, buf);


#define INIT
#ifdef INIT
	/////////////////////////////////////////////////////////////////////////////////////////
	// load sim
	/////////////////////////////////////////////////////////////////////////////////////////

    printf("Loading simulation data...\n");

	SimulationCuda sim;

	ErrorCode error_code;
    error_code = sim.loadMaterial(argv[material_file_index]);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s    %s\n", argv[material_file_index], message, argv[material_file_index]);
		return EXIT_SUCCESS;
	}


    printf("Generating agglomerate...\n");
    SimLib::generateRBDSphere(&sim, argv[material_file_index], agg_mass*1.e-6, seed_agg);

    printf("Saving to file...\n");
    sim.saveToFile(shpere1_file);

    // get timestep size
#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
    double sliding_osc_period = 2.0 * M_PI * sqrt(mass / (k_s * equilibrium_contact_radius) );
#else
    double sliding_osc_period = 2.0 * M_PI * sqrt(mass / k_s);
#endif

    timestep = sliding_osc_period / 30.0;

    printf("Timestep = %e\n", timestep);


        SimLib::centerCMS(&sim);

        double sphere_radius = 0.0;
        double sphere_gradius = 0.0;

        SimLib::getSize(sim, &sphere_gradius, &sphere_radius);

        sim.m_sphere_radius = sphere_radius + particle_radius;

	/////////////////////////////////////////////////////////////////////////////////////////
	// setup cuda
	/////////////////////////////////////////////////////////////////////////////////////////

	printf("Setting up CUDA...\n");
    error_code = sim.initCuda(GPU_id);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nWhile setting up CUDA the following error occurred:\n%s\n", message);
		return EXIT_SUCCESS;
	}

	error_code = sim.toggleGPUMode(true);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to enable GPU mode!\n%s\n", message);
		return EXIT_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	// start simulation
	/////////////////////////////////////////////////////////////////////////////////////////

    error_code = sim.startSimulation(sim.m_sphere_radius/sim.m_sphere_compaction_speed, timestep);

	if(error_code != EC_OK)
	{
		char message[200];
		sim.getErrorMessage(error_code, message);
		printf("ERROR:\nFailed to start simulation!\n%s\n", message);
		return EXIT_SUCCESS;
        }


	printf("Running simulation with %i particles on GPU...\n", sim.number_of_particles);
            double inc = 10.0;
            double dest_radius = pow(sim.number_of_particles*particle_radius*particle_radius*particle_radius / destination_filling_factor, 1.0/3.0);
            double start_radius = sim.m_sphere_radius;
            double d_rad = (start_radius - dest_radius)/inc;


            double compress_later = (start_radius - dest_radius)*0.0;

            dest_radius += compress_later;



            double speed1 = 60.0; // in cm
            double speed2 = 60.0;

            const double damping_factor = 0.90;

            sim.m_sphere_compaction_speed = speed1; // in cm




            double time = ((start_radius - dest_radius)/(0.5*(speed1 + speed2)));
            time = time/280.1;


            printf("current filling factor = %e destination filling factor = %e number of timesteps: %d\n", sim.number_of_particles*particle_radius*particle_radius*particle_radius/(sim.m_sphere_radius*sim.m_sphere_radius*sim.m_sphere_radius), destination_filling_factor, int((start_radius - dest_radius)/(sim.m_sphere_compaction_speed * sim.timestep)));


            start_radius -= d_rad;
            double percent = 0.0;

            double switch_time = time;
            sim.setDampingFactor(0.0);
	while(!sim.stop_simulation)
        {

            if(sim.current_time > switch_time)
            {



                sphere_radius = sim.m_sphere_radius;

                sim.copySimDataFromGPU();
                sim.cleanUpCuda();
                SimLib::disturbParticles(&sim, 0.01*particle_radius);
                error_code = sim.initCuda(GPU_id);
                error_code = sim.toggleGPUMode(true);

                sim.m_sphere_radius = sphere_radius + 0.01*particle_radius;




                // damp away energy
                sim.setDampingFactor(damping_factor);
                for(int i = 0; i < 30; ++i)
                {
                    sim.m_sphere_compaction_speed *= damping_factor;
                    sim.update();
                }




                // accelerate back
                sim.setDampingFactor(0.999);
                sim.m_sphere_compaction_speed = speed1;
                switch_time = sim.current_time + time;



                // alternate compaction speed

                speed1 = speed2;
                speed2 = sim.m_sphere_compaction_speed;
                sim.m_sphere_compaction_speed = speed1;
                //printf("speed = %e\n", sim.m_sphere_compaction_speed);



            }

            double fill_factor = sim.number_of_particles*particle_radius*particle_radius*particle_radius/(sim.m_sphere_radius*sim.m_sphere_radius*sim.m_sphere_radius);

            if(sim.m_sphere_radius < start_radius)
            {
                percent+= 100.0/inc;
                start_radius -= d_rad;
                printf("Compacting %.2f percent done    coord_n = %.3f\n", percent, 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
            }

            if(sim.m_sphere_radius < dest_radius)
            {
                sim.stop_simulation = true;
            }

                sim.update();
        }


    sim.stop_simulation = false;
    sim.m_sphere_compaction_speed = 0.0;


    sim.setDampingFactor(0.85);

    for(int i = 0; i < 500; ++i)
        sim.update(); // damp away energies


    sim.copySimDataFromGPU();
    sim.saveToFile(shpere1_file);

    sim.copySimDataFromGPU();
    SimLib::printEnergy(&sim, NULL);
    sim.cleanUpCuda();
    //SimLib::disturbParticles(&sim, 0.01*particle_radius);
    SimLib::resetContacts(&sim);
    error_code = sim.initCuda(GPU_id);
    error_code = sim.toggleGPUMode(true);
    sim.setDampingFactor(0.90);



    printf("\n\n\n");

    for(int i = 0; i < 2000; ++i)
    {
        if(i % 500 == 0)
        {

            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            sim.cleanUpCuda();
            //SimLib::disturbParticles(&sim, 0.02*particle_radius);
            SimLib::resetContacts(&sim);
            error_code = sim.initCuda(GPU_id);
            error_code = sim.toggleGPUMode(true);
            sim.setDampingFactor(0.90);

            printf("relaxing %.2f\ percent done coord_number = %.2f\n", 100.0*double(i) / double(2000), 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }
        sim.update(); // relax
    }

    sim.setDampingFactor(0.99);

    printf("\n\n\n");

    /*
    for(int i = 0; i < 2000; ++i)
    {
        if(i % 500 == 0)
        {
            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            sim.cleanUpCuda();
            SimLib::resetContacts(&sim);
            error_code = sim.initCuda(GPU_id);
            error_code = sim.toggleGPUMode(true);
            sim.setDampingFactor(0.99);

            printf("relaxing %.2f\ percent done coord_number = %.2f\n", 100.0*double(i) / double(2000), 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }
        sim.update(); // relax
    }
    */

    /*
    printf("\n\n\n");

    dest_radius -= compress_later;
    sim.setDampingFactor(0.95);
    sim.m_sphere_compaction_speed = 100; // in cm
    while(!sim.stop_simulation)
    {
        double fill_factor = sim.number_of_particles*particle_radius*particle_radius*particle_radius/(sim.m_sphere_radius*sim.m_sphere_radius*sim.m_sphere_radius);

        if(sim.m_sphere_radius < start_radius)
        {
            percent+= 100.0/inc;
            start_radius -= d_rad;
            printf("Compacting %.2f percent done    coord_n = %.3f\n", percent, 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }

        if(destination_filling_factor < fill_factor)
        {
            sim.stop_simulation = true;
        }

            sim.update();
    }
    sim.stop_simulation = false;
    sim.m_sphere_compaction_speed = 0.0;


    sim.setDampingFactor(0.85);

    for(int i = 0; i < 500; ++i)
        sim.update(); // damp away energies
    */
    printf("\n\n\n");


    sim.setDampingFactor(0.99);
    for(int i = 0; i < 15000; ++i)
    {
        if(i % 500 == 0)
        {

            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            sim.cleanUpCuda();
            SimLib::resetContacts(&sim);
            error_code = sim.initCuda(GPU_id);
            error_code = sim.toggleGPUMode(true);
            sim.setDampingFactor(0.99);

            printf("relaxing %.2f\ percent done coord_number = %.2f\n", 100.0*double(i) / double(15000), 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }
        sim.update(); // relax
    }

#endif





#ifdef RELAX_MORE

    // more relaxing
    ErrorCode error_code;
    SimulationCuda sim;
    error_code = sim.loadMaterial(argv[material_file_index]);


    if(error_code != EC_OK)
    {
        char message[200];
        sim.getErrorMessage(error_code, message);
        printf("ERROR:\nWhile trying to load material from %s, the following error occurred:\n%s    %s\n", argv[material_file_index], message, argv[material_file_index]);
        return EXIT_SUCCESS;
    }




    sim.loadFromFile(shpere1_file);

    // get timestep size
#ifdef USE_VARYING_POTENTIAL_COEFFICIENTS
    double sliding_osc_period = 2.0 * M_PI * sqrt(mass / (k_s * equilibrium_contact_radius) );
#else
    double sliding_osc_period = 2.0 * M_PI * sqrt(mass / k_s);
#endif

    timestep = sliding_osc_period / 50.0;

    printf("Timestep = %e\n", timestep);

    SimLib::centerCMS(&sim);

    double sphere_radius = 0.0;
    double sphere_gradius = 0.0;

    SimLib::getSize(sim, &sphere_gradius, &sphere_radius);

    sim.m_sphere_radius = sphere_radius + particle_radius;
    sim.m_sphere_compaction_speed = 0.0;

    error_code = sim.startSimulation(1.0, timestep);
    error_code = sim.initCuda(GPU_id);
    error_code = sim.toggleGPUMode(true);





    sim.setDampingFactor(0.99);
    for(int i = 0; i < 6000; ++i)
    {
        if(i % 500 == 0)
        {

            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            sim.cleanUpCuda();
            SimLib::resetContacts(&sim);
            error_code = sim.initCuda(GPU_id);
            error_code = sim.toggleGPUMode(true);
            sim.setDampingFactor(0.99);

            printf("relaxing %.2f\ percent done coord_number = %.2f\n", 100.0*double(i) / double(4000), 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }
        sim.update(); // relax
    }
#endif



    //sim.m_sphere_radius = 0.0;
    sim.setDampingFactor(0.95);


    printf("\n\n\n");

    for(int i = 0; i < 1000; ++i)
    {
        if(i % 250 == 0)
        {

            sim.copySimDataFromGPU();
            SimLib::printEnergy(&sim, NULL);
            printf("relaxing %.2f\ percent done coord_number = %.2f\n", 100.0*double(i) / double(1000), 2.0*(double)sim.getGPUNumContacts()/(double)sim.number_of_particles);
        }
        sim.update(); // relax
    }


    if(sim.get_use_gpu())
        sim.copySimDataFromGPU();

    sim.saveToFile(shpere2_file);

    printf(" ... finished\n");

    return EXIT_SUCCESS;
}
