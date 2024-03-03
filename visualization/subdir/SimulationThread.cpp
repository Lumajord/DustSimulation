
#include "MainWindow.h"
#include "SimulationThread.h"
#include "SimulationLib.h"
#include "SimulationVisualization.h"

extern DisplayInfo sim_info;
extern double F_c;
extern double ENERGY_UNIT;

extern DisplayMode vis_display_mode;

extern Qt::CheckState display_text_info;
extern Qt::CheckState display_pressure;
extern Qt::CheckState display_force;
extern Qt::CheckState display_depth;
extern Qt::CheckState track_fragments;

extern QWaitCondition wait_for_screenshot;
extern QMutex mutex;

SimulationThread::SimulationThread()
{
	last_visualization_update.start();
	pause = false;
	update_counter = 0;

	number_of_particles = 0;
	particle_positions = NULL;
	particle_colors = NULL;

	number_of_walls = 0;
	alpha_values = NULL;
	wall_positions = NULL;
	wall_normals = NULL;

	draw_gyration_radius = false;
}

SimulationThread::~SimulationThread()
{
	pause = true;
	wait();

	delete [] particle_positions;
	delete [] particle_colors;

	delete [] alpha_values;
	delete [] wall_positions;
	delete [] wall_normals;
}

void SimulationThread::run()
{
	while(sim.current_time < sim.end_time && !pause && !sim.stop_simulation)
	{

		sim.update();
		++update_counter;

		if(display_text_info)
			updateInfo();

		// check if visualization needs to be updated
		if(last_visualization_update.elapsed() > 30)
		{
			updateDrawnObjects();

			last_visualization_update.restart();
        }

		// check if screenshot should be taken
		if(video_mode && update_counter == screenshot_interval)
		{
			mutex.lock();
			emit signalTakeScreenshot();
			wait_for_screenshot.wait(&mutex);
			mutex.unlock();

			update_counter = 0;
		}
		/*else
		{
			// slow down if sim is running too fast
			if(update_counter > updates_per_frame)
			{
				while(last_sim_update.elapsed() < desired_update_time)
				{
				}

				update_counter = 0;

				// restart timer
				last_sim_update.restart();
			}
		}*/
	}
	
    updateDrawnObjects();
}

void SimulationThread::updateDrawnObjects()
{
	can_draw = false;

	// update particle positions
	sim.storeParticlePos(particle_positions);

	// update particle colors
	switch(vis_display_mode)
	{
	case DISPLAY_DISSIPATED_ENERGY:
		sim.storeDissipatedContactEnergy(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_PARTICLE_SPEED:
		sim.storeParticleSpeed(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_PARTICLE_CONTACTS:
		sim.storeColorParticleContacts(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_WALL_CONTACTS:
		sim.storeColorWallContacts(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_WALL_ANGLE:
		sim.storeWallAngle(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_DENSITY:
		sim.storeDensity(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_NOTHING:
		sim.storeDepth(particle_colors, (bool)display_depth);
		break;
	case DISPLAY_FRAGMENTS:
		if(track_fragments)
			sim.fragment_id_offset = SimLib::detectFragments(sim, sim.fragment_ids);

		sim.storeFragments(particle_colors, (bool)display_depth);
		break;		
	case DISPLAY_DISLOCATION:
		sim.storeDislocation(particle_colors, (bool)display_depth);
		break;
    case DISPLAY_50_50:
        sim.store50(particle_colors, (bool)display_depth);
        break;
	}

	// update walls
	sim.storeWalls(wall_positions, wall_normals, alpha_values);

	// update gyration radius
	updateGyrationRadius();

	can_draw = true;

	emit signalUpdateGLView();
}

void SimulationThread::resizeArrays(int new_number_of_particles, int new_number_of_walls)
{
	sim_info.particles = new_number_of_particles;
	sim_info.broken_contacts = 0;
	sim_info.created_contacts = 0;
	sim_info.time = 0;
	sim_info.wall_speed = 0;
	sim_info.collision_speed = 0;
	sim_info.filling_factor = 0;
	sim_info.pressure = -1;
	sim_info.force = -1;
	sim_info.coordination_number = 0;

	number_of_particles = new_number_of_particles;
	number_of_walls = new_number_of_walls;

	delete [] particle_positions;
	delete [] particle_colors;
	delete [] alpha_values;
	delete [] wall_positions;
	delete [] wall_normals;

	particle_positions = new double[3 * number_of_particles];
	particle_colors = new float[3 * number_of_particles];

	alpha_values = new float[number_of_walls];
	wall_positions = new double[24*number_of_walls];
	wall_normals = new double[3*number_of_walls];

	emit signalNumberOfParticlesChanged(number_of_particles);
	emit signalNumberOfWallsChanged(number_of_walls);
}

void SimulationThread::loadFromFile(char* filename)
{
	pause = true;
	wait();

	ErrorCode error_code = sim.loadFromFile(filename);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	// update particle data for visualization
	resizeArrays(sim.number_of_particles, sim.number_of_walls);

	// determine wall speed
	if(sim.sim_info.info_storage[2] != 0)
		sim_info.collision_speed = sim.sim_info.info_storage[2];

	updateInfo();
	updateDrawnObjects();

    pause = false;
}


void SimulationThread::importFromFile(char* filename)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::importParticlePositions(&sim, filename);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	// update particle data for visualization
	resizeArrays(sim.number_of_particles, (int)sim.walls.size());

	updateInfo();
	updateDrawnObjects();

	pause = false;
}


void SimulationThread::saveToFile(char* filename)
{
	if(!sim.saveToFile(filename))
		emit signalSaveFailed();
}

void SimulationThread::savePositionsToFile(char* filename)
{
	SimLib::printPositions(sim, filename, false);
}

void SimulationThread::loadMaterial(char* filename)
{
	pause = true;
	wait();

	ErrorCode error_code = sim.loadMaterial(filename);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	emit signalMaterialLoaded();

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initChain(int number_of_particles, int number_of_impact_particles, int target_id, double angular_irregularity, double impact_speed)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initChain(&sim, number_of_particles, number_of_impact_particles, target_id, angular_irregularity, impact_speed);

	if( error_code != EC_OK )
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initChainBox(int number_of_particles, double filling_factor, double angular_irregularity)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initChainBox(&sim, number_of_particles, filling_factor, angular_irregularity);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::impactChainOnAgglomerate(int number_of_impact_particles, double angular_irregularity, double impact_speed)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initImpactChainOnAgglomerate(&sim, number_of_impact_particles, angular_irregularity, impact_speed, true);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initCluster(int number_of_particles, int number_of_impact_particles, double filling_factor, double angular_irregularity, double impact_speed)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initCluster(&sim, number_of_particles, number_of_impact_particles, filling_factor, angular_irregularity, impact_speed);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initBAMAggregate(unsigned int number_of_particles, double migration_probability1, double migration_probability2, BAMSelectionMethod bam_selection_method)
{
	pause = true;
    wait();

	sim.cleanUp();
	ErrorCode error_code = SimLib::initBAMAggregate(&sim, NULL, number_of_particles, migration_probability1, migration_probability2, bam_selection_method);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initFractalAggregate(unsigned int number_of_particles, double migration_probability1, double migration_probability2, double chain_ratio, unsigned int min_chain_length, unsigned int max_chain_length)
{
	pause = true;
    wait();

	sim.cleanUp();
	ErrorCode error_code = SimLib::initFractalAggregate(&sim, number_of_particles, migration_probability1, migration_probability2, chain_ratio, min_chain_length, max_chain_length);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initTreeAggregate(unsigned int number_of_particles, double migration_rate1, double migration_rate2, double contacts0, double contacts1, double contacts2, double contacts3, double contacts4, double contacts5)
{
	pause = true;
    wait();

	sim.cleanUp();

	double contact_distribution[6];
	contact_distribution[0] = contacts0;
	contact_distribution[1] = contacts1;
	contact_distribution[2] = contacts2;
	contact_distribution[3] = contacts3;
	contact_distribution[4] = contacts4;
	contact_distribution[5] = contacts5;

	ErrorCode error_code = SimLib::initTreeAggregate(&sim, number_of_particles, contact_distribution, migration_rate1, migration_rate2);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initSimpleCubicLattice(unsigned int x_particles, unsigned int y_particles, unsigned int z_particles, double filling_factor)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::buildSimpleCubicGrid(&sim, (double)x_particles * 2.0 * particle_radius, (double)y_particles * 2.0 * particle_radius, (double)z_particles * 2.0 * particle_radius, filling_factor);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initHexagonalLattice(unsigned int x_particles, unsigned int y_particles, unsigned int z_particles, double filling_factor)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::buildHexagonalGrid(&sim, x_particles, y_particles, z_particles, filling_factor);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initRBDCake(unsigned int number_of_particles, double x_size, double y_size, double top_slice_factor, double bottom_slice_factor)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initRandomBallisticDeposition(&sim, number_of_particles, x_size, y_size);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	// slice if necessary
	if(top_slice_factor > 0  || bottom_slice_factor > 0)
	{
		vec3 lower_pos, upper_pos;
		sim.getEnclosingBox(&lower_pos, &upper_pos);

		upper_pos[1] -= (upper_pos[1] - lower_pos[1]) * top_slice_factor;
		lower_pos[1] += (upper_pos[1] - lower_pos[1]) * bottom_slice_factor;
		SimLib::sliceBox(&sim, lower_pos, upper_pos);
	}

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initBAMCake(unsigned int number_of_particles, double x_size, double y_size, double migration_probability1, double migration_probability2, double top_slice_factor, double bottom_slice_factor, BAMSelectionMethod bam_selection_method)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initBAMCake(&sim, number_of_particles, x_size, y_size, migration_probability1, migration_probability2);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	// slice if necessary
	if(top_slice_factor > 0  || bottom_slice_factor > 0)
	{
		vec3 lower_pos, upper_pos;
		sim.getEnclosingBox(&lower_pos, &upper_pos);

		upper_pos[1] -= (upper_pos[1] - lower_pos[1]) * top_slice_factor;
		lower_pos[1] += (upper_pos[1] - lower_pos[1]) * bottom_slice_factor;
		SimLib::sliceBox(&sim, lower_pos, upper_pos);
	}

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initNeedleDeposition(unsigned int number_of_particles)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::initNeedleDeposition(&sim, number_of_particles);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::collideAgglomerateWithWall(const char* filename, double impact_speed, int impact_angle, double impact_distance, bool random_orientation)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::collideAgglomerateWithWall(&sim, filename, impact_speed, impact_angle, impact_distance, random_orientation);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	sim_info.collision_speed = impact_speed;

	pause = false;
}

void SimulationThread::initCompression(const char *filename, bool random_orientation, bool side_walls, bool moving_bottom_wall, double wall_speed, double stop_filling_factor, double side_wall_mod)
{
	pause = true;
	wait();

	if(random_orientation)
		SimLib::rotateSimRandomly(&sim);

	ErrorCode error_code = SimLib::initCompressionBox(&sim, filename, side_walls, moving_bottom_wall, wall_speed, stop_filling_factor, 1.0, side_wall_mod, side_wall_mod);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	if(display_pressure == Qt::Checked)
		sim_info.pressure = 0.0;

	if(display_force == Qt::Checked)
		sim_info.force = 0.0;

	pause = false;
}

void SimulationThread::initCompressionRelaxation(const char* filename, bool random_orientation, double wall_speed, double stop_filling_factor, double stop_dissipation_factor)
{
	pause = true;
	wait();

	if(random_orientation)
		SimLib::rotateSimRandomly(&sim);

	ErrorCode error_code = SimLib::initCompressionRelaxationBox(&sim, filename, wall_speed, stop_filling_factor, stop_dissipation_factor);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	if(display_pressure == Qt::Checked)
		sim_info.pressure = 0.0;

	if(display_force == Qt::Checked)
		sim_info.force = 0.0;

	pause = false;
}

void SimulationThread::initShockwave(const char* filename, bool random_orientation, double wall_speed, double stop_pressure, double side_wall_mod)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::initShockwaveBox(&sim, filename, wall_speed, 2.0 * particle_radius, stop_pressure, side_wall_mod, side_wall_mod);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	if(display_pressure == Qt::Checked)
		sim_info.pressure = 0.0;

	if(display_force == Qt::Checked)
		sim_info.force = 0.0;

	pause = false;
}

void SimulationThread::initDynamicCompression(const char* filename, bool random_orientation, double wall_speed, double dynamic_pressure, double side_wall_mod)
{
	pause = true;
	wait();

	if(random_orientation)
		 SimLib::rotateSimRandomly(&sim);

	ErrorCode error_code = SimLib::initDynamicCompressionBox(&sim, filename, wall_speed, dynamic_pressure, 1.0, side_wall_mod, side_wall_mod);
	
	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	if(display_pressure == Qt::Checked)
		sim_info.pressure = 0.0;

	if(display_force == Qt::Checked)
		sim_info.force = 0.0;

	pause = false;
}

void SimulationThread::initOpenBox(const char* filename, bool random_orientation, bool side_walls, double side_wall_mod)
{
	pause = true;
	wait();

	if(random_orientation)
		SimLib::rotateSimRandomly(&sim);

	ErrorCode error_code = SimLib::initOpenBox(&sim, filename, side_walls, 1.0, side_wall_mod, side_wall_mod);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	if(display_pressure == Qt::Checked)
		sim_info.pressure = 0.0;

	if(display_force == Qt::Checked)
		sim_info.force = 0.0;

	pause = false;
}

void SimulationThread::initPullStrengthTest(double pull_speed)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::initPullStrengthTest(&sim, pull_speed);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::initShearStrengthTest(double pull_speed)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::initShearStrengthTest(&sim, pull_speed);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::addBAMParticles(unsigned int number_of_particles, double migration_probability1, double migration_probability2, const char* filename, BAMSelectionMethod bam_selection_method)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::initBAMAggregate(&sim, filename, number_of_particles, migration_probability1, migration_probability2, bam_selection_method);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::addFractalChains(const char* filename, unsigned int number_of_chains, unsigned int min_chain_length, unsigned int max_chain_length)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::addFractalChainsToAggregate(&sim, filename, number_of_chains, min_chain_length, max_chain_length, 0);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::rotateAggregate(double axis_x, double axis_y, double axis_z, double angle)
{
	pause = true;
	wait();

	vec3 axis;
	axis[0] = axis_x;
	axis[1] = axis_y;
	axis[2] = axis_z;
	normalize(&axis);

	SimLib::centerCMS(&sim);
	SimLib::rotateParticles(&sim, axis, M_PI / 180.0 * angle, 0, sim.number_of_particles-1);

	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::sliceBox(const char* filename, double x_size, double y_size, double z_size)
{
	pause = true;
	wait();

	// load file if specified
	if(filename != NULL)
		sim.loadFromFile(filename);

	SimLib::centerCMS(&sim);

	vec3 pos_lower, pos_upper;
	pos_lower[0] = -0.5*x_size;
	pos_lower[1] = -0.5*y_size;
	pos_lower[2] = -0.5*z_size;
	pos_upper[0] = 0.5*x_size;
	pos_upper[1] = 0.5*y_size;
	pos_upper[2] = 0.5*z_size;

	SimLib::sliceBox(&sim, pos_lower, pos_upper);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::sliceSphere(const char* filename, double slice_radius, bool random_orientation)
{
	pause = true;
	wait();

	// load file if specified
	if(filename != NULL)
		sim.loadFromFile(filename);

	if(random_orientation)
		SimLib::rotateSimRandomly(&sim);

	vec3 cms;
	SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);
	SimLib::sliceSphere(&sim, cms, slice_radius);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
    updateInfo();
    sim.sim_info.info_storage[5] = sim_info.coordination_number;
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::sliceCylinder(const char* filename, double slice_radius)
{
	pause = true;
	wait();

	// load file if specified
	if(filename != NULL)
		sim.loadFromFile(filename);

	vec3 cms;
	SimLib::getCenterOfMass(&cms, sim, 0, sim.number_of_particles-1);

	SimLib::sliceCylinder(&sim, cms, slice_radius);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::sliceTop(const char* filename, double top_slice_factor)
{
	pause = true;
	wait();

	// load file if specified
	if(filename != NULL)
		sim.loadFromFile(filename);

	SimLib::sliceTop(&sim, top_slice_factor);
	SimLib::centerCMS(&sim);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::sliceBottom(const char* filename, double bottom_slice_factor)
{
	pause = true;
	wait();

	// load file if specified
	if(filename != NULL)
		sim.loadFromFile(filename);

	SimLib::sliceBottom(&sim, bottom_slice_factor);
	SimLib::centerCMS(&sim);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::duplicateAggregate(int x_duplications, int y_duplications, int z_duplications, bool mirror,  bool random_orientation)
{
	pause = true;
	wait();

	ErrorCode error_code = SimLib::duplicateSlices(&sim, x_duplications, y_duplications, z_duplications, mirror, random_orientation);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::collideAgglomerates(const char *filename1, const char *filename2, bool random_orientation1, bool random_orientation2, double impact_speed, double impact_parameter, double impact_distance, bool minimize_distance)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::collideTwoAgglomerates(&sim, filename1, filename2, impact_speed, impact_parameter, impact_distance, random_orientation1, random_orientation2, minimize_distance, true);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	sim_info.collision_speed = impact_speed;

	pause = false;
}

void SimulationThread::hitAndStick(const char *filename1, const char *filename2, bool random_orientation1, bool random_orientation2, double impact_parameter)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::hitAndStick(&sim, filename1, filename2, impact_parameter, random_orientation1, random_orientation2);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::impactMultiProjectiles(const char* target_filename, const char* projectile_filename, unsigned int projectile_count, double impact_velocity, double impact_parameter, bool use_random_impact_parameter, double initial_distance, bool smart_distance)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::impactMultiProjectiles(&sim, target_filename, true, projectile_filename, projectile_count, 0, impact_velocity, impact_parameter, use_random_impact_parameter, initial_distance, initial_distance, smart_distance);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::impactMultiPCAProjectiles(const char* target_filename, unsigned int projectile_count, unsigned int projectile_size, double impact_velocity, double impact_parameter, bool use_random_impact_parameter, double initial_distance, bool smart_distance)
{
	pause = true;
    wait();

	ErrorCode error_code = SimLib::impactMultiProjectiles(&sim, target_filename, true, NULL, projectile_count, projectile_size, impact_velocity, impact_parameter, use_random_impact_parameter, initial_distance, initial_distance, smart_distance);

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::filterFragments()
{
	pause = true;
    wait();

	SimLib::filterFragments(&sim);

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::setWallProperties(int wall_id, double wall_speed, double dynamic_pressure, double alpha, double compression_modifier, double rolling_modifier, double sliding_modifier)
{
	pause = true;
    wait();

	if(wall_id < sim.number_of_walls)
	{
		sim.walls[wall_id].velocity[0] = wall_speed * sim.walls[wall_id].normal[0];
		sim.walls[wall_id].velocity[1] = wall_speed * sim.walls[wall_id].normal[1];
		sim.walls[wall_id].velocity[2] = wall_speed * sim.walls[wall_id].normal[2];
		sim.walls[wall_id].mass = dynamic_pressure;
		sim.walls[wall_id].alpha = alpha;
		sim.walls[wall_id].compression_modifier = compression_modifier;
		sim.walls[wall_id].rolling_modifier = rolling_modifier;
		sim.walls[wall_id].sliding_modifier = sliding_modifier;
	}
	else
	{
		for(int w = 0; w < sim.number_of_walls; ++w)
		{
			sim.walls[w].velocity[0] = wall_speed * sim.walls[wall_id].normal[0];
			sim.walls[w].velocity[1] = wall_speed * sim.walls[wall_id].normal[1];
			sim.walls[w].velocity[2] = wall_speed * sim.walls[wall_id].normal[2];
			sim.walls[w].mass = dynamic_pressure;
			sim.walls[w].alpha = alpha;
			sim.walls[w].compression_modifier = compression_modifier;
			sim.walls[w].rolling_modifier = rolling_modifier;
			sim.walls[w].sliding_modifier = sliding_modifier;
		}
	}

	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::removeWalls()
{
	pause = true;
    wait();

	sim.removeWalls();

	resizeArrays(sim.number_of_particles, (int)sim.walls.size());
	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::resetContacts()
{
	pause = true;
    wait();

	SimLib::resetContacts(&sim);

	updateInfo();
	updateDrawnObjects();

	pause = false;
}

void SimulationThread::disturbParticles(double max_dist)
{
	pause = true;
    wait();

	SimLib::disturbParticles(&sim, max_dist);

	updateInfo();
	updateDrawnObjects();

	pause = false;
}

#ifdef DRAW_CONTACT_POINTERS
void SimulationThread::StoreContactPointers(multimap<int,vec3> *dest)
{
	sim.storeContactPointers(dest);

	emit ContactPointerUpdateDone();
}
#endif

void SimulationThread::startSimulation(double duration, double update_interval, int print_energy_interval, int print_positions_interval, int pressure_averaging_steps, int screenshot_interval, bool video_mode, bool use_gravity)
{
	this->video_mode = video_mode;
	this->screenshot_interval = screenshot_interval; 

	sim.initPressureAveraging(pressure_averaging_steps);

	update_counter = 0;

#ifdef _DEBUG
	ErrorCode error_code = sim.startSimulation(duration, update_interval, print_energy_interval, print_positions_interval, "energies_debug.txt", "positions\\", use_gravity);
#else

#ifdef _WIN32
	ErrorCode error_code = sim.startSimulation(duration, update_interval, print_energy_interval, print_positions_interval, "energies.txt", "positions\\", use_gravity);
#else
	ErrorCode error_code = sim.startSimulation(duration, update_interval, print_energy_interval, print_positions_interval, "energies.txt", "positions/", use_gravity);
#endif

#endif

	if(error_code != EC_OK)
		emit signalInitFailed(error_code);

	last_visualization_update.restart();

	if(!isRunning()) 
        start(NormalPriority);
}

void SimulationThread::pauseSimulation()
{
	pause = true;
}

void SimulationThread::continueSimulation()
{
	pause = false;
	last_visualization_update.restart();

	if(!isRunning()) 
        start(NormalPriority);
}

void SimulationThread::moveToOrigin()
{
	pause = true;
    wait();

	SimLib::centerCMS(&sim);

	updateDrawnObjects();
	emit signalParticleUpdateDone(sim.number_of_particles);

	pause = false;
}

void SimulationThread::getCrossSection(double *cross_section, double *sigma_cross_section)
{
	pause = true;
    wait();

    //SimLib::getCrossSection(sim, 0.2 * particle_radius, 8, cross_section, sigma_cross_section);
    SimLib::getCrossSectionNoRotation(sim, 0.2*particle_radius, *cross_section);
    *sigma_cross_section = -1.0;

	pause = false;
}

void SimulationThread::getCMS(float *cms)
{
	cms[0] = 0;
	cms[1] = 0;
	cms[2] = 0;

    if(number_of_particles == 0)
            return;

	for(int p = 0; p < number_of_particles; ++p)
	{
		cms[0] += particle_positions[X_COORD_(p)];
		cms[1] += particle_positions[Y_COORD_(p)];
		cms[2] += particle_positions[Z_COORD_(p)];
	}

	cms[0] /= (float)number_of_particles;
	cms[1] /= (float)number_of_particles;
	cms[2] /= (float)number_of_particles;

    return;
}

float SimulationThread::getGyrationRadius(float *cms)
{
	// calc gyration radius
	float result = 0;
	float dist;

	for(int p = 0; p < number_of_particles; ++p)
	{
		dist = (particle_positions[X_COORD_(p)] - cms[0]) * (particle_positions[X_COORD_(p)] - cms[0]) 
			+ (particle_positions[Y_COORD_(p)] - cms[1]) * (particle_positions[Y_COORD_(p)] - cms[1]) 
			+ (particle_positions[Z_COORD_(p)] - cms[2]) * (particle_positions[Z_COORD_(p)] - cms[2]);
	
		result += dist; 
	}

	result /= (float)number_of_particles;

	return sqrt(result);
}

void SimulationThread::updateGyrationRadius()
{
	gyration_radii.clear();

	if(draw_gyration_radius)
	{
		GyrationRadiusInfo temp;

		// get CMS
		getCMS(temp.pos);

		temp.radius = getGyrationRadius(temp.pos);

		temp.color[0] = 1.0f;
		temp.color[1] = 0.0f;
		temp.color[2] = 0.0f;
		temp.color[3] = 0.3f;

		gyration_radii.push_back(temp);
	}
}

void SimulationThread::toggleGyrationRadius()
{
	if(draw_gyration_radius)
		draw_gyration_radius = false;
	else
		draw_gyration_radius = true;

	updateDrawnObjects();
}

void SimulationThread::updateInfo()
{
	if(sim.number_of_particles > 0)
		sim_info.coordination_number = 2.0 * (double)sim.number_of_contacts / (double)sim.number_of_particles;

	if(sim.current_time > 0)
		sim_info.time = sim.current_time;

	sim_info.created_contacts = sim.created_contacts;
    sim_info.broken_contacts = sim.broken_contacts;

	if(sim.box)
	{
		if(sim.sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
			sim_info.filling_factor = sim.getCrossSectionFillingFactor();
		else
            sim_info.filling_factor = sim.getBoxFillingFactor();
		
		if(display_pressure == Qt::Checked)
			sim_info.pressure = sim.getAvgWallPressure();

		if(display_force == Qt::Checked)
			sim_info.force = sim.getAvgWallForce();

		if(sim.box->top_wall_id >= 0)
			sim_info.wall_speed = norm(sim.walls[sim.box->top_wall_id].velocity);
	}
}

void SimulationThread::updateAgglomerateInfo(AgglomerateInfo *agg_info, bool only_inner_particles)
{
	if(sim.number_of_particles > 0)
	{
		SimLib::getSize(sim, &(agg_info->gyration_radius), &(agg_info->outer_radius));

		agg_info->filling_factor_box = sim.getBoxFillingFactor();
		agg_info->filling_factor_sphere = (double)sim.number_of_particles * particle_radius*particle_radius*particle_radius / (agg_info->outer_radius*agg_info->outer_radius*agg_info->outer_radius);

		std::vector<int> fragment_ids;
		std::vector<int> size_of_fragment;

		SimLib::detectFragments(sim, &fragment_ids, &size_of_fragment, NULL);
		agg_info->fragments = size_of_fragment.size();

		if(sim.box)
		{
			agg_info->box_height = sim.box->height;

			if(sim.sim_info.sim_type == SIM_TYPE_COMPRESSION_NO_SIDE_WALLS)
				agg_info->box_base = sim.getCrossSection();
			else
				agg_info->box_base = sim.box->base;
		}
		else
		{
			agg_info->box_height = 0.0;
			agg_info->box_base = 0.0;
		}

		SimLib::getContactHistogram(sim, agg_info->contact_histogram, only_inner_particles);
	}
	else
	{
		agg_info->filling_factor_box = 0.0;
		agg_info->filling_factor_sphere = 0.0;
		agg_info->fragments = 0;
		agg_info->gyration_radius = 0.0;
		agg_info->outer_radius = 0.0;
		agg_info->box_height = 0.0;
		agg_info->box_base = 0.0;
	}
}

void SimulationThread::printContactHistogram(const char *filename, bool only_inner_particles)
{
	SimLib::printContactHistogram(sim, filename, only_inner_particles);
}

void SimulationThread::printElasticCharge(const char *filename)
{
	FILE *file = fopen(filename, "w+");

	if(file)
	{
		fprintf(file, "# number of contacts    charged contacts    elastic_charge    kinetic_energy\n");
		SimLib::printElasticCharge(sim, file);

		fclose(file);
	}
}

void SimulationThread::printFragmentVelocities(const char *filename)
{
	SimLib::printFragmentVelocities(sim, filename);
}

void SimulationThread::printFillingFactorProfile(const char *filename)
{
	SimLib::printFillingFactorProfile(sim, filename, 2.0 * particle_radius);
}

void SimulationThread::printFillingFactorProfileSphere(const char *filename)
{
    SimLib::printFillingFactorProfileSphere(sim, filename, 2.0 * particle_radius);
}

void SimulationThread::printFractalDimension(const char *filename)
{
	SimLib::printFractalDimension(sim, filename, 4.0 * particle_radius);
}
