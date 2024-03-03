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

#include "MainWindow.h"

#include <QtGui>
#include <string>
#include <locale.h>

extern char material_name[200];
extern double particle_radius;
extern double density;
extern double surface_energy;
extern double nu;
extern double young_mod;
extern double crit_rolling_displacement;
extern double osc_damping_factor;
extern double T_vis;
extern double viscoelastic_damping_constant;
extern double rolling_modifier;
extern double sliding_modifier;
extern double twisting_modifier;
extern double gravity_modifier;
extern double wall_inertia;
extern double delta_0;

const char *CFG_FILE_VERSION = "CFG_FILE_VERSION_1_88";
#define GUI_VERSION "1.88"

// background/border color
float clear_color[3] = {1.0f, 1.0f, 1.0f};
float border_color[3] = {0.0f, 0.0f, 0.0f};
float wall_color[3] = {0.58f, 0.67f, 0.8f};
float text_color[3] = {0.0f, 0.0f, 0.0f};

// parameters for border shader
float border_distance_factor = 20.0f;
float border_min_depth = 6.0f;
Qt::CheckState draw_borders = Qt::Unchecked;
Qt::CheckState draw_particles = Qt::Checked;

// depth/density/speed visualization
int vis_neighbour_search_dist = 2;
int vis_neighbour_max_particles = 150;
int vis_neighbour_min_offset = 40;
int vis_neighbour_max_offset = 20;
double vis_density_max_filling_factor = 0.75;
double vis_dislocation_min_value = 2;
double vis_dislocation_max_value = 5;
int vis_velocity_averaging_steps = 10;
bool vis_density_consider_walls = true;
DisplayMode vis_display_mode = DISPLAY_NOTHING;

Qt::CheckState display_key = Qt::Checked;
Qt::CheckState display_text_info = Qt::Checked;
Qt::CheckState display_pressure = Qt::Checked;
Qt::CheckState display_force = Qt::Unchecked;
Qt::CheckState display_depth = Qt::Unchecked;
Qt::CheckState display_changed_contacts = Qt::Unchecked;
int font_size = 16;
Qt::CheckState track_fragments = Qt::Unchecked;

// settings for StartSimWidget
double sim_timestep = 0.1;
double sim_duration = 1000.0;
int sim_time_display_mode = 0;
int sim_averaging_steps = 10;
int print_energy_interval = 0;
int print_positions_interval = 0;
int take_screenshot_interval = 50;
Qt::CheckState follow_cms = Qt::Unchecked;
Qt::CheckState sim_use_gravity = Qt::Unchecked;
double sim_azimuthal_acceleration = 1e9;
Qt::CheckState sim_use_sim_azimuthal_acceleration = Qt::Unchecked;

// init chains widget
int init_chain_particles = 10;
int init_chain_impact_particles = 1;
int init_chain_target_particle= 5;
double init_chain_impact_speed = 100;
double init_chain_angular_irregularity = 0.0;
Qt::CheckState init_chain_impact_current_agglomerate = Qt::Unchecked;

// init aggregate widget
unsigned int init_aggregate_particles = 1000;
double init_aggregate_migration_probability1 = 0.0;
double init_aggregate_migration_probability2 = 0.0;
BAMSelectionMethod init_aggregate_bam_selection_method = BAM_SELECT_CLOSEST;
double init_aggregate_chain_ratio = 0.5;
unsigned int init_aggregate_min_chain_length = 2;
unsigned int init_aggregate_max_chain_length = 6;
double init_aggregate_contact_distribution[6] = {0.5, 0.25, 0.25, 0.0, 0.0, 0.0};

// init cake widget
unsigned int init_cake_x_particles = 10;
unsigned int init_cake_y_particles = 10;
unsigned int init_cake_z_particles = 10;
double init_cake_filling_factor = 0.5;
unsigned int init_cake_particles = 10000;
double init_cake_migration_probability1 = 0.0;
double init_cake_migration_probability2 = 0.0;
BAMSelectionMethod init_cake_bam_selection_method = BAM_SELECT_INTERIOR;
double init_cake_x_size = 30.0;
double init_cake_y_size = 30.0;
double init_cake_top_slice_factor = 0.0;
double init_cake_bottom_slice_factor = 0.0;

// agglomerates collision widget
double agglomerates_collision_impact_speed = 100.0;
double agglomerates_collision_impact_parameter = 0.0;
double agglomerates_collision_initial_separation = 2.0;
Qt::CheckState agglomerates_collision_minimize_distance = Qt::Checked;
Qt::CheckState agglomerates_collision_random_orientation1 = Qt::Checked;
Qt::CheckState agglomerates_collision_random_orientation2 = Qt::Checked;
Qt::CheckState agglomerates_collision_slice1 = Qt::Unchecked;
Qt::CheckState agglomerates_collision_slice2 = Qt::Unchecked;
unsigned int agglomerates_collision_projectile_count = 10;
unsigned int agglomerates_collision_projectile_size = 1;
Qt::CheckState agglomerates_collision_random_impact_parameter = Qt::Unchecked;
Qt::CheckState agglomerates_collision_use_PCA_projectiles = Qt::Unchecked;

// wall collision widget
double wall_collision_impact_speed = 2000.0;
int wall_collision_impact_angle = 0;
double wall_collision_impact_distance = 4.0;
Qt::CheckState wall_collision_random_orientation = Qt::Checked;

// box widget
double box_wall_speed = 100.0;
double box_pressure = 1000.0;
double box_stop_filling_factor = 0.6;
double box_stop_dissipation_factor = 0.01;
double box_stop_pressure = 1.0;
double box_dynamic_pressure = 1.0;
double box_side_wall_mod = 0.001;
int box_mode = 0;
Qt::CheckState box_random_orientation = Qt::Unchecked;
Qt::CheckState box_moving_bottom_wall = Qt::Unchecked;
Qt::CheckState box_side_walls = Qt::Checked;

// add particles widget
unsigned int add_particles_number_of_particles = 100;
double add_particles_migration_probability1 = 0;
double add_particles_migration_probability2 = 0;
BAMSelectionMethod add_particles_bam_selection_method = BAM_SELECT_CLOSEST;
unsigned int add_particles_number_of_chains = 10;
unsigned int add_particles_min_chain_length = 10;
unsigned int add_particles_max_chain_length = 20;

// rotation widget
double rotation_axis_x = 1.0;
double rotation_axis_y = 0.0;
double rotation_axis_z = 0.0;
double rotation_angle = 90.0;

// slice widget
Qt::CheckState slice_random_orientation = Qt::Checked;
double slice_sphere_radius = 1.0;
double slice_box_x_size = 10.0;
double slice_box_y_size = 10.0;
double slice_box_z_size = 10.0;
double slice_top_slice_factor = 0.1;

// duplicate widget
int duplicate_x_duplications = 1;
int duplicate_y_duplications = 1;
int duplicate_z_duplications = 1;
Qt::CheckState duplicate_random_orientation = Qt::Checked;
Qt::CheckState duplicate_slice_agglomerate = Qt::Checked;
Qt::CheckState duplicate_mirror = Qt::Unchecked;

// last browsed directory
QString file_path;

// properties of current sim
DisplayInfo sim_info;
AgglomerateInfo agg_info;

// if true only the interior particles are taken into account when generating a contact histogram
Qt::CheckState agg_info_only_inner_particles = Qt::Unchecked;

// thread synchronization
extern QWaitCondition wait_for_screenshot;

// additional description text
QString info_text = "";

MainWindow::MainWindow(SimulationThread *sim_thread)
{
    setlocale(LC_ALL, "C"); // for some reasong local was set to comma as digit separator, causing the files to be read incorrectly


	this->sim_thread = sim_thread;

    const int width = 1920;     // 1024
    const int height = 1080;    // 768
    resize(width, height);
    setMinimumWidth(width);
    setMinimumHeight(height);

	file_path = QDir::currentPath();

	// load settings
	FILE *file = fopen("settings.cfg", "r");

	if(file)
	{
		// check file version
		char buffer[200];
		fscanf(file, "%200s\n", buffer);

		if( !strcmp(buffer, CFG_FILE_VERSION) )
		{
			// load material
			double T_vis;
			fscanf(file, "%200s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", material_name, &particle_radius, &density, &surface_energy, &nu, &young_mod, &crit_rolling_displacement, &osc_damping_factor, &T_vis,
																	&rolling_modifier, &sliding_modifier, &twisting_modifier);
			sim_thread->sim.setMaterialConstants(particle_radius, density, surface_energy, nu, young_mod, crit_rolling_displacement, osc_damping_factor, T_vis, rolling_modifier, sliding_modifier, twisting_modifier);

            // load background/border colors
			fscanf(file, "%f %f %f", &clear_color[0], &clear_color[1], &clear_color[2]);
			fscanf(file, "%f %f %f", &border_color[0], &border_color[1], &border_color[2]);
			fscanf(file, "%f %f %f", &wall_color[0], &wall_color[1], &wall_color[2]);
			fscanf(file, "%f %f %f", &text_color[0], &text_color[1], &text_color[2]);

			// load shader settings
            fscanf(file, "%f %f %i %i %i %i %i %i %i %lf %lf %lf %i %i", &border_distance_factor, &border_min_depth, (int*)&draw_borders, (int*)&draw_particles, (int*)&track_fragments,
													&vis_neighbour_search_dist, &vis_neighbour_max_particles, &vis_neighbour_min_offset, &vis_neighbour_max_offset, 
                                                    &vis_density_max_filling_factor, &vis_dislocation_min_value, &vis_dislocation_max_value, &vis_velocity_averaging_steps, (int*)&vis_density_consider_walls);


            fscanf(file, "%i %i %i %i %i %i %i", (int*)&display_text_info, (int*)&display_pressure, (int*)&display_force, (int*)&display_changed_contacts, (int*)&display_depth, (int*)&display_key, (int*)&vis_display_mode);

            fscanf(file, "%i", (int*)&agg_info_only_inner_particles);

			// load settings for SimStartWidget
            fscanf(file, "%lf %lf %i %i %i %i %i %i %i %lf %lf %lf %i", &sim_timestep, &sim_duration, &print_energy_interval, &print_positions_interval, &take_screenshot_interval, (int*)&follow_cms,
                                                            &sim_time_display_mode, &sim_averaging_steps, (int*)&sim_use_gravity, &gravity_modifier, &wall_inertia, &sim_azimuthal_acceleration, (int*)&sim_use_sim_azimuthal_acceleration);

			// load settings for init chain widget
			fscanf(file, "%i %i %i %lf %lf %i", &init_chain_particles, &init_chain_impact_particles, &init_chain_target_particle, &init_chain_impact_speed, 
                                                &init_chain_angular_irregularity, (int*)&init_chain_impact_current_agglomerate);

			// load settings for init aggregate widget
            fscanf(file, "%u %lf %lf %i %lf %u %u %lf %lf %lf %lf %lf %lf", &init_aggregate_particles, &init_aggregate_migration_probability1,  &init_aggregate_migration_probability2, (int*)&init_aggregate_bam_selection_method,
												&init_aggregate_chain_ratio, &init_aggregate_min_chain_length, &init_aggregate_max_chain_length, &(init_aggregate_contact_distribution[0]), &(init_aggregate_contact_distribution[1]),
												&(init_aggregate_contact_distribution[2]), &(init_aggregate_contact_distribution[3]), &(init_aggregate_contact_distribution[4]), &(init_aggregate_contact_distribution[5]));

			// load settings for init cake widget
			fscanf(file, "%u %u %u %lf %u %lf %lf %i %lf %lf %lf %lf",  &init_cake_x_particles, &init_cake_y_particles, &init_cake_z_particles, &init_cake_filling_factor, &init_cake_particles,
                                                                        &init_cake_migration_probability1, &init_cake_migration_probability2, (int*)&init_cake_bam_selection_method,
																		&init_cake_x_size, &init_cake_y_size, &init_cake_top_slice_factor, &init_cake_bottom_slice_factor);

			// load settings for agglomerates collision widget
			fscanf(file, "%lf %lf %lf %i %i %i %i %u %u %i %i", &agglomerates_collision_impact_speed, &agglomerates_collision_impact_parameter, &agglomerates_collision_initial_separation,
                (int*)&agglomerates_collision_random_orientation1, (int*)&agglomerates_collision_random_orientation2, (int*)&agglomerates_collision_slice1, (int*)&agglomerates_collision_slice2,
                &agglomerates_collision_projectile_count, &agglomerates_collision_projectile_size, (int*)&agglomerates_collision_use_PCA_projectiles, (int*)&agglomerates_collision_random_impact_parameter);


			// load settings for wall collision widget
			fscanf(file, "%lf %i %lf %i", &wall_collision_impact_speed, &wall_collision_impact_angle, &wall_collision_impact_distance, 
                                            (int*)&wall_collision_random_orientation);

			// load settings box widget
			fscanf(file, "%lf %lf %lf %lf %lf %lf %i %i %i %i", &box_wall_speed, &box_pressure, &box_stop_filling_factor, &box_stop_dissipation_factor, &box_stop_pressure,
                                                            &box_dynamic_pressure, &box_mode, (int*)&box_random_orientation, (int*)&box_moving_bottom_wall, (int*)&box_side_walls);

			// load settings for add particles widget
            fscanf(file, "%u %lf %lf %i %u %u %u", &add_particles_number_of_particles, &add_particles_migration_probability1, &add_particles_migration_probability2, (int*)&add_particles_bam_selection_method, &add_particles_number_of_chains, &add_particles_min_chain_length, &add_particles_max_chain_length);

			// load settings for rotation widget
			fscanf(file, "%lf %lf %lf %lf", &rotation_axis_x, &rotation_axis_y, &rotation_axis_z, &rotation_angle);

			// load settings for slice widget
            fscanf(file, "%lf %lf %lf %lf %lf %i", &slice_sphere_radius, &slice_box_x_size, &slice_box_y_size, &slice_box_z_size, &slice_top_slice_factor, (int*)&slice_random_orientation);

			// load settings for duplicate widget
            fscanf(file, "%i %i %i %i %i %i", &duplicate_x_duplications, &duplicate_y_duplications, &duplicate_z_duplications, (int*)&duplicate_slice_agglomerate,
                                            (int*)&duplicate_random_orientation, (int*)&duplicate_mirror);

			// load load/save path
			char buffer[300];
			fscanf(file, "%s", buffer);

			file_path = buffer;

			fclose(file);
		}
		else
			QMessageBox::about(this, tr("About ParticleTest"), tr("Settings file outdated - using default settings") );
	}

	setWindowTitle( tr("Particle Viewer v").append(GUI_VERSION) );

	// Init window for visualization
#ifdef DRAW_CONTACT_POINTERS
	glView = new OpenGLWidget(this, &current_pos, &current_particle_colors, &current_walls, &current_contact_pointers);
#else
	glView = new OpenGLWidget(this, sim_thread);
#endif

	// create child widgets
	startSimWidget = new StartSimWidget();
	initChainWidget = new InitChainWidget();
	initAggregateWidget = new InitAggregateWidget();
	initCakeWidget = new InitCakeWidget();
	collideAgglomeratesWidget = new CollideAgglomeratesWidget();
	wallCollisionWidget = new WallCollisionWidget();
	initBoxWidget = new InitBoxWidget();
	addParticlesWidget = new AddParticlesWidget();
	sliceWidget = new SliceWidget();
	rotationWidget = new RotationWidget();
	duplicationWidget = new DuplicationWidget();
	modifyWallWidget = new ModifyWallWidget(sim_thread);
	agglomerateInfoWidget = new AgglomerateInfoWidget();
	materialInfoWidget = new MaterialInfoWidget();
	visualizationOptionsWidget = new VisualizationOptionsWidget();
	displayOptionsWidget = new DisplayOptionsWidget();

#ifdef DRAW_CONTACT_POINTERS
	current_contact_pointers = &contact_pointers1;
	new_contact_pointers = &contact_pointers2;

	updating_contact_pointers = false;
#endif

	pause = false;

	// menu
	runSimAct = new QAction(tr("Run Simulation"), this);
	runSimAct->setStatusTip(tr("Runs simulation"));
	runSimAct->setShortcut(tr("Ctrl+R"));
	connect(runSimAct, SIGNAL(triggered()), this, SLOT(showRunSimWidget()));

	loadSimAct = new QAction(tr("Load Data"), this);
	loadSimAct->setStatusTip(tr("Loads positions and velocities of particles from a file"));
	loadSimAct->setShortcut(tr("Ctrl+L"));
	connect(loadSimAct, SIGNAL(triggered()), this, SLOT(showLoadWidget()));

	importAct = new QAction(tr("Import"), this);
	importAct->setStatusTip(tr("Import positions from other formats"));
	connect(importAct, SIGNAL(triggered()), this, SLOT(showImportWidget()));

	saveSimAct = new QAction(tr("Save Data"), this);
	saveSimAct->setStatusTip(tr("Stores current positions and velocities of particles to a file"));
	saveSimAct->setShortcut(tr("Ctrl+S"));
	connect(saveSimAct, SIGNAL(triggered()), this, SLOT(showSaveWidget()));

	savePositionsAct = new QAction(tr("Save Positions"), this);
	savePositionsAct->setStatusTip(tr("Stores current positions of particles to a file"));
	connect(savePositionsAct, SIGNAL(triggered()), this, SLOT(showSavePositionsWidget()));

	screenshotAct = new QAction(tr("Save Image"), this);
	screenshotAct->setStatusTip(tr("Saves the current visualization as a png file"));
	screenshotAct->setShortcut(tr("Ctrl+I"));
	connect(screenshotAct, SIGNAL(triggered()), this, SLOT(customScreenshot()));

    folderToScreenshotsAct = new QAction(tr("Folder To Screenshots"), this);
    folderToScreenshotsAct->setStatusTip(tr("Turns all files in folder into screenshots"));
    connect(folderToScreenshotsAct, SIGNAL(triggered()), this, SLOT(showFolderToScreenshotsWidget()));

    posFolderToScreenshotsAct = new QAction(tr("Folder with position data To Screenshots"), this);
    posFolderToScreenshotsAct->setStatusTip(tr("Turns all files in folder into screenshots"));
    connect(posFolderToScreenshotsAct, SIGNAL(triggered()), this, SLOT(showPosFolderToScreenshotsWidget()));


	enterSeedAct = new QAction(tr("Enter Seed"), this);
	enterSeedAct->setStatusTip(tr("initializes random number generator with specified seed"));
    connect(enterSeedAct, SIGNAL(triggered()), this, SLOT(enterSeed()));

	exitAct = new QAction(tr("Exit"), this);
	exitAct->setShortcut(tr("Ctrl+Q"));
	exitAct->setStatusTip(tr("Exit the application"));
	connect(exitAct, SIGNAL(triggered()), this, SLOT(close()));

	initChainAct = new QAction(tr("Init chain"), this);
	initChainAct->setStatusTip(tr("Setup collision of chains"));
	connect(initChainAct, SIGNAL(triggered()), this, SLOT(showInitChainWidget()));

	initAggregateAct = new QAction(tr("Init Aggregate"), this);
	initAggregateAct->setStatusTip(tr("Provides different aggregation building methods"));
	connect(initAggregateAct, SIGNAL(triggered()), this, SLOT(showInitAggregateWidget()));

	initCakeAct = new QAction(tr("Init Cake"), this);
	initCakeAct->setStatusTip(tr("Provides different cake building methods"));
	connect(initCakeAct, SIGNAL(triggered()), this, SLOT(showInitCakeWidget()));

	collideAgglomeratesAct = new QAction(tr("Setup Collision"), this);
	collideAgglomeratesAct->setStatusTip(tr("Setup collision of two agglomerates"));
	connect(collideAgglomeratesAct, SIGNAL(triggered()), this, SLOT(showCollideAgglomeratesWidget()));

	collideWallAct = new QAction(tr("Setup Wall Collision"), this);
	collideWallAct->setStatusTip(tr("Setup collision with wall"));
	connect(collideWallAct, SIGNAL(triggered()), this, SLOT(showInitWallWidget()));

	initBoxAct = new QAction(tr("Setup Box"), this);
	initBoxAct->setStatusTip(tr("Initializes confining box"));
	connect(initBoxAct, SIGNAL(triggered()), this, SLOT(showInitBoxWidget()));

	initPullStrengthTestAct = new QAction(tr("Pull Strength Test"), this);
	initPullStrengthTestAct->setStatusTip(tr("Initializes testing of tensile strength"));
	connect(initPullStrengthTestAct, SIGNAL(triggered()), this, SLOT(initPullStrengthTest()));

	initShearStrengthTestAct = new QAction(tr("Shear Strength Test"), this);
	initShearStrengthTestAct->setStatusTip(tr("Initializes testing of shear strength"));
	connect(initShearStrengthTestAct, SIGNAL(triggered()), this, SLOT(initShearStrengthTest()));

	addParticlesAct = new QAction(tr("Add Particles"), this);
	addParticlesAct->setStatusTip(tr("Adds more particles to current aggregate"));
	connect(addParticlesAct, SIGNAL(triggered()), this, SLOT(showAddParticlesWidget()));

	rotateAct = new QAction(tr("Rotate Aggregate"), this);
	rotateAct->setStatusTip(tr("Rotates aggregate around desired axis"));
	connect(rotateAct, SIGNAL(triggered()), this, SLOT(showRotationWidget()));

	sliceAct = new QAction(tr("Slice Aggregate"), this);
	sliceAct->setStatusTip(tr("Slices aggregate in various ways"));
	connect(sliceAct, SIGNAL(triggered()), this, SLOT(showSliceWidget()));

	duplicateAct = new QAction(tr("Duplicate Aggregate"), this);
	duplicateAct->setStatusTip(tr("Duplicates aggregates in various ways"));
	connect(duplicateAct, SIGNAL(triggered()), this, SLOT(showDuplicationWidget()));

	filterFragmentsAct = new QAction(tr("&Filter fragments"), this);
	filterFragmentsAct->setStatusTip(tr("Remove all fragments"));
	filterFragmentsAct->setShortcut(tr("Ctrl+F"));
	connect(filterFragmentsAct, SIGNAL(triggered()), sim_thread, SLOT(filterFragments()));

	removeWallsAct = new QAction(tr("Remove Walls"), this);
	removeWallsAct->setStatusTip(tr("Removes all walls from the simulation"));
	connect(removeWallsAct, SIGNAL(triggered()), sim_thread, SLOT(removeWalls()));

	modifyWallsAct = new QAction(tr("Modify Walls"), this);
	modifyWallsAct->setStatusTip(tr("Allows to change the properties of existing walls"));
	connect(modifyWallsAct, SIGNAL(triggered()), this, SLOT(showModifyWallsWidget()));

	resetContactsAct = new QAction(tr("Reset Contacts"), this);
	resetContactsAct->setStatusTip(tr("All contacts will be set to the initial state"));
	connect(resetContactsAct, SIGNAL(triggered()), sim_thread, SLOT(resetContacts()));

	disturbParticlesAct = new QAction(tr("Disturb Particles"), this);
	disturbParticlesAct->setStatusTip(tr("Positions of particles will be randomly disturbed"));
	connect(disturbParticlesAct, SIGNAL(triggered()), this, SLOT(disturbParticles()));

	displayCrossSectionAct = new QAction(tr("Determine Cross Section"), this);
	displayCrossSectionAct->setStatusTip(tr("Determines the projected cross section of the current aggregate"));
	connect(displayCrossSectionAct, SIGNAL(triggered()), this, SLOT(displayCrossSection()));

	moveToOriginAct = new QAction(tr("Center agglomerate"), this);
	moveToOriginAct->setStatusTip(tr("Translates the agglomerate to move the center of mass to the origin"));
	connect(moveToOriginAct, SIGNAL(triggered()), sim_thread, SLOT(moveToOrigin()));

	showAgglomerateInfoAct = new QAction(tr("Properties"), this);
	showAgglomerateInfoAct->setStatusTip(tr("Displays various properties about the current agglomerate"));
	connect(showAgglomerateInfoAct, SIGNAL(triggered()), this, SLOT(showAgglomerateInfo()));

	printContactHistogramAct = new QAction(tr("Export Contact Histogramm"), this);
	printContactHistogramAct->setStatusTip(tr("Prints the histogramm of the number of contacts per particle to an external file"));
	connect(printContactHistogramAct, SIGNAL(triggered()), this, SLOT(printContactHistogram()));

	printFragmentVelocitiesAct = new QAction(tr("Print Fragment Velocities"), this);
	printFragmentVelocitiesAct->setStatusTip(tr("Prints the velocities of the individual fragments to an external file"));
	connect(printFragmentVelocitiesAct, SIGNAL(triggered()), this, SLOT(printFragmentVelocities()));

	printElasticChargeAct = new QAction(tr("Print Elastic Charge"), this);
	printElasticChargeAct->setStatusTip(tr("Prints the current elastic charge to an external file"));
	connect(printElasticChargeAct, SIGNAL(triggered()), this, SLOT(printElasticCharge()));

    printFillingFactorProfileAct = new QAction(tr("Export Filling Factor Profile Box"), this);
	printFillingFactorProfileAct->setStatusTip(tr("Prints a vertical profile of the local filling factor to an external file"));
	connect(printFillingFactorProfileAct, SIGNAL(triggered()), this, SLOT(printFillingFactorProfile()));

    printFillingFactorProfileSphereAct = new QAction(tr("Export Filling Factor Profile Sphere"), this);
    printFillingFactorProfileSphereAct->setStatusTip(tr("Prints a profile of the local filling factor of the spherical layer to an external file"));
    connect(printFillingFactorProfileSphereAct, SIGNAL(triggered()), this, SLOT(printFillingFactorProfileSphere()));

	printFractalDimensionAct = new QAction(tr("Export Fractal Dimension"), this);
	printFractalDimensionAct->setStatusTip(tr("Prints information required to determine the fractal dimension of the current aggregate to an external file"));
	connect(printFractalDimensionAct, SIGNAL(triggered()), this, SLOT(printFractalDimension()));

	loadMaterialAct = new QAction(tr("Load material"), this);
	loadMaterialAct->setStatusTip(tr("Load material properties from file"));
	connect(loadMaterialAct, SIGNAL(triggered()), this, SLOT(loadMaterial()));

	showMaterialInfoAct = new QAction(tr("Show Material Info"), this);
	showMaterialInfoAct->setStatusTip(tr("Displays additional information about the current material"));
	connect(showMaterialInfoAct, SIGNAL(triggered()), this, SLOT(showMaterialInfo()));

	enterInfoTextAct = new QAction(tr("Display Description"), this);
	enterInfoTextAct->setStatusTip(tr("Adds a short description text to the visualization"));
	connect(enterInfoTextAct, SIGNAL(triggered()), this, SLOT(enterInfoText()));

	centerViewAct = new QAction(tr("Center View"), this);
	centerViewAct->setStatusTip(tr("Adjusts the view on the center of mass"));
	connect(centerViewAct, SIGNAL(triggered()), this, SLOT(centerView()));

	toggleGyrationRadiusAct = new QAction(tr("Toggle gyration radius"), this);
	toggleGyrationRadiusAct->setStatusTip(tr("Toggles the display of the gyration radius"));
	connect(toggleGyrationRadiusAct, SIGNAL(triggered()), sim_thread, SLOT(toggleGyrationRadius()));

	showVisualizationOptionsAct = new QAction(tr("Configure Visualization"), this);
	showVisualizationOptionsAct->setStatusTip(tr("Shows various options for particle coloring"));
	connect(showVisualizationOptionsAct, SIGNAL(triggered()), this, SLOT(showVisualizationOptions()));

	showDisplayOptionsAct = new QAction(tr("Display Options"), this);
	showDisplayOptionsAct->setStatusTip(tr("Allows to select various display options"));
	connect(showDisplayOptionsAct, SIGNAL(triggered()), this, SLOT(showDisplayOptions()));

	aboutAct = new QAction(tr("About"), this);
	aboutAct->setStatusTip(tr("Show the application's About box"));
	connect(aboutAct, SIGNAL(triggered()), this, SLOT(showAboutMessage()));

	// set up menu bar
	fileMenu = menuBar()->addMenu(tr("&Simulation"));
	fileMenu->addAction(runSimAct);
	fileMenu->addAction(loadSimAct);
	fileMenu->addAction(importAct);
	fileMenu->addAction(saveSimAct);
	fileMenu->addAction(savePositionsAct);
	fileMenu->addAction(screenshotAct);
    fileMenu->addAction(folderToScreenshotsAct);
    fileMenu->addAction(posFolderToScreenshotsAct);
	fileMenu->addAction(enterSeedAct);
	fileMenu->addAction(exitAct);

	initMenu = menuBar()->addMenu(tr("&Init"));
	initMenu->addAction(initChainAct);
	initMenu->addAction(initAggregateAct);
	initMenu->addAction(initCakeAct);
	initMenu->addAction(collideAgglomeratesAct);
	initMenu->addAction(collideWallAct);
	initMenu->addAction(initBoxAct);
	initMenu->addAction(initPullStrengthTestAct);
	initMenu->addAction(initShearStrengthTestAct);

	modifyMenu = menuBar()->addMenu(tr("&Modify"));
	modifyMenu->addAction(addParticlesAct);
	modifyMenu->addAction(rotateAct);
	modifyMenu->addAction(sliceAct);
	modifyMenu->addAction(duplicateAct);
	modifyMenu->addAction(filterFragmentsAct);
	modifyMenu->addAction(modifyWallsAct);
	modifyMenu->addAction(removeWallsAct);
	modifyMenu->addAction(resetContactsAct);
	modifyMenu->addAction(disturbParticlesAct);

	aggregateMenu = menuBar()->addMenu(tr("&Aggregate"));
	aggregateMenu->addAction(showAgglomerateInfoAct);
	aggregateMenu->addAction(displayCrossSectionAct);
	aggregateMenu->addAction(printContactHistogramAct);
	aggregateMenu->addAction(printFragmentVelocitiesAct);
	aggregateMenu->addAction(printElasticChargeAct);
	aggregateMenu->addAction(printFillingFactorProfileAct);
    aggregateMenu->addAction(printFillingFactorProfileSphereAct);
	aggregateMenu->addAction(printFractalDimensionAct);
	aggregateMenu->addAction(moveToOriginAct);

	materialMenu = menuBar()->addMenu(tr("&Material"));
	materialMenu->addAction(loadMaterialAct);
	materialMenu->addAction(showMaterialInfoAct);

	viewMenu = menuBar()->addMenu(tr("&View"));
	viewMenu->addAction(showVisualizationOptionsAct);
	viewMenu->addAction(showDisplayOptionsAct);
	viewMenu->addAction(enterInfoTextAct);
	viewMenu->addAction(toggleGyrationRadiusAct);
	viewMenu->addAction(centerViewAct);

	helpMenu = menuBar()->addMenu(tr("&Help"));
	helpMenu->addAction(aboutAct);

	// disable when no data present
	centerViewAct->setDisabled(true);

	// Communication GUI -> simulation thread
	connect(this, SIGNAL(signalLoadDataFromFile(char*)), sim_thread, SLOT(loadFromFile(char*)));
    connect(this, SIGNAL(signalImportDataFromFile(char*)), sim_thread, SLOT(importFromFile(char*)));
	connect(this, SIGNAL(signalSaveDataToFile(char*)), sim_thread, SLOT(saveToFile(char*)));
	connect(this, SIGNAL(signalSavePositionsToFile(char*)), sim_thread, SLOT(savePositionsToFile(char*)));
	connect(this, SIGNAL(signalLoadMaterial(char*)), sim_thread, SLOT(loadMaterial(char*)));
	connect(this, SIGNAL(signalInitSimChain(int, int, int, double, double)), sim_thread, SLOT(initChain(int, int, int, double, double)));
	connect(this, SIGNAL(signalImpactChainOnAgglomerate(int, double, double)), sim_thread, SLOT(impactChainOnAgglomerate(int, double, double)));
	connect(this, SIGNAL(signalInitPullStrengthTest(double)), sim_thread, SLOT(initPullStrengthTest(double)));
	connect(this, SIGNAL(signalInitShearStrengthTest(double)), sim_thread, SLOT(initShearStrengthTest(double)));
	connect(this, SIGNAL(signalDisturbParticles(double)), sim_thread, SLOT(disturbParticles(double)));

	connect(this, SIGNAL(signalPauseSimulation()), sim_thread, SLOT(pauseSimulation()));
	connect(this, SIGNAL(signalContinueSimulation()), sim_thread, SLOT(continueSimulation()));

	// Communication simulation thread -> GUI
	connect(sim_thread, SIGNAL(signalInitFailed(ErrorCode)), this , SLOT(showErrorMessage(ErrorCode)));
	connect(sim_thread, SIGNAL(signalNumberOfParticlesChanged(int)), this, SLOT(numberOfParticlesChanged(int)));
	connect(sim_thread, SIGNAL(signalNumberOfWallsChanged(int)), modifyWallWidget, SLOT(setupWallSelectionBox(int)));
	connect(sim_thread, SIGNAL(signalTakeScreenshot()), this, SLOT(autoScreenshot()));
	connect(agglomerateInfoWidget->refreshButton, SIGNAL(clicked(bool)), this, SLOT(updateAgglomerateInfo()));
	connect(agglomerateInfoWidget, SIGNAL(signalUpdateAgglomerateInfo()), this, SLOT(updateAgglomerateInfo()));

	// Communication  simulation thread -> visualization
	connect(sim_thread, SIGNAL(signalMaterialLoaded()), glView, SLOT(updateGL()));
	connect(sim_thread, SIGNAL(signalUpdateGLView()), this, SLOT(updateVisualization()));
	connect(initChainWidget->initChainButton, SIGNAL(clicked(bool)), this, SLOT(initChain()));

	// Communication widgets -> simulation thread
	connect(startSimWidget, SIGNAL(signalStartSim(double, double, int, int, int, int, bool, bool)), sim_thread, SLOT(startSimulation(double, double, int, int, int, int, bool, bool)));

	connect(initCakeWidget, SIGNAL(signalInitSimpleCubicLattice(unsigned int, unsigned int, unsigned int, double)), sim_thread, SLOT(initSimpleCubicLattice(unsigned int, unsigned int, unsigned int, double)));
	connect(initCakeWidget, SIGNAL(signalInitHexagonalLattice(unsigned int, unsigned int, unsigned int, double)), sim_thread, SLOT(initHexagonalLattice(unsigned int, unsigned int, unsigned int, double)));
	connect(initCakeWidget, SIGNAL(signalInitRBDCake(unsigned int, double, double, double, double)), sim_thread, SLOT(initRBDCake(unsigned int, double, double, double, double)));
	connect(initCakeWidget, SIGNAL(signalInitBAMCake(unsigned int, double, double, double, double, double, double, BAMSelectionMethod)), sim_thread, SLOT(initBAMCake(unsigned int, double, double, double, double, double, double, BAMSelectionMethod)));

	connect(initAggregateWidget, SIGNAL(signalInitBAMAggregate(unsigned int, double, double, BAMSelectionMethod)), sim_thread, SLOT(initBAMAggregate(unsigned int, double, double, BAMSelectionMethod)));
	connect(initAggregateWidget, SIGNAL(signalInitFractalAggregate(unsigned int, double, double, double, unsigned int, unsigned int)), sim_thread, SLOT(initFractalAggregate(unsigned int, double, double, double, unsigned int, unsigned int)));
	connect(initAggregateWidget, SIGNAL(signalInitTreeAggregate(unsigned int, double, double, double, double, double, double, double, double)), sim_thread, SLOT(initTreeAggregate(unsigned int, double, double, double, double, double, double, double, double)));
	connect(initAggregateWidget, SIGNAL(signalInitNeedleDeposition(unsigned int)), sim_thread, SLOT(initNeedleDeposition(unsigned int)));

	connect(collideAgglomeratesWidget, SIGNAL(signalCollideAgglomerates(const char*, const char*, bool, bool, double, double, double, bool)), sim_thread, SLOT(collideAgglomerates(const char*, const char*, bool, bool, double, double, double, bool)));
	connect(collideAgglomeratesWidget, SIGNAL(signalImpactMultiProjectiles(const char*, const char*, unsigned int, double, double, bool, double, bool)), sim_thread, SLOT(impactMultiProjectiles(const char*, const char*, unsigned int, double, double, bool, double, bool)));
	connect(collideAgglomeratesWidget, SIGNAL(signalImpactMultiPCAProjectiles(const char*, unsigned int, unsigned int, double, double, bool, double, bool)), sim_thread, SLOT(impactMultiPCAProjectiles(const char*, unsigned int, unsigned int, double, double, bool, double, bool)));

	connect(collideAgglomeratesWidget, SIGNAL(signalHitAndStick(const char*, const char*, bool, bool, double)), sim_thread, SLOT(hitAndStick(const char*, const char*, bool, bool, double)));
	connect(wallCollisionWidget, SIGNAL(signalCollideAgglomerateWithWall(const char*, double, int, double, bool)), sim_thread, SLOT(collideAgglomerateWithWall(const char*, double, int, double, bool)));
	connect(initBoxWidget, SIGNAL(signalInitCompression(const char*, bool, bool, bool, double, double, double)), sim_thread, SLOT(initCompression(const char*, bool, bool, bool, double, double, double)));
	connect(initBoxWidget, SIGNAL(signalInitCompressionRelaxation(const char*, bool, double, double, double)), sim_thread, SLOT(initCompressionRelaxation(const char*, bool, double, double, double)));
	connect(initBoxWidget, SIGNAL(signalInitShockwave(const char*, bool, double, double, double)), sim_thread, SLOT(initShockwave(const char*, bool, double, double, double)));
	connect(initBoxWidget, SIGNAL(signalInitDynamicCompression(const char*, bool, double, double, double)), sim_thread, SLOT(initDynamicCompression(const char*, bool, double, double, double)));
	connect(initBoxWidget, SIGNAL(signalInitOpenBox(const char*, bool, bool, double)), sim_thread, SLOT(initOpenBox(const char*, bool, bool, double)));

	connect(addParticlesWidget, SIGNAL(signalAddBAMParticles(unsigned int, double, double, const char*, BAMSelectionMethod)), sim_thread, SLOT(addBAMParticles(unsigned int, double, double, const char*, BAMSelectionMethod)));
	connect(addParticlesWidget, SIGNAL(signalAddFractalChains(const char*, unsigned int, unsigned int, unsigned int)), sim_thread, SLOT(addFractalChains(const char*, unsigned int, unsigned int, unsigned int)));
	connect(rotationWidget, SIGNAL(signalRotateAggregate(double, double, double, double)), sim_thread, SLOT(rotateAggregate(double, double, double, double)));
	connect(sliceWidget, SIGNAL(signalSliceBox(const char*, double, double, double)), sim_thread, SLOT(sliceBox(const char*, double, double, double)));
	connect(sliceWidget, SIGNAL(signalSliceSphere(const char*, double, bool)), sim_thread, SLOT(sliceSphere(const char*, double, bool)));
	connect(sliceWidget, SIGNAL(signalSliceCylinder(const char*, double)), sim_thread, SLOT(sliceCylinder(const char*, double)));
	connect(sliceWidget, SIGNAL(signalSliceTop(const char*, double)), sim_thread, SLOT(sliceTop(const char*, double)));
	connect(sliceWidget, SIGNAL(signalSliceBottom(const char*, double)), sim_thread, SLOT(sliceBottom(const char*, double)));
	connect(duplicationWidget, SIGNAL(signalDuplicate(int, int, int, bool, bool)), sim_thread, SLOT(duplicateAggregate(int, int, int, bool, bool)));
	connect(modifyWallWidget, SIGNAL(signalSetWallProperties(int, double, double, double, double, double, double)), sim_thread, SLOT(setWallProperties(int, double, double, double, double, double, double)));

	// Communication DisplayOptions -> OpenGLView
	connect(startSimWidget, SIGNAL(signalUpdateView()), glView, SLOT(updateGL()));
	connect(displayOptionsWidget, SIGNAL(signalUpdateView()), glView, SLOT(updateGL()));
	connect(displayOptionsWidget, SIGNAL(signalUpdateShaders()), glView, SLOT(selectShader()));
	connect(displayOptionsWidget, SIGNAL(signalUpdateColors()), sim_thread, SLOT(updateDrawnObjects()));
	connect(visualizationOptionsWidget, SIGNAL(signalUpdateColors()), sim_thread, SLOT(updateDrawnObjects()));
	connect(visualizationOptionsWidget, SIGNAL(signalUpdateColors()), glView, SLOT(updateKey()));

#ifdef DRAW_CONTACT_POINTERS
	connect(this, SIGNAL(signalGetContactPointers(multimap<int, vec3>*)), sim_thread, SLOT(StoreContactPointers(multimap<int, vec3>*)));
	connect(sim_thread, SIGNAL(ContactPointerUpdateDone()), this, SLOT(contactPointerUpdateDone()));
#endif

	cms[0] = 0;
	cms[1] = 0;
	cms[2] = 0;

	agg_info.filling_factor_box = 0.0;
	agg_info.filling_factor_sphere = 0.0;
	agg_info.fragments = 0;
	agg_info.gyration_radius = 0.0;
	agg_info.outer_radius = 0.0;
	agg_info.box_height = 0.0;
	agg_info.box_base = 0.0;
	memset(agg_info.contact_histogram, 0, 13 * sizeof(int));

	sim_info.particles = 0;
	sim_info.broken_contacts = 0;
	sim_info.created_contacts = 0;
	sim_info.time = 0;
	sim_info.filling_factor = 0;
	sim_info.wall_speed = 0;
	sim_info.collision_speed = 0;
	sim_info.pressure = -1;
	sim_info.force = -1;
	sim_info.coordination_number = 0;

	// init key
	glView->updateKey();

	// screenshot
	screenshot_counter = 0;

#if defined (_WIN32) || (_WIN64)
	QString path = "images\\";
#else
	QString path = QDir::currentPath();
	path.append("/images/");
#endif

    strcpy(screenshot_path, path.toLatin1().data() );
}

MainWindow::~MainWindow()
{
	// save settings
	FILE *file = fopen("settings.cfg", "w+");

	if(file)
	{
		fprintf(file, "%s\n", CFG_FILE_VERSION);

		// save material
		fprintf(file, "%s %g %g %g %g %g %g %g %g %g %g %g\n", material_name, particle_radius, density, surface_energy, nu, young_mod, crit_rolling_displacement, osc_damping_factor, T_vis,
														rolling_modifier, sliding_modifier, twisting_modifier);
		
		// save background/border colors
		fprintf(file, "%f %f %f\n", clear_color[0], clear_color[1], clear_color[2]);
		fprintf(file, "%f %f %f\n", border_color[0], border_color[1], border_color[2]);
		fprintf(file, "%f %f %f\n", wall_color[0], wall_color[1], wall_color[2]);
		fprintf(file, "%f %f %f\n", text_color[0], text_color[1], text_color[2]);

		// save shader settings
		fprintf(file, "%f %f %i %i %i %i %i %i %i %g %g %g %i %i\n", border_distance_factor, border_min_depth, (int)draw_borders, (int)draw_particles, (int)track_fragments, 
													vis_neighbour_search_dist, vis_neighbour_max_particles, vis_neighbour_min_offset, vis_neighbour_max_offset,
													vis_density_max_filling_factor, vis_dislocation_min_value, vis_dislocation_max_value, vis_velocity_averaging_steps, (int)vis_density_consider_walls);

		fprintf(file, "%i %i %i %i %i %i %i\n", (int)display_text_info, (int)display_pressure, (int)display_force, (int)display_changed_contacts, (int)display_depth, (int)display_key, (int)vis_display_mode);

		fprintf(file, "%i\n", (int)agg_info_only_inner_particles);

		// save settings for SimStartWidget
		fprintf(file, "%g %g %i %i %i %i %i %i %i %g %g %g %i\n", sim_timestep, sim_duration, print_energy_interval, print_positions_interval, take_screenshot_interval, (int)follow_cms, sim_time_display_mode, sim_averaging_steps,
														(int)sim_use_gravity, gravity_modifier, wall_inertia, sim_azimuthal_acceleration, (int)sim_use_sim_azimuthal_acceleration);

		// save settings for init chains widget
		fprintf(file, "%i %i %i %g %g %i\n", init_chain_particles, init_chain_impact_particles, init_chain_target_particle, init_chain_impact_speed, init_chain_angular_irregularity, (int)init_chain_impact_current_agglomerate);	

		// save settings for init aggregate widget
		fprintf(file, "%u %lf %lf %i %lf %u %u %lf %lf %lf %lf %lf %lf\n", init_aggregate_particles, init_aggregate_migration_probability1,  init_aggregate_migration_probability2, (int)init_aggregate_bam_selection_method,
												init_aggregate_chain_ratio, init_aggregate_min_chain_length, init_aggregate_max_chain_length,  init_aggregate_contact_distribution[0], init_aggregate_contact_distribution[1],
												init_aggregate_contact_distribution[2], init_aggregate_contact_distribution[3], init_aggregate_contact_distribution[4], init_aggregate_contact_distribution[5]);

		// save settings for init cake widget
		fprintf(file, "%u %u %u %lf %u %lf %lf %i %lf %lf %lf %lf\n", init_cake_x_particles, init_cake_y_particles, init_cake_z_particles, init_cake_filling_factor, init_cake_particles,
																	init_cake_migration_probability1, init_cake_migration_probability2, (int)init_cake_bam_selection_method, 
																	init_cake_x_size, init_cake_y_size, init_cake_top_slice_factor, init_cake_bottom_slice_factor);

		// save settings for agglomerates collision widget
		fprintf(file, "%g %g %g %i %i %i %i %u %u %i %i\n", agglomerates_collision_impact_speed, agglomerates_collision_impact_parameter, agglomerates_collision_initial_separation,
						(int)agglomerates_collision_random_orientation1, (int)agglomerates_collision_random_orientation2, (int)agglomerates_collision_slice1, (int)agglomerates_collision_slice2,
						agglomerates_collision_projectile_count, agglomerates_collision_projectile_size, (int)agglomerates_collision_use_PCA_projectiles, (int)agglomerates_collision_random_impact_parameter);

		// save settings for wall collision widget
		fprintf(file, "%g %i %g %i\n", wall_collision_impact_speed, wall_collision_impact_angle, wall_collision_impact_distance, (int)wall_collision_random_orientation);

		// save settings for box widget
		fprintf(file, "%g %g %g %g %g %g %i %i %i %i\n", box_wall_speed, box_pressure, box_stop_filling_factor, box_stop_dissipation_factor, box_stop_pressure,
														box_dynamic_pressure, box_mode, (int)box_random_orientation, (int)box_moving_bottom_wall, (int)box_side_walls);

		// load settings for add particles widget
		fprintf(file, "%u %lf %lf %i %u %u %u\n", add_particles_number_of_particles, add_particles_migration_probability1, add_particles_migration_probability2, (int)add_particles_bam_selection_method, add_particles_number_of_chains, add_particles_min_chain_length, add_particles_max_chain_length);

		// save settings for rotation widget
		fprintf(file, "%g %g %g %g\n", rotation_axis_x, rotation_axis_y, rotation_axis_z, rotation_angle);

		// save settings for slice widget
		fprintf(file, "%lf %lf %lf %lf %lf %i\n", slice_sphere_radius, slice_box_x_size, slice_box_y_size, slice_box_z_size, slice_top_slice_factor, (int)slice_random_orientation);

		// save settings for duplicate widget
		fprintf(file, "%i %i %i %i %i %i\n", duplicate_x_duplications, duplicate_y_duplications, duplicate_z_duplications, (int)duplicate_slice_agglomerate, (int)duplicate_random_orientation, (int)duplicate_mirror);

		// save load/save file path
		QByteArray temp = file_path.toLocal8Bit();
		fprintf(file, "%s\n", temp.data());

		fclose(file);
	}

	// clean up
	delete startSimWidget;
	delete initChainWidget;
	delete initAggregateWidget;
	delete initCakeWidget;
	delete collideAgglomeratesWidget;
	delete wallCollisionWidget;
	delete initBoxWidget;
	delete addParticlesWidget;
	delete rotationWidget;
	delete sliceWidget;
	delete duplicationWidget;
	delete modifyWallWidget;
	delete agglomerateInfoWidget;
	delete materialInfoWidget;
	delete visualizationOptionsWidget;
	delete displayOptionsWidget;
}

void MainWindow::closeEvent(QCloseEvent *event)
 {
    startSimWidget->close();
	initChainWidget->close();
	initAggregateWidget->close();
	initCakeWidget->close();
	collideAgglomeratesWidget->close();
	wallCollisionWidget->close();
	initBoxWidget->close();
	addParticlesWidget->close();
	rotationWidget->close();
	sliceWidget->close();
	duplicationWidget->close();
	modifyWallWidget->close();
	agglomerateInfoWidget->close();
	materialInfoWidget->close();
	visualizationOptionsWidget->close();
	displayOptionsWidget->close();
 }

void MainWindow::resizeEvent(QResizeEvent *event)
{
	glView->setGeometry(0, 20, width(), height());

	QMainWindow::resizeEvent(event);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Slots
//////////////////////////////////////////////////////////////////////////////////////////////////////

void MainWindow::showAboutMessage()
{
	char s[100];
	sprintf(s, "Core version: %s\nGUI version: %s\n\nWritten by Alexander Seizinger, 2009-2012", CORE_VERSION, GUI_VERSION);
	QMessageBox::about(this, tr("About ParticleTest"), s);
}

void MainWindow::showRunSimWidget()
{
	startSimWidget->show();
	startSimWidget->activateWindow();
}

void MainWindow::showInitChainWidget()
{
	initChainWidget->show();
	initChainWidget->activateWindow();
}

void MainWindow::showInitWallWidget()
{
	wallCollisionWidget->show();
	wallCollisionWidget->activateWindow();
}

void MainWindow::showInitBoxWidget()
{
	initBoxWidget->show();
	initBoxWidget->activateWindow();
}

void MainWindow::initPullStrengthTest()
{
	bool ok;
	double pull_speed = QInputDialog::getDouble(this, tr("Pull Speed"), tr("Pull speed:"), 10, 0, 10000, 2, &ok);
	
	if(ok)
		emit signalInitPullStrengthTest(pull_speed);
}

void MainWindow::initShearStrengthTest()
{
	bool ok;
	double pull_speed = QInputDialog::getDouble(this, tr("Pull Speed"), tr("Pull speed:"), 10, 0, 10000, 2, &ok);
	
	if(ok)
		emit signalInitShearStrengthTest(pull_speed);
}

void MainWindow::showCollideAgglomeratesWidget()
{
	collideAgglomeratesWidget->show();
	collideAgglomeratesWidget->activateWindow();
}

void MainWindow::showInitAggregateWidget()
{
	initAggregateWidget->show();
	initAggregateWidget->activateWindow();
}

void MainWindow::showInitCakeWidget()
{
	initCakeWidget->show();
	initCakeWidget->activateWindow();
}

void MainWindow::showAddParticlesWidget()
{
	addParticlesWidget->show();
	addParticlesWidget->activateWindow();
}

void MainWindow::showRotationWidget()
{
	rotationWidget->show();
	rotationWidget->activateWindow();
}

void MainWindow::showSliceWidget()
{
	sliceWidget->show();
	sliceWidget->activateWindow();
}

void MainWindow::showDuplicationWidget()
{
	duplicationWidget->show();
	duplicationWidget->activateWindow();
}

void MainWindow::showModifyWallsWidget()
{
	modifyWallWidget->show();
	modifyWallWidget->activateWindow();
}

void MainWindow::disturbParticles()
{
	bool ok;
	double max_dist = QInputDialog::getDouble(this, tr("Enter disturbance"), tr("Maximum dislocation in particle radii:"), 0.01, 0, 1, 4, &ok) * particle_radius;
	
	if(ok)
	{
		emit signalDisturbParticles(max_dist);
	}
}

void MainWindow::showVisualizationOptions()
{
	visualizationOptionsWidget->show();
	visualizationOptionsWidget->activateWindow();
}

void MainWindow::showDisplayOptions()
{
	displayOptionsWidget->show();
	displayOptionsWidget->activateWindow();
}

void MainWindow::showAgglomerateInfo()
{
	sim_thread->updateAgglomerateInfo(&agg_info, (bool)agg_info_only_inner_particles);

	agglomerateInfoWidget->show();
	agglomerateInfoWidget->activateWindow();
}

void MainWindow::printContactHistogram()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Select File"), "", tr("*.txt"));

	if(!filename.isEmpty())
        sim_thread->printContactHistogram( filename.toLatin1().data(), (bool) agg_info_only_inner_particles );
}

void MainWindow::printFragmentVelocities()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Select File"), "", tr("*.txt"));

	if(!filename.isEmpty())
        sim_thread->printFragmentVelocities( filename.toLatin1().data() );
}

void MainWindow::printElasticCharge()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Select File"), "", tr("*.txt"));

	if(!filename.isEmpty())
        sim_thread->printElasticCharge( filename.toLatin1().data() );
}

void MainWindow::printFillingFactorProfile()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Select File"), "", tr("*.txt"));

	if(!filename.isEmpty())
        sim_thread->printFillingFactorProfile( filename.toLatin1().data() );
}

void MainWindow::printFillingFactorProfileSphere()
{
    QString filename = QFileDialog::getSaveFileName(this, tr("Select File"), "", tr("*.txt"));

    if(!filename.isEmpty())
        sim_thread->printFillingFactorProfileSphere( filename.toLatin1().data() );
}

void MainWindow::printFractalDimension()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Select File"), "", tr("*.txt"));

	if(!filename.isEmpty())
        sim_thread->printFractalDimension( filename.toLatin1().data() );
}

void MainWindow::showMaterialInfo()
{
	materialInfoWidget->show();
	materialInfoWidget->activateWindow();
}

void MainWindow::enterInfoText()
{
	bool ok;
	QString text = QInputDialog::getText(this, tr("Display Description"),  tr("Enter short description:"), QLineEdit::Normal, tr(""), &ok);

	if(ok)
	{
		info_text = text;
		glView->update();
	}
}

void MainWindow::showLoadWidget()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Load Particle File"), file_path, tr("Particle Files (*.dat)"));

	if(!filename.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(filename).path();

        emit signalLoadDataFromFile(  filename.toLatin1().data() );
	}
}



void MainWindow::showFolderToScreenshotsWidget()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Load Particle File"), file_path, tr("Particle Files (*.dat)"));

    if(!filename.isEmpty())
    {
        // store path for next time
        file_path = QFileInfo(filename).path();



        // find the index of the filename before the '.dat'
        char* file = filename.toLatin1().data();

        int i = 0;
        int filename_end_index = 0;
        while(file[i] != '\0')
        {
            if(filename[i] == '.')
                filename_end_index = i;

            ++i;
        }

        file[filename_end_index-1] = '\0';

        strcpy(screenshot_path, file);
        screenshot_counter = 0;
    }

    char newFile[1024];
    char newImg[1024];

    while(!filename.isEmpty())
    {

        sprintf(newFile, "%s%04i.dat", screenshot_path, screenshot_counter);
        sprintf(newImg, "%s%04i.jpg", screenshot_path, screenshot_counter);

        screenshot_counter++;

        printf("filename = %s\n", newFile);
        printf("imagename = %s", newImg);

        // load data
        centerView();

        emit signalLoadDataFromFile(newFile);

        centerView();

        printf("    number of particles %d\n", sim_thread->sim.number_of_particles);

        // get image
        glView->paintGL();
        QImage image = glView->grabFrameBuffer();
        image.save(newImg, 0, 100);

        wait_for_screenshot.wakeAll();

    }


}





void MainWindow::showPosFolderToScreenshotsWidget()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Load Particle File"), file_path, tr("Particle Files (*.dat)"));

    if(!filename.isEmpty())
    {
        // store path for next time
        file_path = QFileInfo(filename).path();



        // find the index of the filename before the '.dat'
        char* file = filename.toLatin1().data();

        int i = 0;
        int filename_end_index = 0;
        while(file[i] != '\0')
        {
            if(filename[i] == '.')
                filename_end_index = i;

            ++i;
        }

        file[filename_end_index-1] = '\0';

        strcpy(screenshot_path, file);
        screenshot_counter = 0;
    }

    char newFile[1024];
    char newImg[1024];

    while(!filename.isEmpty())
    {

        sprintf(newFile, "%s%i.dat", screenshot_path, screenshot_counter);
        sprintf(newImg, "%s%i.jpg", screenshot_path, screenshot_counter);

        screenshot_counter++;

        printf("filename = %s\n", newFile);
        printf("imagename = %s", newImg);

        centerView();

        // load data
        emit signalImportDataFromFile(newFile);

        //if(screenshot_counter == 0)
        centerView();

        printf("    number of particles %d\n", sim_thread->sim.number_of_particles);


        // get image
        glView->paintGL();
        QImage image = glView->grabFrameBuffer();
        image.save(newImg, 0, 100);

        wait_for_screenshot.wakeAll();

    }


}


void MainWindow::showImportWidget()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Import Particle File"), file_path, tr("Data Files (*.dat)"));

	if(!filename.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(filename).path();

        emit signalImportDataFromFile(  filename.toLatin1().data() );
	}
}

void MainWindow::showSaveWidget()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Particle File"), file_path, tr("Particle Files (*.dat)"));

	if(!filename.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(filename).path();
        emit signalSaveDataToFile(  filename.toLatin1().data() );
	}
}

void MainWindow::showSavePositionsWidget()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Position File"), file_path, tr("Position Files (*.pos)"));

	if(!filename.isEmpty())
	{
		// store path for next time
		file_path = QFileInfo(filename).path();
        emit signalSavePositionsToFile(  filename.toLatin1().data() );
	}
}

void MainWindow::enterSeed()
{
	bool ok;
    int seed = QInputDialog::getInt(this, tr("Seed"), tr("New Seed:"), 1337, 0, 2147483647, 1, &ok);
	
	if(ok)
	{
        sim_thread->sim.rand_generator.seed((unsigned)seed);
	}
}

void MainWindow::loadMaterial()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Load Material File"), "", tr("Material description (*.mat)"));

	if(!filename.isEmpty())
        emit signalLoadMaterial(  filename.toLatin1().data() );

	centerView();
}

void MainWindow::customScreenshot()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Save Image"), "", tr("JPG (*.jpg)"));

	if(!filename.isEmpty())
	{
		// get image
		QImage image = glView->grabFrameBuffer();
		image.save(filename, 0, 100);
	}
}

void MainWindow::autoScreenshot()
{
	char temp[200];
	sprintf(temp, "%simage_%i.jpg", screenshot_path, screenshot_counter);
	QString filename = temp;

	// get image
	QImage image = glView->grabFrameBuffer();
	image.save(filename, 0, 100);

	wait_for_screenshot.wakeAll();

	++screenshot_counter;
}

void MainWindow::displayCrossSection()
{
	double cross_section, sigma_cross_section;
	sim_thread->getCrossSection(&cross_section, &sigma_cross_section);

	char buffer[200];
    sprintf(buffer, "Average projected cross section: %4.2le m^2\nDeviation of cross section:           %4.2le m^2", cross_section, sigma_cross_section);
	QMessageBox::about(this, tr("Cross Section"), buffer);
}

void MainWindow::initChain()
{
	if(initChainWidget->impactOnCurrentAgglomerateCheckBox->isChecked())
		emit signalImpactChainOnAgglomerate(initChainWidget->particleImpactNumberEdit->value(), initChainWidget->angularIrregularityEdit->value(), initChainWidget->impactSpeedEdit->value());
	else
		emit signalInitSimChain(initChainWidget->particleNumberEdit->value(), initChainWidget->particleImpactNumberEdit->value(), initChainWidget->targetParticleEdit->value()-1, initChainWidget->angularIrregularityEdit->value(), initChainWidget->impactSpeedEdit->value());
}

void MainWindow::showErrorMessage(ErrorCode error_code)
{
	char message[200];
	Simulation::getErrorMessage(error_code, message);

	QMessageBox::about(this, tr("Error"), message);
}

void MainWindow::updateVisualization()
{
	if(follow_cms)
		followCMS();

	glView->updateVertexBuffers();
}

void MainWindow::updateAgglomerateInfo()
{
	sim_thread->updateAgglomerateInfo(&agg_info, (bool) agg_info_only_inner_particles);
	agglomerateInfoWidget->repaint();
}

void MainWindow::centerView()
{
	if(sim_thread->number_of_particles <= 0)
		return;

	// get CMS
	float cms[3];
	sim_thread->getCMS(cms);

	// set cam
    glView->setCameraPosition(Vector<GLdouble, 3>(3, cms[0], cms[1], cms[2] + 350.0f * particle_radius));
	glView->setCameraLookAt(Vector<GLdouble, 3>(3, cms[0], cms[1], cms[2]));
	glView->setCameraUp(Vector<GLdouble, 3>(3, 0.0, 1.0, 0.0));

	glView->update();
}

void MainWindow::followCMS()
{
	if(sim_thread->number_of_particles <= 0)
		return;

	// get new cms
	float new_cms[3];
	sim_thread->getCMS(new_cms);

	// calculate new cam pos
	Vector<GLdouble, 3> cam_pos = glView->getCameraPosition();
	cam_pos += Vector<GLdouble, 3>(3, new_cms[0]-cms[0], new_cms[1]-cms[1], new_cms[2]-cms[2]);
	glView->setCameraPosition(cam_pos);

	// calculate new look at pos
	Vector<GLdouble, 3> look_at_pos = glView->getCameraLookAt();
	look_at_pos += Vector<GLdouble, 3>(3, new_cms[0]-cms[0], new_cms[1]-cms[1], new_cms[2]-cms[2]);
	glView->setCameraLookAt(look_at_pos);

	cms[0] = new_cms[0];
	cms[1] = new_cms[1];
	cms[2] = new_cms[2];

	glView->update();
}


#ifdef DRAW_CONTACT_POINTERS
void MainWindow::contactPointerUpdateDone()
{
	if(new_contact_pointers == &contact_pointers1)
	{
		current_contact_pointers = &contact_pointers1;
		new_contact_pointers = &contact_pointers2;
	}
	else
	{
		current_contact_pointers = &contact_pointers2;
		new_contact_pointers = &contact_pointers1;
	}

	updating_contact_pointers = false;

	glView->update();
}
#endif

void MainWindow::numberOfParticlesChanged(int number_of_particles)
{
	// en-/disable action dependent on present data
	if(number_of_particles > 0)
		centerViewAct->setDisabled(false);
	else
		centerViewAct->setDisabled(true);

	sim_thread->getCMS(cms);
}

void MainWindow::keyPressEvent(QKeyEvent *event) 
{
	if(event->key() == Qt::Key_P)
	{
		if(pause)
		{
			pause = false;
			emit signalContinueSimulation();
		}
		else
		{
			pause = true;
			emit signalPauseSimulation();
		}
	}
	

	QMainWindow::keyPressEvent(event);
}
