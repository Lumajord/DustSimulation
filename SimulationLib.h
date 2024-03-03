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
 *   Free Software Foundation, Inc.,                                        *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              *
 ***************************************************************************/

#ifndef SIMULATIONLIB_H
#define SIMULATIONLIB_H

#include <list>
#include <vector>

#include "Constants.h"
#include "Simulation.h"

namespace SimLib
{
ErrorCode disturbParticles(Simulation *sim, double max_dist);

ErrorCode duplicateSlices(Simulation *sim, int x_duplications, int y_duplications, int z_duplications, bool mirror, bool random_orientation);

ErrorCode particleClusterAggregation(Simulation *sim, int final_cluster_size, double fractal_prefactor, double fractal_dimension);

bool getNeighbourPos(vec3 *pos, Simulation &sim, int particle_id);

bool getMigrationPos1(vec3 *migration_pos, Simulation &sim, int particle_id, vec3 &default_pos, BAMSelectionMethod bam_selection_method);

bool getMigrationPos2(vec3 *migration_pos, Simulation &sim, int particle_id, vec3 &default_pos, BAMSelectionMethod bam_selection_method);

ErrorCode initBAMAggregate(Simulation *sim, const char *filename, unsigned int number_of_particles, double migration_probability1, double migration_probability2, BAMSelectionMethod bam_selection_method);

ErrorCode initBAMCake(Simulation *sim, unsigned int number_of_particles, double x_size, double y_size, double migration_probability1, double migration_probability2);

ErrorCode addFractalChainsToAggregate(Simulation *sim, const char *filename, unsigned int number_of_chains, unsigned int min_chain_length, unsigned int max_chain_length, double migration_probability);

ErrorCode initFractalAggregate(Simulation *sim, unsigned int number_of_particles, double migration_probability1, double migration_probability2, double chain_ratio, unsigned int min_chain_length, unsigned int max_chain_length);

ErrorCode initTreeAggregate(Simulation *sim, unsigned int number_of_particles, double contact_distribution[6], double migration_rate1, double migration_rate2);

ErrorCode buildSimpleCubicGrid(Simulation *sim, double x_size, double y_size, double z_size, double filling_factor);

ErrorCode buildHexagonalGrid(Simulation *sim, int x_particles, int y_particles, int z_particles, double filling_factor);

double getLowerXPos(Simulation *sim, double pos_y, double pos_z, double x_max, double x_min);

double getYPos(Simulation *sim, double pos_x, double pos_z, double y_max, double y_min);

ErrorCode initRandomBallisticDeposition(Simulation *sim, int number_of_particles, double x_size, double y_size, double dest_filling_factor = 0);

ErrorCode reduceBoxFillingFactor(Simulation *sim, double dest_filling_factor, bool side_walls);

ErrorCode initNeedleDeposition(Simulation *sim, int number_of_particles);

ErrorCode initCylindricalRBD(Simulation *sim, int number_of_particles, double radius, double slice_factor);

ErrorCode initChain(Simulation *sim, int number_of_chain_particles, int number_of_impact_particles, int target_id, double angular_irregularity, double impact_speed);

ErrorCode initImpactChainOnAgglomerate(Simulation *sim, int number_of_chain_particles, double angular_irregularity, double impact_speed, bool random_orientation);

ErrorCode initCluster(Simulation *sim, int number_of_cluster_particles, int number_of_impact_particles, double filling_factor, double angular_irregularity, double impact_speed);

ErrorCode initChainBox(Simulation *sim, int number_of_particles, double filling_factor, double angular_irregularity);

ErrorCode collideAgglomerateWithWall(Simulation *sim, const char *filename, double impact_speed, int impact_angle, double impact_distance, bool random_orientation, double wall_rolling_modifier = 1.0, double wall_sliding_modifier = 1.0, double wall_twisting_modifier = 1.0);

ErrorCode initBox(Simulation *sim, const char *filename, bool side_walls, double side_wall_height_modifier, bool top_wall, double top_wall_speed, double dynamic_pressure);

ErrorCode setWallInteractionModifiers(Simulation *sim, double side_wall_compression_modifier, double side_wall_rolling_modifier, double side_wall_sliding_modifier, double top_wall_compression_modifier, double top_wall_rolling_modifier, double top_wall_sliding_modifier);

ErrorCode initCompressionBox(Simulation *sim, const char *filename, bool side_walls, bool moving_bottom_wall, double wall_speed, double stop_filling_factor, double side_wall_compression_modifier = 1.0, double side_wall_rolling_modifier = 1.0, double side_wall_sliding_modifier = 1.0, double top_wall_compression_modifier = 1.0, double top_wall_rolling_modifier = 1.0, double top_wall_sliding_modifier = 1.0);

ErrorCode initCompressionRelaxationBox(Simulation *sim, const char *filename, double wall_speed, double stop_filling_factor, double stop_dissipation_factor, double side_wall_compression_modifier = 1.0, double side_wall_rolling_modifier = 1.0, double side_wall_sliding_modifier = 1.0);

ErrorCode initShockwaveBox(Simulation *sim, const char *filename, double compression_speed, double perturbation_layer_thickness, double stop_pressure, double side_wall_rolling_modifier = 1.0, double side_wall_sliding_modifier = 1.0);

ErrorCode initDynamicCompressionBox(Simulation *sim, const char *filename, double wall_speed, double wall_thickness, double side_wall_compression_modifier = 1.0, double side_wall_rolling_modifier = 1.0, double side_wall_sliding_modifier = 1.0);

ErrorCode initOpenBox(Simulation *sim, const char* filename, bool side_walls, double side_wall_height_modifier = 1.0, double side_wall_compression_modifier = 1.0, double side_wall_rolling_modifier = 1.0, double side_wall_sliding_modifier = 1.0);

ErrorCode initCompactionBox(Simulation *sim, const char *filename, double wall_speed, double stop_filling_factor);

ErrorCode initPullStrengthTest(Simulation *sim, double pull_speed);

ErrorCode initShearStrengthTest(Simulation *sim, double pull_speed);

ErrorCode collideTwoAgglomerates(Simulation *sim, const char *filename1, const char *filename2, double impact_speed, double impact_parameter, double impact_distance, bool random_orientation1, bool random_orientation2, bool smart_distance, bool cms, int seed = -1);

ErrorCode collideTwoAgglomerates(Simulation *sim,
        const char *material_name,
        const char *file1,
        const char *file2,
        const unsigned int seed1,
        const unsigned int seed2,
        const double impact_speed,
        const double impact_parameter,
        const double impact_distance = 0.0,
        const bool smart_distance = true,
        const bool cms = true
        );

ErrorCode generateBAMsphereAgglomerate(Simulation *sim,
        const char *material_name,
        const char *file1,
        double agg_size,
        const double agg_mass,
        const unsigned int seed1);

ErrorCode generateRBDSphere(
        Simulation *sim,
        const char *material_name,
        const double agg_mass,
        const unsigned int seed1
        );

ErrorCode generateBAMSphere(
        Simulation *sim,
        const char *material_name,
        const double agg_mass,
        const unsigned int seed1
        );

ErrorCode impactMultiProjectiles(Simulation *sim, const char *target_filename, bool random_target_orientation, const char *projectile_filename, unsigned int projectile_count, unsigned int projectile_size, double impact_velocity, double impact_parameter, bool random_impact_parameter, double min_distance, double max_distance, bool smart_distance);

ErrorCode hitAndStick(Simulation *sim, const char *filename1, const char *filename2, double impact_parameter, bool random_orientation1, bool random_orientation2);

ErrorCode hitAndStick(Simulation *sim, Simulation *sim2, double impact_parameter, bool random_orientation1, bool random_orientation2);

ErrorCode sandblastAggregate(Simulation *sim, const char *filename, int number_of_projectiles, double impact_velocity, double impact_distance);

void getCenterOfMass(vec3 *cms, const Simulation &sim, int first_particle, int last_particle);

void centerCMS(Simulation *sim);

void rotateSimRandomly(Simulation *sim);

void rotateParticles(Simulation *sim, vec3 &axis, double angle, int first_particle, int last_particle);

void rotateParticles(Simulation *sim, double angle1, double angle2, double angle3, int first_particle, int last_particle);

void detectFragments(Simulation &sim, std::vector<int> *fragment_id, std::vector<int> *size_of_fragment, std::vector< std::list<int> > *particles_of_fragment = NULL);

int detectFragments(Simulation &sim, int *fragment_ids);

void filterFragments(Simulation *sim);

void removeFragments(Simulation *sim, int remaining_fragment_id, std::vector<int> &fragment_ids, std::vector<int> &size_of_fragment);

void sliceBox(Simulation *sim, vec3 &pos_lower, vec3 &pos_upper);

void sliceSphere(Simulation *sim, vec3 &center, double slice_radius);

void sliceCylinder(Simulation *sim, vec3 &center, double slice_radius);

ErrorCode sliceTop(Simulation *sim, double top_slice_factor);

ErrorCode sliceBottom(Simulation *sim, double bottom_slice_factor);

void getSize(const Simulation &sim, double *gyration_radius, double *outer_radius);
//void accelerate_inward(const Simulation &sim, double speed);

void getCrossSection(Simulation &sim, double cell_size, unsigned int rotations, double *cross_section, double *sigma_cross_section, double radius_cutoff_factor = 1.01);
void getCrossSectionNoRotation(const Simulation &sim, const double cell_size, double &cross_section);

bool inContactWithWall(Simulation *sim, int particle_id);

void printContactList(Simulation &sim, const char* filename);

void printForces(Simulation &sim, const char *filename, bool new_forces);

void printContactPointers(Simulation &sim, const char *filename);

void printContactNormals(Simulation &sim, const char *filename);

void printSimData(Simulation &sim, const char *filename);

void printPositions(Simulation &sim, const char *filename, bool new_positions = false, bool orientation = true);

void getContactHistogram(Simulation &sim, int *particles, bool only_inner_particles = false);

ErrorCode printContactHistogram(Simulation &sim, const char *filename, bool only_inner_particles = false);

ErrorCode printFragmentVelocities(Simulation &sim, const char *filename);

ErrorCode printElasticCharge(Simulation &sim, FILE *file);

double getElasticCharge(Simulation &sim, int *charged_contacts);

double getKineticEnergy(Simulation &sim);

// prints various kinds of energies to the specified file
void printEnergy(Simulation *sim, const char *filename);

double getTotalEnergy(Simulation* sim);

#if defined(TRACK_DISSIPATED_ENERGY_PER_PARTICLE) && defined(TRACK_CONTACTS_PER_PARTICLE)
void saveParticleEnergyDissipation(Simulation &sim, char* filename);
#endif

void getOrthoVector(Simulation *sim, vec3 *ortho_vec, vec3 &vec);

double getMaxSpeed(Simulation *sim);

double getFillingFactorOfBox(Simulation &sim, vec3 &lower_pos, vec3 &upper_pos);

double getFillingFactorOfSlice(Simulation &sim, double y_lower, double y_upper, vec3 &lower_pos, vec3 &upper_pos);

void getFillingFactorOfSphere(Simulation &sim, double r_start, double r_end, double &filling_fact, double &coord_number);

void getParticlesInSlice(Simulation &sim, vec3 &lower_pos, vec3 &upper_pos, int *particles = NULL, int *contacts = NULL);

void printFillingFactorProfile(Simulation &sim, const char *filename, double slice_height);

void printFillingFactorProfileSphere(Simulation &sim, const char *filename, double layer_width);


ErrorCode printFractalDimension(Simulation &sim, const char *filename, double radius_increment, double inner_cutoff = 0.0, double outer_cutoff = 1.0);

void checkEquiblibriumDistances(Simulation &sim, const char *filename);

CollisionResult getCollisionResult(int size_largest_fragment, int size_second_largest_fragment, int number_of_particles, double sticking_threshhold = 0.6);

ErrorCode importParticlePositions(Simulation *sim, const char *filename);

#ifdef ENABLE_FIXED_PARTICLES
ErrorCode fixateParticles(Simulation *sim, int first_id, int last_id);
#endif

ErrorCode resetContacts(Simulation *sim);
}

#endif
