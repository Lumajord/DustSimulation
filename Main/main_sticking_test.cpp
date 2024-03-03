#include <stdlib.h>਍⌀椀渀挀氀甀搀攀 㰀猀琀搀椀漀⸀栀㸀ഀഀ
਍⼀⼀⌀搀攀昀椀渀攀 匀唀刀䘀䄀䌀䔀开䄀匀倀䔀刀䤀吀䤀䔀匀ഀഀ
#ifdef SURFACE_ASPERITIES਍ ⌀甀渀搀攀昀 匀唀刀䘀䄀䌀䔀开䄀匀倀䔀刀䤀吀䤀䔀匀ഀഀ
#endif਍ഀഀ
#include "Simulation.h"਍⌀椀渀挀氀甀搀攀 ∀匀椀洀甀氀愀琀椀漀渀䰀椀戀⸀栀∀ഀഀ
਍攀砀琀攀爀渀 搀漀甀戀氀攀 瀀愀爀琀椀挀氀攀开爀愀搀椀甀猀㬀ഀഀ
extern double delta_0;਍攀砀琀攀爀渀 搀漀甀戀氀攀 搀攀渀猀椀琀礀㬀ഀഀ
extern double surface_energy;਍攀砀琀攀爀渀 搀漀甀戀氀攀 渀甀㬀ഀഀ
extern double young_mod;਍攀砀琀攀爀渀 搀漀甀戀氀攀 挀爀椀琀开爀漀氀氀椀渀最开搀椀猀瀀氀愀挀攀洀攀渀琀㬀ഀഀ
extern double yield_strength;਍攀砀琀攀爀渀 搀漀甀戀氀攀 爀漀氀氀椀渀最开洀漀搀椀昀椀攀爀㬀ഀഀ
extern double sliding_modifier;਍攀砀琀攀爀渀 搀漀甀戀氀攀 琀眀椀猀琀椀渀最开洀漀搀椀昀椀攀爀㬀ഀഀ
਍戀漀漀氀 挀栀攀挀欀倀愀爀琀椀挀氀攀倀愀爀琀椀挀氀攀䈀漀甀渀挀椀渀最⠀匀椀洀甀氀愀琀椀漀渀 ⨀猀椀洀Ⰰ 搀漀甀戀氀攀 琀椀洀攀猀琀攀瀀Ⰰ 搀漀甀戀氀攀 挀漀氀氀椀猀椀漀渀开猀瀀攀攀搀⤀ഀഀ
{਍ऀ⼀⼀ 爀攀猀攀琀 猀椀洀ഀഀ
	sim->resizeArrays(2, 0);਍ഀഀ
	memset(sim->pos_old, 0, 8 * sizeof(double));਍ऀ洀攀洀猀攀琀⠀猀椀洀ⴀ㸀瘀攀氀Ⰰ 　Ⰰ 㠀 ⨀ 猀椀稀攀漀昀⠀搀漀甀戀氀攀⤀⤀㬀ഀഀ
	memset(sim->vel_angular, 0, 8 * sizeof(double));਍ऀ猀椀洀ⴀ㸀瀀漀猀开漀氀搀嬀㐀崀 㴀 ㈀⸀　㄀ ⨀ 瀀愀爀琀椀挀氀攀开爀愀搀椀甀猀㬀ഀഀ
	sim->vel[0] = collision_speed;਍ഀഀ
	sim->initSimData(SIM_TYPE_GENERAL);਍ऀ猀椀洀ⴀ㸀猀琀愀爀琀匀椀洀甀氀愀琀椀漀渀⠀　⸀㄀ ⨀ 瀀愀爀琀椀挀氀攀开爀愀搀椀甀猀 ⼀ 挀漀氀氀椀猀椀漀渀开猀瀀攀攀搀Ⰰ 琀椀洀攀猀琀攀瀀⤀㬀ഀഀ
਍ऀ⼀⼀ 爀甀渀 猀椀洀ഀഀ
	while(!sim->stop_simulation)਍ऀऀ猀椀洀ⴀ㸀甀瀀搀愀琀攀⠀⤀㬀ഀഀ
਍ऀ⼀⼀ 挀栀攀挀欀 椀昀 戀漀甀渀挀椀渀最 漀挀挀甀爀攀搀ഀഀ
	return (SimLib::getNumberOfContacts(sim) == 0);਍紀ഀഀ
਍搀漀甀戀氀攀 最攀琀倀愀爀琀椀挀氀攀倀愀爀琀椀挀氀攀䈀漀甀渀挀椀渀最匀瀀攀攀搀⠀挀漀渀猀琀 挀栀愀爀⨀ 洀愀琀攀爀椀愀氀开昀椀氀攀渀愀洀攀Ⰰ 搀漀甀戀氀攀 爀愀搀椀甀猀Ⰰ 搀漀甀戀氀攀 琀椀洀攀猀琀攀瀀Ⰰ 搀漀甀戀氀攀 猀琀愀爀琀开猀瀀攀攀搀Ⰰ 搀漀甀戀氀攀 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀Ⰰ 椀渀琀 洀愀砀开爀甀渀猀⤀ഀഀ
{਍ऀ椀昀⠀猀琀愀爀琀开猀瀀攀攀搀 㰀㴀 　⤀ഀഀ
		return 0;਍ഀഀ
	Simulation sim;਍ऀ猀椀洀⸀氀漀愀搀䴀愀琀攀爀椀愀氀⠀洀愀琀攀爀椀愀氀开昀椀氀攀渀愀洀攀⤀㬀ഀഀ
	sim.setMaterialConstants(radius, density, surface_energy, nu, young_mod, crit_rolling_displacement, yield_strength, rolling_modifier, sliding_modifier, twisting_modifier);਍ഀഀ
	//printf("Particle<->Particle bouncing:\nCurrent speed:\n");਍ഀഀ
	// check if bouncing occurs at max impact speed (if not, subsequent test at lower impact speeds are useless)਍ऀ椀昀⠀℀挀栀攀挀欀倀愀爀琀椀挀氀攀倀愀爀琀椀挀氀攀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 猀琀愀爀琀开猀瀀攀攀搀 ⬀ ⠀搀漀甀戀氀攀⤀⠀洀愀砀开爀甀渀猀ⴀ㄀⤀ ⨀ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀⤀⤀ഀഀ
		return 0;਍ऀഀഀ
	// determine lower bound for search਍ऀ椀渀琀 爀甀渀 㴀 　㬀ഀഀ
਍ऀ椀昀⠀℀挀栀攀挀欀倀愀爀琀椀挀氀攀倀愀爀琀椀挀氀攀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 猀琀愀爀琀开猀瀀攀攀搀 ⬀ ⠀搀漀甀戀氀攀⤀⠀洀愀砀开爀甀渀猀⼀㐀⤀ ⨀ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀⤀⤀ഀഀ
		run = max_runs/4;਍ഀഀ
	if(!checkParticleParticleBouncing(&sim, timestep, start_speed + (double)(max_runs/3) * speed_interval))਍ऀऀ爀甀渀 㴀 洀愀砀开爀甀渀猀⼀㌀㬀ഀഀ
਍ऀ椀昀⠀℀挀栀攀挀欀倀愀爀琀椀挀氀攀倀愀爀琀椀挀氀攀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 猀琀愀爀琀开猀瀀攀攀搀 ⬀ ⠀搀漀甀戀氀攀⤀⠀洀愀砀开爀甀渀猀⤀ ⨀ 　⸀㜀㔀 ⨀ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀⤀⤀ഀഀ
		run = (int)((double)max_runs * 0.75)  - 1;਍ഀഀ
	while(run < max_runs - 1)਍ऀ笀ഀഀ
		double current_speed = start_speed + (double)run * speed_interval;਍ऀऀ⼀⼀瀀爀椀渀琀昀⠀∀尀戀尀戀尀戀尀戀尀戀尀戀─⸀㈀氀昀∀Ⰰ 挀甀爀爀攀渀琀开猀瀀攀攀搀⤀㬀ഀഀ
਍ऀऀ椀昀⠀挀栀攀挀欀倀愀爀琀椀挀氀攀倀愀爀琀椀挀氀攀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 挀甀爀爀攀渀琀开猀瀀攀攀搀⤀⤀ഀഀ
			return current_speed;਍ഀഀ
		++run;਍ऀ紀ഀഀ
਍ऀ爀攀琀甀爀渀 　㬀ഀഀ
}਍ഀഀ
bool checkParticleWallBouncing(Simulation *sim, double timestep, double collision_speed)਍笀ഀഀ
	sim->resizeArrays(1, 1);਍ऀ洀攀洀猀攀琀⠀猀椀洀ⴀ㸀瀀漀猀开漀氀搀Ⰰ 　Ⰰ 㐀 ⨀ 猀椀稀攀漀昀⠀搀漀甀戀氀攀⤀⤀㬀ഀഀ
	memset(sim->vel, 0, 4 * sizeof(double));਍ऀ洀攀洀猀攀琀⠀猀椀洀ⴀ㸀瘀攀氀开愀渀最甀氀愀爀Ⰰ 　Ⰰ 㐀 ⨀ 猀椀稀攀漀昀⠀搀漀甀戀氀攀⤀⤀㬀ഀഀ
	sim->pos_old[0] = -2.01 * particle_radius;਍ऀ猀椀洀ⴀ㸀瘀攀氀嬀　崀 㴀 挀漀氀氀椀猀椀漀渀开猀瀀攀攀搀㬀ഀഀ
਍ऀ猀椀洀ⴀ㸀眀愀氀氀猀嬀　崀⸀渀漀爀洀愀氀嬀　崀 㴀 ⴀ㄀⸀　㬀ഀഀ
	sim->walls[0].normal[1] = 0.0;਍ऀ猀椀洀ⴀ㸀眀愀氀氀猀嬀　崀⸀渀漀爀洀愀氀嬀㈀崀 㴀 　⸀　㬀ഀഀ
	sim->walls[0].pos[0] = 0.0;਍ऀ猀椀洀ⴀ㸀眀愀氀氀猀嬀　崀⸀瀀漀猀嬀㄀崀 㴀 　⸀　㬀 ഀഀ
	sim->walls[0].pos[2] = 0.0; ਍ഀഀ
	sim->initSimData(SIM_TYPE_GENERAL);਍ऀ猀椀洀ⴀ㸀猀琀愀爀琀匀椀洀甀氀愀琀椀漀渀⠀　⸀㄀ ⨀ 瀀愀爀琀椀挀氀攀开爀愀搀椀甀猀 ⼀ 挀漀氀氀椀猀椀漀渀开猀瀀攀攀搀Ⰰ 琀椀洀攀猀琀攀瀀⤀㬀ഀഀ
਍ऀ⼀⼀ 爀甀渀 猀椀洀ഀഀ
	while(!sim->stop_simulation)਍ऀऀ猀椀洀ⴀ㸀甀瀀搀愀琀攀⠀⤀㬀ഀഀ
਍ऀ⼀⼀ 挀栀攀挀欀 椀昀 戀漀甀渀挀椀渀最 漀挀挀甀爀攀搀ഀഀ
	return (SimLib::getNumberOfContacts(sim) == 0);਍紀ഀഀ
਍搀漀甀戀氀攀 最攀琀倀愀爀琀椀挀氀攀圀愀氀氀䈀漀甀渀挀椀渀最匀瀀攀攀搀⠀挀漀渀猀琀 挀栀愀爀⨀ 洀愀琀攀爀椀愀氀开昀椀氀攀渀愀洀攀Ⰰ 搀漀甀戀氀攀 爀愀搀椀甀猀Ⰰ 搀漀甀戀氀攀 琀椀洀攀猀琀攀瀀Ⰰ 搀漀甀戀氀攀 猀琀愀爀琀开猀瀀攀攀搀Ⰰ 搀漀甀戀氀攀 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀Ⰰ 椀渀琀 洀愀砀开爀甀渀猀⤀ഀഀ
{਍ऀ椀昀⠀猀琀愀爀琀开猀瀀攀攀搀 㰀㴀 　⤀ഀഀ
		return 0;਍ഀഀ
	Simulation sim;਍ऀ猀椀洀⸀氀漀愀搀䴀愀琀攀爀椀愀氀⠀洀愀琀攀爀椀愀氀开昀椀氀攀渀愀洀攀⤀㬀ഀഀ
	sim.setMaterialConstants(radius, density, surface_energy, nu, young_mod, crit_rolling_displacement, yield_strength, rolling_modifier, sliding_modifier, twisting_modifier);਍ഀഀ
	//printf("Particle<->Particle bouncing:\nCurrent speed:\n");਍ഀഀ
	// check if bouncing occurs at max impact speed (if not, subsequent test at lower impact speeds are useless)਍ऀ椀昀⠀℀挀栀攀挀欀倀愀爀琀椀挀氀攀圀愀氀氀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 猀琀愀爀琀开猀瀀攀攀搀 ⬀ ⠀搀漀甀戀氀攀⤀⠀洀愀砀开爀甀渀猀ⴀ㄀⤀ ⨀ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀⤀⤀ഀഀ
		return 0;਍ऀഀഀ
	// determine lower bound for search਍ऀ椀渀琀 爀甀渀 㴀 　㬀ഀഀ
਍ऀ椀昀⠀℀挀栀攀挀欀倀愀爀琀椀挀氀攀圀愀氀氀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 猀琀愀爀琀开猀瀀攀攀搀 ⬀ ⠀搀漀甀戀氀攀⤀⠀洀愀砀开爀甀渀猀⼀㐀⤀ ⨀ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀⤀⤀ഀഀ
		run = max_runs/4;਍ഀഀ
	if(!checkParticleWallBouncing(&sim, timestep, start_speed + (double)(max_runs/3) * speed_interval))਍ऀऀ爀甀渀 㴀 洀愀砀开爀甀渀猀⼀㌀㬀ഀഀ
਍ऀ椀昀⠀℀挀栀攀挀欀倀愀爀琀椀挀氀攀圀愀氀氀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 猀琀愀爀琀开猀瀀攀攀搀 ⬀ ⠀搀漀甀戀氀攀⤀⠀洀愀砀开爀甀渀猀⤀ ⨀ 　⸀㜀㔀 ⨀ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀⤀⤀ഀഀ
		run = (int)((double)max_runs * 0.75)  - 1;਍ഀഀ
	while(run < max_runs)਍ऀ笀ഀഀ
		double current_speed = start_speed + (double)run * speed_interval;਍ऀऀ⼀⼀瀀爀椀渀琀昀⠀∀尀戀尀戀尀戀尀戀尀戀尀戀─⸀㈀氀昀∀Ⰰ 挀甀爀爀攀渀琀开猀瀀攀攀搀⤀㬀ഀഀ
		਍ऀऀ椀昀⠀挀栀攀挀欀倀愀爀琀椀挀氀攀圀愀氀氀䈀漀甀渀挀椀渀最⠀☀猀椀洀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 挀甀爀爀攀渀琀开猀瀀攀攀搀⤀⤀ഀഀ
			return current_speed;਍ഀഀ
		++run;਍ऀ紀ഀഀ
਍ऀ爀攀琀甀爀渀 　㬀ഀഀ
}਍ഀഀ
int main(int argc, char **argv)਍笀ഀഀ
	double timestep;਍ऀ搀漀甀戀氀攀 猀琀愀爀琀开猀瀀攀攀搀㬀ഀഀ
	double wall_start_speed;਍ऀ搀漀甀戀氀攀 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀㬀ഀഀ
	int max_runs;਍ऀ搀漀甀戀氀攀 洀椀渀开爀愀搀椀甀猀㬀ഀഀ
	double radius_increment;਍ऀ椀渀琀 爀愀搀椀甀猀开椀渀琀攀爀瘀愀氀猀㬀ഀഀ
	int run = 0;਍ഀഀ
	if(argc == 11)਍ऀ笀ഀഀ
		timestep = atof(argv[1]);਍ऀऀ猀琀愀爀琀开猀瀀攀攀搀 㴀 愀琀漀昀⠀愀爀最瘀嬀㈀崀⤀㬀ഀഀ
		wall_start_speed = atof(argv[3]);਍ऀऀ猀瀀攀攀搀开椀渀琀攀爀瘀愀氀 㴀 愀琀漀昀⠀愀爀最瘀嬀㐀崀⤀㬀ഀഀ
		max_runs = atoi(argv[5]);਍ऀऀ洀椀渀开爀愀搀椀甀猀 㴀 愀琀漀昀⠀愀爀最瘀嬀㜀崀⤀㬀ഀഀ
		radius_intervals = atoi(argv[9]);਍ऀऀ爀愀搀椀甀猀开椀渀挀爀攀洀攀渀琀 㴀 ⠀愀琀漀昀⠀愀爀最瘀嬀㠀崀⤀ⴀ愀琀漀昀⠀愀爀最瘀嬀㜀崀⤀⤀ ⼀ ⠀搀漀甀戀氀攀⤀爀愀搀椀甀猀开椀渀琀攀爀瘀愀氀猀㬀ഀഀ
	}਍ऀ攀氀猀攀ഀഀ
	{਍ऀऀ瀀爀椀渀琀昀⠀∀䤀渀猀甀昀昀椀挀椀攀渀琀 愀爀最甀洀攀渀琀猀 ⴀ 甀猀攀㨀尀渀∀⤀㬀ഀഀ
		printf("-timestep -start_speed -wall_start_speed -speed_interval -max_runs -material_filename -radius_min -radius_max -radius_intervals -output_file\n");਍ऀऀ爀攀琀甀爀渀 　㬀ഀഀ
	}਍ഀഀ
	// prepare log file਍ऀ䘀䤀䰀䔀 ⨀氀漀最开昀椀氀攀 㴀 昀漀瀀攀渀⠀愀爀最瘀嬀㄀　崀Ⰰ ∀眀⬀∀⤀㬀ഀഀ
਍ऀ椀昀⠀氀漀最开昀椀氀攀⤀ഀഀ
	{਍⌀椀昀搀攀昀  匀唀刀䘀䄀䌀䔀开䄀匀倀䔀刀䤀吀䤀䔀匀ഀഀ
		fprintf(log_file, "# Surface asperities enabled\n#\n# particle radius    particle<->particle sticking velocity    particle<->wall sticking velocity\n");਍⌀攀氀猀攀ഀഀ
		fprintf(log_file, "# Surface asperities disabled\n#\n# particle radius    particle<->particle sticking velocity    particle<->wall sticking velocity\n");਍⌀攀渀搀椀昀ഀഀ
	}਍ऀ攀氀猀攀ഀഀ
		printf("ERROR: Could not open %s!\n", argv[6]);਍ഀഀ
	/////////////////////////////////////////////////////////////////////////////////////////////////////// ਍ऀ⼀⼀ 瀀攀爀昀漀爀洀 戀漀甀渀挀椀渀最 琀攀猀琀猀ഀഀ
	///////////////////////////////////////////////////////////////////////////////////////////////////////਍ഀഀ
	for(int i = 0; i < radius_intervals; ++i)਍ऀ笀ഀഀ
		printf("Radius: %g\n", min_radius + (double)i * radius_increment);਍ഀഀ
		double critical_speed = getParticleParticleBouncingSpeed(argv[6], min_radius + (double)i * radius_increment, timestep, start_speed, speed_interval, max_runs);਍ऀऀ搀漀甀戀氀攀 眀愀氀氀开挀爀椀琀椀挀愀氀开猀瀀攀攀搀 㴀 最攀琀倀愀爀琀椀挀氀攀圀愀氀氀䈀漀甀渀挀椀渀最匀瀀攀攀搀⠀愀爀最瘀嬀㘀崀Ⰰ 洀椀渀开爀愀搀椀甀猀 ⬀ ⠀搀漀甀戀氀攀⤀椀 ⨀ 爀愀搀椀甀猀开椀渀挀爀攀洀攀渀琀Ⰰ 琀椀洀攀猀琀攀瀀Ⰰ 眀愀氀氀开猀琀愀爀琀开猀瀀攀攀搀Ⰰ 猀瀀攀攀搀开椀渀琀攀爀瘀愀氀Ⰰ 洀愀砀开爀甀渀猀⤀㬀ഀഀ
਍ऀऀ椀昀⠀猀琀愀爀琀开猀瀀攀攀搀 㸀 　⤀ഀഀ
		{਍ऀऀऀ椀昀⠀挀爀椀琀椀挀愀氀开猀瀀攀攀搀 㸀 　⤀ഀഀ
				printf("Bouncing observed at %lf m/s\n\n", critical_speed);਍ऀऀऀ攀氀猀攀ഀഀ
				printf("No bouncing observed between %lf and %lf m/s\n\n", start_speed, start_speed + (double)(max_runs-1) * speed_interval);਍ऀऀ紀ഀഀ
਍ऀऀ椀昀⠀眀愀氀氀开猀琀愀爀琀开猀瀀攀攀搀 㸀 　⤀ഀഀ
		{਍ऀऀऀ椀昀⠀挀爀椀琀椀挀愀氀开猀瀀攀攀搀 㸀 　⤀ഀഀ
				printf("\nWall bouncing observed at %lf m/s\n", wall_critical_speed);਍ऀऀऀ攀氀猀攀ഀഀ
				printf("\nNo wall bouncing observed between %lf and %lf m/s\n", start_speed, wall_start_speed + (double)(max_runs-1) * speed_interval);਍ऀऀ紀ഀഀ
਍ऀऀ椀昀⠀氀漀最开昀椀氀攀⤀ഀഀ
			fprintf(log_file, "%g %g %g\n", particle_radius, critical_speed, wall_critical_speed);਍ऀ紀ഀഀ
਍ऀ椀昀⠀氀漀最开昀椀氀攀⤀ഀഀ
		fclose(log_file);਍ഀഀ
    return 0;਍紀ഀഀ
