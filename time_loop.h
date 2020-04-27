//!  @file time_loop.h
#ifndef TIME_LOOP_H
#define TIME_LOOP_H

#include "constants.h"
#include "rate_of_change.h"
#include "time_integration.h"
#include "output.h"
#include <string.h>

/** 	
*		@brief Compute dt for the next iteration using CFL criterion
*	
*		@param all_particle pointer to an array containing information of all the particles
*		@return dt
*/
double ComputeTimeStep (Particle* all_particle) {
	//~ double max = 40 * sqrt(2*gravity*dam_height);
	//~ double v;
	//~ for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
		//~ if (all_particle[i].tag == interior) {
			//~ v = sqrt(pow(all_particle[i].velocity.first, 2) + pow(all_particle[i].velocity.second, 2));
			//~ max = (v > max) ? v : max;
		//~ }
	//~ }
	//~ return 0.4 * H / max;
    return 0.00005;
}

double TimeLoop() {
  double dt, t = 0;

  // Initialization
  Particle *all_particle = Init();
  Particle *initial_configuration = Init();
  printf("init completed.\n");

  // choose which time integration method to use. By default using Explicit
  // Euler
  Set_Integration_Method(EXPLICIT_EULER);

  int N = NUMBER_OF_PARTICLE; // get the number of particles

  // Initial steps without moving the boundary
  for (int step = 0; step < 20000; step++) {
    dt = ComputeTimeStep(all_particle);

     // if using Heun or Midpoint method
    if (integration == HEUN || integration == MIDPOINT)
      Time_Integration_Half(all_particle, dt);

    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
      SearchNeighbors(all_particle, i);
    }

    ComputeGlobalKernel(all_particle);
    ComputeGlobalDensity(all_particle);
    DensityAndBCVelocityCorrection(all_particle);
    ComputeGlobalPressure(all_particle, t);
    ComputeInteriorLaminarAcceleration(all_particle, t);
    AddRepulsiveForce(all_particle, t);

    Time_Integration(all_particle, dt);
  }

  // Write output of Initialization
  WriteData(all_particle, t);

  for (int step = 0; step < 100000; step++) {
    dt = ComputeTimeStep(all_particle);

     // if using Heun or Midpoint method
    if (integration == HEUN || integration == MIDPOINT)
      Time_Integration_Half(all_particle, dt);

    DisplaceBoundaries(all_particle, initial_configuration, t);
    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
      SearchNeighbors(all_particle, i);
    }

    ComputeGlobalKernel(all_particle);
    ComputeGlobalDensity(all_particle);
    DensityAndBCVelocityCorrection(all_particle);
    ComputeGlobalPressure(all_particle, t);
    ComputeInteriorLaminarAcceleration(all_particle, t);
    AddRepulsiveForce(all_particle, t);

    Time_Integration(all_particle, dt);
    t += dt;

    // output data to file
    if ((step + 1) % 2000 == 0) {
      WriteData(all_particle, t);
    }
    printf("time t = %f\n", t);
  }
  return t;
}
#endif // TIME_LOOP_H
