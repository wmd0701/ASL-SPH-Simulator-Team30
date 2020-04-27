//!  @file time_loop.h
#ifndef TIME_LOOP_H
#define TIME_LOOP_H

#include "constants.h"
#include "rate_of_change.h"
#include "time_integration.h"
#include "output.h"
#include <string.h>


double TimeLoop() {
  double t = 0;
  double dt = 0.00005; //Warning: this timestep is valid only for number of interior particles < 4000! Else use dt = 0.000005

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

  //-------------------------------------------------------------------
  // MEASURE FROM HERE
  //-------------------------------------------------------------------
  for (int step = 0; step < 1000; step++) {
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
    //~ if ((step + 1) % 2000 == 0) {
      //~ WriteData(all_particle, t);
    //~ }
    //~ printf("time t = %f\n", t);
  }
  return t;
}
#endif // TIME_LOOP_H
