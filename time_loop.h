//!  @file time_loop.h
#ifndef TIME_LOOP_H
#define TIME_LOOP_H

#include "constants.h"
#include "rate_of_change.h"
#include "time_integration.h"
#include "output.h"
#include "tsc_x86.h"
#include <string.h>

double cycles_DispBoundary   = 0;
double cycles_SearchNeighbor = 0;
double cycles_CompGlbKernel  = 0;
double cycles_CompGlbDensity = 0;
double cycles_DensityCorr    = 0;
double cycles_CompPressure   = 0;
double cycles_CompAccelerat  = 0;
double cycles_RepulsiveForce = 0;
double cycles_TimeIntegral   = 0;


double TimeLoop() {
  myInt64 start;
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
  /*
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
  */
  // Write output of Initialization
  // WriteData(all_particle, t);

  //-------------------------------------------------------------------
  // MEASURE FROM HERE
  //-------------------------------------------------------------------
  int overall_step = 1000;
  for (int step = 0; step < overall_step; step++) {
    // ------------------------
    start = start_tsc();
    DisplaceBoundaries(all_particle, initial_configuration, t);
    cycles_DispBoundary += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
      SearchNeighbors(all_particle, i);
    }
    cycles_SearchNeighbor += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    ComputeGlobalKernel(all_particle);
    cycles_CompGlbKernel  += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    ComputeGlobalDensity(all_particle);
    cycles_CompGlbDensity += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    DensityAndBCVelocityCorrection(all_particle);
    cycles_DensityCorr += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    ComputeGlobalPressure(all_particle, t);
    cycles_CompPressure += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    ComputeInteriorLaminarAcceleration(all_particle, t);
    cycles_CompAccelerat += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    AddRepulsiveForce(all_particle, t);
    cycles_RepulsiveForce += (double)stop_tsc(start);
    // ------------------------
    start = start_tsc();
    Time_Integration(all_particle, dt);
    cycles_TimeIntegral += (double)stop_tsc(start);

    t += dt;

    // output data to file
    //~ if ((step + 1) % 2000 == 0) {
      //~ WriteData(all_particle, t);
    //~ }
    //~ printf("time t = %f\n", t);
  }

  cycles_DispBoundary   /= overall_step;
  cycles_SearchNeighbor /= overall_step;
  cycles_CompGlbKernel  /= overall_step;
  cycles_CompGlbDensity /= overall_step;
  cycles_DensityCorr    /= overall_step;
  cycles_CompPressure   /= overall_step;
  cycles_CompAccelerat  /= overall_step;
  cycles_RepulsiveForce /= overall_step;
  cycles_TimeIntegral   /= overall_step;
  return t;
}
#endif // TIME_LOOP_H
