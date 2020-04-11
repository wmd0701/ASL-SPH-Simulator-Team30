//!  @file time_loop.h
#ifndef TIME_LOOP_H
#define TIME_LOOP_H

#include "rate_of_change.h"
#include "time_integration.h"
#include "output.h"

/** 	
*		@brief Compute dt for the next iteration
*	
*			   Compute the maximum dt that won't cause every particle to move more than range of a certain radius.
*			   This is needed to validly update the nearby relation in all_particle
*
*		@param all_particle pointer to an array containing information of all the particles
*		@return dt
*/
double ComputeTimeStep (Particle* all_particle) {
	return 0.01;
}

/** 	
*		@brief Iterate until t = t_end
*		@param all_particle pointer to an array containing information of all the particles
*		@param t_end end time
*/
void TimeLoop (double t_end) {
	
	double dt, t = 0;
	
	Particle* all_particle = Init();
	WriteData               (all_particle, t);

	// choose which time integration method to use, details in time_integration.h
	Set_Integration_Method(EXPLICIT_EULER);
    
  int N = NUMBER_OF_PARTICLE;   // get the number of particles
	
	while (t < t_end) {
		dt = ComputeTimeStep     (all_particle);
		
		//----------------------------------
		// Tianwei
		ComputeGlobalKernel      (all_particle);
		ComputeGlobalDensity     (all_particle);
        
		// DensityCorrection        (all_particle);
        
		//KernelGradientCorrection (all_particle);
		//----------------------------------
		
		ComputeGhostAndRepulsiveVelocity     (all_particle);
		
		//----------------------------------
		// Valerie
		ComputeInteriorLaminarAcceleration (all_particle);
		//AddTurbulentModel        (all_particle);
		AddRepulsiveForce	     (all_particle);
		//----------------------------------

		//----------------------------------
		// Silvia
		ComputeGlobalPressure    (all_particle);
		//----------------------------------
		
		//----------------------------------
		// Mengdi
		Time_Integration		(all_particle, dt);
    //----------------------------------
		
		t += dt;
		
		for(int i = 0; i < NUMBER_OF_PARTICLE; i++)
			SearchNeighbors(all_particle, i);

		//output data to file
		printf("write at time t = %f\n",t);
		WriteData               (all_particle, t);
	}

	free(all_particle);
	
}

#endif // TIME_LOOP_H
