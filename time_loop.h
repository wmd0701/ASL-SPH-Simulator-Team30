//!  @file time_loop.h
#ifndef TIME_LOOP_H
#define TIME_LOOP_H

#include "rate_of_change.h"
#include "time_integration.h"
#include "output.h"

/** 	
*		@brief Compute dt for the next iteration using CFL criterion
*	
*		@param all_particle pointer to an array containing information of all the particles
*		@return dt
*/
double ComputeTimeStep (Particle* all_particle) {
	double max = 40 * sqrt(2*gravity*dam_height);
	double v;
	for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
		if (all_particle[i].tag == interior) {
			v = sqrt(pow(all_particle[i].velocity.first, 2) + pow(all_particle[i].velocity.second, 2));
			max = (v > max) ? v : max;
		}
	}
	return 0.4 * H / max;
}

/** 	
*		@brief Iterate the simulation
*
*/
double TimeLoop () {
	double dt, t = 0;
	
	Particle* all_particle = Init4();
	printf("init completed.\n");
	WriteData(all_particle, t);

	FILE *fp = fopen("data/wave_height.csv", "w");
	fprintf(fp, "t, height\n");

	// choose which time integration method to use, details in time_integration.h
	Set_Integration_Method(EXPLICIT_EULER);
    
  	int N = NUMBER_OF_PARTICLE;   // get the number of particles
	
	for (int step = 0; step < 80000; step ++) {
		dt = ComputeTimeStep                 (all_particle);

		ComputeGlobalKernel                  (all_particle);
		ComputeGlobalDensity                 (all_particle);
		
        ComputeGhostAndRepulsiveVelocity     (all_particle);
		DensityAndBCVelocityCorrection       (all_particle);
    	//KernelGradientCorrection           (all_particle);
		ComputeGlobalPressure2               (all_particle);
		ComputeInteriorLaminarAcceleration   (all_particle);
		//AddTurbulentModel                  (all_particle);
		AddRepulsiveForce	                 (all_particle);
		AddInertialForce		             (all_particle, t);   
		
		Time_Integration		             (all_particle, dt);
		
		t += dt;
		
		for(int i = 0; i < NUMBER_OF_PARTICLE; i++){
			SearchNeighbors(all_particle, i);
		}

		//output data to file
		if (step % 200 == 0) {
			WriteData(all_particle, t);
			RecordWaveHeight(all_particle, fp, t);
		}
		printf("time t = %f\n",t);
	
	}

	WriteData(all_particle, t);
	RecordWaveHeight(all_particle, fp, t);
	fclose(fp);
	free(all_particle);

	return t;
}

#endif // TIME_LOOP_H
