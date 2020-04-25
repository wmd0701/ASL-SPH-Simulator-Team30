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

/** 	
*		@brief Iterate the simulation
*
*/
double TimeLoop(){
    double dt, t = 0;

	// Initialization
	Particle* all_particle = Init4();
	Particle* initial_configuration = Init4();
	printf("init completed.\n");
    
    // Define a file for the wave height
	char output_path_wave[40];
	strcpy(output_path_wave, folder_name);
	strcat(output_path_wave, "/wave_height.csv");
	// open file for wave height
	FILE *fp = fopen(output_path_wave, "w");
	fprintf(fp, "t, height\n");
	
	// Write output of Initialization
	WriteData(all_particle, t);

	// Define a file for the wave height
	char output_path_wave[40];
	strcpy(output_path_wave, folder_name);
	strcat(output_path_wave, "/wave_height.csv");
	// open file for wave height
	FILE *fp = fopen(output_path_wave, "w");
	fprintf(fp, "t, height\n");

	// choose which time integration method to use. By default using Explicit Euler
	Set_Integration_Method(EXPLICIT_EULER);
    int N = NUMBER_OF_PARTICLE;   // get the number of particles

	// loop over time steps
	for (int step = 0; step < 80000; step ++) {
        dt = ComputeTimeStep                 (all_particle);

		// if using Heun or Midpoint method
		if(integration == HEUN || integration == MIDPOINT)
			Time_Integration_Half(all_particle, dt);

		ComputeGlobalKernel                  (all_particle);
		ComputeGlobalDensity                 (all_particle);

		ComputeGhostAndRepulsiveVelocity     (all_particle);
		DensityAndBCVelocityCorrection       (all_particle);
        //KernelGradientCorrection           (all_particle);
		
		ComputeGlobalPressure2               (all_particle);
		ComputeInteriorLaminarAcceleration   (all_particle);
		//AddTurbulentModel                  (all_particle);

		AddRepulsiveForce	                 (all_particle);
		// AddInertialForce		             (all_particle, t);   

		Time_Integration		             (all_particle, dt);

		t += dt;
	
		// Search the neighbors for the next time step
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
    RecordWaveHeight(all_particle, fp, t);
	fclose(fp);
	free(all_particle);
    return t;
}

double TimeLoop2 () {
	double dt, t = 0;

	// Initialization
	Particle* all_particle = Init4();
	Particle* initial_configuration = Init4();
	printf("init completed.\n");
	
	// Write output of Initialization
	WriteData(all_particle, t);

	// choose which time integration method to use. By default using Explicit Euler
	Set_Integration_Method(EXPLICIT_EULER);
    
  	int N = NUMBER_OF_PARTICLE;   // get the number of particles
    
    //Initial steps without moving the boundary
    for (int step = 0; step < 20000; step ++) {
		dt = ComputeTimeStep                 (all_particle);

		ComputeGlobalKernel                  (all_particle);
		ComputeGlobalDensity                 (all_particle);

		ComputeGhostAndRepulsiveVelocity     (all_particle);
		DensityAndBCVelocityCorrection       (all_particle);
		ComputeGlobalPressure                (all_particle, t);
		ComputeInteriorLaminarAcceleration   (all_particle, t);
		AddRepulsiveForce2	                 (all_particle, t);
		
		Time_Integration		             (all_particle, dt);
		
		for(int i = 0; i < NUMBER_OF_PARTICLE; i++){
			SearchNeighbors(all_particle, i);
		}	
	}
	
	for (int step = 0; step < 100000; step ++) {
		dt = ComputeTimeStep                 (all_particle);

		ComputeGlobalKernel                  (all_particle);
		ComputeGlobalDensity                 (all_particle);
		
        ComputeGhostAndRepulsiveVelocity     (all_particle);
		DensityAndBCVelocityCorrection       (all_particle);
		ComputeGlobalPressure                (all_particle, t);
		ComputeInteriorLaminarAcceleration   (all_particle, t);
		AddRepulsiveForce2	                 (all_particle, t);
		
		Time_Integration		             (all_particle, dt);
		
		t += dt;
        //output data to file
		if (step % 100 == 0) {
			WriteData(all_particle, t);
			RecordWaveHeight(all_particle, fp, t);
		}

		printf("time t = %f\n",t);
		
        DisplaceBoundaries(all_particle, initial_configuration, t);
		for(int i = 0; i < NUMBER_OF_PARTICLE; i++){
			SearchNeighbors(all_particle, i);
		}	
	}
	return t;
}

#endif // TIME_LOOP_H
