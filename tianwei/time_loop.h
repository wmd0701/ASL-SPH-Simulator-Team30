//!     @file time_loop.h


/** 	
*		@brief Compute dt for the next iteration
*	
*			   Compute the maximum dt that won't cause every particle to move more than range of a certain radius.
*			   This is needed to validly update the nearby relation in all_particle
*
*		@param all_particle pointer to an array containing information of all the particles
*		@return dt
*/
double ComputeTimeStep (Particle*& all_particle) {
	
}

/** 	
*		@brief Iterate until t = t_end
*		@param all_particle pointer to an array containing information of all the particles
*		@param t_end end time
*/
void TimeLoop (Particle*& all_particle, double t_end) {
	double dt, t = 0;
	Particle* all_particle = Init();
	
	while (t < t_end) {
		ComputeGlobalDensity  (all_particle);
		DensityCorrection     (all_particle);
		ComputeGlobalPressure (all_particle);
		ComputeGhostVelocity  (all_particle);
		ComputeInteriorLaminarAcceleration (all_particle);
		AddTurbulentModel     (all_particle);
		AddRepulsiveForce	  (all_particle);
		dt = ComputeTimeStep  (all_particle);
		TimeUpdate			  (all_particle, dt);
		t += dt;
	}
	
}