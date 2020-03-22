//!   @file rate_of_change.h

/**   
*     @brief Compute the value of density field at a certain particle, namely rho_i
*     @param all_particle pointer to an array containing information of all the particles
*	  @param par_idx index of the particle that is being considered
*     @return density value
*/
double ComputeLocalDensity (particle*& all_particle, int par_idx) {
	
}


/**!  
*     @brief Compute the value of density field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [density] attribute in all_particle
*/
void ComputeGlobalDensity (particle*& all_particle) {
	
}

/**   
*     @brief Compute the "corrected" value of density field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [density] attribute in all_particle
*/
void DensityCorrection (particle*& all_particle) {
	
}

/**   
*     @brief Compute the value of pressure field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [pressure] attribute in all_particle
*/
void ComputeGlobalPressure (particle*& all_particle){
	
}


/**   
*     @brief Compute the velocity of every ghost particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [velocity] attribute in all_particle
*/
void ComputeGhostVelocity (particle*& all_particle){
	
}

/**   
*     @brief Compute dvdt for every particle without considering the influence of boundary points and turbulence
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void ComputeInteriorLaminarAcceleration (particle*& all_particle) {
	
}

/**   
*     @brief Add the turbulent term into the acceleration of every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddTurbulentModel(particle*& all_particle){
	
}

/**   
*     @brief Add the repulsive force into the acceleration of involved particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddRepulsiveForce(particle*& all_particle){
	
}






