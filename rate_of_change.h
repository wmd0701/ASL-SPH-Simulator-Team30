//!  @file rate_of_change.h
#ifndef RATE_OF_CHANGE_H
#define RATE_OF_CHANGE_H

#include "data_set.h"
#include "kernel.h"
#include "constants.h"

/**   
*     @brief Compute the value of density field at a certain particle, namely rho_i
*     @param all_particle pointer to an array containing information of all the particles
*	  @param ptc_idx index of the particle that is being considered
*     @return density value
*/
double ComputeLocalDensity (Particle *all_particle, index ptc_idx) {
	double sum = 0;
    for (Neighbor_p p = all_particle[ptc_idx].neighbors; p != NULL; p = p->next) {
        sum += p->Wij * all_particle[ptc_idx].mass;
    }
    return sum;
}


/**!  
*     @brief Compute the value of density field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [density] attribute in all_particle
*/
void ComputeGlobalDensity (Particle *all_particle) {
	int N = sizeof(all_particle) / sizeof(all_particle[0]);
    for (index i = 0; i < N; i++) {
        ComputeLocalDensity(all_particle, i);
    }
}

/**   
*     @brief Compute the "corrected" value of density field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [density] attribute in all_particle
*/
void DensityCorrection (Particle *all_particle) {
    double sum;
	int N = sizeof(all_particle) / sizeof(all_particle[0]);   // get the number of particles
    for (index i = 0; i < N; i++) {     // traverse particles
        sum = 0;
        for (Neighbor_p p = all_particle[i].neighbors; p != NULL; p = p->next) {    // traverse neighbors
            sum += p->Wij * all_particle[i].mass / all_particle[i].density;
        }
        all_particle[i].density /= sum;
    }
}

/**   
*     @brief Compute the velocity of sound
*     @param all_particle pointer to an array containing information of all the particles
*     @return sound speed value squared
*/
double ComputeSoundSpeedSquared(Particle *all_particle){
	//Bulk velocity for a dam break
    double v = 2*dam_height*gravity;
    
    //Computing delta
    int N = sizeof(all_particle) / sizeof(all_particle[0]);
    double delta = abs(all_particle[0].density - initial_density)/initial_density;
    double temp = delta;
	for (index i = 1; i < N; i++) {
        temp = abs(all_particle[i].density - initial_density)/initial_density;
        if (temp < delta){
            delta = temp;
        }
    }
    
    //This is only the first term of the three shown in the paper
    return v*v/delta;
}

/**   
*     @brief Compute the value of pressure field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [pressure] attribute in all_particle
*/
void ComputeGlobalPressure (Particle *all_particle){
	double c2 = ComputeSoundSpeedSquared(all_particle);
    
    int N = sizeof(all_particle) / sizeof(all_particle[0]);
	for (index i = 0; i < N; i++) {
        all_particle[i].pressure = c2 * all_particle[i].density;
    }
    return;
}

/**   
*     @brief Compute the velocity of every ghost particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [velocity] attribute in all_particle
*/
void ComputeGhostVelocity (Particle *all_particle){
	
}

/**   
*     @brief Compute dvdt for every particle without considering the influence of boundary points and turbulence
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void ComputeInteriorLaminarAcceleration (Particle *all_particle) {
	
}

/**   
*     @brief Add the turbulent term into the acceleration of every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddTurbulentModel(Particle *all_particle){
	
}

/**   
*     @brief Add the repulsive force into the acceleration of involved particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddRepulsiveForce(Particle *all_particle){
	
}

#endif // RATE_OF_CHANGE_H






