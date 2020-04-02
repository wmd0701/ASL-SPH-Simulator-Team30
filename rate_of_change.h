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
	int N = NUMBER_OF_PARTICLE;
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
	int N = NUMBER_OF_PARTICLE;   // get the number of particles
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
    int N = NUMBER_OF_PARTICLE;
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
    
    int N = NUMBER_OF_PARTICLE;
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
void ComputeInteriorLaminarAcceleration(Particle *all_particle) {
  for (index i = 0; i < NUMBER_OF_PARTICLE; i++) {
    Neighbor_p n = all_particle[i].neighbors;
    Particle *pi = &all_particle[i];

    while (n != NULL) {
      Particle *pj = &all_particle[n->idx];
      vector gradient = KernelGradient(pi->position, pj->position);
      double constant1 =
          pj->mass * (pi->pressure / (pi->density * pi->density) +
                      pj->pressure / (pj->density * pj->density));
      all_particle[i].accelerat.first -= constant1 * gradient.first;
      all_particle[i].accelerat.second -= constant1 * gradient.second;

      double xij = sqrt(pow((pj->position.first - pi->position.first), 2) +
                        pow((pj->position.second - pi->position.second), 2));
      double vij = sqrt(pow((pj->velocity.first - pi->velocity.first), 2) +
                        pow((pj->velocity.second - pi->velocity.second), 2));

      double constant2 =
          (4.0 * pj->mass * (2.0 * dynamic_viscosity) * xij * vij) /
          (pow((pi->density + pj->density), 2) *
           (pow(xij, 2) + 0.01 * pow(H, 2)));

      all_particle[i].accelerat.first += constant2 * gradient.first;
      all_particle[i].accelerat.second += constant2 * gradient.second - gravity;

      n = n->next;
    }
  }
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






