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
double ComputeLocalDensity (Particle *all_particle, int ptc_idx) {
	double sum = 0;
    for (Neighbor_p *p = all_particle[ptc_idx].neighbors; p != NULL; p = p->next) {
        sum += p->Wij * all_particle[p->idx].mass;
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
    for (int i = 0; i < N; i++) {
        all_particle[i].density = ComputeLocalDensity(all_particle, i);
    }
}

/**   
*     @brief Compute the "corrected" value of density field at every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [density] attribute in all_particle
*/
void DensityCorrection (Particle *all_particle) {
	int N = NUMBER_OF_PARTICLE;   // get the number of particles
    double* sum = (double*)malloc(sizeof(double)*N);
    for (int i = 0; i < N; i++) {     // traverse particles
        double sum_Wij = 0;
        for (Neighbor_p *nk = all_particle[i].neighbors; nk != NULL; nk = nk->next) {    // traverse neighbors
                Particle * pk = &all_particle[nk->idx];
                    sum_Wij += nk->Wij * pk->mass / pk->density;
                    printf("Wij = %f \t mass = %f \t density = %f\n", nk->Wij, pk->mass, pk->density);
        }
        sum[i] = sum_Wij;
    }
    for (int i = 0; i < N; i++) {
        all_particle[i].density /= sum[i];
        printf("%u: density = %f\n", i, all_particle[i].density);
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
	for (int i = 1; i < N; i++) {
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
	for (int i = 0; i < N; i++) {
        all_particle[i].pressure = c2 * all_particle[i].density;
    }
    return;
}

/**   
*     @brief Compute the velocity of every ghost particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [velocity] attribute in all_particle
*/
void ComputeGhostAndRepulsiveVelocity (Particle *all_particle){
	int N = NUMBER_OF_PARTICLE;
	for (int i = 0; i < N; i++) {
        if(all_particle[i].tag != 0){
            vector sum;
            sum.first = 0;
            sum.second = 0;
            for (Neighbor_p *p = all_particle[i].neighbors; p != NULL; p = p->next) {    // traverse neighbors
                sum.first -= all_particle[p->idx].velocity.first * p->Wij * all_particle[p->idx].mass * all_particle[p->idx].density;
                sum.second -= all_particle[p->idx].velocity.second * p->Wij * all_particle[p->idx].mass * all_particle[p->idx].density;
            }
            all_particle[i].velocity.first = 1.0;//sum.first;
            all_particle[i].velocity.second = 1.0; //sum.second;
        }
    }
}

/**   
*     @brief Compute dvdt for every particle without considering the influence of boundary points and turbulence
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void ComputeInteriorLaminarAcceleration(Particle *all_particle) {
  for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
    if (all_particle[i].tag == interior) {
        Neighbor_p *n = all_particle[i].neighbors;
        Particle *pi = &all_particle[i];
        pi->accelerat.first  = 0;
        pi->accelerat.second = 0;
    
        while (n != NULL) {
          Particle *pj = &all_particle[n->idx];
          vector gradient = n->Wij_grad_i;
          double constant1 =
              pj->mass * (pi->pressure / (pi->density * pi->density) +
                          pj->pressure / (pj->density * pj->density));
          pi->accelerat.first -= constant1 * gradient.first;
          pi->accelerat.second -= constant1 * gradient.second;
    
          vector xij = vec_sub_vec(pj->position, pi->position);
          vector vij = vec_sub_vec(pj->velocity, pi->velocity);
    
          double constant2 =
              (4.0 * pj->mass * (2.0 * dynamic_viscosity) * vec_dot_vec(xij, gradient)) /
              (pow((pi->density + pj->density), 2) *
               (vec_dot_vec(xij, xij) + 0.01 * pow(H, 2)));
    
          pi->accelerat.first += constant2 * vij.first;
          pi->accelerat.second += constant2 * vij.second;
         	
    	  n = n->next;
        }
        
        pi->accelerat.second -= gravity;
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






