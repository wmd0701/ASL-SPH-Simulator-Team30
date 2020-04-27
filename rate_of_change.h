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
*     @brief Compute - the "corrected" value of density field at every particle using new Wij
*                    - the "corrected" value of velocity of boundary particles  using new Wij
*
*     @param all_particle pointer to an array containing information of all the particles
*/
void DensityAndBCVelocityCorrection (Particle *all_particle) {
	int N = NUMBER_OF_PARTICLE;   // get the number of particles
	double* sum = (double*)malloc(sizeof(double)*N);
	for (int i = 0; i < N; i++) {     // traverse particles
		double sum_Wij = 0;
		for (Neighbor_p *nk = all_particle[i].neighbors; nk != NULL; nk = nk->next) {    // traverse neighbors
			Particle * pk = &all_particle[nk->idx];
			sum_Wij += nk->Wij * pk->mass / pk->density;
		}
		sum[i] = sum_Wij;
	}
	for (int i = 0; i < N; i++) {
		all_particle[i].density /= sum[i];
		if (all_particle[i].tag != interior) {
			all_particle[i].velocity.first /= (-sum[i]);
			all_particle[i].velocity.second /= (-sum[i]);
		}
	}
	free(sum);
}

/**   
*     @brief Compute the velocity of sound
*            For weakly compressible fluid, it is assigned with 10 * maximum velocity.
*     
*     @note  It is dependent on concrete cases.
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return sound speed value squared
*/
double ComputeSoundSpeedSquared(Particle *all_particle, double t){
    double bulk_velocity = sqrt(2*dam_height*gravity);
    double l = 1.7;
    double delta = 0.01;
    double nu = dynamic_viscosity/initial_density;
    double force_x = 0;
    if (t > 0){
        force_x = - 0.032* (2*M_PI/1.5)*(2*M_PI/1.5)*cos(2*M_PI*t/1.5)*all_particle[0].mass;
    }
    double force_y = -gravity*all_particle[0].mass;
    double force = sqrt((force_x*force_x) + (force_y*force_y));
    
    double term1 = bulk_velocity*bulk_velocity/delta;
    double term2 = nu*bulk_velocity/(delta*l);
    double term3 = force*l/delta;
    
    //Max computation
    double temp;
    if(term1 >= term2){
        temp = term1;
    }
    else{
        temp = term2;
    }
    
    if(term3 >= temp){
        temp = term3;
    }
    return temp;
}

/**   
*     @brief Compute the value of pressure field at every particle
*           
*              p = c2 * max(rho - rho0, 0)
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [pressure] attribute in all_particle
*/

void ComputeGlobalPressure (Particle *all_particle, double t){
	double c2 = ComputeSoundSpeedSquared(all_particle, t);
    
    int N = NUMBER_OF_PARTICLE;
	for (int i = 0; i < N; i++) {
		all_particle[i].pressure = c2 * (all_particle[i].density - initial_density);
		if (all_particle[i].pressure < 0) {
			all_particle[i].pressure = 0;
		}
	}
	return;
}

/**   
*     @brief Compute dvdt for every particle without considering the influence of boundary points and turbulence
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void ComputeInteriorLaminarAcceleration(Particle *all_particle, double t) {
    double c = sqrt(ComputeSoundSpeedSquared(all_particle, t));
    double alpha = 0.2;
    double mu_ij, PI_ij;

    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
        if (all_particle[i].tag == interior) {
            Neighbor_p *n = all_particle[i].neighbors;
            Particle *pi = &all_particle[i];
            
            pi->accelerat.first = 0;
            pi->accelerat.second = 0;
            
            while (n != NULL) {
                Particle *pj = &all_particle[n->idx];
                vector gradient = n->Wij_grad_i;

                // Pressure force
                double constant1 =
                    pj->mass * (pi->pressure / (pi->density * pi->density) +
                                pj->pressure / (pj->density * pj->density));

                
                pi->accelerat.first -= constant1 * gradient.first;
                pi->accelerat.second -= constant1 * gradient.second;

                // Viscosity force
                vector xij = vec_sub_vec(pj->position, pi->position);
                vector vij = vec_sub_vec(pj->velocity, pi->velocity);
             
                if (vec_dot_vec(xij, vij) < 0) {
                    mu_ij = H * vec_dot_vec(xij, vij) / (vec_dot_vec(xij, xij) + 0.01 * H * H);
                    PI_ij = - alpha *c * mu_ij / (pi->density + pj->density);
                    pi->accelerat.first  -= pj->mass * PI_ij * gradient.first;
                    pi->accelerat.second -= pj->mass * PI_ij * gradient.second;
                }

                n = n->next;
            }

            // Gravity
            pi->accelerat.second -= gravity;
        }
    }
}

/**   
*     @brief Implementation of function f, formula 28 in the paper
*/
double f(double eta){
	const double c = 2.0 / 3.0;
	if(0.0 < eta && eta <= c)
		return c;
	else if(c < eta && eta <= 1.0)
		return 2.0 * eta - 1.5 * eta * eta;
	else if(1.0 < eta && eta < 2.0)
		return 0.5 * (2.0 - eta) * (2.0 - eta);
	else
		return 0.0;
}

/**   
*     @brief Add the repulsive force into the acceleration of involved particle
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddRepulsiveForce(Particle *all_particle, double t){
    double d = H;
    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
        // interior particles
        Neighbor_p *n = all_particle[i].neighbors;
        Particle *pi = &all_particle[i];
        
        if (pi->tag == interior) {
            while (n != NULL) {
                Particle *pj = &all_particle[n->idx];
                if (pj->tag == repulsive) {
                    vector xij = vec_sub_vec(pj->position, pi->position);
                    double r2 = vec_dot_vec(xij, xij);
                    double r = sqrt(r2);
                    double c2 = ComputeSoundSpeedSquared(all_particle, t);
                    double eta = r / (0.75 * H);
                    if (0 < r && r < d) {
                        double chi = 1 - r / d;
                        double constant = 0.01 * c2 * chi * f(eta) / r2;
                        
                        pi->accelerat.first -= constant * xij.first;
                        pi->accelerat.second -= constant * xij.second;
                    }
                }
                n = n->next;
            }
        }
    }
}

/**   
*     @brief Displace the boundaries and update their velocity according to x(t) = A*cos(2*M_PI*t/T)
*
*     @param all_particle: pointer to an array containing information of all the particles
*     @param initial_configuration: pointer to an array containing the initial position of all particles
*     @return no returns. Update the attributes in all_particle
*/
void DisplaceBoundaries(Particle* all_particle, Particle* initial_configuration, double t){
	double A = amplitude;
	double T = period;
	for(int i = 0; i < NUMBER_OF_PARTICLE; ++i){
		if(all_particle[i].tag != interior){
			all_particle[i].position.first = initial_configuration[i].position.first + (A*cos(2*M_PI*t/T) - A); 
			all_particle[i].velocity.first = - 2 * M_PI * A * sin(2 * M_PI * t / T) / T;
		}
	}
}

#endif // RATE_OF_CHANGE_H







