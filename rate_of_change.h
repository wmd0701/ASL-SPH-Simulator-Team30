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
            all_particle[i].velocity.first /= sum[i];
            all_particle[i].velocity.second /= sum[i];
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
double ComputeSoundSpeedSquared(Particle *all_particle){
	//Bulk velocity for a dam break
    double v2 = 2*dam_height*gravity;
    return v2*10;
}

/**   
*     @brief Compute the value of pressure field at every particle
*           
*              p = c2 * max(rho - rho0, 0)
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [pressure] attribute in all_particle
*/
void ComputeGlobalPressure (Particle *all_particle){
	double c2 = ComputeSoundSpeedSquared(all_particle);
    
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
*     @brief Compute the value of pressure field at every particle
*           
*              p = c2 * rho0 / 7 * max((rho/rho0)^7 - 1, 0)
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [pressure] attribute in all_particle
*/
void ComputeGlobalPressure2 (Particle *all_particle){
	double c2 = ComputeSoundSpeedSquared(all_particle);
    
    int N = NUMBER_OF_PARTICLE;
	for (int i = 0; i < N; i++) {
        all_particle[i].pressure = c2 * initial_density / 7 * (pow(all_particle[i].density / initial_density, 7) - 1);
        if (all_particle[i].pressure < 0) {
            all_particle[i].pressure = 0;
        }
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

                sum.first  -= all_particle[p->idx].velocity.first  * p->Wij * all_particle[p->idx].mass / all_particle[p->idx].density;
                sum.second -= all_particle[p->idx].velocity.second * p->Wij * all_particle[p->idx].mass / all_particle[p->idx].density;

            }
            all_particle[i].velocity.first = sum.first;//sum.first;
            all_particle[i].velocity.second = sum.second; //sum.second;
        }
    }
}

/**   
*     @brief Compute dvdt for every particle without considering the influence of boundary points and turbulence
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void ComputeInteriorLaminarAcceleration(Particle *all_particle) {
    double c = sqrt(ComputeSoundSpeedSquared(all_particle));
    double alpha = 0.2;
    double mu_ij, PI_ij;

    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
        if (all_particle[i].tag == interior) {
            Neighbor_p *n = all_particle[i].neighbors;
            Particle *pi = &all_particle[i];
            pi->accelerat.first  = 0;
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
                    PI_ij = - alpha * c * mu_ij / (pi->density + pj->density);
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
*     @brief Add the turbulent term into the acceleration of every particle
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddTurbulentModel(Particle *all_particle){
}

/**   
*     @brief Add the repulsive force into the acceleration of involved particle
*
*     @note  This is dependent on concrete cases.
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddRepulsiveForce(Particle *all_particle){
    int N = NUMBER_OF_PARTICLE;
    int i;
    Particle *pi;
    double c2 = ComputeSoundSpeedSquared(all_particle);
    double beta, q;
    vector pos;
    for (i = 0; i < N; i++) {
        pi = &all_particle[i];
        if (pi->tag == interior) {                     
            pos = pi->position;

            // for bottom boundary
            q = pos.second / (0.75 * H) ;
            if (q > 0) {
                if (q <= 2 / 3) {
                    beta = 0.02 * c2 / pos.second;
                    pi->accelerat.second += 2 / 3 * beta / 2; 
                }
                else if (q <= 1) {
                    beta = 0.02 * c2 / pos.second;
                    pi->accelerat.second += beta * (2 * q - 3/2 * q * q) / 2;
                }
                else if (q <= 2) {
                    beta = 0.02 * c2 / pos.second;
                    pi->accelerat.second += 0.5 * beta * (2 - q) * (2 - q) / 2;
                }
            }

            // for left boundary
            q = pos.first / (0.75 * H);
            if (q > 0) {
                if (q <= 2 / 3) {
                    beta = 0.02 * c2 / pos.first;
                    pi->accelerat.first += 2 / 3 * beta / 2;
                }
                else if (q <= 1) {
                    beta = 0.02 * c2 / pos.first;
                    pi->accelerat.first += beta * (2 * q - 3/2 * q * q) / 2;
                }
                else if (q <= 2) {
                    beta = 0.02 * c2 / pos.first;
                    pi->accelerat.first += 0.5 * beta * (2 - q) * (2 - q) / 2;
                }
            }

            // for right boundary
            q = (160 * H - pos.first) / (0.75 * H);
            if (q > 0) {
                if (q <= 2 / 3) {
                    beta = 0.02 * c2 / (160 * H - pos.first);
                    pi->accelerat.first -= 2 / 3 * beta / 2;
                }
                else if (q <= 1) {
                    beta = 0.02 * c2 / (160 * H - pos.first);
                    pi->accelerat.first -= beta * (2 * q - 3/2 * q * q) / 2;
                }
                else if (q <= 2) {
                    beta = 0.02 * c2 / (160 * H - pos.first);
                    pi->accelerat.first -= 0.5 * beta * (2 - q) * (2 - q) / 2;
                }
            }
        }
    }	
}


/**   
*     @brief Add the repulsive force into the acceleration of involved particle
*
*     @note  This is dependent on concrete cases.
*
*     @param all_particle pointer to an array containing information of all the particles
*     @return no returns. Update the [accelerat] attribute in all_particle
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

void AddRepulsiveForce2(Particle *all_particle){
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
                    double c2 = ComputeSoundSpeedSquared(all_particle);
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
 *      @brief use an extra inertial force to deal with moving boundary
 */
void AddInertialForce (Particle* all_particle, double t_now) {
    int N = NUMBER_OF_PARTICLE;
    double amplitude = 0.32, T = 1.5;
    
    for (int i = 0; i < N; i++) {
        if (all_particle[i].tag == interior) {

            all_particle[i].accelerat.first += - amplitude * 2 * M_PI / T * sin(2 * M_PI * t_now / T);

        }
    }
    return;
}

void DisplaceBoundaries(Particle* all_particle, Particle* initial_configuration, double t){
    double A = amplitude;
    double T = period;
    for(int i = 0; i < NUMBER_OF_PARTICLE; ++i){
        if(all_particle[i].tag != interior){
            all_particle[i].position.first = initial_configuration[i].position.first + (A*cos(2*M_PI*t/T) - A); 
        }
    }
}

#endif // RATE_OF_CHANGE_H







