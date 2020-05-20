//! @file rate_of_change.h
//! @file rate_of_change.h
#ifndef RATE_OF_CHANGE_H
#define RATE_OF_CHANGE_H

#include "data_set.h"
#include "constants.h"

/**   
*     @brief Compute the value of density field at a certain particle, namely rho_i
*	  @param ptc_idx index of the particle that is being considered
*     @return density value
*/
double ComputeLocalDensity (int ptc_idx) {
	double sum = 0;
    int counts = neighbor_counts[ptc_idx];
    for (int i = 0 ; i < counts ; i++)
        sum += Wijs[ptc_idx][i] * mass;
	return sum;
}

/**!  
*     @brief Compute the value of density field at every particle
*     @return no returns. Update the [density] attribute in all_particle
*/
void ComputeGlobalDensity () {
	for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
		densities[i] = ComputeLocalDensity(i);
	}
}

/**   
*     @brief Compute - the "corrected" value of density field at every particle using new Wij
*                    - the "corrected" value of velocity of boundary particles  using new Wij
*/
void DensityAndBCVelocityCorrection () {
	double* sum = (double*)malloc(sizeof(double) * NUMBER_OF_PARTICLE);
	for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {     // traverse particles
		double sum_Wij = 0;
		for (int j = 0; j < neighbor_counts[i]; j++) {    // traverse neighbors
			double density = densities[neighbor_indices[i][j]]; 
            sum_Wij += Wijs[i][j] * mass / density;
		}
		sum[i] = sum_Wij;
	}

	for (int i = 0; i < NUMBER_OF_PARTICLE; i++)
		densities[i] /= sum[i];
    
    // if not interior
    for (int i = N_interior ; i < NUMBER_OF_PARTICLE ; i++){
        x_velocities[i] /= (-sum[i]);
        y_velocities[i] /= (-sum[i]);
    }

	free(sum);
}

/**   
*     @brief Compute the velocity of sound
*            For weakly compressible fluid, it is assigned with 10 * maximum velocity.
*     
*     @note  It is dependent on concrete cases.
*
*     @return sound speed value squared
*/
double ComputeSoundSpeedSquared(double t){
    double bulk_velocity = sqrt(2*dam_height*gravity);
    double l = 1.7;
    double delta = 0.01;
    double nu = dynamic_viscosity/initial_density;
    double force_x = 0;
    if (t > 0){
        force_x = - 0.032* (2*M_PI/1.5)*(2*M_PI/1.5)*cos(2*M_PI*t/1.5)*mass;
    }
    double force_y = -gravity*mass;
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
*     @return no returns. Update the [pressure] attribute in all_particle
*/

void ComputeGlobalPressure (double t){
	double c2 = ComputeSoundSpeedSquared(t);
    
	for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
		pressures[i] = c2 * (densities[i] - initial_density);
		if (pressures[i] < 0) {
			pressures[i] = 0;
		}
	}
	return;
}

/**   
*     @brief Compute dvdt for every particle without considering the influence of boundary points and turbulence
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void ComputeInteriorLaminarAcceleration(double t) {
    double c = sqrt(ComputeSoundSpeedSquared(t));
    double alpha = 0.2;
    double mu_ij, PI_ij;

    // only for interior particles
    for (int i = 0; i < N_interior; i++) {
        x_accelerats[i] = 0.;
        y_accelerats[i] = 0.;
                
        for (int j = 0; j < neighbor_counts[i]; j++) {
            int neighbor    = neighbor_indices[i][j];      // index of neighbor
            double x_gradient = x_Wij_grads[i][j];
            double y_gradient = y_Wij_grads[i][j];
            
            // Pressure force
            double constant1 =
                mass * (pressures[i] / (densities[i] * densities[i]) +
                        pressures[neighbor] / (densities[neighbor] * densities[neighbor]));
            
            x_accelerats[i] -= constant1 * x_gradient;
            y_accelerats[i] -= constant1 * y_gradient;

            // Viscosity force
            double xij_first  = x_positions[neighbor] - x_positions[i];
            double xij_second = y_positions[neighbor] - y_positions[i];
            double vij_first  = x_velocities[neighbor] - x_velocities[i];
            double vij_second = y_velocities[neighbor] - y_velocities[i];
            double xij_dot_vij = xij_first * vij_first + xij_second * vij_second;
            
            if (xij_dot_vij < 0) {
                mu_ij = H * xij_dot_vij / ((xij_first*xij_first + xij_second*xij_second) + 0.01 * H * H);
                PI_ij = - alpha *c * mu_ij / (densities[i] + densities[neighbor]);
                x_accelerats[i] -= mass * PI_ij * x_gradient;
                y_accelerats[i] -= mass * PI_ij * y_gradient;
            }
        }

        // Gravity
        y_accelerats[i] -= gravity;
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
*     @return no returns. Update the [accelerat] attribute in all_particle
*/
void AddRepulsiveForce(double t){
    double d = H;
    for (int i = 0; i < N_interior; i++) {
        for (int j = 0; j < neighbor_counts[i]; j++) {
            int neighbor = neighbor_indices[i][j];  // index of neighbor
            // if neighbor is repulsive 
            if (neighbor >= N_interior && neighbor < N_interior + N_repulsive) {
                double xij_first  = x_positions[neighbor] - x_positions[i];
                double xij_second = y_positions[neighbor] - y_positions[i];
                double r2 = xij_first*xij_first + xij_second*xij_second;
                double r = sqrt(r2);
                double c2 = ComputeSoundSpeedSquared(t);
                double eta = r / (0.75 * H);
                if (0 < r && r < d) {
                    double chi = 1 - r / d;
                    double constant = 0.01 * c2 * chi * f(eta) / r2;
                    
                    x_accelerats[i] -= constant * xij_first;
                    y_accelerats[i] -= constant * xij_second;
                }
            }
        }

    }
}

/**   
*     @brief Displace the boundaries and update their velocity according to x(t) = A*cos(2*M_PI*t/T)
*
*     @return no returns. Update the attributes in all_particle
*/
void DisplaceBoundaries(double t){
	double A = amplitude;
	double T = period;
	for(int i = 0; i < N_boundary; ++i){
        x_positions[i + N_interior] = x_init_positions[i] + (A*cos(2*M_PI*t/T) - A); 
        x_velocities[i + N_interior] = - 2 * M_PI * A * sin(2 * M_PI * t / T) / T;
	}
}

#endif // RATE_OF_CHANGE_H







