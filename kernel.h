//!  @file kernel.h
#ifndef KERNEL_H
#define KERNEL_H

#include "data_set.h"
#include <math.h>

#define M_PI 3.141592653579
#define H 1   //!< smoothing length

/**  
*    @brief Compute the kernel function
*    @param xi coordinate of the particle which is being considered now
*    @param xj coordinate of the particle nearby
*	 @return value of kernel function
*/
double Kernel (const vector xi,const vector xj) {
	
	double r = sqrt(pow((xi.first - xj.first),2) + pow((xi.second - xj.second),2));
	double epsilon = r / H;
	if (r <= 1) {
		return (1 - 3/2*epsilon*epsilon + 3/4*epsilon*epsilon*epsilon) / M_PI / (H*H*H);
	}
	else if(r <= 2) {
		return 0.25 * (2 - epsilon) * (2 - epsilon) * (2 - epsilon) / M_PI / (H*H*H); 
	}
}

/** 
*   @brief Compute the gradient of the kernel function
*   @param xi coordinate of the particle which is being considered now
*   @param xj coordinate of the particle nearby
*	@return gradient of kernel function
*/   
vector KernelGradient (const vector xi,const vector xj) {
	double r = sqrt(pow((xi.first - xj.first),2) + pow((xi.second - xj.second),2));
	double epsilon = r / H;
	double tmp = 1 / (M_PI * pow(H, 4) * r);
	vector grad;
	if (r <= 1) {

		tmp *= (-3) * epsilon + 2.25 * epsilon * epsilon;
		grad.first  = tmp * (xi.first  - xj.first);
		grad.second = tmp * (xi.second - xj.second);
		return grad;

	}
	else if (r <= 2) {

		tmp *= (-0.75) * (2 - epsilon) * (2 - epsilon);
		grad.first  = tmp * (xi.first  - xj.first);
		grad.second = tmp * (xi.second - xj.second);
		return grad;

	}
}

/** 
*   @brief Compute the [kernel] & [kernel gradient] for every particle w.r.t. its nearby particles
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [neighbors] attribute in all_particle
*/   
void ComputeGlobalKernel (Particle *all_particle) {
	vector xi, xj;
	int N = sizeof(all_particle) / sizeof(all_particle[0]);
	for (index i = 0; i < N; i++) {

		xi = all_particle[i].position;
		for (Neighbor_p p = all_particle[i].neighbors; p != NULL; p = p->next) {

			xj = all_particle[p->idx].position;
			
			p->Wij 		  = Kernel(xi, xj);
			p->Wij_grad_i = KernelGradient(xi, xj);

		}

	}	
}


/** 
*   @brief Compute the "corrected" gradient of the kernel function
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the L attribute in all_particle, then update the Wij_grad attribute.
*/   
void KernelGradientCorrection (Particle *all_particle) {
	double a00, a01, a10, a11, V, xji, yji;
	double determinant;
	vector new_grad;
	int N = sizeof(all_particle) / sizeof(all_particle[0]);

	for (index i = 0; i < N; i++) {

		a00 = 0;
		a01 = 0;
		a10 = 0;
		a11 = 0;
		for (Neighbor_p p = all_particle[i].neighbors; p != NULL; p = p->next) {
			
			V = all_particle[p->idx].mass / all_particle[p->idx].density;
			xji = all_particle[p->idx].position.first  - all_particle[i].position.first;
			yji = all_particle[p->idx].position.second - all_particle[i].position.second;
			a00 += xji * p->Wij_grad_i.first  * V;
			a01 += yji * p->Wij_grad_i.first  * V;
			a10 += xji * p->Wij_grad_i.second * V;
			a11 += xji * p->Wij_grad_i.second * V;

		}

		determinant = a00 * a11 - a10 * a01;

		for (Neighbor_p p = all_particle[i].neighbors; p != NULL; p = p->next) {
			
			new_grad.first  = (  a00 * p->Wij_grad_i.first - a01 * p->Wij_grad_i.second) / determinant;
			new_grad.second = (- a10 * p->Wij_grad_i.first + a11 * p->Wij_grad_i.second) / determinant;
			p->Wij_grad_i = new_grad;
			
		}
		
	}
}

#endif // KERNEL_H

