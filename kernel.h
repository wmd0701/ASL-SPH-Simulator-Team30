//!  @file kernel.h
#ifndef KERNEL_H
#define KERNEL_H

#define _USE_MATH_DEFINES

#include <math.h>
#include "constants.h"
#include "data_set.h"

/**  
*    @brief Compute the kernel function
*    @param xi coordinate of the particle which is being considered now
*    @param xj coordinate of the particle nearby
*	 @return value of kernel function
*/
double Kernel (const vector xi,const vector xj) {

	double r = sqrt(pow((xi.first - xj.first),2) + pow((xi.second - xj.second),2));
	double q = r / H;

	double factor = 10 / 7 / M_PI / H / H;
	if (q <= 1) {
		return factor * ( 1 - 1.5 * q * q * (1 - 0.5*q));
	}
	else if (q <= 2) {
		return factor * 0.25 * pow((2 - q), 3);
	}
	else{
		printf("something wrong with searchneibors.");
		return 0;
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
	double q = r / H;
	double factor = 10 / 7 / M_PI / H / H;
	vector grad;

	if(1e-12 <= q && q <= 1){
		double temp = factor * (-3 * q + 9/4 * q * q);
		grad.first = temp*(xi.first - xj.first)/(H * r);
		grad.second = temp*(xi.second - xj.second)/(H * r);
	}
	else if(1 <= q && q <= 2){
		double temp = - factor * 3 / 4 * (2 - q) * (2 - q);
		grad.first = temp*(xi.first - xj.first)/(H * r);
		grad.second = temp*(xi.second - xj.second)/(H * r); 
	}
	else{
		grad.first = 0;
		grad.second = 0;
	}

	return grad;

}

/** 
*   @brief Compute the [kernel] & [kernel gradient] for every particle w.r.t. its nearby particles
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [neighbors] attribute in all_particle
*/   
void ComputeGlobalKernel (Particle *all_particle) {
	vector xi, xj;
	int N = NUMBER_OF_PARTICLE;
	for (int i = 0; i < N; i++) {
		xi = all_particle[i].position;
		for (Neighbor_p *p = all_particle[i].neighbors; p != NULL; p = p->next) {
			xj = all_particle[p->idx].position;
			p->Wij 		  = Kernel(xi, xj);
			p->Wij_grad_i = KernelGradient(xi, xj);

		}
	}	
}
#endif // KERNEL_H

