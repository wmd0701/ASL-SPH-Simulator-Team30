//!  @file kernel.h
#ifndef KERNEL_H
#define KERNEL_H

#define _USE_MATH_DEFINES

#include "data_set.h"
#include "constants.h"
#include <math.h>

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
    //double prefactor = 40./(7*M_PI*H*H);
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
	for (Index i = 0; i < N; i++) {

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
*	@param all_particle all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [Wij_grad_i] in [neighbors] of [all_particle]
*/   
void KernelGradientCorrection (Particle *all_particle) {
	double a00, a01, a10, a11, V, xji, yji;
	double determinant;
	vector new_grad;
	int N = NUMBER_OF_PARTICLE;

	for (Index i = 0; i < N; i++) {

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
			a11 += yji * p->Wij_grad_i.second * V;

		}

		determinant = a00 * a11 - a10 * a01;

		for (Neighbor_p p = all_particle[i].neighbors; p != NULL; p = p->next) {
			
			new_grad.first  = (  a11 * p->Wij_grad_i.first - a01 * p->Wij_grad_i.second) / determinant;
			new_grad.second = (- a10 * p->Wij_grad_i.first + a00 * p->Wij_grad_i.second) / determinant;
			p->Wij_grad_i = new_grad;
			
		}		
	}
}

#endif // KERNEL_H

