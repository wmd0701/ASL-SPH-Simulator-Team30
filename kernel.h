//!  @file kernel.h

#define M_PI 3.141592653579

#include "data_set.h"

#include <math.h>

/**  
*    @brief Compute the kernel function
*    @param xi coordinate of the particle which is being considered now
*    @param xj coordinate of the particle nearby
*    @param h smoothing length 
*	 @return value of kernel function
*/
double Kernel (const vector xi,const vector xj, double h) {
	
	double r = sqrt(pow((xi.first - xj.first),2) + pow((xi.second - xj.second),2));
	double epsilon = r / h;
	if (r <= 1) {
		return (1 - 3/2*epsilon*epsilon + 3/4*epsilon*epsilon*epsilon) / M_PI / (h*h*h);
	}
	else {
		return 0.25 * (2 - epsilon) * (2 - epsilon) * (2 - epsilon) / M_PI / (h*h*h); 
	}
}

/** 
*   @brief Compute the kernel for every particle w.r.t. its nearby particles
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [Wij] attribute in all_particle
*/   
vector ComputeGlobalKernel (Particle* all_particle) {
	
}


/** 
*   @brief Compute the gradient of the kernel function
*   @param xi coordinate of the particle which is being considered now
*   @param xj coordinate of the particle nearby
*   @param h smoothing length 
*	@return gradient of kernel function
*/   
vector KernelGradient (const vector xi,const vector xj, double h) {
	
}

/** 
*   @brief Compute the kernel gradient for every particle w.r.t. its nearby particles
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [Wij_grad] attribute in all_particle
*/   
vector ComputeGlobalKernelGradient (Particle* all_particle) {
	
}


/** 
*   @brief Compute the "corrected" gradient of the kernel function
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the L attribute in all_particle, then update the Wij_grad attribute.
*/   
void KernelGradientCorrection (Particle* all_particle) {
	
}

