//!  @file kernel.h

/**  
*    @brief Compute the kernel function
*    @param xi coordinate of the particle which is being considered now
*    @param xj coordinate of the particle nearby
*    @param h smoothing length 
*	 @return value of kernel function
*/
double Kernel (const double*& xi,const double*& xj, double h) {
	
}

/** 
*   @brief Compute the kernel for every particle w.r.t. its nearby particles
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [Wij] attribute in all_particle
*/   
vector ComputeGlobalKernel (particle*& all_particle) {
	
}


/** 
*   @brief Compute the gradient of the kernel function
*   @param xi coordinate of the particle which is being considered now
*   @param xj coordinate of the particle nearby
*   @param h smoothing length 
*	@return gradient of kernel function
*/   
vector KernelGradient (const double*& xi,const double*& xj, double h) {
	
}

/** 
*   @brief Compute the kernel gradient for every particle w.r.t. its nearby particles
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the [Wij_grad] attribute in all_particle
*/   
vector ComputeGlobalKernelGradient (particle*& all_particle) {
	
}


/** 
*   @brief Compute the "corrected" gradient of the kernel function
*	@param all_particle pointer to an array containing information of all the particles
*   @return no returns. Update the L attribute in all_particle, then update the Wij_grad attribute.
*/   
void KernelGradientCorrection (particle*& all_particle) {
	
}

