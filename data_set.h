//! @file data_set.h
#ifndef DATA_SET_H
#define DATA_SET_H

#include <stdlib.h>

typedef int index;

typedef struct{
	double first;
	double second;
} vector;


/**  
*	 @brief A struct containing some variables of a particle
*/
typedef struct {
	
	vector   position;      //!< 2d coordinate
	double   mass;          //!< mass
	vector   velocity;      //!< 2d velocity
	double   density;       //!< the value of density field 
	double   pressure;      //!< the value of pressure field
	double 	 accelerat;     //!< acceleration of the particle, namely dv/dt
	
	index    tag;           //!< whether it's an interior particle, repulsive particle or ghost particle
		
	index*   ptc_nearby;    //!< array containing the index array of nearby particles 
	double*  Wij;           /*!< array containing the value of kernel function w.r.t. the nearby particles. 
						         The sequence is corresponding to the sequence in ptc_nearby */
	vector*  Wij_grad;      /*!< 2 x N matrix, N is the total number of nearby particles.  
								 The sequence is corresponding to the sequence in ptc_nearby */
	double   L;              //!< used for kernel gradient correction : ▽W_new = L * ▽W
	
								
} Particle;


/**  
*    @brief search for the nearby particles 
*    @param all_particle pointer to an array containing information of all the particles
*    @param par_idx index of the particle that is being considered
*    @param h searching radius (smoothing length)
*    @return pointer to an array containing the index of nearby particles
*/
index* SearchNearby (Particle* all_particle, index par_idx, double h) {
	
}

#endif DATA_SET_H


	


