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
*	 @brief A struct containing some variables of a neighbor particle
*/
struct Neighbor{

	index  idx;				//!< global index of this neighbor particle
	double Wij;				//!< the value of kernel function
	vector Wij_grad_i;		//!< the gradient of kernel function w.r.t. the position of [particle i]

	struct Neighbor *next;			//!< pointer to the next neighbor particle.

};

typedef struct Neighbor *Neighbor_p;
typedef Neighbor_p NeighborList;


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
	
	int    	 tag;           //!< whether it's an interior particle (0), repulsive particle (1) or ghost particle (2)
		
	NeighborList   neighbors;    /*!< List containing the index array of nearby particles
									  This is the pointer to the first neighbor. */

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

#endif // DATA_SET_H


	


