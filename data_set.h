//!  @file data_set.h
#ifndef DATA_SET_H
#define DATA_SET_H

#include <stdlib.h>
#include <math.h>

// index for the particles
typedef int index;  

// struct for position, velocity, grad
typedef struct vec{
	double first;
	double second;
} vector;

// define mul for vector
vector vec_mul_scalar(const vector v, const double d){
	vector vv = {.first = v.first*d, .second = v.second*d};
	return vv; 
}

// define div for vector
vector vec_div_scalar(const vector v, const double d){
	vector vv = {.first = v.first/d, .second = v.second/d};
	return vv;
}

// define add for vector
vector vec_add_vec(const vector v1, const vector v2){
	vector vv = {.first = v1.first+v2.first, .second = v1.second+v2.second};
	return vv;
}

// tag used to tell different types of particles 
enum Particle_Type {interior, repulsive, ghost};
typedef enum Particle_Type ParticleType;

/**  
*	 @brief A struct containing some variables of a neighbor particle
*/
struct Neighbor{

	index  idx;				//!< global index of this neighbor particle
	double Wij;				//!< the value of kernel function
	vector Wij_grad_i;		//!< the gradient of kernel function w.r.t. the position of [particle i]

	struct Neighbor *next;	//!< pointer to the next neighbor particle.

};
typedef struct Neighbor *Neighbor_p;
typedef struct Neighbor *NeighborList;


/**  
*	 @brief A struct containing some variables of a particle
*/
typedef struct {
	
	vector   position;      //!< 2d coordinate
	double   mass;          //!< mass
	vector   velocity;      //!< 2d velocity
	double   density;       //!< the value of density field 
	double   pressure;      //!< the value of pressure field
	vector	 pressure_force;//!< pressure force
	vector 	 accelerat;     //!< acceleration of the particle, namely dv/dt
	
	ParticleType   tag;          //!< whether it's an interior particle (0), repulsive particle (1) or ghost particle (2)
		
	NeighborList   neighbors;    /*!< List containing the index array of nearby particles
									  This is the pointer to the first neighbor. */
								
} Particle;

/**     
*		@brief Initiate the data structure as well as initial condition
*		
*			   Step 1: Memory allocation
*			   Step 2: Assign tag, position, mass and velocity to every particles
*			   Step 3: Establish the nearby relaiton 
*
*		@return	pointer to an array containing information of all the particles
*/
Particle* Init(){
	// TODO: initialization
	Particle* particles = (Particle*)malloc(sizeof(Particle)*NUMBER_OF_PARTICLE);
	return particles;
}

/**  
*    @brief search for the neighbor particles and  allocate memory for [neighbors]
*    @param all_particle pointer to an array containing information of all the particles
*    @param ptc_idx index of the particle that is being considered
*	 @note   - search radius = 2h !!
*			 - only be called at initiation step
*			 - repulsive particles only need information of interior particles
*			 - ghost particles need information of both interior and repulsive particles
*/
void SearchNeighbors (Particle* all_particle, index ptc_idx) {
	vector xi = all_particle[ptc_idx].position;
	double r2;   // distance of two particles
	Neighbor_p p, tmp;
	int N = sizeof(all_particle) / sizeof(all_particle[0]);

	if (all_particle[ptc_idx].tag == interior) {
		for (index j = 0; j < N; j++) {
			if (j != ptc_idx) {   // check if itself
				r2 = sqrt(pow((all_particle[j].position.first  - xi.first ), 2) +     \
						  pow((all_particle[j].position.second - xi.second), 2));
				if (r2 < 2*H) {	  // check if neighbor
					if (all_particle[ptc_idx].neighbors == NULL) {	// if it's the first pointer of list
						p = (Neighbor_p)malloc(sizeof(struct Neighbor));
						p->idx = j;
						p->next = NULL;
					}
					else {	 // if it's not the first pointer of list
						tmp = p;
						p = (Neighbor_p)malloc(sizeof(struct Neighbor));
						p->idx = j;
						p->next = NULL;
						tmp->next = p;
					}
				}
			}
		}
	}
	else if (all_particle[ptc_idx].tag == repulsive) {
		for (index j = 0; j < N; j++) {
			if (all_particle[ptc_idx].tag = interior) {   // check it's interior particle
				r2 = sqrt(pow((all_particle[j].position.first  - xi.first ), 2) +     \
						  pow((all_particle[j].position.second - xi.second), 2));
				if (r2 < 2*H) {	  // check if neighbor
					if (all_particle[ptc_idx].neighbors == NULL) {	// if it's the first pointer of list
						p = (Neighbor_p)malloc(sizeof(struct Neighbor));
						p->idx = j;
						p->next = NULL;
					}
					else {	 // if it's not the first pointer of list
						tmp = p;
						p = (Neighbor_p)malloc(sizeof(struct Neighbor));
						p->idx = j;
						p->next = NULL;
						tmp->next = p;
					}
				}
			}
		}
	}
	else {	// it's a ghost particle
		for (index j = 0; j < N; j++) {
			if (all_particle[ptc_idx].tag != ghost) {   // check it's not ghost particle
				r2 = sqrt(pow((all_particle[j].position.first  - xi.first ), 2) +     \
						  pow((all_particle[j].position.second - xi.second), 2));
				if (r2 < 2*H) {	  // check if neighbor
					if (all_particle[ptc_idx].neighbors == NULL) {	// if it's the first pointer of list
						p = (Neighbor_p)malloc(sizeof(struct Neighbor));
						p->idx = j;
						p->next = NULL;
					}
					else {	 // if it's not the first pointer of list
						tmp = p;
						p = (Neighbor_p)malloc(sizeof(struct Neighbor));
						p->idx = j;
						p->next = NULL;
						tmp->next = p;
					}
				}
			}
		}
	}
}



#endif // DATA_SET_H


	


