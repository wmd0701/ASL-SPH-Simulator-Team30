//!  @file data_set.h
#ifndef DATA_SET_H
#define DATA_SET_H

#include <stdlib.h>
#include <math.h>
#include "constants.h"


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

// define inner product for vector
double vec_dot_vec(const vector v1, const vector v2){
	return v1.first * v2.first + v1.second * v2.second;
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

// define sub for vector
vector vec_sub_vec(const vector v1, const vector v2){
	vector vv = {.first = v1.first-v2.first, .second = v1.second-v2.second};
	return vv;
}

// define Euler L2 norm for vector
double vec_distance_vec(const vector v1, const vector v2){
	return sqrt(pow((v1.first  - v2.first ), 2) + pow((v1.second - v2.second), 2));
}

// define square Euler L2 norm for vector, should be used for optimization
double vec_distance_vec_square(const vector v1, const vector v2){
	return pow((v1.first  - v2.first ), 2) + pow((v1.second - v2.second), 2);
}

// tag used to tell different types of particles 
enum Particle_Type {interior, repulsive, ghost};
typedef enum Particle_Type ParticleType;

/**  
*	 @brief A struct containing some variables of a neighbor particle
*/
struct Neighbor{

	int  idx;				//!< global index of this neighbor particle
	double Wij;				//!< the value of kernel function
	vector Wij_grad_i;		//!< the gradient of kernel function w.r.t. the position of [particle i]

	struct Neighbor *next;	//!< pointer to the next neighbor particle.

};
typedef struct Neighbor Neighbor_p;


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
		
	Neighbor_p   *neighbors;    /*!< List containing the index array of nearby particles
									  This is the pointer to the first neighbor. */
								
} Particle;

/**  
*    @brief search for the neighbor particles and  allocate memory for [neighbors]
*    @param all_particle pointer to an array containing information of all the particles
*    @param ptc_idx index of the particle that is being considered
*	 @note   - search radius = 2h !!
*			 - only be called at initiation step
*			 - repulsive particles only need information of interior particles
*			 - ghost particles need information of both interior and repulsive particles
*/
void SearchNeighbors (Particle* all_particle, int ptc_idx, const int N) {
	vector xi = all_particle[ptc_idx].position;
	double r2;   // distance between two particles
	// double H_square_4 = 4 * H * H;
	Neighbor_p *p;
	all_particle[ptc_idx].neighbors = NULL;

	if (all_particle[ptc_idx].tag == interior) {
		for (int j = 0; j < N; j++) {
			if (j != ptc_idx) {   // the particle should not be neighbor of itself
				r2 = vec_distance_vec(all_particle[j].position, xi);
				if (r2 < 2*H) {
				/* 
					TODO: potential optimization:
					r2 = vec_distance_vec_square(all_particle[j].position, xi);
					if (r2 < H_square_4)
					...
				*/
					if (all_particle[ptc_idx].neighbors == NULL) {	// if it's the first pointer of list
						Neighbor_p new_p = {.idx = j, .next = NULL}; 
						p = &new_p;
                        all_particle[ptc_idx].neighbors = p;
					}
					else {	 // if it's not the first pointer of list
						Neighbor_p new_p = {.idx = j, .next = NULL};
						p->next = &new_p;
						p = p->next;
					}
				}
			}
		}
	}
	else if (all_particle[ptc_idx].tag == repulsive) {
		for (int j = 0; j < N; j++) {
			if (all_particle[ptc_idx].tag = interior) {
				r2 = vec_distance_vec(all_particle[j].position, xi);
				if (r2 < 2*H) {	  // check if neighbor
					if (all_particle[ptc_idx].neighbors == NULL) {	// if it's the first pointer of list
						Neighbor_p new_p = {.idx = j, .next = NULL}; 
						p = &new_p;
                        all_particle[ptc_idx].neighbors = p;
					}
					else {	 // if it's not the first pointer of list
						Neighbor_p new_p = {.idx = j, .next = NULL};
						p->next = &new_p;
						p = p->next;
					}
				}
			}
		}
	}
	else {	// it's a ghost particle
		for (int j = 0; j < N; j++) {
			if (all_particle[ptc_idx].tag != ghost) {   // check it's not ghost particle
				r2 = vec_distance_vec(all_particle[j].position, xi);
				if (r2 < 2*H) {
					if (all_particle[ptc_idx].neighbors == NULL) {	// if it's the first pointer of list
						Neighbor_p new_p = {.idx = j, .next = NULL}; 
						p = &new_p;
                        all_particle[ptc_idx].neighbors = p;
					}
					else {	 // if it's not the first pointer of list
						Neighbor_p new_p = {.idx = j, .next = NULL};
						p->next = &new_p;
						p = p->next;
					}
				}
			}
		}
	}
}

/**     
*		@brief Initiate the data structure as well as initial condition
*		
*			   Step 1: Memory allocation
*			   Step 2: Assign tag, position, mass and velocity to every particles
*			   Step 3: Establish the nearby relaiton 
*
*		@return	pointer to an array containing information of all the particles
*/
Particle *Init() {
  // TODO: initialization
  Particle *particles =
      (Particle *)malloc(sizeof(Particle) * NUMBER_OF_PARTICLE);
  for (int i = 0; i < 16; ++i) {
    for (int j = 0; j < 16; ++j) {
      particles[i * 16 + j].position.first = 0.5 * H * i;
      particles[i * 16 + j].position.second = 0.5 * H * j;

      if ((i > 2 && i < 13) && (j > 2 && j < 13))
        particles[i * 16 + j].tag = interior; // interior
      else if (i < 2 || i > 13 || j < 2 || j > 13)
        particles[i * 16 + j].tag = ghost; // ghost
      else {
        particles[i * 16 + j].tag = 1; // repulsive
      }
    }
  }
  int N = NUMBER_OF_PARTICLE;     // get the number of particles
  for (int i = 0; i < N; i++) { // traverse particles
    particles[i].velocity.first = 0.;
    particles[i].velocity.second = 0.;
    particles[i].mass = 10.;
    particles[i].density = initial_density;
    particles[i].pressure = 1.;
    particles[i].pressure_force.first = 0.;
    particles[i].pressure_force.second = 0.;
    particles[i].accelerat.first = 0.;
    particles[i].accelerat.second = 0.;

    SearchNeighbors(particles, i, N);
  }
  return particles;
}



#endif // DATA_SET_H


	


