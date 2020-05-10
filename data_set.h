//!  @file data_set.h
#ifndef DATA_SET_H
#define DATA_SET_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"

// struct for position, velocity, grad
typedef struct vec {
  double first;
  double second;
} vector;

// define mul for vector
vector vec_mul_scalar(const vector v, const double d) {
  vector vv = {.first = v.first * d, .second = v.second * d};
  return vv;
}

// define inner product between vectors
double vec_dot_vec(const vector v1, const vector v2) {
  return v1.first * v2.first + v1.second * v2.second;
}

// define div for vector
vector vec_div_scalar(const vector v, const double d) {
  vector vv = {.first = v.first / d, .second = v.second / d};
  return vv;
}

// define add for vector
vector vec_add_vec(const vector v1, const vector v2) {
  vector vv = {.first = v1.first + v2.first, .second = v1.second + v2.second};
  return vv;
}

// define sub for vector
vector vec_sub_vec(const vector v1, const vector v2) {
  vector vv = {.first = v1.first - v2.first, .second = v1.second - v2.second};
  return vv;
}

// define Euler distance between vectors
double vec_distance_vec(const vector v1, const vector v2) {
  return sqrt(pow((v1.first - v2.first), 2) + pow((v1.second - v2.second), 2));
}

// define square of Euler distance between vectors, should be used for
// optimization
double vec_distance_vec_square(const vector v1, const vector v2) {
  return (v1.first - v2.first) * (v1.first - v2.first) +
         (v1.second - v2.second) * (v1.second - v2.second);
}

// tag used to tell different types of particles
enum Particle_Type { interior, repulsive, ghost };
typedef enum Particle_Type ParticleType;

/**  
*	 @brief A struct containing some variables of a neighbor particle
*/
struct Neighbor {

  int idx;           //!< global index of this neighbor particle
  double Wij;        //!< the value of kernel function
  vector Wij_grad_i; //!< the gradient of kernel function w.r.t. the position of
                     //!< [particle i]
	/* 
		TODO: 
		possible improvement: use array with fixed length, instead of linked list
		e.g. int neightbors[100]
		pros: faster access, no need for delete
		cons: fixed size, cannot handle more neighbors
	*/
	struct Neighbor *next;	//!< pointer to the next neighbor particle.

};

typedef struct Neighbor Neighbor_p;
// delete a whole linked list
void deleteNeighbors(Neighbor_p **p) {
  Neighbor_p *prev = *p;
  while (*p) {
    *p = (*p)->next;
    free(prev);
    prev = *p;
  }
}

/**  
*	 @brief A struct containing some variables of a particle
*/
typedef struct {
  vector position;  //!< 2d coordinate
  double mass;      //!< mass
  vector velocity;  //!< 2d velocity
  double density;   //!< the value of density field
  double pressure;  //!< the value of pressure field
  vector accelerat; //!< acceleration of the particle, namely dv/dt

  vector position_help;  //!< help variable for Heun and Midpoint methods
  vector accelerat_help; //!< help variable for Heun and Midpoint methods

  ParticleType tag; //!< whether it's an interior particle (0), repulsive
                    //!< particle (1) or ghost particle (2)

  Neighbor_p * neighbors; /*!< List containing the index array of nearby particles
                               This is the pointer to the first neighbor. */

} Particle;

    double factor = 10 / 7 / M_PI / H / H;
    double q = r / H;
	
		kernel = factor * ( 1 - 1.5 * q * q * (1 - 0.5*q));
        if(1e-12 <= q){
            double temp = factor * (-3 * q + 9/4 * q * q);
            grad.first = temp*diff.first/(H * r);
            grad.second = temp*diff.second/(H * r);
        }
        else{
            grad.first = 0.;
            grad.second = 0.;
        }
	}
		kernel = factor * 0.25 * pow((2 - q), 3);
        double temp = - factor * 3 / 4 * (2 - q) * (2 - q);
		grad.first = temp*diff.first/(H * r);
		grad.second = temp*diff.second/(H * r); 
	}
void KernelAndGradient(vector diff, Neighbor_p *neighbor1,
                       Neighbor_p *neighbor2, double r) {
  double kernel;
  vector grad;
  if (q <= 1.) {
    }
  } else if (1. <= q && q <= 2.) {
  } else {
    printf("something wrong with searchneibors.");
    kernel = 0;
    grad.first = 0;
    grad.second = 0;
  }

  neighbor1->Wij = kernel;
  neighbor1->Wij_grad_i = grad;
  if (neighbor2 != NULL) {
    neighbor2->Wij = kernel;
    grad.first = -grad.first;
    grad.second = -grad.second;
    neighbor2->Wij_grad_i = grad;
  }
}

void AddNeighbor(Particle *all_particle, Neighbor_p **LastNeighbor, int i,
                 int j) {
  if (all_particle[i].neighbors == NULL) { // if it's the first pointer of list
    Neighbor_p *new_p = (Neighbor_p *)malloc(sizeof(Neighbor_p));
    new_p->idx = j;
    new_p->next = NULL;
    all_particle[i].neighbors = new_p;
    LastNeighbor[i] = new_p;
  } else { // if it's not the first pointer of list
    Neighbor_p *new_p = (Neighbor_p *)malloc(sizeof(Neighbor_p));
    new_p->idx = j;
    new_p->next = NULL;
    LastNeighbor[i]->next = new_p;
    LastNeighbor[i] = new_p;
  }
}



/**  
*    @brief search for the neighbor particles and  allocate memory for [neighbors]
*    @param all_particle pointer to an array containing information of all the particles
*    @param ptc_idx index of the particle that is being considered
*	 @note   - search radius = 2H
*			 - only be called at initiation step
*			 - repulsive particles only need information of interior particles
*			 - ghost particles need information of both interior and repulsive particles
*/
void SearchNeighbors(Particle *all_particle) {
    for(int i = 0; i < NUMBER_OF_PARTICLE; ++i){
        deleteNeighbors(&(all_particle[i].neighbors));
    }  
    
    double r;
    Neighbor_p *LastNeighbor[NUMBER_OF_PARTICLE];    
    
    for(int i = 0; i < NUMBER_OF_PARTICLE; ++i){
        vector xi = all_particle[i].position;
        for(int j = i; j < NUMBER_OF_PARTICLE; ++j){
            vector xj = all_particle[j].position;
            if(i == j){
                //Add particle j to neighbors of particle i
                AddNeighbor(all_particle, LastNeighbor, i, j);
                vector diff = vec_sub_vec(xi, xj);
                KernelAndGradient(diff, LastNeighbor[i], NULL, 0.);
            }
            else if(all_particle[i].tag == interior){
                vector diff = vec_sub_vec(xi, xj);
                r = sqrt(pow(diff.first, 2) + pow(diff.second, 2));
                if (r < 2 * H) {
                    //Add particle j to neighbors of particle i
                    AddNeighbor(all_particle, LastNeighbor, i, j);
                    //Add particle i to neighbors of particle j
                    AddNeighbor(all_particle, LastNeighbor, j, i);
                    //Compute kernel and its gradient
                    KernelAndGradient(diff, LastNeighbor[i], LastNeighbor[j], r);
                }
            }
            else if(all_particle[j].tag == interior){
                vector diff = vec_sub_vec(xj, xi);
                r = sqrt(pow(diff.first, 2) + pow(diff.second, 2));
                if (r < 2 * H) {
                    //Add particle i to neighbors of particle j
                    AddNeighbor(all_particle, LastNeighbor, j, i);
                    //Add particle j to neighbors of particle i
                    AddNeighbor(all_particle, LastNeighbor, i, j);
                    //Compute kernel and its gradient
                    KernelAndGradient(diff, LastNeighbor[j], LastNeighbor[i], r);
                }
            }
            else if(all_particle[i].tag == repulsive && all_particle[j].tag == ghost){
                vector diff = vec_sub_vec(xj, xi);
                r = sqrt(pow(diff.first, 2) + pow(diff.second, 2));
                if (r < 2 * H) {
                    //Add particle i to neighbors of particle j
                    AddNeighbor(all_particle, LastNeighbor, j, i);
                    //Compute kernel and its gradient
                    KernelAndGradient(diff, LastNeighbor[j], NULL, r);
                }
            }
            else if(all_particle[i].tag == ghost && all_particle[j].tag == repulsive){
                vector diff = vec_sub_vec(xi, xj);
                r = sqrt(pow(diff.first, 2) + pow(diff.second, 2));
                if (r < 2 * H) {
                    //Add particle j to neighbors of particle i
                    AddNeighbor(all_particle, LastNeighbor, i, j);
                    //Compute kernel and its gradient
                    KernelAndGradient(diff, LastNeighbor[i], NULL, r);
                }
            }
        }
    }
}


/**
 * 		@brief initialize a tank with water
 * 			   This case is corresponding to that in SHAO_2012
 * 
 * 				|	        |
 * 				|	        |
 * 				|■ ■ ■ ■ ■|
 * 				|■_■_■_■_■|
 * 
 *		@return pointer to an array containing information of all the particles
 */
Particle *Init() {
  // TODO: initialization
  Particle *particles = (Particle *)malloc(sizeof(Particle) * NUMBER_OF_PARTICLE);
  int now = 0;
  
  // Set interior particles
  for (int i = 0; i < Nx_interior; ++i) {
    for (int j = 0; j < Ny_interior; ++j) {
      particles[now].position.first = (i + 1) * H;
      particles[now].position.second = (j + 1) * H;
      particles[now].tag = interior;
      ++now;
    }
  }
	
  // Set repulsive particles
  for (int i = 0; i < Nx_boundary; i++) {
    particles[now].position.first = i * H / 2;
    particles[now].position.second = 0;
    particles[now].tag = repulsive;
    ++now;
  }
	
  for (int j = 0; j < Ny_boundary; j++) {
    particles[now].position.first = 0;
    particles[now].position.second = j * H / 2;
    particles[now].tag = repulsive;
    ++now;
  }
  for (int j = 0; j < Ny_boundary; j++) {
    particles[now].position.first = (Nx_interior + 1) * H;
    particles[now].position.second = j * H / 2;
    particles[now].tag = repulsive;
    ++now;
  }
		
  // Set ghost particles
  for (int i = -2; i < Nx_boundary + 2; i++) {
    particles[now].position.first = i * H / 2;
    particles[now].position.second = -H / 2;
    particles[now].tag = ghost;
    ++now;

    particles[now].position.first = i * H / 2;
    particles[now].position.second = -H;
    particles[now].tag = ghost;
    ++now;
  }
		
  for (int j = 0; j < Ny_boundary; j++) {
    particles[now].position.first = -H;
    particles[now].position.second = j * H / 2;
    particles[now].tag = ghost;
    ++now;

    particles[now].position.first = -H / 2;
    particles[now].position.second = j * H / 2;
    particles[now].tag = ghost;
    ++now;

    particles[now].position.first = (Nx_interior + 1.5) * H;
    particles[now].position.second = j * H / 2;
    particles[now].tag = ghost;
    ++now;

    particles[now].position.first = (Nx_interior + 2) * H;
    particles[now].position.second = j * H / 2;
    particles[now].tag = ghost;
		++now;
	}
	
  //printf("%i \n", now);
  if (NUMBER_OF_PARTICLE != now)
    printf("number of particles doesn't match with init,\n");  
	
  int N = NUMBER_OF_PARTICLE;   // get the number of particles
  for (int i = 0; i < N; i++) { // traverse particles
    particles[i].position.first += amplitude;
    particles[i].velocity.first = 0.;
    particles[i].velocity.second = 0.;
    particles[i].mass = 7 * M_PI * H * H * initial_density / 40 / 384 * 997;
    particles[i].density = initial_density;
    particles[i].pressure = 1.;
    particles[i].accelerat.first = 0.;
    particles[i].accelerat.second = 0.;
    particles[i].neighbors = NULL;
  }
  return particles;
}

/**
 * 		@brief initialize from existing data file
 * 		@return pointer to all particles
 */
Particle *Read_Init(char filename[]) {
  FILE *fp = fopen(filename, "r");
  if (!fp)
    printf("fail to read the file.\n");
  Particle *all_particle =
      (Particle *)malloc(sizeof(Particle) * NUMBER_OF_PARTICLE);
  char *str;
  double x1, x2, v1, v2, m;
  int t;
  fgets(str, 99, fp);
  for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
    fgets(str, 99, fp);
    sscanf(str, "%lf,%lf,%d,%lf,%lf,%lf", &x1, &x2, &t, &v1, &v2, &m);
    all_particle[i].position.first = x1;
    all_particle[i].position.second = x2;
    all_particle[i].tag = t;
    all_particle[i].velocity.first = v1;
    all_particle[i].velocity.second = v2;
    all_particle[i].mass = m;
    all_particle[i].neighbors = NULL;
  }
  SearchNeighbors(all_particle);
  return all_particle;
}

#endif // DATA_SET_H


	


 
