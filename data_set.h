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

typedef int bool;

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

  int    neighbor_index[40];
  double Wij           [40];
  vector Wij_grad      [40];
  int    num_of_neighbors;

} Particle;

void AddNeighborWithKernel(Particle* all_particle, int idx_1, int idx_2, vector diff, double r, bool flag) {
  double kernel;
  vector grad;
  double q = r * Hinv;
	double q2 = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q2 * (1. - 0.5 * q));
    if (1e-12 <= q) {
      temp = factor * (-3. * q + 2.25 * q2) * Hinv / r;
      grad.first = temp * diff.first;
      grad.second = temp * diff.second;
    } else {
      grad.first = 0.;
      grad.second = 0.;
    }
  } else if (1. <= q && q <= 2.) {
		c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    grad.first = temp * diff.first;
    grad.second = temp * diff.second;
  } else {
    printf("Something wrong with SearchNeighbors.");
    kernel = 0;
    grad.first = 0;
    grad.second = 0;
  }

  all_particle[idx_1].Wij           [all_particle[idx_1].num_of_neighbors] = kernel;
  all_particle[idx_1].Wij_grad      [all_particle[idx_1].num_of_neighbors] = grad;
  all_particle[idx_1].neighbor_index[all_particle[idx_1].num_of_neighbors] = idx_2;
  all_particle[idx_1].num_of_neighbors++;

  if (flag == 1) {
    grad.first  = -grad.first;
    grad.second = -grad.second;
    
    all_particle[idx_2].Wij           [all_particle[idx_2].num_of_neighbors] = kernel; 
    all_particle[idx_2].Wij_grad      [all_particle[idx_2].num_of_neighbors] = grad;
    all_particle[idx_2].neighbor_index[all_particle[idx_2].num_of_neighbors] = idx_1;
    all_particle[idx_2].num_of_neighbors++;
  }
}

/**  
*    @brief search for the neighbor particles and  allocate memory for [neighbors]
*    @param all_particle pointer to an array containing information of all the particles
*    @param ptc_idx index of the particle that is being considered
*	   @note   - search radius = 2H
*			 - only be called at initiation step
*			 - repulsive particles only need information of interior particles
*			 - ghost particles need information of both interior and repulsive particles
*/
void SearchNeighbors(Particle *all_particle) {
  int i, j;
  double r;
  vector xi, xj, diff;

  // - Delete neighbors list for all particles
  // - Add itself to neighbor list
  for (i = 0; i < NUMBER_OF_PARTICLE; i++) {
    all_particle[i].num_of_neighbors = 0;
    AddNeighborWithKernel(all_particle, i, i, (vector){0.0, 0.0}, 0.0, -1);
  }

  // Search neighbors of [interior] particles
  for (i = 0; i < N_interior; i++) {
    xi = all_particle[i].position;
    for (j = i + 1; j < NUMBER_OF_PARTICLE; j++) {
      xj = all_particle[j].position;
      r = vec_distance_vec(xi, xj);
      if (r < Hradius) {
        diff = vec_sub_vec(xi, xj);
        AddNeighborWithKernel(all_particle, i, j, diff, r, 1);
      }
    }  
  }
  // Search [impulsive] neighbors of [ghost] particles
  for (i = N_interior + N_replusive; i < NUMBER_OF_PARTICLE; i++) {
    xi = all_particle[i].position;
    for (j = N_interior; j < N_interior + N_replusive; j++) {
      xj = all_particle[j].position;
      r = vec_distance_vec(xi, xj);
      if (r < Hradius) {
        diff = vec_sub_vec(xi, xj);
        AddNeighborWithKernel(all_particle, i, j, diff, r, -1);
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
  }
  SearchNeighbors(all_particle);
  return all_particle;
}

#endif // DATA_SET_H


	


 
