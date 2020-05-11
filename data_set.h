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
*	 @brief A struct containing some variables of a particle
*/
typedef struct {
  vector position;  //!< 2d coordinate
  vector velocity;  //!< 2d velocity
  double density;   //!< the value of density field
  double pressure;  //!< the value of pressure field
  vector accelerat; //!< acceleration of the particle, namely dv/dt

  vector position_help;  //!< help variable for Heun and Midpoint methods
  vector accelerat_help; //!< help variable for Heun and Midpoint methods

  ParticleType tag; //!< whether it's an interior particle (0), repulsive
                    //!< particle (1) or ghost particle (2)

  double Wij[40];        //!< store the kernel
  vector Wij_grad[40];   //!< store the kernel gradient
  int    index[40];      //!< store the indices of neighbors
  int    neighbor_num;   //!< number of neighbors

} Particle;

/**
 *      @brief Add neighbors for one/two particle(s), meanwhile compute the kernel and gradient
 *      @param par_idx_1 index of particle_1
 *      @param par_idx_2 index of particle_2
 *      @param diff  position(1) - position(2)
 *      @param r distance of two particles
 *      @param flag  1: particle_1 and particle_2 are neighbors of each other
 *                  -1: particle_2 is a neighbor of particle_1 but the opposite is false.
 */
void KernelAndGradient(Particle *all_particle, vector diff, int par_idx_1,
                       int par_idx_2, double r, int flag) {
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
    printf("Something wrong with SearchnNeighbors.");
    kernel = 0;
    grad.first = 0;
    grad.second = 0;
  }

  all_particle[par_idx_1].Wij[all_particle[par_idx_1].neighbor_num] = kernel;
  all_particle[par_idx_1].Wij_grad[all_particle[par_idx_1].neighbor_num] = grad;
  all_particle[par_idx_1].index[all_particle[par_idx_1].neighbor_num] =
      par_idx_2;
  all_particle[par_idx_1].neighbor_num += 1;
  if (flag == 1) {

    grad.first = -grad.first;
    grad.second = -grad.second;

    all_particle[par_idx_2].Wij[all_particle[par_idx_2].neighbor_num] = kernel;
    all_particle[par_idx_2].Wij_grad[all_particle[par_idx_2].neighbor_num] =
        grad;
    all_particle[par_idx_2].index[all_particle[par_idx_2].neighbor_num] =
        par_idx_1;
    all_particle[par_idx_2].neighbor_num += 1;
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
  vector xi, xj, diff;
  double r;

  vector zero = {0.0, 0.0};
  for (int i = 0; i < NUMBER_OF_PARTICLE; ++i) {
    all_particle[i].neighbor_num = 0;
    KernelAndGradient(all_particle, zero, i, i, 0.0, -1);
  }

  for (int i = 0; i < N_interior; i++) {
    xi = all_particle[i].position;
    for (int j = i + 1; j < NUMBER_OF_PARTICLE; j++) {
      xj = all_particle[j].position;
      r = vec_distance_vec(xi, xj);
      if (r < Hradius) {
        diff = vec_sub_vec(xi, xj);
        KernelAndGradient(all_particle, diff, i, j, r, 1);
      }
    }
  }
  for (int i = N_interior + N_repulsive; i < NUMBER_OF_PARTICLE; i++) {
    xi = all_particle[i].position;
    for (int j = N_interior; j < N_interior + N_repulsive; j++) {
      xj = all_particle[j].position;
      r = vec_distance_vec(xi, xj);
      if (r < Hradius) {
        diff = vec_sub_vec(xi, xj);
        KernelAndGradient(all_particle, diff, i, j, r, -1);
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
  Particle *particles =
      (Particle *)malloc(sizeof(Particle) * NUMBER_OF_PARTICLE);
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

  // printf("%i \n", now);
  if (NUMBER_OF_PARTICLE != now)
    printf("number of particles doesn't match with init,\n");

  int N = NUMBER_OF_PARTICLE;   // get the number of particles
  for (int i = 0; i < N; i++) { // traverse particles
    particles[i].position.first += amplitude;
    particles[i].velocity.first = 0.;
    particles[i].velocity.second = 0.;
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
  }
  SearchNeighbors(all_particle);
  return all_particle;
}

#endif // DATA_SET_H


	


 
