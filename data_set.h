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




// instead of using array of structs, use now multiple arrays
vector* positions;
vector* velocities;
double* densities;
double* pressures;
vector* accelerats;
double** Wijs;
vector** Wij_grads;
int** neighbor_indices;
int* neighbor_counts;

// initial positions and velocities of boundary particles, used for DisplaceBoundaries in rate_of_change.h
vector* init_positions;
vector* init_velocities;







/**
 *      @brief Add neighbors for one/two particle(s), meanwhile compute the kernel and gradient
 *      @param par_idx_1 index of particle_1
 *      @param par_idx_2 index of particle_2
 *      @param diff  position(1) - position(2)
 *      @param r distance of two particles
 *      @param bi_directional: 1 if two particles are neighbor to each other, otherwise 0
 */
void KernelAndGradient(vector diff, int par_idx_1, int par_idx_2, double r, int bi_directional) {
  double kernel;
  vector grad;
  double q = r * Hinv;
  double q2 = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q2 * (1. - 0.5 * q));
    if (1.0e-12 <= q) {
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

  int count = neighbor_counts[par_idx_1]++;
  // TODO: remove following 2 lines
  if(count >= max_num_neighbors - 1)
    printf("too much neighbors!!!\n");

  Wijs            [par_idx_1][count] = kernel;
  Wij_grads       [par_idx_1][count] = grad;
  neighbor_indices[par_idx_1][count] = par_idx_2;

  if (bi_directional == 1) {
    grad.first = -grad.first;
    grad.second = -grad.second;
    
    count = neighbor_counts[par_idx_2]++;
    Wijs            [par_idx_2][count] = kernel;
    Wij_grads       [par_idx_2][count] = grad;
    neighbor_indices[par_idx_2][count] = par_idx_1;
  }

}

void ClearNeighbors(){
  for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++)
    neighbor_counts[i] = 0;
}

void SearchNeighbors() {
  vector xi, xj, diff;
  double r;
  vector zero = {0.0, 0.0};
  
  // firstly, each particle is neighbor to itself
  for (int i = 0; i < NUMBER_OF_PARTICLE; ++i)
    KernelAndGradient(zero, i, i, 0.0, 0);

  // secondly, check interior-other pair
  for (int i = 0; i < N_interior; i++) {
    xi = positions[i];
    for (int j = i + 1; j < NUMBER_OF_PARTICLE; j++) {
      xj = positions[j];
      r = vec_distance_vec(xi, xj);
      if (r < Hradius) {
        diff = vec_sub_vec(xi, xj);
        KernelAndGradient(diff, i, j, r, 1);
      }
    }
  }

  // lastly, check repulsive-ghost pair
  for (int i = N_interior + N_repulsive; i < NUMBER_OF_PARTICLE; i++) {
    xi = positions[i];
    for (int j = N_interior; j < N_interior + N_repulsive; j++) {
      xj = positions[j];
      r = vec_distance_vec(xi, xj);
      if (r < Hradius) {
        diff = vec_sub_vec(xi, xj);
        KernelAndGradient(diff, i, j, r, 0);
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
void *Init() {
  // use calloc instead of malloc, so that initial values are set to 0 by default
  positions  = (vector*)calloc(NUMBER_OF_PARTICLE, sizeof(vector));
  velocities = (vector*)calloc(NUMBER_OF_PARTICLE, sizeof(vector));
  densities  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  pressures  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  accelerats = (vector*)calloc(NUMBER_OF_PARTICLE, sizeof(vector));
  
  Wijs             = (double**)calloc(NUMBER_OF_PARTICLE, sizeof(double*));
  Wij_grads        = (vector**)calloc(NUMBER_OF_PARTICLE, sizeof(vector*));
  neighbor_indices = (int**   )calloc(NUMBER_OF_PARTICLE, sizeof(int*));
  for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
    Wijs[i]             = (double*)calloc(max_num_neighbors, sizeof(double));
    Wij_grads[i]        = (vector*)calloc(max_num_neighbors, sizeof(vector));
    neighbor_indices[i] = (int*   )calloc(max_num_neighbors, sizeof(int)); 
  }
  neighbor_counts = (int*)calloc(NUMBER_OF_PARTICLE, sizeof(int));

  init_positions  = (vector*)calloc(N_boundary, sizeof(vector));
  init_velocities = (vector*)calloc(N_boundary, sizeof(vector));

  int now = 0;

  // Set interior particles
  for (int i = 0; i < Nx_interior; ++i)
    for (int j = 0; j < Ny_interior; ++j)
      positions[now++] = (vector){(i + 1) * H, (j + 1) * H};

  // Set repulsive particles
  for (int i = 0; i < Nx_repulsive; i++)
    positions[now++] = (vector){i * H / 2., 0};

  for (int j = 0; j < Ny_repulsive; j++)
    positions[now++] = (vector){0, j * H / 2.};

  for (int j = 0; j < Ny_repulsive; j++)
    positions[now++] = (vector){(Nx_interior + 1) * H, j * H / 2.};

  // Set ghost particles
  for (int i = -2; i < Nx_repulsive + 2; i++) {
    positions[now++] = (vector){i * H / 2., -H / 2.};
    positions[now++] = (vector){i * H / 2., -H};
  }

  for (int j = 0; j < Ny_repulsive; j++) {
    positions[now++] = (vector){-H, j * H / 2.};
    positions[now++] = (vector){-H / 2., j * H / 2.};
    positions[now++] = (vector){(Nx_interior + 1.5) * H, j * H / 2.};
    positions[now++] = (vector){(Nx_interior + 2.) * H, j * H / 2.};
  }

  if (NUMBER_OF_PARTICLE != now)
    printf("number of particles doesn't match with init,\n");

  for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
    positions[i].first += amplitude;
    densities[i] = initial_density;
    pressures[i] = 1.;
  }

  // copy initial values
  for(int i = 0 ; i < N_boundary ; i++){
    init_positions [i] = positions [i + N_interior];
    init_velocities[i] = velocities[i + N_interior];
  }

}

void Destroy(){
  free(positions);
  free(velocities);
  free(densities);
  free(pressures);
  free(accelerats);

  for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
    free(Wijs[i]);
    free(Wij_grads[i]);
    free(neighbor_indices[i]);
  }
  free(Wijs);
  free(Wij_grads);
  free(neighbor_indices);
  free(neighbor_counts);

  free(init_positions);
  free(init_velocities);
}

/**
 * 		@brief initialize from existing data file
 * 		@return pointer to all particles
 */

/*
  !!! comment out Read_Init to simplify work

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
*/

#endif // DATA_SET_H


	


 
