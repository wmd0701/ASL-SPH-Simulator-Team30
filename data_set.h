//! @file data_set.h
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

// zero vector
vector zero = {0.0, 0.0};

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
 */
void KernelAndGradient_unidirectional(vector diff, int par_idx_1, int par_idx_2, double r) {
  double kernel;
  vector grad;
  double q = r * Hinv;
  double q_square = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q_square * (1. - 0.5 * q));
    if (1.0e-12 <= q) {
      temp = factor * (-3. * q + 2.25 * q_square) * Hinv / r;
      grad.first = temp * diff.first;
      grad.second = temp * diff.second;
    } else {
      grad.first = 0.;
      grad.second = 0.;
    }
  } 
  else {
    c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    grad.first = temp * diff.first;
    grad.second = temp * diff.second;
  }

  int count = neighbor_counts[par_idx_1]++;
  
  Wijs            [par_idx_1][count] = kernel;
  Wij_grads       [par_idx_1][count] = grad;
  neighbor_indices[par_idx_1][count] = par_idx_2;

}

void KernelAndGradient_bidirectional(vector diff, int par_idx_1, int par_idx_2, double r) {
  double kernel;
  vector grad;
  double q = r * Hinv;
  double q_square = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q_square * (1. - 0.5 * q));
    if (1.0e-12 <= q) {
      temp = factor * (-3. * q + 2.25 * q_square) * Hinv / r;
      grad.first = temp * diff.first;
      grad.second = temp * diff.second;
    } else {
      grad.first = 0.;
      grad.second = 0.;
    }
  } 
  else {
    c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    grad.first = temp * diff.first;
    grad.second = temp * diff.second;
  }

  int count = neighbor_counts[par_idx_1]++;
  
  Wijs            [par_idx_1][count] = kernel;
  Wij_grads       [par_idx_1][count] = grad;
  neighbor_indices[par_idx_1][count] = par_idx_2;

  grad.first = -grad.first;
  grad.second = -grad.second;
  count = neighbor_counts[par_idx_2]++;
  Wijs            [par_idx_2][count] = kernel;
  Wij_grads       [par_idx_2][count] = grad;
  neighbor_indices[par_idx_2][count] = par_idx_1;

}

void KernelAndGradient_zero(int par_idx) {
  double kernel = factor;
  vector grad   = zero;
  
  int count = neighbor_counts[par_idx]++;
  
  Wijs            [par_idx][count] = kernel;
  Wij_grads       [par_idx][count] = grad;
  neighbor_indices[par_idx][count] = par_idx;
}

void ClearNeighbors(){
  for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++)
    neighbor_counts[i] = 0;
}

void SearchNeighbors() {
  int block_size = 64;          // size of block
  int unrolling_factor1 = 4;    // unrolling factor 1
  int unrolling_factor2 = 4;    // unrolling factor 2
  int block_end_i, block_end_j; // end of array that fits into block, e.g. if array size is 230, block size is 50, then block end is 200, 
                                // the array can be divided into 4 block, with 30 elements left that do not fit into a block
  int block_i, block_j;         // begin of block
  int limit_i, limit_j;         //  end  of block
  
  vector xi, xj1, xj2, xj3, xj4, diff1, diff2, diff3, diff4;
  double r1, r2, r3, r4;
  int unrolling_limit;          // unrolling limit
  int i, j;
  

  //---------------------------------------------------------------------------
  // firstly, check interior-interior pair
  block_end_i = N_interior - N_interior % block_size;
  block_end_j = block_end_i;
  for (block_i = 0 ; block_i < block_end_i ; block_i += block_size){
    // search neighbor inside a block
    limit_i = block_i + block_size;
    limit_j = limit_i;
    unrolling_limit = limit_j - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      // each particle is neighbor to itself
      KernelAndGradient_zero(i);

      xi = positions[i];
      
      // unrolling with factor 2
      for (j = i + 1 ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = positions[j];
        xj2 = positions[j + 1];
        xj3 = positions[j + 2];
        xj4 = positions[j + 3];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
        diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
        diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
        r1 = diff1.first * diff1.first + diff1.second * diff1.second;
        r2 = diff2.first * diff2.first + diff2.second * diff2.second;
        r3 = diff3.first * diff3.first + diff3.second * diff3.second;
        r4 = diff4.first * diff4.first + diff4.second * diff4.second;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
        
      }

      // the particles left which do not fit into unrolling factor 2
      for (; j < limit_j; j++){
        xj1 = positions[j];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(diff1, i, j, r1);
      }

    }

    // search neighbor between two blocks
    for (block_j = block_i + block_size ; block_j < block_end_j ; block_j += block_size){
      limit_j = block_j + block_size;
      for (i = block_i ; i < limit_i ; i++){
        xi = positions[i];
        
        // unrolling with factor 1
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1  = positions[j];
          xj2  = positions[j + 1];
          xj3  = positions[j + 2];
          xj4  = positions[j + 3];
          diff1  = (vector){xi.first - xj1.first, xi.second - xj1.second};
          diff2  = (vector){xi.first - xj2.first, xi.second - xj2.second};
          diff3  = (vector){xi.first - xj3.first, xi.second - xj3.second};
          diff4  = (vector){xi.first - xj4.first, xi.second - xj4.second};
          r1  = diff1.first  * diff1.first  + diff1.second  * diff1.second;
          r2  = diff2.first  * diff2.first  + diff2.second  * diff2.second;
          r3  = diff3.first  * diff3.first  + diff3.second  * diff3.second;
          r4  = diff4.first  * diff4.first  + diff4.second  * diff4.second;
          r1  = sqrt(r1);
          r2  = sqrt(r2);
          r3  = sqrt(r3);
          r4  = sqrt(r4);
          if(r1  < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi = positions[i];
      
      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = positions[j];
        xj2 = positions[j + 1];
        xj3 = positions[j + 2];
        xj4 = positions[j + 3];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
        diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
        diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
        r1 = diff1.first * diff1.first + diff1.second * diff1.second;
        r2 = diff2.first * diff2.first + diff2.second * diff2.second;
        r3 = diff3.first * diff3.first + diff3.second * diff3.second;
        r4 = diff4.first * diff4.first + diff4.second * diff4.second;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
      }

      // particles left that do not fit into unrolling factor 2
      for (; j < N_interior; j++){
        xj1 = positions[j];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(diff1, i, j, r1);
      }

    }

  }
  // rest part that does not fit into block_i
  unrolling_limit = N_interior - unrolling_factor2;
  for(i = block_end_i ; i < N_interior ; i++){
    // each particle is neighbor to itself
    KernelAndGradient_zero(i);

    xi = positions[i];

    // unrolling with factor 2
    for (j = i + 1 ; j <= unrolling_limit ; j += unrolling_factor2){
      xj1 = positions[j];
      xj2 = positions[j + 1];
      xj3 = positions[j + 2];
      xj4 = positions[j + 3];
      diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
      diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
      diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
      diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
      r1 = diff1.first * diff1.first + diff1.second * diff1.second;
      r2 = diff2.first * diff2.first + diff2.second * diff2.second;
      r3 = diff3.first * diff3.first + diff3.second * diff3.second;
      r4 = diff4.first * diff4.first + diff4.second * diff4.second;
      r1 = sqrt(r1);
      r2 = sqrt(r2);
      r3 = sqrt(r3);
      r4 = sqrt(r4);
      if(r1 < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
      if(r2 < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
      if(r3 < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
      if(r4 < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
    }

    // particles left that do not fit into unrolling factor 2
    for (; j < N_interior; j++){
        xj1 = positions[j];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(diff1, i, j, r1);
    }

  }
  //---------------------------------------------------------------------------
  // secondly, check boundary-interior pair, here boundary = repulsive + ghost
  block_end_i = NUMBER_OF_PARTICLE - N_boundary % block_size;
  block_end_j = N_interior - N_interior % block_size;
  for (block_i = N_interior ; block_i < block_end_i ; block_i += block_size){
    limit_i = block_i + block_size;
    
    // each particle is neighbor to itself
    for (i = block_i ; i < limit_i ; i++)
      KernelAndGradient_zero(i);

    // search neighbor between two blocks
    for (block_j = 0 ; block_j < block_end_j ; block_j += block_size){
      limit_j = block_j + block_size;
      for (i = block_i ; i < limit_i ; i++){
        xi = positions[i];
        
        // unrolling with factor 1
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1  = positions[j];
          xj2  = positions[j + 1];
          xj3  = positions[j + 2];
          xj4  = positions[j + 3];
          diff1  = (vector){xi.first - xj1.first, xi.second - xj1.second};
          diff2  = (vector){xi.first - xj2.first, xi.second - xj2.second};
          diff3  = (vector){xi.first - xj3.first, xi.second - xj3.second};
          diff4  = (vector){xi.first - xj4.first, xi.second - xj4.second};
          r1  = diff1.first  * diff1.first  + diff1.second  * diff1.second;
          r2  = diff2.first  * diff2.first  + diff2.second  * diff2.second;
          r3  = diff3.first  * diff3.first  + diff3.second  * diff3.second;
          r4  = diff4.first  * diff4.first  + diff4.second  * diff4.second;
          r1  = sqrt(r1);
          r2  = sqrt(r2);
          r3  = sqrt(r3);
          r4  = sqrt(r4);
          if(r1  < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi = positions[i];

      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = positions[j];
        xj2 = positions[j + 1];
        xj3 = positions[j + 2];
        xj4 = positions[j + 3];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
        diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
        diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
        r1 = diff1.first * diff1.first + diff1.second * diff1.second;
        r2 = diff2.first * diff2.first + diff2.second * diff2.second;
        r3 = diff3.first * diff3.first + diff3.second * diff3.second;
        r4 = diff4.first * diff4.first + diff4.second * diff4.second;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
      }

      // particles left that do not fit into unrolling factor 2
      for (; j < N_interior; j++){
        xj1 = positions[j];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(diff1, i, j, r1);
      }
    }

  }
  // rest part that does not fit into block_i
  unrolling_limit = N_interior - unrolling_factor2;
  for(i = block_end_i ; i < NUMBER_OF_PARTICLE ; i++){
    // each particle is neighbor to itself
    KernelAndGradient_zero(i);

    xi = positions[i];

    // unrolling with factor 2
    for (j = 0 ; j <= unrolling_limit ; j += unrolling_factor2){
      xj1 = positions[j];
      xj2 = positions[j + 1];
      xj3 = positions[j + 2];
      xj4 = positions[j + 3];
      diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
      diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
      diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
      diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
      r1 = diff1.first * diff1.first + diff1.second * diff1.second;
      r2 = diff2.first * diff2.first + diff2.second * diff2.second;
      r3 = diff3.first * diff3.first + diff3.second * diff3.second;
      r4 = diff4.first * diff4.first + diff4.second * diff4.second;
      r1 = sqrt(r1);
      r2 = sqrt(r2);
      r3 = sqrt(r3);
      r4 = sqrt(r4);
      if(r1 < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
      if(r2 < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
      if(r3 < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
      if(r4 < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
    }

    // particles left that do not fit into unrolling factor 2
    for (; j < N_interior; j++){
      xj1 = positions[j];
      diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
      r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
      if(r1 < Hradius)
        KernelAndGradient_bidirectional(diff1, i, j, r1);
    }
    
  }
  //---------------------------------------------------------------------------
  // thirdly, check ghost-repulsive pair
  block_end_i = NUMBER_OF_PARTICLE - N_ghost % block_size;
  block_end_j = N_interior + N_repulsive - N_repulsive % block_size;
  for (block_i = N_interior + N_repulsive ; block_i < block_end_i ; block_i += block_size){
    limit_i = block_i + block_size;
    
    // search neighbor between two blocks
    for (block_j = N_interior ; block_j < block_end_j ; block_j += block_size){
      limit_j = block_j + block_size;
      for (i = block_i ; i < limit_i ; i++){
        xi = positions[i];

        // unrolling
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1  = positions[j];
          xj2  = positions[j + 1];
          xj3  = positions[j + 2];
          xj4  = positions[j + 3];
          diff1  = (vector){xi.first - xj1.first, xi.second - xj1.second};
          diff2  = (vector){xi.first - xj2.first, xi.second - xj2.second};
          diff3  = (vector){xi.first - xj3.first, xi.second - xj3.second};
          diff4  = (vector){xi.first - xj4.first, xi.second - xj4.second};
          r1  = diff1.first  * diff1.first  + diff1.second  * diff1.second;
          r2  = diff2.first  * diff2.first  + diff2.second  * diff2.second;
          r3  = diff3.first  * diff3.first  + diff3.second  * diff3.second;
          r4  = diff4.first  * diff4.first  + diff4.second  * diff4.second;
          r1  = sqrt(r1);
          r2  = sqrt(r2);
          r3  = sqrt(r3);
          r4  = sqrt(r4);
          if(r1  < Hradius) KernelAndGradient_unidirectional(diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_unidirectional(diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_unidirectional(diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_unidirectional(diff4, i, j + 3, r4);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior + N_repulsive - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi = positions[i];

      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = positions[j];
        xj2 = positions[j + 1];
        xj3 = positions[j + 2];
        xj4 = positions[j + 3];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
        diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
        diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
        r1 = diff1.first * diff1.first + diff1.second * diff1.second;
        r2 = diff2.first * diff2.first + diff2.second * diff2.second;
        r3 = diff3.first * diff3.first + diff3.second * diff3.second;
        r4 = diff4.first * diff4.first + diff4.second * diff4.second;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_unidirectional(diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_unidirectional(diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_unidirectional(diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_unidirectional(diff4, i, j + 3, r4);
      }

      // particles left that do not fit into unrolling factor 2
      for (; j < N_interior + N_repulsive ; j++){
        xj1 = positions[j];
        diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
        r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
        if(r1 < Hradius)
          KernelAndGradient_unidirectional(diff1, i, j, r1);
      }
    
    } 
  }

  // rest part that does not fit into block_i
  unrolling_limit = N_interior + N_repulsive - unrolling_factor2;
  for(i = block_end_i ; i < NUMBER_OF_PARTICLE ; i++){
    xi = positions[i];

    // unrolling with factor 2
    for (j = N_interior ; j <= unrolling_limit ; j += unrolling_factor2){
      xj1 = positions[j];
      xj2 = positions[j + 1];
      xj3 = positions[j + 2];
      xj4 = positions[j + 3];
      diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
      diff2 = (vector){xi.first - xj2.first, xi.second - xj2.second};
      diff3 = (vector){xi.first - xj3.first, xi.second - xj3.second};
      diff4 = (vector){xi.first - xj4.first, xi.second - xj4.second};
      r1 = diff1.first * diff1.first + diff1.second * diff1.second;
      r2 = diff2.first * diff2.first + diff2.second * diff2.second;
      r3 = diff3.first * diff3.first + diff3.second * diff3.second;
      r4 = diff4.first * diff4.first + diff4.second * diff4.second;
      r1 = sqrt(r1);
      r2 = sqrt(r2);
      r3 = sqrt(r3);
      r4 = sqrt(r4);
      if(r1 < Hradius) KernelAndGradient_unidirectional(diff1, i, j    , r1);
      if(r2 < Hradius) KernelAndGradient_unidirectional(diff2, i, j + 1, r2);
      if(r3 < Hradius) KernelAndGradient_unidirectional(diff3, i, j + 2, r3);
      if(r4 < Hradius) KernelAndGradient_unidirectional(diff4, i, j + 3, r4);
    }

    // particles left that do not fit into unrolling factor 2
    for (; j < N_interior + N_repulsive ; j++){
      xj1 = positions[j];
      diff1 = (vector){xi.first - xj1.first, xi.second - xj1.second};
      r1 = sqrt(diff1.first * diff1.first + diff1.second * diff1.second);
      if(r1 < Hradius)
        KernelAndGradient_unidirectional(diff1, i, j, r1);
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

  for (int j = 1; j < Ny_repulsive; j++)
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

#endif // DATA_SET_H
