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
  int unrolling_factor1 = 16;   // unrolling factor 1
  int unrolling_factor2 = 4;    // unrolling factor 2
  int block_end_i, block_end_j; // end of array that fits into block, e.g. if array size is 230, block size is 50, then block end is 200, 
                                // the array can be divided into 4 block, with 30 elements left that do not fit into a block
  int block_i, block_j;         // begin of block
  int limit_i, limit_j;         //  end  of block
  
  vector xi, xj1, xj2, xj3, xj4, xj5, xj6, xj7, xj8, xj9, xj10, xj11, xj12, xj13, xj14, xj15, xj16,
         diff1, diff2, diff3, diff4, diff5, diff6, diff7, diff8, diff9, diff10, diff11, diff12, 
         diff13, diff14, diff15, diff16;
  double r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16;
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
          xj5  = positions[j + 4];
          xj6  = positions[j + 5];
          xj7  = positions[j + 6];
          xj8  = positions[j + 7];
          xj9  = positions[j + 8];
          xj10 = positions[j + 9];
          xj11 = positions[j + 10];
          xj12 = positions[j + 11];
          xj13 = positions[j + 12];
          xj14 = positions[j + 13];
          xj15 = positions[j + 14];
          xj16 = positions[j + 15];
          diff1  = (vector){xi.first - xj1.first, xi.second - xj1.second};
          diff2  = (vector){xi.first - xj2.first, xi.second - xj2.second};
          diff3  = (vector){xi.first - xj3.first, xi.second - xj3.second};
          diff4  = (vector){xi.first - xj4.first, xi.second - xj4.second};
          diff5  = (vector){xi.first - xj5.first, xi.second - xj5.second};
          diff6  = (vector){xi.first - xj6.first, xi.second - xj6.second};
          diff7  = (vector){xi.first - xj7.first, xi.second - xj7.second};
          diff8  = (vector){xi.first - xj8.first, xi.second - xj8.second};
          diff9  = (vector){xi.first - xj9.first, xi.second - xj9.second};
          diff10 = (vector){xi.first - xj10.first, xi.second - xj10.second};
          diff11 = (vector){xi.first - xj11.first, xi.second - xj11.second};
          diff12 = (vector){xi.first - xj12.first, xi.second - xj12.second};
          diff13 = (vector){xi.first - xj13.first, xi.second - xj13.second};
          diff14 = (vector){xi.first - xj14.first, xi.second - xj14.second};
          diff15 = (vector){xi.first - xj15.first, xi.second - xj15.second};
          diff16 = (vector){xi.first - xj16.first, xi.second - xj16.second};
          r1  = diff1.first  * diff1.first  + diff1.second  * diff1.second;
          r2  = diff2.first  * diff2.first  + diff2.second  * diff2.second;
          r3  = diff3.first  * diff3.first  + diff3.second  * diff3.second;
          r4  = diff4.first  * diff4.first  + diff4.second  * diff4.second;
          r5  = diff5.first  * diff5.first  + diff5.second  * diff5.second;
          r6  = diff6.first  * diff6.first  + diff6.second  * diff6.second;
          r7  = diff7.first  * diff7.first  + diff7.second  * diff7.second;
          r8  = diff8.first  * diff8.first  + diff8.second  * diff8.second;
          r9  = diff9.first  * diff9.first  + diff9.second  * diff9.second;
          r10 = diff10.first * diff10.first + diff10.second * diff10.second;
          r11 = diff11.first * diff11.first + diff11.second * diff11.second;
          r12 = diff12.first * diff12.first + diff12.second * diff12.second;
          r13 = diff13.first * diff13.first + diff13.second * diff13.second;
          r14 = diff14.first * diff14.first + diff14.second * diff14.second;
          r15 = diff15.first * diff15.first + diff15.second * diff15.second;
          r16 = diff16.first * diff16.first + diff16.second * diff16.second;
          r1  = sqrt(r1);
          r2  = sqrt(r2);
          r3  = sqrt(r3);
          r4  = sqrt(r4);
          r5  = sqrt(r5);
          r6  = sqrt(r6);
          r7  = sqrt(r7);
          r8  = sqrt(r8);
          r9  = sqrt(r9);
          r10 = sqrt(r10);
          r11 = sqrt(r11);
          r12 = sqrt(r12);
          r13 = sqrt(r13);
          r14 = sqrt(r14);
          r15 = sqrt(r15);
          r16 = sqrt(r16);
          if(r1  < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
          if(r5  < Hradius) KernelAndGradient_bidirectional(diff5, i, j + 4, r5);
          if(r6  < Hradius) KernelAndGradient_bidirectional(diff6, i, j + 5, r6);
          if(r7  < Hradius) KernelAndGradient_bidirectional(diff7, i, j + 6, r7);
          if(r8  < Hradius) KernelAndGradient_bidirectional(diff8, i, j + 7, r8);
          if(r9  < Hradius) KernelAndGradient_bidirectional(diff9, i, j + 8, r9);
          if(r10 < Hradius) KernelAndGradient_bidirectional(diff10, i, j + 9, r10);
          if(r11 < Hradius) KernelAndGradient_bidirectional(diff11, i, j + 10, r11);
          if(r12 < Hradius) KernelAndGradient_bidirectional(diff12, i, j + 11, r12);
          if(r13 < Hradius) KernelAndGradient_bidirectional(diff13, i, j + 12, r13);
          if(r14 < Hradius) KernelAndGradient_bidirectional(diff14, i, j + 13, r14);
          if(r15 < Hradius) KernelAndGradient_bidirectional(diff15, i, j + 14, r15);
          if(r16 < Hradius) KernelAndGradient_bidirectional(diff16, i, j + 15, r16);
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
          xj5  = positions[j + 4];
          xj6  = positions[j + 5];
          xj7  = positions[j + 6];
          xj8  = positions[j + 7];
          xj9  = positions[j + 8];
          xj10 = positions[j + 9];
          xj11 = positions[j + 10];
          xj12 = positions[j + 11];
          xj13 = positions[j + 12];
          xj14 = positions[j + 13];
          xj15 = positions[j + 14];
          xj16 = positions[j + 15];
          diff1  = (vector){xi.first - xj1.first, xi.second - xj1.second};
          diff2  = (vector){xi.first - xj2.first, xi.second - xj2.second};
          diff3  = (vector){xi.first - xj3.first, xi.second - xj3.second};
          diff4  = (vector){xi.first - xj4.first, xi.second - xj4.second};
          diff5  = (vector){xi.first - xj5.first, xi.second - xj5.second};
          diff6  = (vector){xi.first - xj6.first, xi.second - xj6.second};
          diff7  = (vector){xi.first - xj7.first, xi.second - xj7.second};
          diff8  = (vector){xi.first - xj8.first, xi.second - xj8.second};
          diff9  = (vector){xi.first - xj9.first, xi.second - xj9.second};
          diff10 = (vector){xi.first - xj10.first, xi.second - xj10.second};
          diff11 = (vector){xi.first - xj11.first, xi.second - xj11.second};
          diff12 = (vector){xi.first - xj12.first, xi.second - xj12.second};
          diff13 = (vector){xi.first - xj13.first, xi.second - xj13.second};
          diff14 = (vector){xi.first - xj14.first, xi.second - xj14.second};
          diff15 = (vector){xi.first - xj15.first, xi.second - xj15.second};
          diff16 = (vector){xi.first - xj16.first, xi.second - xj16.second};
          r1  = diff1.first  * diff1.first  + diff1.second  * diff1.second;
          r2  = diff2.first  * diff2.first  + diff2.second  * diff2.second;
          r3  = diff3.first  * diff3.first  + diff3.second  * diff3.second;
          r4  = diff4.first  * diff4.first  + diff4.second  * diff4.second;
          r5  = diff5.first  * diff5.first  + diff5.second  * diff5.second;
          r6  = diff6.first  * diff6.first  + diff6.second  * diff6.second;
          r7  = diff7.first  * diff7.first  + diff7.second  * diff7.second;
          r8  = diff8.first  * diff8.first  + diff8.second  * diff8.second;
          r9  = diff9.first  * diff9.first  + diff9.second  * diff9.second;
          r10 = diff10.first * diff10.first + diff10.second * diff10.second;
          r11 = diff11.first * diff11.first + diff11.second * diff11.second;
          r12 = diff12.first * diff12.first + diff12.second * diff12.second;
          r13 = diff13.first * diff13.first + diff13.second * diff13.second;
          r14 = diff14.first * diff14.first + diff14.second * diff14.second;
          r15 = diff15.first * diff15.first + diff15.second * diff15.second;
          r16 = diff16.first * diff16.first + diff16.second * diff16.second;
          r1  = sqrt(r1);
          r2  = sqrt(r2);
          r3  = sqrt(r3);
          r4  = sqrt(r4);
          r5  = sqrt(r5);
          r6  = sqrt(r6);
          r7  = sqrt(r7);
          r8  = sqrt(r8);
          r9  = sqrt(r9);
          r10 = sqrt(r10);
          r11 = sqrt(r11);
          r12 = sqrt(r12);
          r13 = sqrt(r13);
          r14 = sqrt(r14);
          r15 = sqrt(r15);
          r16 = sqrt(r16);
          if(r1  < Hradius) KernelAndGradient_bidirectional(diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_bidirectional(diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_bidirectional(diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_bidirectional(diff4, i, j + 3, r4);
          if(r5  < Hradius) KernelAndGradient_bidirectional(diff5, i, j + 4, r5);
          if(r6  < Hradius) KernelAndGradient_bidirectional(diff6, i, j + 5, r6);
          if(r7  < Hradius) KernelAndGradient_bidirectional(diff7, i, j + 6, r7);
          if(r8  < Hradius) KernelAndGradient_bidirectional(diff8, i, j + 7, r8);
          if(r9  < Hradius) KernelAndGradient_bidirectional(diff9, i, j + 8, r9);
          if(r10 < Hradius) KernelAndGradient_bidirectional(diff10, i, j + 9, r10);
          if(r11 < Hradius) KernelAndGradient_bidirectional(diff11, i, j + 10, r11);
          if(r12 < Hradius) KernelAndGradient_bidirectional(diff12, i, j + 11, r12);
          if(r13 < Hradius) KernelAndGradient_bidirectional(diff13, i, j + 12, r13);
          if(r14 < Hradius) KernelAndGradient_bidirectional(diff14, i, j + 13, r14);
          if(r15 < Hradius) KernelAndGradient_bidirectional(diff15, i, j + 14, r15);
          if(r16 < Hradius) KernelAndGradient_bidirectional(diff16, i, j + 15, r16);
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
          xj5  = positions[j + 4];
          xj6  = positions[j + 5];
          xj7  = positions[j + 6];
          xj8  = positions[j + 7];
          xj9  = positions[j + 8];
          xj10 = positions[j + 9];
          xj11 = positions[j + 10];
          xj12 = positions[j + 11];
          xj13 = positions[j + 12];
          xj14 = positions[j + 13];
          xj15 = positions[j + 14];
          xj16 = positions[j + 15];
          diff1  = (vector){xi.first - xj1.first, xi.second - xj1.second};
          diff2  = (vector){xi.first - xj2.first, xi.second - xj2.second};
          diff3  = (vector){xi.first - xj3.first, xi.second - xj3.second};
          diff4  = (vector){xi.first - xj4.first, xi.second - xj4.second};
          diff5  = (vector){xi.first - xj5.first, xi.second - xj5.second};
          diff6  = (vector){xi.first - xj6.first, xi.second - xj6.second};
          diff7  = (vector){xi.first - xj7.first, xi.second - xj7.second};
          diff8  = (vector){xi.first - xj8.first, xi.second - xj8.second};
          diff9  = (vector){xi.first - xj9.first, xi.second - xj9.second};
          diff10 = (vector){xi.first - xj10.first, xi.second - xj10.second};
          diff11 = (vector){xi.first - xj11.first, xi.second - xj11.second};
          diff12 = (vector){xi.first - xj12.first, xi.second - xj12.second};
          diff13 = (vector){xi.first - xj13.first, xi.second - xj13.second};
          diff14 = (vector){xi.first - xj14.first, xi.second - xj14.second};
          diff15 = (vector){xi.first - xj15.first, xi.second - xj15.second};
          diff16 = (vector){xi.first - xj16.first, xi.second - xj16.second};
          r1  = diff1.first  * diff1.first  + diff1.second  * diff1.second;
          r2  = diff2.first  * diff2.first  + diff2.second  * diff2.second;
          r3  = diff3.first  * diff3.first  + diff3.second  * diff3.second;
          r4  = diff4.first  * diff4.first  + diff4.second  * diff4.second;
          r5  = diff5.first  * diff5.first  + diff5.second  * diff5.second;
          r6  = diff6.first  * diff6.first  + diff6.second  * diff6.second;
          r7  = diff7.first  * diff7.first  + diff7.second  * diff7.second;
          r8  = diff8.first  * diff8.first  + diff8.second  * diff8.second;
          r9  = diff9.first  * diff9.first  + diff9.second  * diff9.second;
          r10 = diff10.first * diff10.first + diff10.second * diff10.second;
          r11 = diff11.first * diff11.first + diff11.second * diff11.second;
          r12 = diff12.first * diff12.first + diff12.second * diff12.second;
          r13 = diff13.first * diff13.first + diff13.second * diff13.second;
          r14 = diff14.first * diff14.first + diff14.second * diff14.second;
          r15 = diff15.first * diff15.first + diff15.second * diff15.second;
          r16 = diff16.first * diff16.first + diff16.second * diff16.second;
          r1  = sqrt(r1);
          r2  = sqrt(r2);
          r3  = sqrt(r3);
          r4  = sqrt(r4);
          r5  = sqrt(r5);
          r6  = sqrt(r6);
          r7  = sqrt(r7);
          r8  = sqrt(r8);
          r9  = sqrt(r9);
          r10 = sqrt(r10);
          r11 = sqrt(r11);
          r12 = sqrt(r12);
          r13 = sqrt(r13);
          r14 = sqrt(r14);
          r15 = sqrt(r15);
          r16 = sqrt(r16);
          if(r1  < Hradius) KernelAndGradient_unidirectional(diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_unidirectional(diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_unidirectional(diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_unidirectional(diff4, i, j + 3, r4);
          if(r5  < Hradius) KernelAndGradient_unidirectional(diff5, i, j + 4, r5);
          if(r6  < Hradius) KernelAndGradient_unidirectional(diff6, i, j + 5, r6);
          if(r7  < Hradius) KernelAndGradient_unidirectional(diff7, i, j + 6, r7);
          if(r8  < Hradius) KernelAndGradient_unidirectional(diff8, i, j + 7, r8);
          if(r9  < Hradius) KernelAndGradient_unidirectional(diff9, i, j + 8, r9);
          if(r10 < Hradius) KernelAndGradient_unidirectional(diff10, i, j + 9, r10);
          if(r11 < Hradius) KernelAndGradient_unidirectional(diff11, i, j + 10, r11);
          if(r12 < Hradius) KernelAndGradient_unidirectional(diff12, i, j + 11, r12);
          if(r13 < Hradius) KernelAndGradient_unidirectional(diff13, i, j + 12, r13);
          if(r14 < Hradius) KernelAndGradient_unidirectional(diff14, i, j + 13, r14);
          if(r15 < Hradius) KernelAndGradient_unidirectional(diff15, i, j + 14, r15);
          if(r16 < Hradius) KernelAndGradient_unidirectional(diff16, i, j + 15, r16);
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
