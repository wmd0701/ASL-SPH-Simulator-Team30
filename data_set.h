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

void KernelAndGradient_zero(int par_idx) {
  double kernel = factor;
  vector grad   = zero;
  
  int count = neighbor_counts[par_idx]++;
  
  Wijs            [par_idx][count] = kernel;
  Wij_grads       [par_idx][count] = grad;
  neighbor_indices[par_idx][count] = par_idx;
}

/**
 *      @brief Add neighbors for one/two particle(s), meanwhile compute the kernel and gradient
 *      @param par_idx_1 index of particle_1
 *      @param par_idx_2 index of particle_2
 *      @param diff  position(1) - position(2)
 *      @param r distance of two particles
 */
void KernelAndGradient_unidirectional(double x_diff, double y_diff, int par_idx_1, int par_idx_2, double r) {
  double kernel;
  vector grad;
  double q = r * Hinv;
  double q_square = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q_square * (1. - 0.5 * q));
    if (1.0e-12 <= q) {
      temp = factor * (-3. * q + 2.25 * q_square) * Hinv / r;
      grad.first = temp * x_diff;
      grad.second = temp * y_diff;
    } else {
      grad.first = 0.;
      grad.second = 0.;
    }
  } 
  else {
    c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    grad.first = temp * x_diff;
    grad.second = temp * y_diff;
  }

  int count = neighbor_counts[par_idx_1]++;
  
  Wijs            [par_idx_1][count] = kernel;
  Wij_grads       [par_idx_1][count] = grad;
  neighbor_indices[par_idx_1][count] = par_idx_2;

}

void KernelAndGradient_bidirectional(double x_diff, double y_diff, int par_idx_1, int par_idx_2, double r) {
  double kernel;
  vector grad;
  double q = r * Hinv;
  double q_square = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q_square * (1. - 0.5 * q));
    if (1.0e-12 <= q) {
      temp = factor * (-3. * q + 2.25 * q_square) * Hinv / r;
      grad.first = temp * x_diff;
      grad.second = temp * y_diff;
    } else {
      grad.first = 0.;
      grad.second = 0.;
    }
  } 
  else {
    c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    grad.first = temp * x_diff;
    grad.second = temp * y_diff;
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
  
  double x_diff1, x_diff2, x_diff3, x_diff4, x_diff5, x_diff6, x_diff7, x_diff8, x_diff9, x_diff10, x_diff11, x_diff12, 
         x_diff13, x_diff14, x_diff15, x_diff16;
  double y_diff1, y_diff2, y_diff3, y_diff4, y_diff5, y_diff6, y_diff7, y_diff8, y_diff9, y_diff10, y_diff11, y_diff12, 
         y_diff13, y_diff14, y_diff15, y_diff16;
  double r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16;
  double xi, xj1, xj2, xj3, xj4, xj5, xj6, xj7, xj8, xj9, xj10, xj11, xj12, xj13, xj14, xj15, xj16;
  double yi, yj1, yj2, yj3, yj4, yj5, yj6, yj7, yj8, yj9, yj10, yj11, yj12, yj13, yj14, yj15, yj16;
  int unrolling_limit;          // unrolling limit
  int i, j;
  
  //--------------------------------------------------------------------
  // Copy the positions in two arrays: x_positions, y_positions
  double x_positions[NUMBER_OF_PARTICLE];
  double y_positions[NUMBER_OF_PARTICLE];
  for(int i = 0; i < NUMBER_OF_PARTICLE; ++i){
      vector position = positions[i];
      x_positions[i] = position.first;
      y_positions[i] = position.second;
  }

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

      xi = x_positions[i];
      yi = y_positions[i];
      
      // unrolling with factor 2
      for (j = i + 1 ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = x_positions[j];
        xj2 = x_positions[j + 1];
        xj3 = x_positions[j + 2];
        xj4 = x_positions[j + 3];
        yj1 = y_positions[j];
        yj2 = y_positions[j + 1];
        yj3 = y_positions[j + 2];
        yj4 = y_positions[j + 3];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        x_diff2 = xi - xj2;
        y_diff2 = yi - yj2;
        x_diff3 = xi - xj3;
        y_diff3 = yi - yj3;
        x_diff4 = xi - xj4;
        y_diff4 = yi - yj4;
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
        r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
        r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
        
      }

      // the particles left which do not fit into unrolling factor 2
      for (; j < limit_j; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
      }

    }

    // search neighbor between two blocks
    for (block_j = block_i + block_size ; block_j < block_end_j ; block_j += block_size){
      limit_j = block_j + block_size;
      for (i = block_i ; i < limit_i ; i++){
        xi = x_positions[i];
        yi = y_positions[i];
        
        // unrolling with factor 1
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1  = x_positions[j];
          xj2  = x_positions[j + 1];
          xj3  = x_positions[j + 2];
          xj4  = x_positions[j + 3];
          xj5  = x_positions[j + 4];
          xj6  = x_positions[j + 5];
          xj7  = x_positions[j + 6];
          xj8  = x_positions[j + 7];
          xj9  = x_positions[j + 8];
          xj10 = x_positions[j + 9];
          xj11 = x_positions[j + 10];
          xj12 = x_positions[j + 11];
          xj13 = x_positions[j + 12];
          xj14 = x_positions[j + 13];
          xj15 = x_positions[j + 14];
          xj16 = x_positions[j + 15];
          yj1  = y_positions[j];
          yj2  = y_positions[j + 1];
          yj3  = y_positions[j + 2];
          yj4  = y_positions[j + 3];
          yj5  = y_positions[j + 4];
          yj6  = y_positions[j + 5];
          yj7  = y_positions[j + 6];
          yj8  = y_positions[j + 7];
          yj9  = y_positions[j + 8];
          yj10 = y_positions[j + 9];
          yj11 = y_positions[j + 10];
          yj12 = y_positions[j + 11];
          yj13 = y_positions[j + 12];
          yj14 = y_positions[j + 13];
          yj15 = y_positions[j + 14];
          yj16 = y_positions[j + 15];
          x_diff1 = xi - xj1;
          y_diff1 = yi - yj1;
          x_diff2 = xi - xj2;
          y_diff2 = yi - yj2;
          x_diff3 = xi - xj3;
          y_diff3 = yi - yj3;
          x_diff4 = xi - xj4;
          y_diff4 = yi - yj4;
          x_diff5 = xi - xj5;
          y_diff5 = yi - yj5;
          x_diff6 = xi - xj6;
          y_diff6 = yi - yj6;
          x_diff7 = xi - xj7;
          y_diff7 = yi - yj7;
          x_diff8 = xi - xj8;
          y_diff8 = yi - yj8;
          x_diff9 = xi - xj9;
          y_diff9 = yi - yj9;
          x_diff10 = xi - xj10;
          y_diff10 = yi - yj10;
          x_diff11 = xi - xj11;
          y_diff11 = yi - yj11;
          x_diff12 = xi - xj12;
          y_diff12 = yi - yj12;
          x_diff13 = xi - xj13;
          y_diff13 = yi - yj13;
          x_diff14 = xi - xj14;
          y_diff14 = yi - yj14;
          x_diff15 = xi - xj15;
          y_diff15 = yi - yj15;
          x_diff16 = xi - xj16;
          y_diff16 = yi - yj16;
          r1  = x_diff1  * x_diff1  + y_diff1  * y_diff1;
          r2  = x_diff2  * x_diff2  + y_diff2  * y_diff2;
          r3  = x_diff3  * x_diff3  + y_diff3  * y_diff3;
          r4  = x_diff4  * x_diff4  + y_diff4  * y_diff4;
          r5  = x_diff5  * x_diff5  + y_diff5  * y_diff5;
          r6  = x_diff6  * x_diff6  + y_diff6  * y_diff6;
          r7  = x_diff7  * x_diff7  + y_diff7  * y_diff7;
          r8  = x_diff8  * x_diff8  + y_diff8  * y_diff8;
          r9  = x_diff9  * x_diff9  + y_diff9  * y_diff9;
          r10 = x_diff10 * x_diff10 + y_diff10 * y_diff10;
          r11 = x_diff11 * x_diff11 + y_diff11 * y_diff11;
          r12 = x_diff12 * x_diff12 + y_diff12 * y_diff12;
          r13 = x_diff13 * x_diff13 + y_diff13 * y_diff13;
          r14 = x_diff14 * x_diff14 + y_diff14 * y_diff14;
          r15 = x_diff15 * x_diff15 + y_diff15 * y_diff15;
          r16 = x_diff16 * x_diff16 + y_diff16 * y_diff16;
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
          if(r1  < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
          if(r5  < Hradius) KernelAndGradient_bidirectional(x_diff5, y_diff5, i, j + 4, r5);
          if(r6  < Hradius) KernelAndGradient_bidirectional(x_diff6, y_diff6, i, j + 5, r6);
          if(r7  < Hradius) KernelAndGradient_bidirectional(x_diff7, y_diff7, i, j + 6, r7);
          if(r8  < Hradius) KernelAndGradient_bidirectional(x_diff8, y_diff8, i, j + 7, r8);
          if(r9  < Hradius) KernelAndGradient_bidirectional(x_diff9, y_diff9, i, j + 8, r9);
          if(r10 < Hradius) KernelAndGradient_bidirectional(x_diff10, y_diff10, i, j + 9, r10);
          if(r11 < Hradius) KernelAndGradient_bidirectional(x_diff11, y_diff11, i, j + 10, r11);
          if(r12 < Hradius) KernelAndGradient_bidirectional(x_diff12, y_diff12, i, j + 11, r12);
          if(r13 < Hradius) KernelAndGradient_bidirectional(x_diff13, y_diff13, i, j + 12, r13);
          if(r14 < Hradius) KernelAndGradient_bidirectional(x_diff14, y_diff14, i, j + 13, r14);
          if(r15 < Hradius) KernelAndGradient_bidirectional(x_diff15, y_diff15, i, j + 14, r15);
          if(r16 < Hradius) KernelAndGradient_bidirectional(x_diff16, y_diff16, i, j + 15, r16);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi = x_positions[i];
      yi = y_positions[i];
      
      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = x_positions[j];
        xj2 = x_positions[j + 1];
        xj3 = x_positions[j + 2];
        xj4 = x_positions[j + 3];
        yj1 = y_positions[j];
        yj2 = y_positions[j + 1];
        yj3 = y_positions[j + 2];
        yj4 = y_positions[j + 3];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        x_diff2 = xi - xj2;
        y_diff2 = yi - yj2;
        x_diff3 = xi - xj3;
        y_diff3 = yi - yj3;
        x_diff4 = xi - xj4;
        y_diff4 = yi - yj4;
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
        r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
        r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
      }

      // particles left that do not fit into unrolling factor 2
      for (; j < N_interior; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
      }

    }

  }
  // rest part that does not fit into block_i
  unrolling_limit = N_interior - unrolling_factor2;
  for(i = block_end_i ; i < N_interior ; i++){
    // each particle is neighbor to itself
    KernelAndGradient_zero(i);

    xi = x_positions[i];
    yi = y_positions[i];

    // unrolling with factor 2
    for (j = i + 1 ; j <= unrolling_limit ; j += unrolling_factor2){
      xj1 = x_positions[j];
      xj2 = x_positions[j + 1];
      xj3 = x_positions[j + 2];
      xj4 = x_positions[j + 3];
      yj1 = y_positions[j];
      yj2 = y_positions[j + 1];
      yj3 = y_positions[j + 2];
      yj4 = y_positions[j + 3];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      x_diff2 = xi - xj2;
      y_diff2 = yi - yj2;
      x_diff3 = xi - xj3;
      y_diff3 = yi - yj3;
      x_diff4 = xi - xj4;
      y_diff4 = yi - yj4;
      r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
      r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
      r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
      r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
      r1 = sqrt(r1);
      r2 = sqrt(r2);
      r3 = sqrt(r3);
      r4 = sqrt(r4);
      if(r1 < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
      if(r2 < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
      if(r3 < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
      if(r4 < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
    }

    // particles left that do not fit into unrolling factor 2
    for (; j < N_interior; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
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
        xi = x_positions[i];
        yi = y_positions[i];
        
        // unrolling with factor 1
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1  = x_positions[j];
          xj2  = x_positions[j + 1];
          xj3  = x_positions[j + 2];
          xj4  = x_positions[j + 3];
          xj5  = x_positions[j + 4];
          xj6  = x_positions[j + 5];
          xj7  = x_positions[j + 6];
          xj8  = x_positions[j + 7];
          xj9  = x_positions[j + 8];
          xj10 = x_positions[j + 9];
          xj11 = x_positions[j + 10];
          xj12 = x_positions[j + 11];
          xj13 = x_positions[j + 12];
          xj14 = x_positions[j + 13];
          xj15 = x_positions[j + 14];
          xj16 = x_positions[j + 15];
          yj1  = y_positions[j];
          yj2  = y_positions[j + 1];
          yj3  = y_positions[j + 2];
          yj4  = y_positions[j + 3];
          yj5  = y_positions[j + 4];
          yj6  = y_positions[j + 5];
          yj7  = y_positions[j + 6];
          yj8  = y_positions[j + 7];
          yj9  = y_positions[j + 8];
          yj10 = y_positions[j + 9];
          yj11 = y_positions[j + 10];
          yj12 = y_positions[j + 11];
          yj13 = y_positions[j + 12];
          yj14 = y_positions[j + 13];
          yj15 = y_positions[j + 14];
          yj16 = y_positions[j + 15];
          x_diff1 = xi - xj1;
          y_diff1 = yi - yj1;
          x_diff2 = xi - xj2;
          y_diff2 = yi - yj2;
          x_diff3 = xi - xj3;
          y_diff3 = yi - yj3;
          x_diff4 = xi - xj4;
          y_diff4 = yi - yj4;
          x_diff5 = xi - xj5;
          y_diff5 = yi - yj5;
          x_diff6 = xi - xj6;
          y_diff6 = yi - yj6;
          x_diff7 = xi - xj7;
          y_diff7 = yi - yj7;
          x_diff8 = xi - xj8;
          y_diff8 = yi - yj8;
          x_diff9 = xi - xj9;
          y_diff9 = yi - yj9;
          x_diff10 = xi - xj10;
          y_diff10 = yi - yj10;
          x_diff11 = xi - xj11;
          y_diff11 = yi - yj11;
          x_diff12 = xi - xj12;
          y_diff12 = yi - yj12;
          x_diff13 = xi - xj13;
          y_diff13 = yi - yj13;
          x_diff14 = xi - xj14;
          y_diff14 = yi - yj14;
          x_diff15 = xi - xj15;
          y_diff15 = yi - yj15;
          x_diff16 = xi - xj16;
          y_diff16 = yi - yj16;
          r1  = x_diff1  * x_diff1  + y_diff1  * y_diff1;
          r2  = x_diff2  * x_diff2  + y_diff2  * y_diff2;
          r3  = x_diff3  * x_diff3  + y_diff3  * y_diff3;
          r4  = x_diff4  * x_diff4  + y_diff4  * y_diff4;
          r5  = x_diff5  * x_diff5  + y_diff5  * y_diff5;
          r6  = x_diff6  * x_diff6  + y_diff6  * y_diff6;
          r7  = x_diff7  * x_diff7  + y_diff7  * y_diff7;
          r8  = x_diff8  * x_diff8  + y_diff8  * y_diff8;
          r9  = x_diff9  * x_diff9  + y_diff9  * y_diff9;
          r10 = x_diff10 * x_diff10 + y_diff10 * y_diff10;
          r11 = x_diff11 * x_diff11 + y_diff11 * y_diff11;
          r12 = x_diff12 * x_diff12 + y_diff12 * y_diff12;
          r13 = x_diff13 * x_diff13 + y_diff13 * y_diff13;
          r14 = x_diff14 * x_diff14 + y_diff14 * y_diff14;
          r15 = x_diff15 * x_diff15 + y_diff15 * y_diff15;
          r16 = x_diff16 * x_diff16 + y_diff16 * y_diff16;
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
          if(r1  < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
          if(r5  < Hradius) KernelAndGradient_bidirectional(x_diff5, y_diff5, i, j + 4, r5);
          if(r6  < Hradius) KernelAndGradient_bidirectional(x_diff6, y_diff6, i, j + 5, r6);
          if(r7  < Hradius) KernelAndGradient_bidirectional(x_diff7, y_diff7, i, j + 6, r7);
          if(r8  < Hradius) KernelAndGradient_bidirectional(x_diff8, y_diff8, i, j + 7, r8);
          if(r9  < Hradius) KernelAndGradient_bidirectional(x_diff9, y_diff9, i, j + 8, r9);
          if(r10 < Hradius) KernelAndGradient_bidirectional(x_diff10, y_diff10, i, j + 9, r10);
          if(r11 < Hradius) KernelAndGradient_bidirectional(x_diff11, y_diff11, i, j + 10, r11);
          if(r12 < Hradius) KernelAndGradient_bidirectional(x_diff12, y_diff12, i, j + 11, r12);
          if(r13 < Hradius) KernelAndGradient_bidirectional(x_diff13, y_diff13, i, j + 12, r13);
          if(r14 < Hradius) KernelAndGradient_bidirectional(x_diff14, y_diff14, i, j + 13, r14);
          if(r15 < Hradius) KernelAndGradient_bidirectional(x_diff15, y_diff15, i, j + 14, r15);
          if(r16 < Hradius) KernelAndGradient_bidirectional(x_diff16, y_diff16, i, j + 15, r16);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi = x_positions[i];
      yi = y_positions[i];

      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = x_positions[j];
        xj2 = x_positions[j + 1];
        xj3 = x_positions[j + 2];
        xj4 = x_positions[j + 3];
        yj1 = y_positions[j];
        yj2 = y_positions[j + 1];
        yj3 = y_positions[j + 2];
        yj4 = y_positions[j + 3];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        x_diff2 = xi - xj2;
        y_diff2 = yi - yj2;
        x_diff3 = xi - xj3;
        y_diff3 = yi - yj3;
        x_diff4 = xi - xj4;
        y_diff4 = yi - yj4;
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
        r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
        r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
      }

      // particles left that do not fit into unrolling factor 2
      for (; j < N_interior; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        if(r1 < Hradius)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
      }
    }

  }
  // rest part that does not fit into block_i
  unrolling_limit = N_interior - unrolling_factor2;
  for(i = block_end_i ; i < NUMBER_OF_PARTICLE ; i++){
    // each particle is neighbor to itself
    KernelAndGradient_zero(i);

    xi = x_positions[i];
    yi = y_positions[i];

    // unrolling with factor 2
    for (j = 0 ; j <= unrolling_limit ; j += unrolling_factor2){
      xj1 = x_positions[j];
      xj2 = x_positions[j + 1];
      xj3 = x_positions[j + 2];
      xj4 = x_positions[j + 3];
      yj1 = y_positions[j];
      yj2 = y_positions[j + 1];
      yj3 = y_positions[j + 2];
      yj4 = y_positions[j + 3];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      x_diff2 = xi - xj2;
      y_diff2 = yi - yj2;
      x_diff3 = xi - xj3;
      y_diff3 = yi - yj3;
      x_diff4 = xi - xj4;
      y_diff4 = yi - yj4;
      r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
      r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
      r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
      r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
      r1 = sqrt(r1);
      r2 = sqrt(r2);
      r3 = sqrt(r3);
      r4 = sqrt(r4);
      if(r1 < Hradius) KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j    , r1);
      if(r2 < Hradius) KernelAndGradient_bidirectional(x_diff2, y_diff2, i, j + 1, r2);
      if(r3 < Hradius) KernelAndGradient_bidirectional(x_diff3, y_diff3, i, j + 2, r3);
      if(r4 < Hradius) KernelAndGradient_bidirectional(x_diff4, y_diff4, i, j + 3, r4);
    }

    // particles left that do not fit into unrolling factor 2
    for (; j < N_interior; j++){
      xj1 = x_positions[j];
      yj1 = y_positions[j];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
      if(r1 < Hradius)
        KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
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
        xi = x_positions[i];
        yi = y_positions[i];

        // unrolling
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1  = x_positions[j];
          xj2  = x_positions[j + 1];
          xj3  = x_positions[j + 2];
          xj4  = x_positions[j + 3];
          xj5  = x_positions[j + 4];
          xj6  = x_positions[j + 5];
          xj7  = x_positions[j + 6];
          xj8  = x_positions[j + 7];
          xj9  = x_positions[j + 8];
          xj10 = x_positions[j + 9];
          xj11 = x_positions[j + 10];
          xj12 = x_positions[j + 11];
          xj13 = x_positions[j + 12];
          xj14 = x_positions[j + 13];
          xj15 = x_positions[j + 14];
          xj16 = x_positions[j + 15];
          yj1  = y_positions[j];
          yj2  = y_positions[j + 1];
          yj3  = y_positions[j + 2];
          yj4  = y_positions[j + 3];
          yj5  = y_positions[j + 4];
          yj6  = y_positions[j + 5];
          yj7  = y_positions[j + 6];
          yj8  = y_positions[j + 7];
          yj9  = y_positions[j + 8];
          yj10 = y_positions[j + 9];
          yj11 = y_positions[j + 10];
          yj12 = y_positions[j + 11];
          yj13 = y_positions[j + 12];
          yj14 = y_positions[j + 13];
          yj15 = y_positions[j + 14];
          yj16 = y_positions[j + 15];
          x_diff1 = xi - xj1;
          y_diff1 = yi - yj1;
          x_diff2 = xi - xj2;
          y_diff2 = yi - yj2;
          x_diff3 = xi - xj3;
          y_diff3 = yi - yj3;
          x_diff4 = xi - xj4;
          y_diff4 = yi - yj4;
          x_diff5 = xi - xj5;
          y_diff5 = yi - yj5;
          x_diff6 = xi - xj6;
          y_diff6 = yi - yj6;
          x_diff7 = xi - xj7;
          y_diff7 = yi - yj7;
          x_diff8 = xi - xj8;
          y_diff8 = yi - yj8;
          x_diff9 = xi - xj9;
          y_diff9 = yi - yj9;
          x_diff10 = xi - xj10;
          y_diff10 = yi - yj10;
          x_diff11 = xi - xj11;
          y_diff11 = yi - yj11;
          x_diff12 = xi - xj12;
          y_diff12 = yi - yj12;
          x_diff13 = xi - xj13;
          y_diff13 = yi - yj13;
          x_diff14 = xi - xj14;
          y_diff14 = yi - yj14;
          x_diff15 = xi - xj15;
          y_diff15 = yi - yj15;
          x_diff16 = xi - xj16;
          y_diff16 = yi - yj16;
          r1  = x_diff1  * x_diff1  + y_diff1  * y_diff1;
          r2  = x_diff2  * x_diff2  + y_diff2  * y_diff2;
          r3  = x_diff3  * x_diff3  + y_diff3  * y_diff3;
          r4  = x_diff4  * x_diff4  + y_diff4  * y_diff4;
          r5  = x_diff5  * x_diff5  + y_diff5  * y_diff5;
          r6  = x_diff6  * x_diff6  + y_diff6  * y_diff6;
          r7  = x_diff7  * x_diff7  + y_diff7  * y_diff7;
          r8  = x_diff8  * x_diff8  + y_diff8  * y_diff8;
          r9  = x_diff9  * x_diff9  + y_diff9  * y_diff9;
          r10 = x_diff10 * x_diff10 + y_diff10 * y_diff10;
          r11 = x_diff11 * x_diff11 + y_diff11 * y_diff11;
          r12 = x_diff12 * x_diff12 + y_diff12 * y_diff12;
          r13 = x_diff13 * x_diff13 + y_diff13 * y_diff13;
          r14 = x_diff14 * x_diff14 + y_diff14 * y_diff14;
          r15 = x_diff15 * x_diff15 + y_diff15 * y_diff15;
          r16 = x_diff16 * x_diff16 + y_diff16 * y_diff16;
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
          if(r1  < Hradius) KernelAndGradient_unidirectional(x_diff1, y_diff1, i, j    , r1);
          if(r2  < Hradius) KernelAndGradient_unidirectional(x_diff2, y_diff2, i, j + 1, r2);
          if(r3  < Hradius) KernelAndGradient_unidirectional(x_diff3, y_diff3, i, j + 2, r3);
          if(r4  < Hradius) KernelAndGradient_unidirectional(x_diff4, y_diff4, i, j + 3, r4);
          if(r5  < Hradius) KernelAndGradient_unidirectional(x_diff5, y_diff5, i, j + 4, r5);
          if(r6  < Hradius) KernelAndGradient_unidirectional(x_diff6, y_diff6, i, j + 5, r6);
          if(r7  < Hradius) KernelAndGradient_unidirectional(x_diff7, y_diff7, i, j + 6, r7);
          if(r8  < Hradius) KernelAndGradient_unidirectional(x_diff8, y_diff8, i, j + 7, r8);
          if(r9  < Hradius) KernelAndGradient_unidirectional(x_diff9, y_diff9, i, j + 8, r9);
          if(r10 < Hradius) KernelAndGradient_unidirectional(x_diff10, y_diff10, i, j + 9, r10);
          if(r11 < Hradius) KernelAndGradient_unidirectional(x_diff11, y_diff11, i, j + 10, r11);
          if(r12 < Hradius) KernelAndGradient_unidirectional(x_diff12, y_diff12, i, j + 11, r12);
          if(r13 < Hradius) KernelAndGradient_unidirectional(x_diff13, y_diff13, i, j + 12, r13);
          if(r14 < Hradius) KernelAndGradient_unidirectional(x_diff14, y_diff14, i, j + 13, r14);
          if(r15 < Hradius) KernelAndGradient_unidirectional(x_diff15, y_diff15, i, j + 14, r15);
          if(r16 < Hradius) KernelAndGradient_unidirectional(x_diff16, y_diff16, i, j + 15, r16);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior + N_repulsive - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi = x_positions[i];
      yi = y_positions[i];

      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj1 = x_positions[j];
        xj2 = x_positions[j + 1];
        xj3 = x_positions[j + 2];
        xj4 = x_positions[j + 3];
        yj1 = y_positions[j];
        yj2 = y_positions[j + 1];
        yj3 = y_positions[j + 2];
        yj4 = y_positions[j + 3];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        x_diff2 = xi - xj2;
        y_diff2 = yi - yj2;
        x_diff3 = xi - xj3;
        y_diff3 = yi - yj3;
        x_diff4 = xi - xj4;
        y_diff4 = yi - yj4;
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
        r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
        r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
        r1 = sqrt(r1);
        r2 = sqrt(r2);
        r3 = sqrt(r3);
        r4 = sqrt(r4);
        if(r1 < Hradius) KernelAndGradient_unidirectional(x_diff1, y_diff1, i, j    , r1);
        if(r2 < Hradius) KernelAndGradient_unidirectional(x_diff2, y_diff2, i, j + 1, r2);
        if(r3 < Hradius) KernelAndGradient_unidirectional(x_diff3, y_diff3, i, j + 2, r3);
        if(r4 < Hradius) KernelAndGradient_unidirectional(x_diff4, y_diff4, i, j + 3, r4);
      }

      // particles left that do not fit into unrolling factor 2
      for (; j < N_interior + N_repulsive ; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        if(r1 < Hradius)
          KernelAndGradient_unidirectional(x_diff1, y_diff1, i, j, r1);
      }
    
    } 
  }

  // rest part that does not fit into block_i
  unrolling_limit = N_interior + N_repulsive - unrolling_factor2;
  for(i = block_end_i ; i < NUMBER_OF_PARTICLE ; i++){
    xi = x_positions[i];
    yi = y_positions[i];

    // unrolling with factor 2
    for (j = N_interior ; j <= unrolling_limit ; j += unrolling_factor2){
      xj1 = x_positions[j];
      xj2 = x_positions[j + 1];
      xj3 = x_positions[j + 2];
      xj4 = x_positions[j + 3];
      yj1 = y_positions[j];
      yj2 = y_positions[j + 1];
      yj3 = y_positions[j + 2];
      yj4 = y_positions[j + 3];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      x_diff2 = xi - xj2;
      y_diff2 = yi - yj2;
      x_diff3 = xi - xj3;
      y_diff3 = yi - yj3;
      x_diff4 = xi - xj4;
      y_diff4 = yi - yj4;
      r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
      r2 = x_diff2 * x_diff2 + y_diff2 * y_diff2;
      r3 = x_diff3 * x_diff3 + y_diff3 * y_diff3;
      r4 = x_diff4 * x_diff4 + y_diff4 * y_diff4;
      r1 = sqrt(r1);
      r2 = sqrt(r2);
      r3 = sqrt(r3);
      r4 = sqrt(r4);
      if(r1 < Hradius) KernelAndGradient_unidirectional(x_diff1, y_diff1, i, j    , r1);
      if(r2 < Hradius) KernelAndGradient_unidirectional(x_diff2, y_diff2, i, j + 1, r2);
      if(r3 < Hradius) KernelAndGradient_unidirectional(x_diff3, y_diff3, i, j + 2, r3);
      if(r4 < Hradius) KernelAndGradient_unidirectional(x_diff4, y_diff4, i, j + 3, r4);
    }

    // particles left that do not fit into unrolling factor 2
    for (; j < N_interior + N_repulsive ; j++){
      xj1 = x_positions[j];
      yj1 = y_positions[j];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
      if(r1 < Hradius)
        KernelAndGradient_unidirectional(x_diff1, y_diff1, i, j, r1);
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
