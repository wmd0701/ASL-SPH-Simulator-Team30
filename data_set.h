//! @file data_set.h
#ifndef DATA_SET_H
#define DATA_SET_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include <immintrin.h>

// instead of using array of structs, use now multiple arrays
double* x_positions;
double* y_positions;
double* x_velocities;
double* y_velocities;
double* densities;
double* pressures;
double* x_accelerats;
double* y_accelerats;
double** Wijs;
double** x_Wij_grads;
double** y_Wij_grads;
int** neighbor_indices;
int* neighbor_counts;

// initial positions and velocities of boundary particles, used for DisplaceBoundaries in rate_of_change.h
double* x_init_positions;
double* y_init_positions;
double* x_init_velocities;
double* y_init_velocities;

void KernelAndGradient_zero(int par_idx) {
  double kernel = factor;
  
  int count = neighbor_counts[par_idx]++;
  
  Wijs            [par_idx][count] = kernel;
  x_Wij_grads     [par_idx][count] = 0.;
  y_Wij_grads     [par_idx][count] = 0.;
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
  r = sqrt(r);
  double kernel;
  double x_grad;
  double y_grad;
  double q = r * Hinv;
  double q_square = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q_square * (1. - 0.5 * q));
    temp = factor * (-3. * q + 2.25 * q_square) * Hinv / r;
    x_grad = temp * x_diff;
    y_grad = temp * y_diff;
  } 
  else {
    c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    x_grad = temp * x_diff;
    y_grad = temp * y_diff;
  }

  int count = neighbor_counts[par_idx_1]++;
  
  Wijs            [par_idx_1][count] = kernel;
  x_Wij_grads     [par_idx_1][count] = x_grad;
  y_Wij_grads     [par_idx_1][count] = y_grad;
  neighbor_indices[par_idx_1][count] = par_idx_2;

}

void KernelAndGradient_bidirectional(double x_diff, double y_diff, int par_idx_1, int par_idx_2, double r) {
  r = sqrt(r);
  double kernel;
  double x_grad;
  double y_grad;
  double q = r * Hinv;
  double q_square = q * q;
  double temp, c;

  if (q <= 1.) {
    kernel = factor * (1. - 1.5 * q_square * (1. - 0.5 * q));
    temp = factor * (-3. * q + 2.25 * q_square) * Hinv / r;
    x_grad = temp * x_diff;
    y_grad = temp * y_diff;
  } 
  else {
    c = (2.0 - q) * (2.0 - q);
    kernel = factor * 0.25 * c * (2.0 - q);
    temp = -factor * 0.75 * c * Hinv / r;
    x_grad = temp * x_diff;
    y_grad = temp * y_diff;
  }

  int count = neighbor_counts[par_idx_1]++;
  
  Wijs            [par_idx_1][count] = kernel;
  x_Wij_grads     [par_idx_1][count] = x_grad;
  y_Wij_grads     [par_idx_1][count] = y_grad;
  neighbor_indices[par_idx_1][count] = par_idx_2;

  x_grad = -x_grad;
  y_grad = -y_grad;
  count = neighbor_counts[par_idx_2]++;
  Wijs            [par_idx_2][count] = kernel;
  x_Wij_grads     [par_idx_2][count] = x_grad;
  y_Wij_grads     [par_idx_2][count] = y_grad;
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
  
  double xi, xj1, yi, yj1, x_diff1, y_diff1, r1;
  int unrolling_limit;          // unrolling limit
  int i, j;
  __m256d xi_vec, yi_vec, xj_vec, yj_vec, x_diff_vec, y_diff_vec, r_vec, temp_vec;
  __m256d xj1_vec, xj2_vec, xj3_vec, xj4_vec, yj1_vec, yj2_vec, yj3_vec, yj4_vec, x_diff1_vec, 
          x_diff2_vec, x_diff3_vec, x_diff4_vec, y_diff1_vec, y_diff2_vec, y_diff3_vec, y_diff4_vec,
          r1_vec, r2_vec, r3_vec, r4_vec, temp1_vec, temp2_vec, temp3_vec, temp4_vec;
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
      
      xi_vec = _mm256_broadcast_sd(x_positions + i);
      yi_vec = _mm256_broadcast_sd(y_positions + i);
      
      // unrolling with factor 2
      for (j = i + 1 ; j <= unrolling_limit ; j += unrolling_factor2){
        xj_vec = _mm256_load_pd(x_positions + j);
        yj_vec = _mm256_load_pd(y_positions + j);
        
        x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
        y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
        temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
        r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
        //~ r_vec = _mm256_sqrt_pd(r_vec);
        
        if(r_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
        if(r_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
        if(r_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
        if(r_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
        
      }

      // the particles left which do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
      for (; j < limit_j; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        if(r1 < Hradius2)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
      }

    }
    
    // search neighbor between two blocks
    for (block_j = block_i + block_size ; block_j < block_end_j ; block_j += block_size){
      limit_j = block_j + block_size;
      for (i = block_i ; i < limit_i ; i++){
        xi_vec = _mm256_broadcast_sd(x_positions + i);
        yi_vec = _mm256_broadcast_sd(y_positions + i);
        
        // unrolling with factor 1
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1_vec = _mm256_load_pd(x_positions + j);
          xj2_vec = _mm256_load_pd(x_positions + j + 4);
          xj3_vec = _mm256_load_pd(x_positions + j + 8);
          xj4_vec = _mm256_load_pd(x_positions + j + 12);
          yj1_vec = _mm256_load_pd(y_positions + j);
          yj2_vec = _mm256_load_pd(y_positions + j + 4);
          yj3_vec = _mm256_load_pd(y_positions + j + 8);
          yj4_vec = _mm256_load_pd(y_positions + j + 12);
          
          x_diff1_vec = _mm256_sub_pd(xi_vec, xj1_vec);
          y_diff1_vec = _mm256_sub_pd(yi_vec, yj1_vec);
          x_diff2_vec = _mm256_sub_pd(xi_vec, xj2_vec);
          y_diff2_vec = _mm256_sub_pd(yi_vec, yj2_vec);
          x_diff3_vec = _mm256_sub_pd(xi_vec, xj3_vec);
          y_diff3_vec = _mm256_sub_pd(yi_vec, yj3_vec);
          x_diff4_vec = _mm256_sub_pd(xi_vec, xj4_vec);
          y_diff4_vec = _mm256_sub_pd(yi_vec, yj4_vec);
          
          temp1_vec = _mm256_mul_pd(y_diff1_vec, y_diff1_vec);
          r1_vec = _mm256_fmadd_pd(x_diff1_vec, x_diff1_vec, temp1_vec);
          //~ r1_vec = _mm256_sqrt_pd(r1_vec);
          temp2_vec = _mm256_mul_pd(y_diff2_vec, y_diff2_vec);
          r2_vec = _mm256_fmadd_pd(x_diff2_vec, x_diff2_vec, temp2_vec);
          //~ r2_vec = _mm256_sqrt_pd(r2_vec);
          temp3_vec = _mm256_mul_pd(y_diff3_vec, y_diff3_vec);
          r3_vec = _mm256_fmadd_pd(x_diff3_vec, x_diff3_vec, temp3_vec);
          //~ r3_vec = _mm256_sqrt_pd(r3_vec);
          temp4_vec = _mm256_mul_pd(y_diff4_vec, y_diff4_vec);
          r4_vec = _mm256_fmadd_pd(x_diff4_vec, x_diff4_vec, temp4_vec);
          //~ r4_vec = _mm256_sqrt_pd(r4_vec);
          
          if(r1_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[0], y_diff1_vec[0], i, j    , r1_vec[0]);
          if(r1_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[1], y_diff1_vec[1], i, j + 1, r1_vec[1]);
          if(r1_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[2], y_diff1_vec[2], i, j + 2, r1_vec[2]);
          if(r1_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[3], y_diff1_vec[3], i, j + 3, r1_vec[3]);
          if(r2_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[0], y_diff2_vec[0], i, j + 4, r2_vec[0]);
          if(r2_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[1], y_diff2_vec[1], i, j + 5, r2_vec[1]);
          if(r2_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[2], y_diff2_vec[2], i, j + 6, r2_vec[2]);
          if(r2_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[3], y_diff2_vec[3], i, j + 7, r2_vec[3]);
          if(r3_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[0], y_diff3_vec[0], i, j + 8, r3_vec[0]);
          if(r3_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[1], y_diff3_vec[1], i, j + 9, r3_vec[1]);
          if(r3_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[2], y_diff3_vec[2], i, j + 10, r3_vec[2]);
          if(r3_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[3], y_diff3_vec[3], i, j + 11, r3_vec[3]);
          if(r4_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[0], y_diff4_vec[0], i, j + 12, r4_vec[0]);
          if(r4_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[1], y_diff4_vec[1], i, j + 13, r4_vec[1]);
          if(r4_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[2], y_diff4_vec[2], i, j + 14, r4_vec[2]);
          if(r4_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[3], y_diff4_vec[3], i, j + 15, r4_vec[3]);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi_vec = _mm256_broadcast_sd(x_positions + i);
      yi_vec = _mm256_broadcast_sd(y_positions + i);
      
      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj_vec = _mm256_load_pd(x_positions + j);
        yj_vec = _mm256_load_pd(y_positions + j);
        
        x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
        y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
        temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
        r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
        //~ r_vec = _mm256_sqrt_pd(r_vec);
        
        if(r_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
        if(r_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
        if(r_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
        if(r_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
      }

      // particles left that do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
      for (; j < N_interior; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        if(r1 < Hradius2)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
      }

    }

  }
  // rest part that does not fit into block_i
  unrolling_limit = N_interior - unrolling_factor2;
  for(i = block_end_i ; i < N_interior ; i++){
    // each particle is neighbor to itself
    KernelAndGradient_zero(i);

    xi_vec = _mm256_broadcast_sd(x_positions+i);
    yi_vec = _mm256_broadcast_sd(y_positions+i);

    // unrolling with factor 2
    for (j = i + 1 ; j <= unrolling_limit ; j += unrolling_factor2){
      xj_vec = _mm256_load_pd(x_positions + j);
      yj_vec = _mm256_load_pd(y_positions + j);
        
      x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
      y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
      temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
      r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
      //~ r_vec = _mm256_sqrt_pd(r_vec);
      
      if(r_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
      if(r_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
      if(r_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
      if(r_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
    }

    // particles left that do not fit into unrolling factor 2
    xi = x_positions[i];
    yi = y_positions[i];
    for (; j < N_interior; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        if(r1 < Hradius2)
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
        xi_vec = _mm256_broadcast_sd(x_positions + i);
        yi_vec = _mm256_broadcast_sd(y_positions + i);
        
        // unrolling with factor 1
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1_vec = _mm256_load_pd(x_positions + j);
          xj2_vec = _mm256_load_pd(x_positions + j + 4);
          xj3_vec = _mm256_load_pd(x_positions + j + 8);
          xj4_vec = _mm256_load_pd(x_positions + j + 12);
          yj1_vec = _mm256_load_pd(y_positions + j);
          yj2_vec = _mm256_load_pd(y_positions + j + 4);
          yj3_vec = _mm256_load_pd(y_positions + j + 8);
          yj4_vec = _mm256_load_pd(y_positions + j + 12);
          
          x_diff1_vec = _mm256_sub_pd(xi_vec, xj1_vec);
          y_diff1_vec = _mm256_sub_pd(yi_vec, yj1_vec);
          x_diff2_vec = _mm256_sub_pd(xi_vec, xj2_vec);
          y_diff2_vec = _mm256_sub_pd(yi_vec, yj2_vec);
          x_diff3_vec = _mm256_sub_pd(xi_vec, xj3_vec);
          y_diff3_vec = _mm256_sub_pd(yi_vec, yj3_vec);
          x_diff4_vec = _mm256_sub_pd(xi_vec, xj4_vec);
          y_diff4_vec = _mm256_sub_pd(yi_vec, yj4_vec);
          
          temp1_vec = _mm256_mul_pd(y_diff1_vec, y_diff1_vec);
          r1_vec = _mm256_fmadd_pd(x_diff1_vec, x_diff1_vec, temp1_vec);
          //~ r1_vec = _mm256_sqrt_pd(r1_vec);
          temp2_vec = _mm256_mul_pd(y_diff2_vec, y_diff2_vec);
          r2_vec = _mm256_fmadd_pd(x_diff2_vec, x_diff2_vec, temp2_vec);
          //~ r2_vec = _mm256_sqrt_pd(r2_vec);
          temp3_vec = _mm256_mul_pd(y_diff3_vec, y_diff3_vec);
          r3_vec = _mm256_fmadd_pd(x_diff3_vec, x_diff3_vec, temp3_vec);
          //~ r3_vec = _mm256_sqrt_pd(r3_vec);
          temp4_vec = _mm256_mul_pd(y_diff4_vec, y_diff4_vec);
          r4_vec = _mm256_fmadd_pd(x_diff4_vec, x_diff4_vec, temp4_vec);
          //~ r4_vec = _mm256_sqrt_pd(r4_vec);
          
          if(r1_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[0], y_diff1_vec[0], i, j    , r1_vec[0]);
          if(r1_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[1], y_diff1_vec[1], i, j + 1, r1_vec[1]);
          if(r1_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[2], y_diff1_vec[2], i, j + 2, r1_vec[2]);
          if(r1_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff1_vec[3], y_diff1_vec[3], i, j + 3, r1_vec[3]);
          if(r2_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[0], y_diff2_vec[0], i, j + 4, r2_vec[0]);
          if(r2_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[1], y_diff2_vec[1], i, j + 5, r2_vec[1]);
          if(r2_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[2], y_diff2_vec[2], i, j + 6, r2_vec[2]);
          if(r2_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff2_vec[3], y_diff2_vec[3], i, j + 7, r2_vec[3]);
          if(r3_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[0], y_diff3_vec[0], i, j + 8, r3_vec[0]);
          if(r3_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[1], y_diff3_vec[1], i, j + 9, r3_vec[1]);
          if(r3_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[2], y_diff3_vec[2], i, j + 10, r3_vec[2]);
          if(r3_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff3_vec[3], y_diff3_vec[3], i, j + 11, r3_vec[3]);
          if(r4_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[0], y_diff4_vec[0], i, j + 12, r4_vec[0]);
          if(r4_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[1], y_diff4_vec[1], i, j + 13, r4_vec[1]);
          if(r4_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[2], y_diff4_vec[2], i, j + 14, r4_vec[2]);
          if(r4_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff4_vec[3], y_diff4_vec[3], i, j + 15, r4_vec[3]);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi_vec = _mm256_broadcast_sd(x_positions + i);
      yi_vec = _mm256_broadcast_sd(y_positions + i);

      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj_vec = _mm256_load_pd(x_positions + j);
        yj_vec = _mm256_load_pd(y_positions + j);
        
        x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
        y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
        temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
        r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
        //~ r_vec = _mm256_sqrt_pd(r_vec);
        
        if(r_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
        if(r_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
        if(r_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
        if(r_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
      }

      // particles left that do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
      for (; j < N_interior; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        if(r1 < Hradius2)
          KernelAndGradient_bidirectional(x_diff1, y_diff1, i, j, r1);
      }
    }

  }
  // rest part that does not fit into block_i
  unrolling_limit = N_interior - unrolling_factor2;
  for(i = block_end_i ; i < NUMBER_OF_PARTICLE ; i++){
    // each particle is neighbor to itself
    KernelAndGradient_zero(i);

    xi_vec = _mm256_broadcast_sd(x_positions + i);
    yi_vec = _mm256_broadcast_sd(y_positions + i);

    // unrolling with factor 2
    for (j = 0 ; j <= unrolling_limit ; j += unrolling_factor2){
      xj_vec = _mm256_load_pd(x_positions + j);
      yj_vec = _mm256_load_pd(y_positions + j);
        
      x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
      y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
      temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
      r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
      //~ r_vec = _mm256_sqrt_pd(r_vec);
      
      if(r_vec[0] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
      if(r_vec[1] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
      if(r_vec[2] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
      if(r_vec[3] < Hradius2) KernelAndGradient_bidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
    }

    // particles left that do not fit into unrolling factor 2
    xi = x_positions[i];
    yi = y_positions[i];
    for (; j < N_interior; j++){
      xj1 = x_positions[j];
      yj1 = y_positions[j];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
      r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
      if(r1 < Hradius2)
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
        xi_vec = _mm256_broadcast_sd(x_positions + i);
        yi_vec = _mm256_broadcast_sd(y_positions + i);

        // unrolling
        for (j = block_j ; j != limit_j ; j += unrolling_factor1){
          xj1_vec = _mm256_load_pd(x_positions + j);
          xj2_vec = _mm256_load_pd(x_positions + j + 4);
          xj3_vec = _mm256_load_pd(x_positions + j + 8);
          xj4_vec = _mm256_load_pd(x_positions + j + 12);
          yj1_vec = _mm256_load_pd(y_positions + j);
          yj2_vec = _mm256_load_pd(y_positions + j + 4);
          yj3_vec = _mm256_load_pd(y_positions + j + 8);
          yj4_vec = _mm256_load_pd(y_positions + j + 12);
          
          x_diff1_vec = _mm256_sub_pd(xi_vec, xj1_vec);
          y_diff1_vec = _mm256_sub_pd(yi_vec, yj1_vec);
          x_diff2_vec = _mm256_sub_pd(xi_vec, xj2_vec);
          y_diff2_vec = _mm256_sub_pd(yi_vec, yj2_vec);
          x_diff3_vec = _mm256_sub_pd(xi_vec, xj3_vec);
          y_diff3_vec = _mm256_sub_pd(yi_vec, yj3_vec);
          x_diff4_vec = _mm256_sub_pd(xi_vec, xj4_vec);
          y_diff4_vec = _mm256_sub_pd(yi_vec, yj4_vec);
          
          temp1_vec = _mm256_mul_pd(y_diff1_vec, y_diff1_vec);
          r1_vec = _mm256_fmadd_pd(x_diff1_vec, x_diff1_vec, temp1_vec);
          //~ r1_vec = _mm256_sqrt_pd(r1_vec);
          temp2_vec = _mm256_mul_pd(y_diff2_vec, y_diff2_vec);
          r2_vec = _mm256_fmadd_pd(x_diff2_vec, x_diff2_vec, temp2_vec);
          //~ r2_vec = _mm256_sqrt_pd(r2_vec);
          temp3_vec = _mm256_mul_pd(y_diff3_vec, y_diff3_vec);
          r3_vec = _mm256_fmadd_pd(x_diff3_vec, x_diff3_vec, temp3_vec);
          //~ r3_vec = _mm256_sqrt_pd(r3_vec);
          temp4_vec = _mm256_mul_pd(y_diff4_vec, y_diff4_vec);
          r4_vec = _mm256_fmadd_pd(x_diff4_vec, x_diff4_vec, temp4_vec);
          //~ r4_vec = _mm256_sqrt_pd(r4_vec);
          
          if(r1_vec[0] < Hradius2) KernelAndGradient_unidirectional(x_diff1_vec[0], y_diff1_vec[0], i, j    , r1_vec[0]);
          if(r1_vec[1] < Hradius2) KernelAndGradient_unidirectional(x_diff1_vec[1], y_diff1_vec[1], i, j + 1, r1_vec[1]);
          if(r1_vec[2] < Hradius2) KernelAndGradient_unidirectional(x_diff1_vec[2], y_diff1_vec[2], i, j + 2, r1_vec[2]);
          if(r1_vec[3] < Hradius2) KernelAndGradient_unidirectional(x_diff1_vec[3], y_diff1_vec[3], i, j + 3, r1_vec[3]);
          if(r2_vec[0] < Hradius2) KernelAndGradient_unidirectional(x_diff2_vec[0], y_diff2_vec[0], i, j + 4, r2_vec[0]);
          if(r2_vec[1] < Hradius2) KernelAndGradient_unidirectional(x_diff2_vec[1], y_diff2_vec[1], i, j + 5, r2_vec[1]);
          if(r2_vec[2] < Hradius2) KernelAndGradient_unidirectional(x_diff2_vec[2], y_diff2_vec[2], i, j + 6, r2_vec[2]);
          if(r2_vec[3] < Hradius2) KernelAndGradient_unidirectional(x_diff2_vec[3], y_diff2_vec[3], i, j + 7, r2_vec[3]);
          if(r3_vec[0] < Hradius2) KernelAndGradient_unidirectional(x_diff3_vec[0], y_diff3_vec[0], i, j + 8, r3_vec[0]);
          if(r3_vec[1] < Hradius2) KernelAndGradient_unidirectional(x_diff3_vec[1], y_diff3_vec[1], i, j + 9, r3_vec[1]);
          if(r3_vec[2] < Hradius2) KernelAndGradient_unidirectional(x_diff3_vec[2], y_diff3_vec[2], i, j + 10, r3_vec[2]);
          if(r3_vec[3] < Hradius2) KernelAndGradient_unidirectional(x_diff3_vec[3], y_diff3_vec[3], i, j + 11, r3_vec[3]);
          if(r4_vec[0] < Hradius2) KernelAndGradient_unidirectional(x_diff4_vec[0], y_diff4_vec[0], i, j + 12, r4_vec[0]);
          if(r4_vec[1] < Hradius2) KernelAndGradient_unidirectional(x_diff4_vec[1], y_diff4_vec[1], i, j + 13, r4_vec[1]);
          if(r4_vec[2] < Hradius2) KernelAndGradient_unidirectional(x_diff4_vec[2], y_diff4_vec[2], i, j + 14, r4_vec[2]);
          if(r4_vec[3] < Hradius2) KernelAndGradient_unidirectional(x_diff4_vec[3], y_diff4_vec[3], i, j + 15, r4_vec[3]);
        }

      }
    }

    // rest part that does not fit into block_j
    unrolling_limit = N_interior + N_repulsive - unrolling_factor2;
    for (i = block_i ; i < limit_i ; i++){
      xi_vec = _mm256_broadcast_sd(x_positions + i);
      yi_vec = _mm256_broadcast_sd(y_positions + i);

      // unrolling with factor 2
      for (j = block_end_j ; j <= unrolling_limit ; j += unrolling_factor2){
        xj_vec = _mm256_load_pd(x_positions + j);
        yj_vec = _mm256_load_pd(y_positions + j);
        
        x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
        y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
        temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
        r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
        //~ r_vec = _mm256_sqrt_pd(r_vec);
        
        if(r_vec[0] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
        if(r_vec[1] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
        if(r_vec[2] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
        if(r_vec[3] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
      }

      // particles left that do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
      for (; j < N_interior + N_repulsive ; j++){
        xj1 = x_positions[j];
        yj1 = y_positions[j];
        x_diff1 = xi - xj1;
        y_diff1 = yi - yj1;
        //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
        r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
        if(r1 < Hradius2)
          KernelAndGradient_unidirectional(x_diff1, y_diff1, i, j, r1);
      }
    
    } 
  }

  // rest part that does not fit into block_i
  unrolling_limit = N_interior + N_repulsive - unrolling_factor2;
  for(i = block_end_i ; i < NUMBER_OF_PARTICLE ; i++){
    xi_vec = _mm256_broadcast_sd(x_positions + i);
    yi_vec = _mm256_broadcast_sd(y_positions + i);

    // unrolling with factor 2
    for (j = N_interior ; j <= unrolling_limit ; j += unrolling_factor2){
      xj_vec = _mm256_load_pd(x_positions + j);
      yj_vec = _mm256_load_pd(y_positions + j);
        
      x_diff_vec = _mm256_sub_pd(xi_vec, xj_vec);
      y_diff_vec = _mm256_sub_pd(yi_vec, yj_vec);
      temp_vec = _mm256_mul_pd(y_diff_vec, y_diff_vec);
      r_vec = _mm256_fmadd_pd(x_diff_vec, x_diff_vec, temp_vec);
      //~ r_vec = _mm256_sqrt_pd(r_vec);
      
      if(r_vec[0] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[0], y_diff_vec[0], i, j    , r_vec[0]);
      if(r_vec[1] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[1], y_diff_vec[1], i, j + 1, r_vec[1]);
      if(r_vec[2] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[2], y_diff_vec[2], i, j + 2, r_vec[2]);
      if(r_vec[3] < Hradius2) KernelAndGradient_unidirectional(x_diff_vec[3], y_diff_vec[3], i, j + 3, r_vec[3]);
    }

    // particles left that do not fit into unrolling factor 2
    xi = x_positions[i];
    yi = y_positions[i];
    for (; j < N_interior + N_repulsive ; j++){
      xj1 = x_positions[j];
      yj1 = y_positions[j];
      x_diff1 = xi - xj1;
      y_diff1 = yi - yj1;
      //~ r1 = sqrt(x_diff1 * x_diff1 + y_diff1 * y_diff1);
      r1 = x_diff1 * x_diff1 + y_diff1 * y_diff1;
      if(r1 < Hradius2)
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
  x_positions  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  y_positions  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  x_velocities = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  y_velocities = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  densities  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  pressures  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  x_accelerats = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  y_accelerats = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));

  Wijs             = (double**)calloc(NUMBER_OF_PARTICLE, sizeof(double*));
  x_Wij_grads      = (double**)calloc(NUMBER_OF_PARTICLE, sizeof(double*));
  y_Wij_grads      = (double**)calloc(NUMBER_OF_PARTICLE, sizeof(double*));
  neighbor_indices = (int**   )calloc(NUMBER_OF_PARTICLE, sizeof(int*));
  for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
    Wijs[i]             = (double*)calloc(max_num_neighbors, sizeof(double));
    x_Wij_grads[i]      = (double*)calloc(max_num_neighbors, sizeof(double));
    y_Wij_grads[i]      = (double*)calloc(max_num_neighbors, sizeof(double));
    neighbor_indices[i] = (int*   )calloc(max_num_neighbors, sizeof(int)); 
  }
  neighbor_counts = (int*)calloc(NUMBER_OF_PARTICLE, sizeof(int));

  x_init_positions  = (double*)calloc(N_boundary, sizeof(double));
  y_init_positions  = (double*)calloc(N_boundary, sizeof(double));
  x_init_velocities = (double*)calloc(N_boundary, sizeof(double));
  y_init_velocities = (double*)calloc(N_boundary, sizeof(double));

  int now = 0;

  // Set interior particles
  for (int i = 0; i < Nx_interior; ++i){
    for (int j = 0; j < Ny_interior; ++j){
      x_positions[now] = (i + 1) * H;
      y_positions[now++] = (j + 1) * H;
    }
  }

  // Set repulsive particles
  for (int i = 0; i < Nx_repulsive; i++){
    x_positions[now] = i * H / 2.;
    y_positions[now++] = 0;
  }

  for (int j = 1; j < Ny_repulsive; j++){
    x_positions[now] = 0;
    y_positions[now++] = j * H / 2.;
  }

  for (int j = 0; j < Ny_repulsive; j++){
    x_positions[now] = (Nx_interior + 1) * H;
    y_positions[now++] = j * H / 2.;
  }

  // Set ghost particles
  for (int i = -2; i < Nx_repulsive + 2; i++) {
    x_positions[now] = i * H / 2.;
    y_positions[now++] = -H / 2.;
    x_positions[now] = i * H / 2.;
    y_positions[now++] = -H;
  }

  for (int j = 0; j < Ny_repulsive; j++) {
    x_positions[now] = -H;
    y_positions[now++] = j * H / 2.;
    x_positions[now] = -H / 2.;
    y_positions[now++] = j * H / 2.;
    x_positions[now] = (Nx_interior + 1.5) * H;
    y_positions[now++] = j * H / 2.;
    x_positions[now] = (Nx_interior + 2.) * H;
    y_positions[now++] = j * H / 2.;
  }

  if (NUMBER_OF_PARTICLE != now)
    printf("number of particles doesn't match with init,\n");

  for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
    x_positions[i] += amplitude;
    densities[i] = initial_density;
    pressures[i] = 1.;
  }

  // copy initial values
  for(int i = 0 ; i < N_boundary ; i++){
    x_init_positions[i] = x_positions[i + N_interior];
    y_init_positions[i] = y_positions[i + N_interior];
    x_init_velocities[i] = x_velocities[i + N_interior];
    y_init_velocities[i] = y_velocities[i + N_interior];
  }

}

void Destroy(){
  free(x_positions);
  free(y_positions);
  free(x_velocities);
  free(y_velocities);
  free(densities);
  free(pressures);
  free(x_accelerats);
  free(y_accelerats);

  for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
    free(Wijs[i]);
    free(x_Wij_grads[i]);
    free(y_Wij_grads[i]);
    free(neighbor_indices[i]);
  }
  free(Wijs);
  free(x_Wij_grads);
  free(y_Wij_grads);
  free(neighbor_indices);
  free(neighbor_counts);

  free(x_init_positions);
  free(y_init_positions);
  free(x_init_velocities);
  free(y_init_velocities);
}

#endif // DATA_SET_H
