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
  double r_bidirectional[9*N_interior];
  double x_diff_bidirectional[9*N_interior];
  double y_diff_bidirectional[9*N_interior];
  int count_bidirectional;
  double r_unidirectional[9*N_ghost];
  double x_diff_unidirectional[9*N_ghost];
  double y_diff_unidirectional[9*N_ghost];
  int count_unidirectional;
   
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
  __m256d xi_vec, yi_vec, xj_vec, yj_vec, x_diff_vec, y_diff_vec, r_vec, temp_vec;
  __m256d xj1_vec, xj2_vec, xj3_vec, xj4_vec, yj1_vec, yj2_vec, yj3_vec, yj4_vec, x_diff1_vec, 
          x_diff2_vec, x_diff3_vec, x_diff4_vec, y_diff1_vec, y_diff2_vec, y_diff3_vec, y_diff4_vec,
          r1_vec, r2_vec, r3_vec, r4_vec, temp1_vec, temp2_vec, temp3_vec, temp4_vec;
  double r[4], x_diff[4], y_diff[4];
  double radii[16], x_diffe[16], y_diffe[16]; 
  __m256d Hradius_vec = _mm256_broadcast_pd(&Hradius);
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
        r_vec = _mm256_sqrt_pd(r_vec);
        _mm256_store_pd(r, r_vec);
        _mm256_store_pd(x_diff, x_diff_vec);
        _mm256_store_pd(y_diff, y_diff_vec);
        if(r[0] < Hradius) KernelAndGradient_bidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
        if(r[1] < Hradius) KernelAndGradient_bidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
        if(r[2] < Hradius) KernelAndGradient_bidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
        if(r[3] < Hradius) KernelAndGradient_bidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
        
      }

      // the particles left which do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
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
          r1_vec = _mm256_sqrt_pd(r1_vec);
          temp2_vec = _mm256_mul_pd(y_diff2_vec, y_diff2_vec);
          r2_vec = _mm256_fmadd_pd(x_diff2_vec, x_diff2_vec, temp2_vec);
          r2_vec = _mm256_sqrt_pd(r2_vec);
          temp3_vec = _mm256_mul_pd(y_diff3_vec, y_diff3_vec);
          r3_vec = _mm256_fmadd_pd(x_diff3_vec, x_diff3_vec, temp3_vec);
          r3_vec = _mm256_sqrt_pd(r3_vec);
          temp4_vec = _mm256_mul_pd(y_diff4_vec, y_diff4_vec);
          r4_vec = _mm256_fmadd_pd(x_diff4_vec, x_diff4_vec, temp4_vec);
          r4_vec = _mm256_sqrt_pd(r4_vec);
          
          _mm256_store_pd(radii, r1_vec);
          _mm256_store_pd(x_diffe, x_diff1_vec);
          _mm256_store_pd(y_diffe, y_diff1_vec);
          _mm256_store_pd(radii + 4, r2_vec);
          _mm256_store_pd(x_diffe + 4, x_diff2_vec);
          _mm256_store_pd(y_diffe + 4, y_diff2_vec);
          _mm256_store_pd(radii + 8, r3_vec);
          _mm256_store_pd(x_diffe + 8, x_diff3_vec);
          _mm256_store_pd(y_diffe + 8, y_diff3_vec);
          _mm256_store_pd(radii + 12, r4_vec);
          _mm256_store_pd(x_diffe + 12, x_diff4_vec);
          _mm256_store_pd(y_diffe + 12, y_diff4_vec);
          
          
          if(radii[0]  < Hradius) KernelAndGradient_bidirectional(x_diffe[0], y_diffe[0], i, j    , radii[0]);
          if(radii[1]  < Hradius) KernelAndGradient_bidirectional(x_diffe[1], y_diffe[1], i, j + 1, radii[1]);
          if(radii[2]  < Hradius) KernelAndGradient_bidirectional(x_diffe[2], y_diffe[2], i, j + 2, radii[2]);
          if(radii[3]  < Hradius) KernelAndGradient_bidirectional(x_diffe[3], y_diffe[3], i, j + 3, radii[3]);
          if(radii[4]  < Hradius) KernelAndGradient_bidirectional(x_diffe[4], y_diffe[4], i, j + 4, radii[4]);
          if(radii[5]  < Hradius) KernelAndGradient_bidirectional(x_diffe[5], y_diffe[5], i, j + 5, radii[5]);
          if(radii[6]  < Hradius) KernelAndGradient_bidirectional(x_diffe[6], y_diffe[6], i, j + 6, radii[6]);
          if(radii[7]  < Hradius) KernelAndGradient_bidirectional(x_diffe[7], y_diffe[7], i, j + 7, radii[7]);
          if(radii[8]  < Hradius) KernelAndGradient_bidirectional(x_diffe[8], y_diffe[8], i, j + 8, radii[8]);
          if(radii[9]  < Hradius) KernelAndGradient_bidirectional(x_diffe[9], y_diffe[9], i, j + 9, radii[9]);
          if(radii[10] < Hradius) KernelAndGradient_bidirectional(x_diffe[10], y_diffe[10], i, j + 10, radii[10]);
          if(radii[11] < Hradius) KernelAndGradient_bidirectional(x_diffe[11], y_diffe[11], i, j + 11, radii[11]);
          if(radii[12] < Hradius) KernelAndGradient_bidirectional(x_diffe[12], y_diffe[12], i, j + 12, radii[12]);
          if(radii[13] < Hradius) KernelAndGradient_bidirectional(x_diffe[13], y_diffe[13], i, j + 13, radii[13]);
          if(radii[14] < Hradius) KernelAndGradient_bidirectional(x_diffe[14], y_diffe[14], i, j + 14, radii[14]);
          if(radii[15] < Hradius) KernelAndGradient_bidirectional(x_diffe[15], y_diffe[15], i, j + 15, radii[15]);
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
        r_vec = _mm256_sqrt_pd(r_vec);
        _mm256_store_pd(r, r_vec);
        _mm256_store_pd(x_diff, x_diff_vec);
        _mm256_store_pd(y_diff, y_diff_vec);
        if(r[0] < Hradius) KernelAndGradient_bidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
        if(r[1] < Hradius) KernelAndGradient_bidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
        if(r[2] < Hradius) KernelAndGradient_bidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
        if(r[3] < Hradius) KernelAndGradient_bidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
      }

      // particles left that do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
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
      r_vec = _mm256_sqrt_pd(r_vec);
      _mm256_store_pd(r, r_vec);
      _mm256_store_pd(x_diff, x_diff_vec);
      _mm256_store_pd(y_diff, y_diff_vec);
      if(r[0] < Hradius) KernelAndGradient_bidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
      if(r[1] < Hradius) KernelAndGradient_bidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
      if(r[2] < Hradius) KernelAndGradient_bidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
      if(r[3] < Hradius) KernelAndGradient_bidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
    }

    // particles left that do not fit into unrolling factor 2
    xi = x_positions[i];
    yi = y_positions[i];
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
          r1_vec = _mm256_sqrt_pd(r1_vec);
          temp2_vec = _mm256_mul_pd(y_diff2_vec, y_diff2_vec);
          r2_vec = _mm256_fmadd_pd(x_diff2_vec, x_diff2_vec, temp2_vec);
          r2_vec = _mm256_sqrt_pd(r2_vec);
          temp3_vec = _mm256_mul_pd(y_diff3_vec, y_diff3_vec);
          r3_vec = _mm256_fmadd_pd(x_diff3_vec, x_diff3_vec, temp3_vec);
          r3_vec = _mm256_sqrt_pd(r3_vec);
          temp4_vec = _mm256_mul_pd(y_diff4_vec, y_diff4_vec);
          r4_vec = _mm256_fmadd_pd(x_diff4_vec, x_diff4_vec, temp4_vec);
          r4_vec = _mm256_sqrt_pd(r4_vec);
          
          _mm256_store_pd(radii, r1_vec);
          _mm256_store_pd(x_diffe, x_diff1_vec);
          _mm256_store_pd(y_diffe, y_diff1_vec);
          _mm256_store_pd(radii + 4, r2_vec);
          _mm256_store_pd(x_diffe + 4, x_diff2_vec);
          _mm256_store_pd(y_diffe + 4, y_diff2_vec);
          _mm256_store_pd(radii + 8, r3_vec);
          _mm256_store_pd(x_diffe + 8, x_diff3_vec);
          _mm256_store_pd(y_diffe + 8, y_diff3_vec);
          _mm256_store_pd(radii + 12, r4_vec);
          _mm256_store_pd(x_diffe + 12, x_diff4_vec);
          _mm256_store_pd(y_diffe + 12, y_diff4_vec);
          
          
          if(radii[0]  < Hradius) KernelAndGradient_bidirectional(x_diffe[0], y_diffe[0], i, j    , radii[0]);
          if(radii[1]  < Hradius) KernelAndGradient_bidirectional(x_diffe[1], y_diffe[1], i, j + 1, radii[1]);
          if(radii[2]  < Hradius) KernelAndGradient_bidirectional(x_diffe[2], y_diffe[2], i, j + 2, radii[2]);
          if(radii[3]  < Hradius) KernelAndGradient_bidirectional(x_diffe[3], y_diffe[3], i, j + 3, radii[3]);
          if(radii[4]  < Hradius) KernelAndGradient_bidirectional(x_diffe[4], y_diffe[4], i, j + 4, radii[4]);
          if(radii[5]  < Hradius) KernelAndGradient_bidirectional(x_diffe[5], y_diffe[5], i, j + 5, radii[5]);
          if(radii[6]  < Hradius) KernelAndGradient_bidirectional(x_diffe[6], y_diffe[6], i, j + 6, radii[6]);
          if(radii[7]  < Hradius) KernelAndGradient_bidirectional(x_diffe[7], y_diffe[7], i, j + 7, radii[7]);
          if(radii[8]  < Hradius) KernelAndGradient_bidirectional(x_diffe[8], y_diffe[8], i, j + 8, radii[8]);
          if(radii[9]  < Hradius) KernelAndGradient_bidirectional(x_diffe[9], y_diffe[9], i, j + 9, radii[9]);
          if(radii[10] < Hradius) KernelAndGradient_bidirectional(x_diffe[10], y_diffe[10], i, j + 10, radii[10]);
          if(radii[11] < Hradius) KernelAndGradient_bidirectional(x_diffe[11], y_diffe[11], i, j + 11, radii[11]);
          if(radii[12] < Hradius) KernelAndGradient_bidirectional(x_diffe[12], y_diffe[12], i, j + 12, radii[12]);
          if(radii[13] < Hradius) KernelAndGradient_bidirectional(x_diffe[13], y_diffe[13], i, j + 13, radii[13]);
          if(radii[14] < Hradius) KernelAndGradient_bidirectional(x_diffe[14], y_diffe[14], i, j + 14, radii[14]);
          if(radii[15] < Hradius) KernelAndGradient_bidirectional(x_diffe[15], y_diffe[15], i, j + 15, radii[15]);
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
        r_vec = _mm256_sqrt_pd(r_vec);
        _mm256_store_pd(r, r_vec);
        _mm256_store_pd(x_diff, x_diff_vec);
        _mm256_store_pd(y_diff, y_diff_vec);
        if(r[0] < Hradius) KernelAndGradient_bidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
        if(r[1] < Hradius) KernelAndGradient_bidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
        if(r[2] < Hradius) KernelAndGradient_bidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
        if(r[3] < Hradius) KernelAndGradient_bidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
      }

      // particles left that do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
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
      r_vec = _mm256_sqrt_pd(r_vec);
      _mm256_store_pd(r, r_vec);
      _mm256_store_pd(x_diff, x_diff_vec);
      _mm256_store_pd(y_diff, y_diff_vec);
      if(r[0] < Hradius) KernelAndGradient_bidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
      if(r[1] < Hradius) KernelAndGradient_bidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
      if(r[2] < Hradius) KernelAndGradient_bidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
      if(r[3] < Hradius) KernelAndGradient_bidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
    }

    // particles left that do not fit into unrolling factor 2
    xi = x_positions[i];
    yi = y_positions[i];
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
          r1_vec = _mm256_sqrt_pd(r1_vec);
          temp2_vec = _mm256_mul_pd(y_diff2_vec, y_diff2_vec);
          r2_vec = _mm256_fmadd_pd(x_diff2_vec, x_diff2_vec, temp2_vec);
          r2_vec = _mm256_sqrt_pd(r2_vec);
          temp3_vec = _mm256_mul_pd(y_diff3_vec, y_diff3_vec);
          r3_vec = _mm256_fmadd_pd(x_diff3_vec, x_diff3_vec, temp3_vec);
          r3_vec = _mm256_sqrt_pd(r3_vec);
          temp4_vec = _mm256_mul_pd(y_diff4_vec, y_diff4_vec);
          r4_vec = _mm256_fmadd_pd(x_diff4_vec, x_diff4_vec, temp4_vec);
          r4_vec = _mm256_sqrt_pd(r4_vec);
          
          _mm256_store_pd(radii, r1_vec);
          _mm256_store_pd(x_diffe, x_diff1_vec);
          _mm256_store_pd(y_diffe, y_diff1_vec);
          _mm256_store_pd(radii + 4, r2_vec);
          _mm256_store_pd(x_diffe + 4, x_diff2_vec);
          _mm256_store_pd(y_diffe + 4, y_diff2_vec);
          _mm256_store_pd(radii + 8, r3_vec);
          _mm256_store_pd(x_diffe + 8, x_diff3_vec);
          _mm256_store_pd(y_diffe + 8, y_diff3_vec);
          _mm256_store_pd(radii + 12, r4_vec);
          _mm256_store_pd(x_diffe + 12, x_diff4_vec);
          _mm256_store_pd(y_diffe + 12, y_diff4_vec);
          
          
          if(radii[0]  < Hradius) KernelAndGradient_unidirectional(x_diffe[0], y_diffe[0], i, j    , radii[0]);
          if(radii[1]  < Hradius) KernelAndGradient_unidirectional(x_diffe[1], y_diffe[1], i, j + 1, radii[1]);
          if(radii[2]  < Hradius) KernelAndGradient_unidirectional(x_diffe[2], y_diffe[2], i, j + 2, radii[2]);
          if(radii[3]  < Hradius) KernelAndGradient_unidirectional(x_diffe[3], y_diffe[3], i, j + 3, radii[3]);
          if(radii[4]  < Hradius) KernelAndGradient_unidirectional(x_diffe[4], y_diffe[4], i, j + 4, radii[4]);
          if(radii[5]  < Hradius) KernelAndGradient_unidirectional(x_diffe[5], y_diffe[5], i, j + 5, radii[5]);
          if(radii[6]  < Hradius) KernelAndGradient_unidirectional(x_diffe[6], y_diffe[6], i, j + 6, radii[6]);
          if(radii[7]  < Hradius) KernelAndGradient_unidirectional(x_diffe[7], y_diffe[7], i, j + 7, radii[7]);
          if(radii[8]  < Hradius) KernelAndGradient_unidirectional(x_diffe[8], y_diffe[8], i, j + 8, radii[8]);
          if(radii[9]  < Hradius) KernelAndGradient_unidirectional(x_diffe[9], y_diffe[9], i, j + 9, radii[9]);
          if(radii[10] < Hradius) KernelAndGradient_unidirectional(x_diffe[10], y_diffe[10], i, j + 10, radii[10]);
          if(radii[11] < Hradius) KernelAndGradient_unidirectional(x_diffe[11], y_diffe[11], i, j + 11, radii[11]);
          if(radii[12] < Hradius) KernelAndGradient_unidirectional(x_diffe[12], y_diffe[12], i, j + 12, radii[12]);
          if(radii[13] < Hradius) KernelAndGradient_unidirectional(x_diffe[13], y_diffe[13], i, j + 13, radii[13]);
          if(radii[14] < Hradius) KernelAndGradient_unidirectional(x_diffe[14], y_diffe[14], i, j + 14, radii[14]);
          if(radii[15] < Hradius) KernelAndGradient_unidirectional(x_diffe[15], y_diffe[15], i, j + 15, radii[15]);
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
        r_vec = _mm256_sqrt_pd(r_vec);
        _mm256_store_pd(r, r_vec);
        _mm256_store_pd(x_diff, x_diff_vec);
        _mm256_store_pd(y_diff, y_diff_vec);
        if(r[0] < Hradius) KernelAndGradient_unidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
        if(r[1] < Hradius) KernelAndGradient_unidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
        if(r[2] < Hradius) KernelAndGradient_unidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
        if(r[3] < Hradius) KernelAndGradient_unidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
      }

      // particles left that do not fit into unrolling factor 2
      xi = x_positions[i];
      yi = y_positions[i];
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
      r_vec = _mm256_sqrt_pd(r_vec);
      _mm256_store_pd(r, r_vec);
      _mm256_store_pd(x_diff, x_diff_vec);
      _mm256_store_pd(y_diff, y_diff_vec);
      if(r[0] < Hradius) KernelAndGradient_unidirectional(x_diff[0], y_diff[0], i, j    , r[0]);
      if(r[1] < Hradius) KernelAndGradient_unidirectional(x_diff[1], y_diff[1], i, j + 1, r[1]);
      if(r[2] < Hradius) KernelAndGradient_unidirectional(x_diff[2], y_diff[2], i, j + 2, r[2]);
      if(r[3] < Hradius) KernelAndGradient_unidirectional(x_diff[3], y_diff[3], i, j + 3, r[3]);
    }

    // particles left that do not fit into unrolling factor 2
    xi = x_positions[i];
    yi = y_positions[i];
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
  //~ x_positions  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  //~ y_positions  = (double*)calloc(NUMBER_OF_PARTICLE, sizeof(double));
  posix_memalign((void**)&x_positions, 32, NUMBER_OF_PARTICLE*sizeof(double));
  posix_memalign((void**)&y_positions, 32, NUMBER_OF_PARTICLE*sizeof(double));
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
