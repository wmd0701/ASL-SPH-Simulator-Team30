//! @file output.h
#ifndef OUTPUT_H
#define OUTPUT_H

#include "data_set.h"
#include <stdlib.h>
#include <string.h>

/**
 * @brief Write data to a .dat file
 * @param all_particle 
 * @param t_now time now
 */
void WriteData(Particle* all_particle, double t_now) {
    Particle this_p;
    char time[23];
    t_now *= 1000000;
    sprintf(time, "data/data-%08.0f.csv", t_now);
    FILE *fp = NULL;
    fp = fopen(time,"w");
    fprintf(fp, "x coord, y coord, tag, u, v, m, rho, p, a1, a2\n"); 
    for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
        this_p = all_particle[i];
        fprintf(fp, "%lf,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",  
        this_p.position.first,
        this_p.position.second,
				this_p.tag,
        this_p.velocity.first,
        this_p.velocity.second,
        this_p.mass,
        this_p.density,
        this_p.pressure,
        this_p.accelerat.first,
        this_p.accelerat.second);
    }
    fclose(fp);
}


void RecordWaveHeight(Particle* all_particle, FILE *fp, double t_now) {
  int N = NUMBER_OF_PARTICLE;
  double wave_height = 0;
  for (int i = 0; i < N; i++) {
    if (all_particle[i].tag == interior && 
        all_particle[i].position.first < 0.5 && 
        all_particle[i].position.second > wave_height) {
      wave_height = all_particle[i].position.second;
    }
  }
  fprintf(fp, "%f, %f\n", t_now, wave_height);
  return;
}

#endif // OUTPUT_H
