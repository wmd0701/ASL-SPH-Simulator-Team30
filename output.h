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
    char time[22];
    t_now *= 1000;
    sprintf(time, "data/data-%08.0f.csv", t_now);
    FILE *fp = NULL;
    fp = fopen(time,"w");
    fprintf(fp, "x coord, y coord, u, v, m, rho, p\n"); 
    for (Index i = 0; i < NUMBER_OF_PARTICLE; i++) {
        this_p = all_particle[i];
        fprintf(fp, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n",  
        this_p.position.first,
        this_p.position.second,
        this_p.velocity.first,
        this_p.velocity.second,
        this_p.mass,
        this_p.density,
        this_p.pressure);
    }
    fclose(fp);
}

#endif // OUTPUT_H
