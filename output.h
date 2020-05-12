//! @file output.h
#ifndef OUTPUT_H
#define OUTPUT_H

#include "data_set.h"
#include "constants.h"
#include <string.h>
char folder_name[40];

void set_output_path(char* output_path){
	strcpy(folder_name, output_path);
}

int get_tag(int i){
	if(i < N_interior)
		return 0;
	else if (i < N_interior + N_repulsive)
		return 1;
	else
		return 2;
}

/**
 * @brief Write data to a .dat file
 * @param all_particle 
 * @param t_now time now
 */
void WriteData(double t_now) {
	char file_name[23];
	char output_path[60];
	t_now *= 1000000;
	sprintf(file_name, "/data-%08.0f.csv", t_now);
	strcpy(output_path, folder_name);
	strcat(output_path, file_name);
	FILE *fp = NULL;
	fp = fopen(output_path,"w");
	fprintf(fp, "x coord, y coord, tag, u, v, rho, p, a1, a2\n"); 
	for (int i = 0; i < NUMBER_OF_PARTICLE; i++) {
		fprintf(fp, "%lf, %lf, %d, %lf, %lf, %lf, %lf, %lf, %lf\n",  
				positions[i].first,
				positions[i].second,
				get_tag(i),
				velocities[i].first,
				velocities[i].second,
				densities[i],
				pressures[i],
				accelerats[i].first,
				accelerats[i].second);
	}
	fclose(fp);
}

void WritePerformance() {
	char file_name[23];
	sprintf(file_name, "/performance.csv");
	char output_path[60];
	strcpy(output_path, folder_name);
	strcat(output_path, file_name);
	
	FILE *fp = NULL;
	fp = fopen(output_path,"w");
	fprintf(fp, "interior_p, boundary_p, smoothing_length, cycles_DispBoundar, cycles_SearchNeighbor, cycles_CompGlbDensity, cycles_DensityCorr, cycles_CompPressure, cycles_CompAccelerat, cycles_RepulsiveForce, cycles_TimeIntegral, cycles_all\n"); 
	fprintf(fp, "%d, %d, %.5lf, %.0lf, %.0lf, %.0lf %.0lf, %.0lf, %.0lf, %.0lf, %.0lf, %.0lf\n",
			N_interior,
			N_boundary,
			H,
			cycles_DispBoundary,
			cycles_SearchNeighbor,
			cycles_CompGlbDensity,
			cycles_DensityCorr,
			cycles_CompPressure,
			cycles_CompAccelerat,
			cycles_RepulsiveForce,
			cycles_TimeIntegral,
			cycles_all	);
	fclose(fp);
}


void RecordWaveHeight(FILE *fp, double t_now) {
	double wave_height = 0;
	for (int i = 0; i < N_interior; i++) {
		if (positions[i].first < 0.1 && positions[i].second > wave_height) {
			wave_height = positions[i].second;
		}
	}
	fprintf(fp, "%f, %f\n", t_now, wave_height);
	return;
}

#endif // OUTPUT_H
