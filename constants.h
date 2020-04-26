//!  @file constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H


int N_interior; 	// interior particles
int Nx_interior; 	// interior particles in x direction
int Ny_interior; 	// interior particles in y direction
double H;  //!< smoothing length
const int    NUMBER_OF_PARTICLE = 891;  //!< number of particles. Need to be change from cases to cases.
const double dam_height = 0.6;
const double gravity = 9.81;
const double initial_density = 997;
const double dynamic_viscosity = 0.8926e-3;
const double amplitude = 0.032;
const double period = 1.5;

void set_particles_interior(int N){
	H = (sqrt(2595 * N + 225) - 15) / (50 * N);
	Nx_interior = round(1.73 / H - 1);
	Ny_interior = round(0.6 / H);
	N_interior = Nx_interior * Ny_interior;	

	printf("interior particles: %i (x: %i, y: %i)\nsmoothing length: %.5f\n\n", 
					N_interior, Nx_interior, Ny_interior, H);
}

#endif // CONSTANTS_H
