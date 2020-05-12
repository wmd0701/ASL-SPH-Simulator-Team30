//!  @file constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H


int N_interior; 
int Nx_interior;
int Ny_interior;

int N_repulsive;
int Nx_repulsive;
int Ny_repulsive;

int N_ghost;
int Nx_ghost;
int Ny_ghost;

int N_boundary;
int NUMBER_OF_PARTICLE;


double H;  //!< smoothing length
double Hinv; // 1 / H
double Hradius;
double factor;


const double dam_height = 0.6;
const double gravity = 9.81;
const double initial_density = 997.;
const double dynamic_viscosity = 0.8926e-3;
const double amplitude = 0.032;
const double period = 1.5;
double mass;

double cycles_DispBoundary   = 0;
double cycles_SearchNeighbor = 0;
double cycles_CompGlbDensity = 0;
double cycles_DensityCorr    = 0;
double cycles_CompPressure   = 0;
double cycles_CompAccelerat  = 0;
double cycles_RepulsiveForce = 0;
double cycles_TimeIntegral   = 0;
double cycles_all            = 0;

void set_particles_interior(int N) {
  H = (sqrt(2595. * N + 225.) - 15.) / (50. * N);
  Hinv = 1.0 / H;
  Hradius = 2.0 * H;
  factor = 10. / 7. / M_PI / H / H;
  Nx_interior = round(1.73 / H - 1);
  Ny_interior = round(0.6 / H);
  N_interior = Nx_interior * Ny_interior;

  Nx_repulsive = 2 * round(1.73 / H) + 1;
  Ny_repulsive = 2 * round(1.15 / H) + 1;
  N_repulsive = Nx_repulsive + 2 * Ny_repulsive;

  Nx_ghost = Nx_repulsive + 4;
  Ny_ghost = Ny_repulsive;
  N_ghost = 2 * Nx_ghost + 4 * Ny_ghost;

  N_boundary = N_repulsive + N_ghost;
  
  mass = H * H * initial_density;

  NUMBER_OF_PARTICLE = N_interior + N_boundary;
  printf("interior particles: %i (x: %i, y: %i)\nboundary particles: "
         "%i\nsmoothing length: %.5f\n\n",
         N_interior, Nx_interior, Ny_interior, N_boundary, H);
}

#endif // CONSTANTS_H
