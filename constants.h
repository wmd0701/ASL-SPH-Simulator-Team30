//!  @file constants.h
#ifndef CONSTANTS_H
#define CONSTANTS_H


const double H = 0.05;  //!< smoothing length
const int    NUMBER_OF_PARTICLE = 891;  //!< number of particles. Need to be change from cases to cases.
const double dam_height = 0.6;
const double gravity = 9.81;
const double initial_density = 997;
const double dynamic_viscosity = 0.8926e-3;
const double amplitude = 0.032;
const double period = 1.5;

void set_smoothing_length(double smoothing_length){
	H = smoothing_length;
}

#endif // CONSTANTS_H
