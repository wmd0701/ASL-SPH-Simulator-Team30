//!     @file time_integration.h
#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H


#include "data_set.h"

enum Integration_Method {EXPLICIT_EULER, HEUN, MIDPOINT, LEAP_FROG};

enum Integration_Method integration;


/**     
*		@brief different time integration methods
*       @param all_particle: pointer to an array of particles
*       @param dt: time step
*		
*			   One of the following four methods will be called in Time_Integration
*/
void Explicit_Euler (Particle* all_particle, double dt);
void Heun           (Particle* all_particle, double dt);
void Midpoint       (Particle* all_particle, double dt);
void Leap_Frog      (Particle* all_particle, double dt);


/**     
*		@brief Choose the integration method that is going to be used
*       @param m: of type enum Integration_Method, determines which integration method to use
*
*               This function should be called by TimeLoop in time_loop.h at very beginning
*/
void Set_Integration_Method(enum Integration_Method m){
    integration = m;
}


/**     
*		@brief Process time integration on all_particles, with time step dt
*       @param all_particle: pointer to an array of particles
+       @param dt: time step
*		
*			   Logically, this function is called by the TimeLoop function in time_loop.h
*/
void Time_Integration(Particle* all_particle, double dt){
    switch(integration){
        case EXPLICIT_EULER:    Explicit_Euler  (all_particle, dt); break;
        case HEUN:              Heun            (all_particle, dt); break;
        case MIDPOINT:          Midpoint        (all_particle, dt); break;
        case LEAP_FROG:         Leap_Frog       (all_particle, dt); break;
        default: break;
    }
}

void Explicit_Euler (Particle* all_particle, double dt){
    int N = sizeof(all_particle) / sizeof(all_particle[0]);    

    for(int i = 0 ; i < N ; i++){
        // only for fluid particles
        if(all_particle[i].tag != interior)
            continue;

        all_particle[i].accelerat = div_vec(all_particle[i].pressure_force, all_particle[i].mass);

        all_particle[i].position = add_vec(all_particle[i].position, mul_vec(all_particle[i].velocity, dt));

        all_particle[i].velocity = add_vec(all_particle[i].velocity, mul_vec(all_particle[i].accelerat, dt));
    }
}

void Heun (Particle* all_particle, double dt){
    // TODO: Heun is not implemented
    // Reason: Heun method requires pressure at midpoint
}

void Midpoint (Particle* all_particle, double dt){
    // TODO: Midpoint is not implemented
    // Reason: Midpoint method requires pressure at midpoint
}

void Leap_Frog (Particle* all_particle, double dt){
    int N = sizeof(all_particle) / sizeof(all_particle[0]);    

    for(int i = 0 ; i < N ; i++){
        // only for fluid particles
        if(all_particle[i].tag != interior)
            continue;

        all_particle[i].accelerat = div_vec(all_particle[i].pressure_force, all_particle[i].mass);

        all_particle[i].velocity = add_vec(all_particle[i].velocity, mul_vec(all_particle[i].accelerat, dt));
        
        all_particle[i].position = add_vec(all_particle[i].position, mul_vec(all_particle[i].velocity, dt));

        
    }
}

#endif // TIME_INTEGRATION_H