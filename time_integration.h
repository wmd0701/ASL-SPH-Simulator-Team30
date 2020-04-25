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
void Leap_Frog      (Particle* all_particle, double dt);
void Heun           (Particle* all_particle, double dt);
void Heun_Half      (Particle* all_particle, double dt);
void Midpoint       (Particle* all_particle, double dt);
void Midpoint_Half  (Particle* all_particle, double dt);

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
void Time_Integration (Particle* all_particle, double dt){
	switch(integration){
		case EXPLICIT_EULER:    Explicit_Euler  (all_particle, dt); break;
		case LEAP_FROG:         Leap_Frog       (all_particle, dt); break;
		case HEUN:              Heun            (all_particle, dt); break;
		case MIDPOINT:          Midpoint        (all_particle, dt); break;
		default: break;
	}
}

/**     
*		@brief Process half time integration, only for Heun and Midpoint methods
*       @param all_particle: pointer to an array of particles
+       @param dt: time step
*		
*			   Logically, this function is called by the TimeLoop function in time_loop.h
*/
void Time_Integration_Half (Particle* all_particle, double dt){
	switch(integration){
		case HEUN:              Heun_Half       (all_particle, dt); break;
		case MIDPOINT:          Midpoint_Half   (all_particle, dt); break;
		default: break;
	}
}

void Explicit_Euler (Particle* all_particle, double dt){
	for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
		// only for fluid particles
		if(all_particle[i].tag != interior)
			continue;

		all_particle[i].position = vec_add_vec(all_particle[i].position, vec_mul_scalar(all_particle[i].velocity, dt));
		// =================== XSPH correction ======================
		/*
			 for (p = all_particle[i].neighbors; p != NULL; p = p->next) {
			 all_particle[i].position.first  += 0.5 * dt * all_particle[p->idx].mass / all_particle[p->idx].density * p->Wij * (all_particle[p->idx].velocity.first  - all_particle[i].velocity.first);
			 all_particle[i].position.second += 0.5 * dt * all_particle[p->idx].mass / all_particle[p->idx].density * p->Wij * (all_particle[p->idx].velocity.second - all_particle[i].velocity.second);
			 }*/
		// ==========================================================
		all_particle[i].velocity = vec_add_vec(all_particle[i].velocity, vec_mul_scalar(all_particle[i].accelerat, dt));
	}
}

void Heun_Half (Particle* all_particle, double dt){
	for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
		// only for fluid particles
		if(all_particle[i].tag != interior)
			continue;

		vector velocity_Euler = vec_add_vec(all_particle[i].velocity, vec_mul_scalar(all_particle[i].accelerat, dt));

		all_particle[i].position_help = vec_add_vec(all_particle[i].position, vec_mul_scalar(vec_add_vec(velocity_Euler, all_particle[i].velocity), 0.5*dt));

		all_particle[i].position = vec_add_vec(all_particle[i].position, vec_mul_scalar(all_particle[i].velocity, dt));

		all_particle[i].accelerat_help = all_particle[i].accelerat;
	}
}

void Heun (Particle* all_particle, double dt){
	for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
		// only for fluid particles
		if(all_particle[i].tag != interior)
			continue;

		all_particle[i].position = all_particle[i].position_help;

		all_particle[i].velocity = vec_add_vec(all_particle[i].velocity, vec_mul_scalar(vec_add_vec(all_particle[i].accelerat, all_particle[i].accelerat_help), 0.5*dt));
	}
}

void Midpoint_Half (Particle* all_particle, double dt){
	for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
		// only for fluid particles
		if(all_particle[i].tag != interior)
			continue;

		vector velocity_Mid = vec_add_vec(all_particle[i].velocity, vec_mul_scalar(all_particle[i].accelerat, 0.5*dt));

		all_particle[i].position_help = vec_add_vec(all_particle[i].position, vec_mul_scalar(velocity_Mid, dt));

		all_particle[i].position = vec_add_vec(all_particle[i].position, vec_mul_scalar(all_particle[i].velocity, 0.5*dt));
	}
}
void Midpoint (Particle* all_particle, double dt){
	for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
		// only for fluid particles
		if(all_particle[i].tag != interior)
			continue;

		all_particle[i].position = all_particle[i].position_help;

		all_particle[i].velocity = vec_add_vec(all_particle[i].velocity, vec_mul_scalar(all_particle[i].accelerat, dt));
	}
}

void Leap_Frog (Particle* all_particle, double dt){
	for(int i = 0 ; i < NUMBER_OF_PARTICLE ; i++){
		// only for fluid particles
		if(all_particle[i].tag != interior)
			continue;

		all_particle[i].velocity = vec_add_vec(all_particle[i].velocity, vec_mul_scalar(all_particle[i].accelerat, dt));

		all_particle[i].position = vec_add_vec(all_particle[i].position, vec_mul_scalar(all_particle[i].velocity, dt));
	}
}

#endif // TIME_INTEGRATION_H
