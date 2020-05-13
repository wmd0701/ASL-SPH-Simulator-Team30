//! @file time_integration.h
#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "data_set.h"

enum Integration_Method {EXPLICIT_EULER, LEAP_FROG};
enum Integration_Method integration;


/**     
*		@brief different time integration methods
*       @param dt: time step
*		
*			   One of the following four methods will be called in Time_Integration
*/

void Explicit_Euler (double dt){
	for(int i = 0 ; i < N_interior ; i++){
		
		positions[i] = (vector){positions[i].first  + velocities[i].first  * dt, 
								positions[i].second + velocities[i].second * dt};
		
		velocities[i] = (vector){velocities[i].first  + accelerats[i].first  * dt,
								 velocities[i].second + accelerats[i].second * dt};
	}
}

void Leap_Frog (double dt){
	for(int i = 0 ; i < N_interior ; i++){
		
		velocities[i] = (vector){velocities[i].first  + accelerats[i].first  * dt,
								 velocities[i].second + accelerats[i].second * dt};

		positions[i] = (vector){positions[i].first  + velocities[i].first  * dt, 
								positions[i].second + velocities[i].second * dt};
	}
}

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
+       @param dt: time step
*		
*			   Logically, this function is called by the TimeLoop function in time_loop.h
*/
void Time_Integration (double dt){
	switch(integration){
		case EXPLICIT_EULER:    Explicit_Euler  (dt); break;
		case LEAP_FROG:         Leap_Frog       (dt); break;
		default: break;
	}
}

#endif // TIME_INTEGRATION_H
