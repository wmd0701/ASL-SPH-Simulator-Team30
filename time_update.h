//!     @file time_update.h

#include "data_set.h"

/**     
*		@brief Initiate the data structure as well as initial condition
*		
*			   Step 1: Memory allocation
*			   Step 2: Assign tag, position, mass and velocity to every particles
*			   Step 3: Establish the nearby relaiton 
*
*		@return	pointer to an array containing information of all the particles
*/
Particle* Init ();

/**     
*		@brief Compute the new velocity and position of every interior particle
*
*			   Meanwhile, update the nearby relation.
*			   Choose suitable time step scheme.
*
*		@param all_particle pointer to an array containing information of all the particles
*		@param dt time step
*		@return no returns. Update the [velocity] [position] attributes in all_particle
*/
void TimeUpdate (Particle*& all_particle, double dt);