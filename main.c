//! @file main.c

#include "data_set.h"
#include "kernel.h"
#include "rate_of_change.h"
#include "time_integration.h"
#include "time_loop.h"
#include "output.h"
#include "constants.h"

/** 
*   @brief: main function
*   @param argc: number of parameters
*   @param argv[0]: name of executable file
*		@param argv[1]: number of interior particles
*		@param argv[2]: output path
*/
int main(int argc, char* argv[]){
	if(argc != 3){
		printf("Wrong parameters!\n");
		return EXIT_FAILURE;
	}

	int particles_interior = atoi(argv[1]);
	set_particles_interior(particles_interior);	

	printf("output path: %s\n", argv[2]);
	set_output_path(argv[2]);

	double t_end = TimeLoop();	

	printf("cycles_DispBoundar    = %.0lf \n", cycles_DispBoundary  );
	printf("cycles_SearchNeighbor = %.0lf \n", cycles_SearchNeighbor);
	printf("cycles_CompGlbDensity = %.0lf \n", cycles_CompGlbDensity);
	printf("cycles_DensityCorr    = %.0lf \n", cycles_DensityCorr   );
	printf("cycles_CompPressure   = %.0lf \n", cycles_CompPressure  );
	printf("cycles_CompAccelerat  = %.0lf \n", cycles_CompAccelerat );
	printf("cycles_RepulsiveForce = %.0lf \n", cycles_RepulsiveForce);
	printf("cycles_TimeIntegral   = %.0lf \n", cycles_TimeIntegral  );
	printf("cycles_all            = %.0lf \n\n", cycles_all);
	printf("Done. End time is: %f\n", t_end);	
	printf("------------------------------------------------\n");

	WritePerformance();

	return EXIT_SUCCESS;
}
