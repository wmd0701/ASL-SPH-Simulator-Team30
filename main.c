//! @file main.c
#include <stdio.h>
#include <stdlib.h>

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
*		@param argv[1]: smoothing length
*		@param argv[2]: output path
*/
int main(int argc, char* argv[]){
	if(argc != 3){
		printf("Wrong parameters!\n");
		return EXIT_FAILURE;
	}

	double smoothing_length = atof(argv[1]);
	printf("smoothing length: %lf\n", smoothing_length);
	set_smoothing_length(smoothing_length);

	printf("output path: %s", argv[2]);
	set_output_path(argv[2]);

	double t_end = TimeLoop2();	

	printf("Done. End time is: %f\n", t_end);	

	return EXIT_SUCCESS;
}
