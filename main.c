//! @file main.c

#include "data_set.h"
#include "kernel.h"
#include "rate_of_change.h"
#include "time_integration.h"
#include "time_loop.h"
#include "output.h"
#include "constants.h"
#include "validation.h"

/** 
*   @brief: main function
*   @param argc: number of parameters
*   @param argv[0]: name of executable file
*		@param argv[1]: number of interior particles
*		@param argv[2]: output path
*/
int main(int argc, char* argv[]){
	validation_main();
	return EXIT_SUCCESS;
}
