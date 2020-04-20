//! @file main.c
#include <stdio.h>
#include <stdlib.h>

#include "data_set.h"
#include "kernel.h"
#include "rate_of_change.h"
#include "time_integration.h"
#include "time_loop.h"
#include "output.h"

int main(int argc, char* argv[]){

	double t_end = TimeLoop();	

	printf("Finish.\nEnd time is: %f\n", t_end);	
	
	return 0;
}
