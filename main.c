//! @file main.c
#include <stdio.h>
#include <stdlib.h>

#include "data_set.h"
#include "kernel.h"
#include "rate_of_change.h"
#include "time_update.h"
#include "time_loop.h"

int main(int argc, char* argv[]){
	
	double endtime = 1.0;
	TimeLoop(endtime);	

	printf("Finish.\nEnd time is: %f\n", endtime);	
	
	return 0;
}
