
CFLAGS = -O3 -lm

run: main.c time_loop.h time_integration.h data_set.h constants.h
	gcc main.c -o simulation $(CFLAGS)

clean:
	- $(RM) simulation
	- $(RM) -r data_num_particles_*
