
CFLAGS = -O3 -lm

run: main.c time_loop.h time_integration.h kernel.h data_set.h constants.h
	gcc main.c -o simulation $(CFLAGS)

clean:
	- $(RM) -r *.dSYM
	- $(RM) *.o 
