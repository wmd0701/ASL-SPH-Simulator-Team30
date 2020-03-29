
CFLAGS = -O3 -lm

run: main.c time_loop.h time_update.h kernel.h data_set.h 
	gcc main.c -o simulation $(CFLAGS)

clean:
	- $(RM) -r *.dSYM
	- $(RM) *.o 