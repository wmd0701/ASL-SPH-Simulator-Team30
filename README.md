# SPH Simulator

This is the project for the course Advanced System Lab. The algorithm that we choose is an improved SPH simulator. For details of the algorithm please have a look into https://www.sciencedirect.com/science/article/abs/pii/S0045794912000429?via%3Dihub.

According to the requirement of the course, the project is fully implemented in C. Target platform should be Intel processor with turbo boost disabled.

## Files
- Makefile:             makefile
- Makefile_gcc_9:       makefile specialized for gcc 9.3 on Windows subsystem Linux
- kernel.h:             SPH kernel function and its gradient
- data_set.h:           data structure and intialization
- constants.h:          constants
- output.h:             output testing results to file, useful for validation
- rate_of_change.h:     computation
- time_integration.h:   time integration methods
- time_loop.h:          main loop
- tsc_x86.h:            count cycles
- main.c:               main function

## Directories
- data:                 output files

### Compile
```sh
make (or make -f Makefile_gcc_9, depending on platform and gcc version)
```
### Usage
```sh
bash simulation_bash.sh (or ./simulation $num_particles $output_path)
```

### Parameters
```sh
```