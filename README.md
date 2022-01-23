# SPH Simulator

This is the group project for the course Advanced System Lab. The algorithm that we optimize is an improved SPH simulator. For details of the algorithm please have a look into this [paper](https://www.sciencedirect.com/science/article/abs/pii/S0045794912000429?via%3Dihub). **However, the goal of this project is not SPH itself, rather how we can accelerate it using techniques learned in the course.** 

According to the requirement of the course, the project is fully implemented in C. Target platform should be Intel processor with turbo boost disabled.

Please note that this repo on GitHub is imported from [ETHZ INFK GitLab](https://gitlab.inf.ethz.ch/COURSE-ASL/team030). There are four contributors to the project: Silvia Nauer, Mengdi Wang, Tianwei Yu and Valerie Kulka. However, in this imported repo only two of them are shown.


## Visualizations

### Trailer

https://user-images.githubusercontent.com/34072813/150687715-ff3f43a9-565b-4af7-9423-e3ce5eb4331e.mp4

### Work flow
<img src="https://user-images.githubusercontent.com/34072813/150687781-d1d5164a-27b6-4249-9eef-c6a3464d977e.jpeg" width=25% height=25%>

### Flops count
<img src="https://user-images.githubusercontent.com/34072813/150687822-0d72f166-626d-4817-a671-b60a75802cbd.png" width=50% height=50%>


### Flops comparison
<img src="https://user-images.githubusercontent.com/34072813/150687807-f7981351-a8ab-48c1-b55c-9647dfadd462.png" width=50% height=50%>

### Validation
<img src="https://user-images.githubusercontent.com/34072813/150687987-ac81f436-ca6f-47df-8512-6baeb5468c97.png" width=50% height=50%>


## Branches
- master:               baseline implementation
- optimization 1:       + changing data structure, + function inling, + scalar replacement, + reducing branchings
- optimization 2:       + unrolling (unrolling factor 4), + blocking (block size 64)
- optimization 3:       + vectorization (unrolling factor 16)
- optimization 4:       + postponing sqrt
- validation:           specialised for validation

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
- *.sh:                 bash code

### Compile
```sh
make (or make -f Makefile_gcc_9, depending on platform and gcc version)
```

### Usage
```sh
bash simulation_bash.sh (or ./simulation $num_particles $output_path)
```
