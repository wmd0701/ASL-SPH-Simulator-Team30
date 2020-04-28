#/bin/bash

# If you modify bash code on Windows, you may get error message "unexpected token \r" when running bash code.
# The reason lies in character difference between Linux and Windows. To solve the problem, please follow:
# https://stackoverflow.com/questions/27176781/bash-file-returns-unexpected-token-do-r

for num_particles in 400 800 1200 1600 2000 2400 2800 3200 3600 4000
do
        mkdir data_num_particles_$num_particles
        ./simulation $num_particles data_num_particles_$num_particles
done
