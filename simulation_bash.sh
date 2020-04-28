#/bin/bash
for num_particles in 400 800 1200 1600 2000 2400 2800 3200 3600 4000
do
        mkdir data_num_particles_$num_particles
        ./simulation $num_particles data_num_particles_$num_particles
done
