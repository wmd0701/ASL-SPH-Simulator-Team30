#/bin/bash
for number_of_particles in 400
do
        mkdir data_$number_of_particles
        ./simulation $number_of_particles data_$number_of_particles
done
