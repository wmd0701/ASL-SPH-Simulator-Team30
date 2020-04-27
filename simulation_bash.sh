#/bin/bash
for number_of_particles in 400
do
        mkdir data_num_of_part_$number_of_particles
        ./simulation $number_of_particles data_num_of_part_$number_of_particles
done
