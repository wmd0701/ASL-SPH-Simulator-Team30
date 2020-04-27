#/bin/bash
for number_of_particles in 400 800 1200 1600 2000 2400 2800 3200 3600 4000
do
        mkdir data_num_of_part_$number_of_particles
        ./simulation $number_of_particles data_num_of_part_$number_of_particles
done
