#/bin/bash
for smoothing_length in 0.02 0.03 0.04 0.05 0.06 0.07
do
        mkdir data_smoothing_length_$smoothing_length
        ./simulation $smoothing_length data_smoothing_length_$smoothing_length
done
