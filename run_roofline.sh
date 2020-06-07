#!/bin/bash
for j in master opt1 opt2 opt3 opt4
do
  for i in 400 800 1200 1600 2000 2400 2800 3200 3600 4000   
  do
    perf stat -e LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,cycles,instructions -x \, -r 10 -o $j/result_$i.txt ./simulation_$j $i data > $j/output_$i.txt && echo "$j $i: finished" 
    perf stat -e LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,cycles,instructions -x \, -r 10 -o $j/result_without_$i.txt ./simulation_without_$j $i data > $j/output_without_$i.txt && echo "$j $i: finished without" 
  done
done
