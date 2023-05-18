#!/bin/bash 
for i in `seq 1 1 129`
do
echo "
$i
5026
50001
283
0.7
2" > input
cat input
survival_time_residuewise *.xtc < input
sleep 59s
        mv survival-time-water-around_residue survival-time-water-around_residue-$i
        wait
done
