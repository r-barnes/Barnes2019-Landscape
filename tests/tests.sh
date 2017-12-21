#!/bin/bash

#Edge length of a dataset. Number of cells is the square of this value.
sizes=( 100 700 1000 7000 10000 ) 
sizes=( 200 500 1000 )

#Number of repetitions for each dataset size. Statistical significance!
reps=( 10 10 10 10 10 )           
reps=( 2 2 2 2 2 )   

steps=120        

#Programs to run
exe_prefix=../
progs=( fastscape_BW.exe fastscape_BW+P.exe fastscape_RB.exe fastscape_RB+P.exe fastscape_RB+PQ.exe )

arraylength=${#array[@]}

# use for loop to read all values and indexes

for prog in "${progs[@]}"; do
for (( s=0;   s<${#sizes[@]}; s++ )); do
for (( rep=0; rep<${reps[s]}; rep++ )); do
  size=${sizes[s]}
  echo "# Prog  = $prog"
  echo "m Size  = $size"
  echo "m Steps = $steps"
  echo "m Rep   = $rep"

  eval "$exe_prefix$prog $size $steps"
done
done
done