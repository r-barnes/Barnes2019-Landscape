#!/bin/bash

host=$(hostname)

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


if [ ! -f z_serial_comparison.dat ]; then
  progs=( fastscape_BW.exe fastscape_RB.exe )
  sizes=( 100 700 1000 7000 10000 ) 
  reps=( 3 3 3 3 3 )
  for prog in "${progs[@]}"; do
  for (( s=0;   s<${#sizes[@]}; s++ )); do
  for (( rep=0; rep<${reps[s]}; rep++ )); do
    size=${sizes[s]}
    echo "# Prog  = $prog"
    echo "m Size  = $size"
    echo "m Steps = $steps"
    echo "m Rep   = $rep"
    echo "H host  = $host"

    echo "R $exe_prefix$prog $size $steps out_${prog}_${size}_${steps}_${rep}.dem"
    eval "$exe_prefix$prog $size $steps out_${prog}_${size}_${steps}_${rep}.dem"
  done
  done
  done > >(tee -i z_serial_comparison.dat)
fi