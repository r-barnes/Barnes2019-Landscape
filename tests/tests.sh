#!/bin/bash

host=$(hostname)


steps=120        

#Programs to run
exe_prefix=../

if [ ! -f "z_parallel_sep_thread_$TESTSYSTEM.dat" ]; then
  echo "RUNNING PARALLEL IMPROVED TESTS"

  progs=( fastscape_RB+PQ.exe )

  #Edge length of a dataset. Number of cells is the square of this value.
  sizes=( 100 700 1000 7000 10000 ) 

  #Number of repetitions for each dataset size. Statistical significance!
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

    echo "R $exe_prefix$prog $size $steps out_simp_parallel_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
    eval "$exe_prefix$prog $size $steps out_simp_parallel_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
  done
  done
  done > >(tee -i "z_parallel_sep_thread_$TESTSYSTEM.dat")
fi

exit 0

if [ ! -f "z_parallel_improved_$TESTSYSTEM.dat" ]; then
  echo "RUNNING PARALLEL IMPROVED TESTS"

  progs=( fastscape_BW+PI.exe fastscape_RB+PI.exe )

  #Edge length of a dataset. Number of cells is the square of this value.
  sizes=( 100 700 1000 7000 10000 ) 

  #Number of repetitions for each dataset size. Statistical significance!
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

    echo "R $exe_prefix$prog $size $steps out_simp_parallel_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
    eval "$exe_prefix$prog $size $steps out_simp_parallel_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
  done
  done
  done > >(tee -i "z_parallel_improved_$TESTSYSTEM.dat")
fi




if [ ! -f "z_serial_comparison_$TESTSYSTEM.dat" ]; then
  echo "RUNNING SERIAL TESTS"

  progs=( fastscape_BW.exe fastscape_RB.exe )

  #Edge length of a dataset. Number of cells is the square of this value.
  sizes=( 100 700 1000 7000 10000 ) 

  #Number of repetitions for each dataset size. Statistical significance!
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

    echo "R $exe_prefix$prog $size $steps out_serial_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
    eval "$exe_prefix$prog $size $steps out_serial_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
  done
  done
  done > >(tee -i "z_serial_comparison_$TESTSYSTEM.dat")
fi


if [ ! -f "z_simple_parallel_$TESTSYSTEM.dat" ]; then
  echo "RUNNING SIMPLE PARALLEL TESTS"

  progs=( fastscape_BW+P.exe fastscape_RB+P.exe )

  #Edge length of a dataset. Number of cells is the square of this value.
  sizes=( 100 700 1000 7000 10000 ) 

  #Number of repetitions for each dataset size. Statistical significance!
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

    echo "R $exe_prefix$prog $size $steps out_simp_parallel_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
    eval "$exe_prefix$prog $size $steps out_simp_parallel_${prog}_${size}_${steps}_${rep}_${TESTSYSTEM}.dem 123"
  done
  done
  done > >(tee -i "z_simple_parallel_$TESTSYSTEM.dat")
fi