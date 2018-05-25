#!/bin/bash

host=$(hostname)

steps=120        

#Programs to run
exe_prefix=../

if [ "$TESTSYSTEM" == "summitdev" ]; then

  if [ ! -f "z_gpu_scaling_$TESTSYSTEM.dat" ]; then
    echo "RUNNING GPU TESTS"

    progs=( fastscape_RB+GPU.exe )

    #Edge length of a dataset. Number of cells is the square of this value.
    sizes=( 100 200 400 800 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000 9500 10000 10500 11000 11500 12000 12500 13000 13500 14000 14500 15000 15500 16000 16500 17000 17500 18000 18500 19000 19500 20000 ) 

    #Number of repetitions for each dataset size. Statistical significance!
    reps=(    1   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1 )
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
    done > >(tee -i "z_gpu_scaling_$TESTSYSTEM.dat")
  fi

  if [ ! -f "z_gpu_$TESTSYSTEM.dat" ]; then
    echo "RUNNING GPU TESTS"

    progs=( fastscape_RB+GPU.exe )

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
    done > >(tee -i "z_gpu_$TESTSYSTEM.dat")
  fi
fi



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
