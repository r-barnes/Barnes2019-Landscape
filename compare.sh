#!/bin/bash

rm -rf out_*.dem*

./fastscape_BW.exe     501 120 out_BW.dem     123
./fastscape_BW+P.exe   501 120 out_BW+P.dem   123
./fastscape_BW+PI.exe  501 120 out_BW+PI.dem  123
./fastscape_RB.exe     501 120 out_RB.dem     123
./fastscape_RB+P.exe   501 120 out_RB+P.dem   123
./fastscape_RB+PI.exe  501 120 out_RB+PI.dem  123
./fastscape_RB+PQ.exe  501 120 out_RB+PQ.dem  123
./fastscape_RB+GPU.exe 501 120 out_RB+GPU.dem 123

rd_compare authoritative2.dem out_BW.dem    
rd_compare authoritative2.dem out_BW+P.dem  
rd_compare authoritative2.dem out_BW+PI.dem  
rd_compare authoritative2.dem out_RB.dem    
rd_compare authoritative2.dem out_RB+P.dem  
rd_compare authoritative2.dem out_RB+PI.dem  
rd_compare authoritative2.dem out_RB+PQ.dem 
rd_compare authoritative2.dem out_RB+GPU.dem 
