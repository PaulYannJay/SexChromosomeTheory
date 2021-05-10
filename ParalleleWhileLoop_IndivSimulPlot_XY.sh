#!/bin/bash
slimFc(){
Init=`ls ../InitialState/slim_g15000_MidSDR_10MChrom_XY_N=1000_r=1e-06_u=$1_s=$3_h=$2_*`#Define the initiation file to use corresponding to the set of parameter
position=$5 #Mid position of the inversion (position of the inversion marker
size=$6
end=$((position+size/2+1))
start=$((position-size/2+1))
~/Software/SLiM/build/slim -d N=1000 -d mu=$1 -d h=$2 -d s=$3 -d r=1e-6 -d rep=$4 -d start=$start -d end=$end -d "Init='$Init'" -d "SexChrom='$7'" Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_XY.slim
}

#Usage: in a bash terminal
## source ParalleleWhileLoop_IndivSimulPlot_XY.sh
## export -f slimFc
## parallel -j30 slimFc ::: 1e-08 ::: 0.5 ::: -0.005 -0.001 -0.01 -0.05 -0.1 ::: {1..10000} ::: 5000000 15000000 ::: 500000 2000000 1000000 5000000 ::: Y X ### Run on 30 cores, 10000 inversion per parameter combination, with one h and five different s, inversions of 4 different sizes are considered, either on a X-bearing genome or in a Y-bearing genome, and either in the sex-chromosome (Mid position of the inversion at 5000000) or in the autosome (Mid position of the inversion at 15000000)

