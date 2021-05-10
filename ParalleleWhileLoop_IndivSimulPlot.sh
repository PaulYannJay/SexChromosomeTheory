#!/bin/bash

slimFc(){
Init=`ls ../InitialState/slim_g15000_MidSDR_10MChrom_N=1000_r=1e-06_u=$1_s=$3_h=$2_*` #Define the initiation files corresponding to the set of parameter
position=$5 #MId position of the inversion (position of the inversion marker
size=$6
end=$((position+size/2+1)) #Start of the inverion
start=$((position-size/2+1)) #End of the inversion
~/Software/SLiM/build/slim -d N=1000 -d mu=$1 -d h=$2 -d s=$3 -d r=1e-6 -d rep=$4 -d start=$start -d end=$end -d "Init='$Init'" Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom.slim
}



