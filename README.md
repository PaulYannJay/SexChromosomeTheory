# SexChromosomeTheory
## Script for analysing the fate of inversion (Figure 2)
Usage: in a bash terminal:

 source ParalleleWhileLoop_IndivSimulPlot_XY.sh #Source the launching file
 export -f slimFc #Export its function
 parallel -j30 slimFc ::: 1e-08 ::: 0.5 ::: -0.005 -0.001 -0.01 -0.05 -0.1 ::: {1..10000} ::: 5000000 15000000 ::: 500000 2000000 1000000 5000000 ::: Y X ### Run on 30 cores, 10000 inversion per parameter combination, with one h and five possible s, inversions of 4 different sizes are considered, either on a X-bearing genome or in a Y-bearing genome, and either in the sex-chromosome (Mid position of the inversion at 5000000) or in the autosome (Mid position of the inversion at 15000000)


## Script for the formation of sexChromosome (Figure 3)
Usage :
parallel -j20 ~/Software/SLiM/build/slim -d N={1} -d mu=1e-08  -d r=1e-5 -d rep={2} -d s={4} -d MaxSizeInv={3} ScriptFormationXYChromosome_VarGamma.slim  ::: 1000 ::: {1..5} ::: 20000000  :::  -0.0001 -0.001 -0.05 -0.01 ### Run on 20 cores. With a maximum inversion size of 20Mb, population of 1000 individuals, and the gamma parameter (for selection coefficient) taking four possible value. For each gamma value, run five simulations. 
The outputs of the long simulation need to be parsed for plotting. For that, use the Parse* files.
../Util/ParseInvFreqOutput.pl -i N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.txt -o N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.Parsed.txt
../Util/ParseRecombinationOutput.pl -i N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.txt -o N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.Parsed.txt
