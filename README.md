This is the readme file to produce de main figures of the paper. Scripts for producing supplementary figures are available upon request.

# Script for analysing the fate of inversion (Figure 2)
## Figure 2a
 Run the R scripts in:
ModelDeterministic.R

## For figure 2b-c, run the simulation
#in a bash terminal:
First, produce the initiation files (population states after 15,000 generations)
> mkdir ../InitialState #Create a directory for saving files

> parallel -j9 slim -d N=10000 -d mu=1e-8 -d h={1} -d s={2} -d r=1e-6 -d rep=1 ScriptNeutralInversion\_DefineInitialState\_XY.slim :::: 0.5 ::: -0.005 -0.001 -0.01 -0.05 -0.1 #Run on 9 cores, 1 inversion per parameter combination, with various selection coefficient (one h and five possible s) 

Then, for each parameter combination, run 10,000 simulations (10,000 inversions) using the appropriate initial state.

> mkdir ../Output

 >source ParalleleWhileLoop\_IndivSimulPlot\_XY.sh #Source the launching file

 >export -f slimFc #Export its function

 >parallel -j30 slimFc ::: 1e-08 ::: 0.5 ::: -0.005 -0.001 -0.01 -0.05 -0.1 ::: {1..10000} ::: 5000000 15000000 ::: 500000 2000000 1000000 5000000 ::: Y X ### Run on 30 cores, 10000 inversion per parameter combination, with one h and five possible s, inversions of 4 different sizes are considered, either on a X-bearing genome or in a Y-bearing genome, and either in the sex-chromosome (Mid position of the inversion at 5000000) or in the autosome (Mid position of the inversion at 15000000)
## To produce the plots, merge the output file (for instance for 2Mb inversions, Figure 3b)
 >cat N=\*Inv=4000001-6000001\*\_XY.txt >> Linked\_2MbInv\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.txt

 >cat N=\*Inv=14000001-16000001\*\_XY.txt >> Unlinked\_2MbInv\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.txt

 >cat Linked\_2MbInv\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.txt >> LinkedAndUnlinked\_2MbInv\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.txt

 >cat Unlinked\_2MbInv\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.txt >> LinkedAndUnlinked\_2MbInv\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.txt

 ##For figure 3c:

 >cat N=\*\_XY.txt >> AllInv\_LinkedAndUnlinked\_IntroduceInvFromInit\_Nmut\_Freq\_Fit\_IndivPlot\_MidSDR\_XY.tx

# Script for the formation of Sex-Chromosomes (Figure 3)
Usage :
>parallel -j20 ~/Software/SLiM/build/slim -d N={1} -d mu=1e-08  -d r=1e-5 -d rep={2} -d s={4} -d MaxSizeInv={3} ScriptFormationXYChromosome\_VarGamma.slim  ::: 1000 ::: {1..5} ::: 20000000  :::  -0.0001 -0.001 -0.05 -0.01 ### Run on 20 cores. With a maximum inversion size of 20Mb, population of 1000 individuals, and the gamma parameter (for selection coefficient) taking four possible value. For each gamma value, run five simulations. 

The outputs of the long simulation need to be parsed for plotting. For that, use the Parse\* files.
>./ParseInvFreqOutput.pl -i N=1000\_u=1e-08\_r=1e-05\_MaxSizeInv=50000000\_Rep\_6\_Nrecomb\_IndivSimulation\_XY.txt -o N=1000\_u=1e-08\_r=1e-05\_MaxSizeInv=50000000\_Rep\_6\_Nrecomb\_IndivSimulation\_XY.Parsed.txt

>./ParseRecombinationOutput.pl -i N=1000\_u=1e-08\_r=1e-05\_MaxSizeInv=50000000\_Rep\_6\_Nrecomb\_IndivSimulation\_XY.txt -o N=1000\_u=1e-08\_r=1e-05\_MaxSizeInv=50000000\_Rep\_6\_Nrecomb\_IndivSimulation\_XY.Parsed.txt

#Script for producing plots are in: 

PaperPlot.R
