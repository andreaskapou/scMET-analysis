#!/bin/sh

# Run the program
#cellArray=( 20 50 100 200 500 1000 2000 5000 )
cellArray=( 200 )
#featArray=( 500 )
featArray=( 20 50 100 200 1000 2000 5000 )
useMCMC="--mcmc" # to use MCMC set it to "--mcmc"
logPath="logs/"
rep=5 # Number of replications
for (( i=1; i <= $rep; i++ )); do
  for c in "${cellArray[@]}"; do
    for f in "${featArray[@]}"; do
      my_command="qsub -N s02_J_mcmc_${c}_${i}_${f} \
                       -o ${logPath}s02_J_mcmc_${c}_${i}_${f}.out
                       -e ${logPath}s02_J_mcmc_${c}_${i}_${f}.error
                       -v REP=$i,CELLS=$c,FEATURES=$f,MCMC=$useMCMC r02_fit.sh"
      eval $my_command
    done
  done
done
