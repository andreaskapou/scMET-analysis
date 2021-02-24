#!/bin/sh

# Run the program
cellArray=( 20 50 100 200 500)
cpgs=8 #(15 50)
orMu=( 2 3 5 )
useMCMC="" # to use MCMC set it to "--mcmc"
logPath=" logs/"
rep=10 # Number of replications
for (( i=1; i <= $rep; i++ )); do
  for c in "${cellArray[@]}"; do
    for cpg in "${cpgs[@]}"; do
      for or in "${orMu[@]}"; do
        my_command="qsub -N vb_orM_${or}_${c}_${cpg}_${i} \
                        -o ${logPath}/vb_orG_${or}_${c}_${cpg}_${i}.out \
                        -e ${logPath}/vb_orG_${or}_${c}_${cpg}_${i}.error \
                        -v REP=$i,CELLS=$c,CPGS=$cpg,ORMU=$or,MCMC=$useMCMC r02_fit.sh"
        eval $my_command
      done
    done
  done
done
