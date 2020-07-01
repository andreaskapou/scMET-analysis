#!/bin/sh

# Run the program
declare -a annos=("prom_2000_2000")

# get length of annos
annoslength=${#annos[@]}

outdir="/exports/igmm/eddie/ckapoura-XDF/scMET_ms/gastrulation/data/"
useMCMC="" # to use MCMC set it to "--mcmc"
logPath="logs/"
for (( i=0; i<${annoslength}; i++ )); do
  my_command="qsub -N vb_${annos[$i]} \
                   -o ${logPath}all_vb_${annos[$i]}.out \
                   -e ${logPath}all_vb_${annos[$i]}.error \
                   -v ANNO=${annos[$i]},OUTDIR=$outdir,MCMC=$useMCMC r_fit_scmet.sh"
  eval $my_command
done
