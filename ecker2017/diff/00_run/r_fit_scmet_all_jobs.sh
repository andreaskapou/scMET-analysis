#!/bin/sh

# Run the program
declare -a annos=("distal_H3K27ac_cortex"
                  "H3K4me1_cortex"
                  "prom_2000_2000")
declare -a groups=("Excitatory"
                   "Inhibitory")
# get length of annos
annoslength=${#annos[@]}
groupslength=${#groups[@]}

outdir="/exports/igmm/eddie/ckapoura-XDF/scMET_ms/ecker2017/diff/data/"
useMCMC="" # to use MCMC set it to "--mcmc"
logPath="logs/"
for (( i=0; i<${annoslength}; i++ )); do
  for (( j=0; j<${groupslength}; j++ )); do
    my_command="qsub -N vb_${annos[$i]}_${groups[$j]} \
                     -o ${logPath}vb_${annos[$i]}_${groups[$j]}.out \
                     -e ${logPath}vb_${annos[$i]}_${groups[$j]}.error \
                     -v ANNO=${annos[$i]},GROUP=${groups[$j]},OUTDIR=$outdir,MCMC=$useMCMC r_fit_scmet.sh"
    eval $my_command
  done
done
