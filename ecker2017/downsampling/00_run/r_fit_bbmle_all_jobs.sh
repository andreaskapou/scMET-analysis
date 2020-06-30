#!/bin/sh

# Run the program
declare -a annos=("distal_H3K27ac_cortex"
                  "H3K4me1_cortex"
                  "prom_2000_2000")

declare -a cells=(200 100 50 20)
# get length of annos
annoslength=${#annos[@]}
cellslength=${#cells[@]}

outdir="/exports/igmm/eddie/ckapoura-XDF/scMET_ms/ecker2017/downsampling/data/"
logPath="logs/"
rep=5 # Number of replications
for (( r=1; r <= $rep; r++ )); do
  for (( c=0; c<${cellslength}; c++ )); do
    for (( i=0; i<${annoslength}; i++ )); do
      my_command="qsub -N mle_${annos[$i]}_${r}_${cells[$c]} \
                       -o ${logPath}mle_${annos[$i]}_${r}_${cells[$c]}.out \
                       -e ${logPath}mle_${annos[$i]}_${r}_${cells[$c]}.error \
                       -v ANNO=${annos[$i]},CELLS=${cells[$c]},REP=${r},OUTDIR=$outdir r_fit_bbmle.sh"
      eval $my_command
    done
  done
done
