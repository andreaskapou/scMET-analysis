#!/bin/sh

# Run the program
declare -a annos=("window20000_step20000") #("window10000_step10000")
                  #"window20000_step20000")
# get length of annos
annoslength=${#annos[@]}

declare -a chrs=("chr2", "chr7", "chr15", "chrX")
# get length of annos
chrslength=${#chrs[@]}

outdir="/exports/igmm/eddie/ckapoura-XDF/scMET_ms/ecker2017/all_cells/data/"
useMCMC="" # to use MCMC set it to "--mcmc"
logPath="logs/"
for (( i=0; i<${annoslength}; i++ )); do
  for (( j=0; j<${chrslength}; j++ )); do
    my_command="qsub -N vb_${annos[$i]}_${chrs[$j]} \
                     -o ${logPath}vb_${annos[$i]}_${chrs[$j]}.out \
                     -e ${logPath}vb_${annos[$i]}_${chrs[$j]}.error \
                     -v ANNO=${annos[$i]},CHR=${chrs[$j]},OUTDIR=$outdir,MCMC=$useMCMC r_fit_scmet_window.sh"
    eval $my_command
  done
done
