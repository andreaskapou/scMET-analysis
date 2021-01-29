#!/bin/sh

# Run the program
declare -a annos=("window20000_step20000") #("window10000_step10000")
                  #"window20000_step20000")
# get length of annos
annoslength=${#annos[@]}

declare -a chrs=("chr1" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr16" "chr17" "chr18" "chr19")
# get length of annos
chrslength=${#chrs[@]}

outdir="/exports/igmm/eddie/ckapoura-XDF/scMET_ms/ecker2017/all_cells/data/"
useMCMC="" # to use MCMC set it to "--mcmc"
logPath="logs/"
i=0
#for (( i=0; i<${annoslength}; i++ )); do
for (( j=0; j<${chrslength}; j++ )); do
    echo $chrs[$j]
    my_command="qsub -N vb_${chrs[$j]}_${annos[$i]} \
                     -o ${logPath}vb_${annos[$i]}_${chrs[$j]}.out \
                     -e ${logPath}vb_${annos[$i]}_${chrs[$j]}.error \
                     -v ANNO=${annos[$i]},CHR=${chrs[$j]},OUTDIR=$outdir,MCMC=$useMCMC r_fit_scmet_window.sh"
    eval $my_command
done
#done
