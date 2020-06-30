#!/bin/sh
# Grid Engine options (lines prefixed with #$)
# -N fit_mle
#$ -cwd
#$ -l h_rt=5:00:00
#$ -l h_vmem=10G
#$ -pe sharedmem 1
#$ -R y

# Initialise PATH to see local scripts
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1
module load phys/compilers/gcc/9.1.0

# Run the program
echo ${ANNO}
echo ${CELLS}
echo ${REP}
echo ${OUTDIR}
echo ""
my_command="Rscript fit_bbmle.R --anno ${ANNO} --cells ${CELLS} --replicate ${REP} --outdir ${OUTDIR} "
eval $my_command
