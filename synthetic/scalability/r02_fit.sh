#!/bin/sh
# Grid Engine options (lines prefixed with #$)
# -N s02_diff_res_disp_fit_models
#$ -cwd
#$ -l h_rt=40:00:00
#$ -l h_vmem=10G
# -pe sharedmem 4
#$ -R y

# Initialise PATH to see local scripts
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load igmm/apps/R/3.6.1
module load phys/compilers/gcc/9.1.0

# Run the program
echo ${REP}
echo ${CELLS}
echo ${FEATURES}
echo ${MCMC}
echo ""
my_command="Rscript 02_fit.R --replicate ${REP} --cells ${CELLS} --features ${FEATURES} ${MCMC}"
eval $my_command
