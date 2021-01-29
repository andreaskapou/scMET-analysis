#!/bin/sh
# Grid Engine options (lines prefixed with #$)
# -N s02_diff_res_disp_fit_models
#$ -cwd
#$ -l h_rt=5:00:00
#$ -l h_vmem=7G
#$ -pe sharedmem 1
#$ -R y

# Initialise PATH to see local scripts
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.0.3
module load phys/compilers/gcc/9.1.0

# Run the program
echo ${REP}
echo ${CELLS}
echo ${CPGS}
echo ${ORGAMMA}
echo ${MCMC}
echo ""
my_command="Rscript 02_fit.R --replicate ${REP} --cells ${CELLS} --cpgs ${CPGS} --oddsratio ${ORGAMMA} ${MCMC}"
eval $my_command
