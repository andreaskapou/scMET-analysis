#!/bin/sh
# Grid Engine options (lines prefixed with #$)
# -N fit_bayes
#$ -cwd
#$ -l h_rt=30:00:00
#$ -l h_vmem=60G
# -pe sharedmem 1
#$ -R y

# Initialise PATH to see local scripts
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH

# Initialise the environment modules
. /etc/profile.d/modules.sh

module load igmm/apps/R/4.0.3
module load phys/compilers/gcc/9.1.0

# Run the program
echo ${ANNO}
echo ${GROUP}
echo ${OUTDIR}
echo ${MCMC}
echo ${CHR}
echo ""
my_command="Rscript fit_scmet_window.R --anno ${ANNO} --group ${GROUP} --chr ${CHR} --outdir ${OUTDIR} ${MCMC}"
eval $my_command
