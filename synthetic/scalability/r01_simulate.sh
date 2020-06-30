#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N s01_simulate
#$ -cwd
#$ -l h_rt=2:00:00
#$ -l h_vmem=5G
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
my_command="Rscript 01_simulate.R "
eval $my_command
