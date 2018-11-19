#!/bin/bash

#SBATCH --job-name=discretize
#SBATCH --partition=fat
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -n 1
#SBATCH --output=discretize-%j.out
#SBATCH --error=discretize-%j.err

ulimit -s unlimited         
ulimit -l unlimited 

source /usr/share/Modules/init/bash
module purge
module load icc impi

./a.out
