#!/bin/bash

#SBATCH --job-name=spheres
#SBATCH --partition=fat
#SBATCH --mail-type=end
#SBATCH --mail-user=qyrong@sjtu.edu.cn
#SBATCH -n 1
#SBATCH --output=spheres-%j.out
#SBATCH --error=spheres-%j.err

ulimit -s unlimited         
ulimit -l unlimited 

source /usr/share/Modules/init/bash
module purge
module load icc impi

./spheres ./input/input.txt
