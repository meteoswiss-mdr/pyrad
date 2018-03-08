#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
 
module use /apps/escha/UES/RH7.3_PE17.02/sandbox/hvictor/easybuild/modules/all
module load greasy/2.1-cscs-gmvolf-17.02
 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
greasy -f tasks


