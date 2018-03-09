#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=postproc
#SBATCH --account=msrad

# Attention: cpus-per-task refers to processors!
# ntasks-per-core referes to processors!

# KESCH specs:
# nodes: 12
# CPUs per node: 2
# Cores per CPU: 12
# Total cores per node : 24
# Memory per node 256 GB
# Total processors : 288
# Total memory : 3 TB
#
# To check specs of partition use: scontrol show partition
# Partition normal:
#   nodes : 11
#   CPUs (processors): 264 
#   Maximum processing time: 1 day
#
# Partition postproc:
#   nodes : 5
#   CPUs (processors) : 120
#   Maximum processing time : 1 day

module use /apps/escha/UES/RH7.3_PE17.02/sandbox/hvictor/easybuild/modules/all
module load greasy/2.1-cscs-gmvolf-17.02

TASKS_FILE=$1
 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
greasy -f ${TASKS_FILE}
