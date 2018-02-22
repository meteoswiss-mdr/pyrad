#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32g
#SBATCH --partition=normal
#SBATCH --account=msrad

EXECUTABLE=$1

echo "Launched executable $EXECUTABLE"

$EXECUTABLE

