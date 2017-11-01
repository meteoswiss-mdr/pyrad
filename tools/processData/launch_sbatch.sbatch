#!/bin/bash -l
#SBATCH --job-name="rad4alp_hydro"
#SBATCH --time=24:00:00
#SBATCH --output=rad4alp_hydro.out
#SBATCH --error=rad4alp_hydro.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=32g
#SBATCH --partition=normal
#SBATCH --account=msrad

./get_and_process_rad4alp_hydro_cscs.sh

