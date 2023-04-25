#!/bin/bash
#SBATCH --job-name=0
#SBATCH --account=nn9827k
#SBATCH --time=1-0:00:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=128

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

## Do some work:
cd /cluster/projects/nn9827k/Simulation2D/seismic/submerged/Japan/sensitive/1g

module load OpenMPI/4.0.3-GCC-9.3.0

mpirun -np 1536 /cluster/projects/nn9827k/Uintah/optQA/StandAlone/sus Japan_1g_0.ups

exit 0
