#!/bin/bash
#SBATCH --job-name=SandSubBa
#SBATCH --account=nn9827k
#SBATCH --time=4-00:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

## Do some work:
cd /cluster/projects/nn9827k/Simulation2D/SandBox

module load OpenMPI/4.0.3-GCC-9.3.0

mpirun -np 1024 /cluster/projects/nn9827k/Uintah/optQA/StandAlone/sus SandSubBa1.ups

exit 0
