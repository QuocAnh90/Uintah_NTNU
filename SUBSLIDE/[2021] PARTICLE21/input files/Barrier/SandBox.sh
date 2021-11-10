#!/bin/bash
#SBATCH --job-name=k1200
#SBATCH --account=nn9827k
#SBATCH --time=4-0:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

## Do some work:
cd /cluster/projects/nn9827k/SHMPM/SandBox/Barrier

module load OpenMPI/4.0.3-GCC-9.3.0

mpirun -np 512 /cluster/projects/nn9827k/Uintah/optQA/StandAlone/sus SandBox.ups
#tar czf SandBox.tar.gz SandBox.uda.000

exit 0
