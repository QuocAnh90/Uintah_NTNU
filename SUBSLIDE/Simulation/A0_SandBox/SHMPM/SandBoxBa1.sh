#!/bin/bash
#SBATCH --job-name=Bak1
#SBATCH --account=nn9827k
#SBATCH --time=4-0:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

## Do some work:
cd /cluster/projects/nn9827k/SHMPM/SandBox/Barrier

module load OpenMPI/4.0.3-GCC-9.3.0

mpirun -np 1024 /cluster/projects/nn9827k/Uintah/optTest/StandAlone/sus SandBoxBa1.ups
#tar czf SandBoxBa.tar.gz SandBoxBa.uda.000

exit 0
