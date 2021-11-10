#!/bin/bash
#SBATCH --job-name=k1e-6
#SBATCH --account=nn9827k
#SBATCH --time=2-0:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

## Do some work:
cd /cluster/projects/nn9827k/SHMPM/SandBox

module load OpenMPI/4.0.3-GCC-9.3.0

mpirun -np 512 /cluster/projects/nn9827k/Uintah/optQA/StandAlone/sus SandBoxk1e-6.ups
#tar czf SandBox.tar.gz SandBox.uda.000

exit 0
