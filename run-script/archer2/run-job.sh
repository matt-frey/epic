#!/bin/bash --login
#SBATCH --job-name=EPICtest
#SBATCH --output=%x.o%j
  # %x gives job-name (SLURM_JOB_NAME)
  # %j gives jobid (individual SLURM_JOB_ID)
  # %A gives jobid (master     SLURM_ARRAY_JOB_ID)
  # %a gives array task id number
  #  https://slurm.schedmd.com/sbatch.html
##SBATCH --open-mode=append
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:19:00
#SBATCH --account=<project-id>
#SBATCH --partition=standard
#SBATCH --qos=short

# refer to ARCHER2 documentation for more information
#     https://docs.archer2.ac.uk/user-guide/

# set number of threads per process
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5
module load cray-netcdf

module list

# assumes MPI-enabled build

# adjust directory to top-level epic install directory (prefix)
EPIC_DIR="/work/<project-id>/<project-id>/$USER/epic"

# adjust config, assumes the case was prepared (.nc file present and namelists/config edited)
srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config model.config

