# Building and running EPIC on ARCHER2

Follow the general build instructions with a few changes (see `../../docs`).

On ARCHER2 some dependencies are pre-installed and can be loaded as modules. Alternatively the dependencies can be installed separately by the user (make sure to use a location on `/work` instead of `/home` and set `PATH` and `LD_LIBRARY_PATH` environment variables as appropriate).

More information about using ARCHER2 can be found in the documentation: [https://docs.archer2.ac.uk/](https://docs.archer2.ac.uk/)


## Build adjustments

Before running the `./bootstrap` command, the `configure.ac` needs to be adjusted in the following way. Find a test for `MPI_DIR` and cahnge the `LIBS=` line to `-lmpichf90` to link with the mpi library.

Then load existing dependencies (hdf5 and netcdf).

`module load cray-hdf5`
`module load cray-netcdf`

Alternatively, to enable parallel IO `cray-hdf5-parallel` and `cray-netcdf-hdf5parallel` modules can be used instead.

Next run the `configure` script and `make && make install`. Make sure you specify the `--prefix=...` option as users on ARCHER2 have no permission to install into system-wide locations as well as `MPI_DIR`. The configure line may look like the example below:

`CC=cc FC=ftn MPI_DIR=$CRAY_MPICH_DIR ../configure --prefix=/work/<project-id>/<project-id>/$USER/epic3d --enable-openmp --enable-debug --enable-verbose --enable-3d --enable-unit-tests`

Note that currently Python version 2.7 is the default on ARCHER2, thus if you need Python 3, you would need to load the respecive module: `module load cray-python`.


## Running a model using slurm

We assume the model has been prepared use the documentation from the examples directory, including initial/boundary conditions, namelist and config settings.

An example batch run script is provided: see `run-job.sh`, which needs to be edited to set project-id and adjust paths as appropriate. Then the job is submitted from an ARCHER2 login node using `sbatch run-job.sh`. More information on Slurm options e.g. queue names and limits can be found in the [ARCHER2 online documentation](https://docs.archer2.ac.uk/user-guide/).

In particular, the `--nodes=...` setting specifies the number of ARCHER2 nodes to be used (note that each node has 128 cores), whilst `--ntasks=...` specifies the number of MPI processes. Generally the value of the sum `--nodes` plus `--ntasks-per-node` should equal `--ntasks` value.

