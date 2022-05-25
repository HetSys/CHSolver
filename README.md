# CHSolver
HetSys PX915 group project (group B): Cahn-Hilliard Solver

## Requirements
Before building the dependencies and compiling the code, ensure that you have the following available:
- GCC/10.2.0
- CMake/3.20.1
- OpenMPI/4.0.5
- HDF5/1.10.7
- Python/3.8.6
- FFTW/3.3.8

If you using [Environment Modules](http://modules.sourceforge.net/) (eg on Warwick's HetMathSys nodes) then this can be accomplished via
```bash
module purge
module load GCC/10.2.0 CMake/3.20.1 OpenMPI/4.0.5 HDF5/1.10.7 Python/3.8.6 FFTW/3.3.8
```

## Installation Instructions
First clone the repository (and the required submodules) onto your machine:
```bash
git clone --recurse-submodules https://github.com/HetSys/CHSolver
```

## Building
### Dependencies
The code requires a number of third party dependencies. We have included a build script
```bash
./build_deps
```
that will call the correct commands to build these locally. If you also require the ability to run unit tests you will also need to run
```bash
./build_pfunit
```
This may take a while (particularly `./build_pfunit`), but this only needs to be done once.

### Compiling
Once these have been built you can compile the code with
```bash
make
```
which will produce the executable `chsolver`

## Running
To run the code you must fully specify the output times and solver parameters. We have included the a script
```bash
./restore_json
```
that we generate a valid JSON input file, but these can be replaced (or overridden) using command line flags (see `./chsolver -h` for a full list).

Since the code uses MPI the code must be run with
```bash
mpirun -np <num_procs> ./chsolver
```
`<num_procs>` must either be 1 (for the pseudo-spectral solver) or a power of 4 (for the finite-difference solver).

## Output Visualisation
Once the output data has been generated (in the `out` directory), the script
```bash
./run_visualisation
```
will produce an mp4 of the system over time and an image of the energy decay.
