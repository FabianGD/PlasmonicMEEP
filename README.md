# PlasmonicMEEP

## What this is

This repo consists of a set of scripts for calculation of plasmon resonance/electric field enhancement on different structures.
Scripts are powered by FDTD solver [MEEP](https://github.com/NanoComp/meep).
The scripts are originally taken from [Trel725/plasmon-meep](https://github.com/Trel725/plasmon-meep).

## How it works

1. Meep simulates interaction of structure of interest with EM waves, producing arrays
   describing field distribution in space and time. Simulation is done twice, for
   the empty cell (reference run) and for the cell containing structure.
2. The FFT is performed on arrays from previous step, thus tranforming electric fields
   to frequency domain.
3. Complex values of tranformed fields are squared, producing energiy density at
   given frequncies.
4. Densities of normal run are normalized by reference run, which directly gives
   EM field enhancement distribution in space and frequency

## Installation

### Installation using Nix (the more reliable way)

First, if not already done, install Nix, the functional package manager from [here](https://nixos.org/download.html#nix-quick-install).

Then, to drop into a production shell with all the entrypoints and dependencies correctly set:

```bash
nix-shell -A production
```

To install the package to have it outside of a nix shell:

```bash
nix-env -f ./nix/default.nix -i plasmonic-meep
```

### The standard (old) way

To install the package classically, first install [MEEP](https://meep.readthedocs.io) using the [install instructions](https://meep.readthedocs.io/en/latest/Installation/) on their website. For the calculations proposed here it's re recommended to install the parallel version via **conda**.

#### Installing dependencies

Basically, on a Unix or WSL system, all you need to run the following commands in order. For more detailed instructions, see the [Anaconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [MEEP](https://meep.readthedocs.io) homepages

```bash
# To install conda on your machine. In the installation wizard you can change
# things like prefixes for installation etc.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda-install.sh && sh miniconda-install.sh

# Now, with conda installed (you might need to reload your terminal), install
# MEEP and some other libraries
conda create -n pmeep -c conda-forge pymeep=*=mpi_mpich_* joblib pandas matplotlib h5py mpi4py

# Optional: You can also install the single-core variant using the following command
conda create -n meep -c conda-forge pymeep joblib pandas matplotlib h5py mpi4py
```

#### Install the scripts via pip

As a second step, `cd` to the package directory and install the package using pip. This is somewhat inconsistent, I'm working on a fix. This installs all the entrypoints to your current conda environment, so be sure to have the correct environment activated.

```bash
pip install .
```

After installation, you should find three executables, namely:

- `plas-meep`, which does the FDTD calculations,
- `plas-field`, which fourier transforms the time series data and calculates field enhancements,
- `plas-vis`, that visualises the calculated enhancements as a spectrum and as a map.


### Install for development

The preferred way to install dependencies for development is using nix, see above.
To drop into a development shell (that skips setting entrypoints, as a
editable package would not be "pure" (in the FP sense of pure)),just run the
following in the root directory of the repo:

```bash
nix-shell -A develop
```

To update the git repositories used for querying the dependencies, run in the repo:

```bash
nix-shell -p niv --run "niv update"
```

## Example usage

```bash
# perform step 1, parallel run on 2 processors. The default geometry are two spherical
# gold nanoparticles with a radius of 50 nm and a spacing of 5 nm. The '-r 200' part
# gives the resolution of the computational box in px / µm. The -x and -y flags govern the
# size of the computational box. The '-o' parameter specifies an output directory.
mpirun -np 2 plas-meep -r 200 -x 1 -y 1 -o data/

# Calculate the FFT and the field enhancements in batches from the data calculated
# in the first step. The '-f' and '-w' flags are not required any more.
plas-field -s 100 ./data/plas-meep-norm.h5 ./data/plas-meep-ref.h5 ./data/freqs.h5

# Visualize the calculated field enhancement arrays. This will produce a spectrum
# window first and, after closing it, show the field enhancement maps.
plas-vis -k 50 ./data/freqs.h5
```

## Questions

In case of any questions, please, firstly have a look at the issues (including closed). If there are no answer there, feel free to open new one.

## Usage on the compute cluster ARA (FSU Jena)

To make this program easy to use on ARA, I provide a conda environment on the cluster.
To use it, copy (or soft-link) the `.condarc` (conda configuration) [file](./.condarc)
to your home directory by running the following command. **Be careful:** If you have
a `.condarc` already, it will be overwritten.

```bash
cp .condarc ~/.condarc
```

Now, load the proper Python module, activate the conda environment and you should be
ready to calculate away.

```bash
# This is optional, but good practice, especially when you have some other conda
# instance running.
module purge

# Now you load the Python 3.8 installation. It automatically comes with conda.
module load tools/python/3.8

# If it's the first time you use conda on ARA, you should initialise conda like so:
conda init bash

# Now activate the conda environment that gives you MEEP.
# It comes with MPI for parallel usage.
conda activate pmeep
```

Now, you're able to rock away. To test that you actually have the correct python and MEEP
installed, you can run the following:

```bash
$ which python3
/beegfs/lu27wil/apps/conda-envs/pmeep/bin/python3

$ python3 -c "import meep; print(meep.__version__)"
Using MPI version 3.1, 1 processes
1.16.1
```

### Custom installation

In case you want to change the PlasmonicMEEP code yourself, you need your own installation
and your own conda environment. In order to get it to work, you need to do the following:

```bash
# 1. Unload all the modules that may disturb the installation.
module purge

# 2. Load the newest Python3 module from the module list
module load tools/python/3.8

# 3. Create a new environment and install the necessary packages
conda create -n pmeep -c conda-forge pymeep=*=mpi_mpich_* joblib pandas matplotlib h5py mpi4py

# 4. After installation (this may take a while ...), you can activate your new environment ...
conda activate pmeep

# 5. ... and install the package.
# 5. a. Clone the upstream repository
git clone \
   https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmonic-meep.git \
   /beegfs/$USER/plasmonic-meep

# 5. b. Now go to the directory ...
cd /beegfs/$USER/plasmonic-meep

# 5. c. ... and install the package using pip. The "-e" allows you to edit the files directly.
pip install -e .
```

### Running jobs on the cluster

**WARNING: Never ever run jobs on the front nodes! Always use output directories on the parallel file system (`/beegfs/$USER/...`)!**

With that out of the way, in the root directory of the repository, you find two files [`pmeep_submit.sh`](./pmeep_submit.sh) and [`pfield_submit.sh`](./pfield_submit.sh) for submitting simulation and field enhancement jobs, respectively. Those files are to be used with the SLURM queuing system.

To submit a classic PlasmonicMEEP job, you need to use `sbatch`:

```bash
# The first argument is the output directory.
# The rest are parameters you would like to parse to plas-meep
sbatch pmeep_submit.sh /beegfs/$USER/<outputdir> -r 400 -x 0.5 -y 0.5 [...]

# After this calculation has finished, you need to calculate the field
# enhancements using the pfield_submit.sh submit script. You may also specify
# additional arguments after the output directory
sbatch pfield_submit.sh /beegfs/$USER/<outputdir> [-s ...]

```

## Notes

The [multiviewer.py](./plasmonicmeep/multiviewer.py) module is taken and modified
from [Datacamp](https://www.datacamp.com/community/tutorials/matplotlib-3d-volumetric-data).
