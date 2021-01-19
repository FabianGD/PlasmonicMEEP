#!/usr/bin/env bash
#SBATCH --job-name=pyQD
#SBATCH --partition=s_standard
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=%J.out
#SBATCH --mem-per-cpu=3GB

# Load modules
module purge
module load tools/python/3.8

eval "$(conda shell.bash hook)"

# Load the pyQD conda environment
conda activate pmeep

set -x

export TMPDIR=/beegfs/$USER/tmp
mkdir -p $TMPDIR

# Get all the args as an array.
args=( $@ )

# Get the necessary info
dir=${args[@]:0:1}
cargs=${args[@]:1}

mpirun plas-meep -o $dir $cargs
