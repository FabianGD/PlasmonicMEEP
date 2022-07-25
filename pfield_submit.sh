#!/usr/bin/env bash
#SBATCH --job-name=pyQD
#SBATCH --partition=s_standard
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --output=%J.out
#SBATCH --mem-per-cpu=5GB

# Load modules
module purge

eval "$(/beegfs/lu27wil/miniconda3/bin/conda shell.bash hook)"

# Load the pyQD conda environment
conda activate pmp

set -x

export TMPDIR=/beegfs/$USER/tmp
mkdir -p $TMPDIR

# Get all the args as an array.
args=( $@ )

# Get the necessary info
dir=$(realpath ${args[@]:0:1})
echo $dir

cargs=${args[@]:1}
echo $cargs

plas-field $dir/plas-meep-norm.h5 $dir/plas-meep-ref.h5 $dir/pfield.h5 -n 12 $cargs
