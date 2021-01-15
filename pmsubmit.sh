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

# Load the pyQD conda environment
conda activate pmeep

export TMPDIR=/beegfs/<your-user-name-here>/tmp
# INPUTDIR=/beegfs/lu27wil/scientific-pyqd/science/inputs/2020-12-11b

mkdir -p $TMPDIR
mkdir -p ./$SLURM_JOB_ID

srun -o %j/plmeep-%j-%2s.out mpirun -np 16 plas-meep -r $1 &

wait