#!/bin/bash
#SBATCH --job-name=PICLtest
#SBATCH --output=ppiclF_test.txt
#SBATCH --error=run.err
#SBATCH --mail-type=none
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8gb                     # Per node
#SBATCH --time=00-00:02:00             # Walltime in hh:mm:ss or d-hh:mm:ss
#SBATCH --account=bala1s
#SBATCH --qos=bala1s
##SBATCH --contiguous
##SBATCH --exclusive
date

module purge 
module load gcc/14.2.0
module load openmpi/5.0.7

EXEC=PICL_Test

ntasks=$SLURM_NTASKS
nodes=$SLURM_JOB_NUM_NODES

echo "Current Directory =" `pwd`

srun -N 1 -n 8 --mpi=${HPC_PMIX} $EXEC -v 2
srun -N 1 -n 4 --mpi=${HPC_PMIX} $EXEC -v 2
