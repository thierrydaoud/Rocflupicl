#!/bin/bash
#SBATCH --job-name=PICLtest
#SBATCH --output=ppiclF_test.txt
#SBATCH --error=run.err
#SBATCH --mail-type=none
#SBATCH --ntasks=27
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=18gb                     # Per node
#SBATCH --time=00-00:02:00             # Walltime in hh:mm:ss or d-hh:mm:ss
#SBATCH --account=bala1s
#SBATCH --qos=bala1s
#SBATCH --partition=hpg-dev
date

module purge 
module load gcc/14.2.0
module load openmpi/5.0.7

EXEC=PICL_Test

ntasks=$SLURM_NTASKS
nodes=$SLURM_JOB_NUM_NODES

echo "Current Directory =" `pwd`

srun -N 1 -n 2 --mpi=${HPC_PMIX} $EXEC -v 2
srun -N 1 -n 27 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 3 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 4 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 5 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 6 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 7 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 8 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 9 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 10 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 11 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 12 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 13 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 14 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 15 --mpi=${HPC_PMIX} $EXEC -v 2
# srun -N 1 -n 16 --mpi=${HPC_PMIX} $EXEC -v 2
