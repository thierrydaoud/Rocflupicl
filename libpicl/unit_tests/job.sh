#!/bin/bash
#SBATCH --job-name=codetest
#SBATCH --output=ppiclfTestResults.txt
#SBATCH --mail-type=none
#SBATCH --ntasks=36
#SBATCH --nodes=2
#SBATCH --mem-per-cpu=20gb                     # Per node
#SBATCH --time=00-00:30:00             # Walltime in hh:mm:ss or d-hh:mm:ss
#SBATCH --account=bala1s
#SBATCH --qos=bala1s
date

module purge 
module load gcc/14.2.0
module load openmpi/5.0.7

EXEC=PICL_Test

ntasks=$SLURM_NTASKS
nodes=$SLURM_JOB_NUM_NODES

echo "Current Directory =" `pwd`

rm PICL_Test
cd ..
make clean
make TEST=1
cd PICL_Test/
make
#srun -N 1 -n 1  --mpi=${HPC_PMIX} $EXEC -v 2
srun -N 1 -n 8  --mpi=${HPC_PMIX} $EXEC -v 2
#srun -N 1 -n 16  --mpi=${HPC_PMIX} $EXEC -v 2
#srun -N 2 -n 27  --mpi=${HPC_PMIX} $EXEC -v 2
srun -N 2 -n 32  --mpi=${HPC_PMIX} $EXEC -v 2
#srun -N 2 -n 36  --mpi=${HPC_PMIX} $EXEC -v 2
#srun -N 2 -n 48  --mpi=${HPC_PMIX} $EXEC -v 2
#srun -N 2 -n 64 --mpi=${HPC_PMIX} $EXEC -v 2
#srun -N 4 -n 128 --mpi=${HPC_PMIX} $EXEC -v 2
