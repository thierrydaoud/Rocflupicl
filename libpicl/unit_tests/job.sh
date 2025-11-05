#!/bin/bash
#SBATCH --job-name=UnitTest
#SBATCH --output=ppiclfTestResults.txt
#SBATCH --error=run.err
##SBATCH --mail-type=none
#SBATCH --ntasks=64
#SBATCH --nodes=4
#SBATCH --mem-per-cpu=20gb                     # Per node
#SBATCH --time=00-02:30:00             # Walltime in hh:mm:ss or d-hh:mm:ss
#SBATCH --account=bala1s
#SBATCH --qos=bala1s
date

module purge 
module load gcc/14.2.0
module load openmpi/5.0.7

EXEC1=CreateBin_UnitTest
EXEC2=Particle2Particle_UnitTest
EXEC3=InterpProj_UnitTest
EXEC4=CodePerformance_UnitTest

ntasks=$SLURM_NTASKS

echo "Current Directory =" `pwd`

cd ..
make clean
make TEST=1
cd unit_tests/
make clean
make

> ppiclfTestResults.txt
> run.err
date
echo "Current Directory =" `pwd`
echo ""
echo "*** Number of Processors Used: 1 ***"

srun -N 1 -n 1  --mpi=${HPC_PMIX} $EXEC1 -v 2
srun -N 1 -n 1  --mpi=${HPC_PMIX} $EXEC2 -v 2
srun -N 1 -n 1  --mpi=${HPC_PMIX} $EXEC3 -v 2
srun -N 1 -n 1  --mpi=${HPC_PMIX} $EXEC4 -v 2

echo ""
for i in {2..64}; do
echo "*** ***  Number of Processors Used: $i  *** ***"

  srun -N 4 -n $i  --mpi=${HPC_PMIX} $EXEC2 -v 2
  srun -N 4 -n $i  --mpi=${HPC_PMIX} $EXEC3 -v 2
  srun -N 4 -n $i  --mpi=${HPC_PMIX} $EXEC4 -v 2

  echo ""
done
