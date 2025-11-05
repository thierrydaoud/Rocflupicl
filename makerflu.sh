module purge
module load gcc/14.2.0 openmpi/5.0.7
rm build_lib/*.f90 
rm build_lib/*.d   
rm build_lib/*.o   
make clean
make RFLU=1 PICL=1 FOLDER=1 SPEC=1 -j16
ls --color=auto
