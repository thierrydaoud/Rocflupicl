#!/bin/bash

# this causes bash to actually exit the script when a command results in an error, instead of just continuing for no reason
set -e

# module purge
# module load python/3.12 gcc/14.2.0 openmpi/5.0.7
# module list
cd libpicl
rm -f ppiclF/source/ppiclf.f
rm -f ppiclF/source/PPICLF_USER.h
rm -f ppiclF/source/PPICLF_USER_COMMON.h
cd ppiclF/source
ln -s ../../user_files/PPICLF_USER.h .
ln -s ../../user_files/PPICLF_USER_COMMON.h .
cd ../../
make clean
make
cd .. 
make clean 

rm -f build_lib/*.f90
rm -f build_lib/*.d
rm -f build_lib/*.o

rm build_util/*/*.f90
rm build_util/*/*.d
rm build_util/*/*.o

make RFLU=1 PICL=1 SPEC=1 FOLDER=1 -j16
ls --color=auto
