FC               = mpifort
FFLAGS          += 	-cpp           \
 	           		-Wall          \
             	   	-fbacktrace    

CC               = mpicc  
CPFLAGS         += -E

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF= ../ppiclF/source
SOURCE_ROOT_USER=$(SOURCE_ROOT_PPICLF)/../../user_files
BUILD_ROOT_PPICLF = ../ppiclf/build
USR_LFLAGS+=-L$(BUILD_ROOT_PPICLF) -lppiclF

PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)/include -I$(SOURCE_ROOT_USER) -I$(SOURCE_ROOT_USER)/include 
PPICLF_FFLAGS+=$(PPICLF_IFLAGS) -cpp -Wall -ffree-line-length-none -Wno-missing-include-dirs -I$(BUILD_ROOT_PPICLF) -fmax-errors=20

FFLAGS+= $(PPICLF_IFLAGS)

SRC1 = unitTestSrc/UT_setup.f90
SRC2 = unitTestSrc/UT_CreateBin.f
SRC3 = unitTestSrc/UT_InterpProj.f
SRC4 = unitTestSrc/UT_P2PNN.f
SRC5 = unitTestSrc/UT_CodePerformance.f90

OBJ0 = unitTestSrc/UT_data.o
OBJ1 = unitTestSrc/UT_setup.o 
OBJ2 = unitTestSrc/UT_CreateBin.o
OBJ3 = unitTestSrc/UT_InterpProj.o
OBJ4 = unitTestSrc/UT_P2PNN.o
OBJ5 = unitTestSrc/UT_CodePerformance.o

ENAME1 = CreateBin_UnitTest
ENAME2 = InterpProj_UnitTest
ENAME3 = Particle2Particle_UnitTest
ENAME4 = CodePerformance_UnitTest

# Make commands
default: getObjs linkObjs

linkObjs:  $(OBJ1)
# 	$(FC) $(FFLAGS) -o $(ENAME1) $(OBJ1) $(OBJ2) $(USR_LFLAGS)
# 	$(FC) $(FFLAGS) -o $(ENAME2) $(OBJ1) $(OBJ3) $(USR_LFLAGS)
# 	$(FC) $(FFLAGS) -o $(ENAME3) $(OBJ1) $(OBJ4) $(USR_LFLAGS)
	$(FC) $(FFLAGS) -o $(ENAME4) $(OBJ0) $(OBJ1) $(OBJ5) $(USR_LFLAGS)
	@echo "                         "
	@echo "*************************"
	@echo "*** USER LINK SUCCESS ***"
	@echo "*************************"
	@echo "                         "

getObjs: $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(SRC5)
	$(FC) $(PPICLF_FFLAGS) -c unitTestSrc/UT_data.f90
	$(FC) $(PPICLF_FFLAGS) -c $(SRC1) 
# 	$(FC) $(FFLAGS) -c $(SRC2) 
# 	$(FC) $(FFLAGS) -c $(SRC3) 
# 	$(FC) $(FFLAGS) -c $(SRC4) 
	$(FC) $(PPICLF_FFLAGS) -c $(SRC5)
	@echo "                           "
	@echo "***************************"
	@echo "*** USER OBJECT SUCCESS ***"
	@echo "***************************"
	@echo "                           "
	@mv *.o unitTestSrc/

# Cleaning stuff
clean:
	@rm $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5)
	@rm $(ENAME1) $(ENAME2) $(ENAME3) $(ENAME4)
	@rm $(POS)
