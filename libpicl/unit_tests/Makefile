FC               = mpif77
FFLAGS          += -cpp           \
 		   -z muldefs     \
 	           -Wall          \
             	   -fbacktrace    \
                   -fallow-argument-mismatch  
CC               = mpicc  
CPFLAGS         += -E

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF= ../ppiclF/source
USR_LFLAGS+=-L $(SOURCE_ROOT_PPICLF) -lppiclF
PRE = unitTestSrc/UNIT_TEST_COMMON.h             \
      $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h         \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h    \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h
POS = unitTestSrc/PPICLF_UNIT_TEST
      

SRC1 = unitTestSrc/UT_setup.f 
SRC2 = unitTestSrc/UT_CreateBin.f
SRC3 = unitTestSrc/UT_InterpProj.f
SRC4 = unitTestSrc/UT_P2PNN.f
SRC5 = unitTestSrc/UT_CodePerformance.f

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
default: preProcess getObjs linkObjs

linkObjs:  $(OBJ1)
	$(FC) $(FFLAGS) -o $(ENAME1) $(OBJ1) $(OBJ2) $(USR_LFLAGS)
	$(FC) $(FFLAGS) -o $(ENAME2) $(OBJ1) $(OBJ3) $(USR_LFLAGS)
	$(FC) $(FFLAGS) -o $(ENAME3) $(OBJ1) $(OBJ4) $(USR_LFLAGS)
	$(FC) $(FFLAGS) -o $(ENAME4) $(OBJ1) $(OBJ5) $(USR_LFLAGS)
	@echo "                         "
	@echo "*************************"
	@echo "*** USER LINK SUCCESS ***"
	@echo "*************************"
	@echo "                         "

getObjs: $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(SRC5)
	$(FC) $(FFLAGS) -c $(SRC1) 
	$(FC) $(FFLAGS) -c $(SRC2) 
	$(FC) $(FFLAGS) -c $(SRC3) 
	$(FC) $(FFLAGS) -c $(SRC4) 
	$(FC) $(FFLAGS) -c $(SRC5)
	@echo "                           "
	@echo "***************************"
	@echo "*** USER OBJECT SUCCESS ***"
	@echo "***************************"
	@echo "                           "
	@mv *.o unitTestSrc/

preProcess: $(PRE)
	$(CC) $(CPFLAGS) $(PRE) > $(POS)_tmp
	grep "^[^#]" $(POS)_tmp > $(POS)
	rm $(POS)_tmp
	@echo "                       "
	@echo "********************"
	@echo "*** PREPROCESSED ***"
	@echo "********************"
	@echo "                       "

# Cleaning stuff
clean:
	@rm $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5)
	@rm $(ENAME1) $(ENAME2) $(ENAME3) $(ENAME4)
	@rm $(POS)
