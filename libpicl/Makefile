FC               = mpif90 # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -w -cpp -O3 -std=legacy -g -fcheck=bounds -fbacktrace -fno-math-errno -funroll-loops -mcmodel=medium
# Add for Turin
# FFLAGS += -march=znver5
# Add for Milan 
FFLAGS += -march=znver3 

CPFLAGS         += -E
ifeq ($(TEST),1)
    FFLAGS += -DTEST
endif

ifeq ($(PERF),1)
    FFLAGS += -DPERF
endif
 	           
INSTALL_LOCATION = .

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(INSTALL_LOCATION)/ppiclF/source
SOURCE_ROOT_USER=$(INSTALL_LOCATION)/user_files
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)

# GSLIB LIBRARY
SOURCE_ROOT_GSLIB=$(INSTALL_LOCATION)/ppiclF/3rd_party/gslib
SOURCE_ROOT_GSLIB_OBJ=$(SOURCE_ROOT_GSLIB)/gslib/src
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include


DYN_SRC = $(SOURCE_ROOT_PPICLF)/ppiclf_DynamicAllocation.F90
DYN_OBJ = $(SOURCE_ROOT_PPICLF)/ppiclf_DynamicAllocation.o

SRC = $(SOURCE_ROOT_USER)/ppiclf_user.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.f       \
      $(SOURCE_ROOT_USER)/ppiclf_user_unitproblems.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_QuasiSteady.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Fluctuations.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_AddedMass.f     \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.f   \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_subroutines.f   \
      $(SOURCE_ROOT_USER)/ppiclf_user_HeatTransfer.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Torque.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Lift.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_ViscousUnsteady.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_InitZero.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_debug.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.f   \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_PBloading.f     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.f       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.f    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.f      
OBJ = $(SOURCE_ROOT_USER)/ppiclf_user.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.o       \
      $(SOURCE_ROOT_USER)/ppiclf_user_unitproblems.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_QuasiSteady.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Fluctuations.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_AddedMass.o     \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.o   \
      $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_subroutines.o   \
      $(SOURCE_ROOT_USER)/ppiclf_user_HeatTransfer.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Torque.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_Lift.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_ViscousUnsteady.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_InitZero.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_debug.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.o   \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_PBloading.o     \
     $(SOURCE_ROOT_PPICLF)/ppiclf_op.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.o    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.o      
PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h         \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h    \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PERFORMANCE.h
POS = $(SOURCE_ROOT_PPICLF)/PPICLF
      

# Make commands
default: makeThird preProcess $(DYN_OBJ) getObjs libObjs 

$(DYN_OBJ): $(DYN_SRC)
	$(FC) -c $(DYN_SRC) -J$(SOURCE_ROOT_PPICLF) -o $(DYN_OBJ) 

libObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.o $(DYN_OBJ)
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(DYN_OBJ) $(SOURCE_ROOT_PPICLF)/ppiclf.o $(SOURCE_ROOT_GSLIB_OBJ)/*.o
	@echo "                       "
	@echo "***********************"
	@echo "*** LIBRARY SUCCESS ***"
	@echo "***********************"
	@echo "                       "

getObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.f $(DYN_OBJ)
	$(FC) $(FFLAGS) -c $(SOURCE_ROOT_PPICLF)/ppiclf.f $(GSLIB_IFLAGS) 
	mv *.o $(SOURCE_ROOT_PPICLF)
	@echo "                              "
	@echo "******************************"
	@echo "*** LIBRARY OBJECT SUCCESS ***"
	@echo "******************************"
	@echo "                              "

preProcess: $(PRE)
	$(CC) $(CPFLAGS) $(PRE) > $(POS)_tmp
	grep "^[^#]" $(POS)_tmp > $(POS)
	rm $(POS)_tmp
	echo '#include "PPICLF_USER.h"' > $(SOURCE_ROOT_PPICLF)/ppiclf.f
	echo '#include "PPICLF_STD.h"' >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	cat $(SRC)  >> $(SOURCE_ROOT_PPICLF)/ppiclf.f
	@echo "                       "
	@echo "********************"
	@echo "*** PREPROCESSED ***"
	@echo "********************"
	@echo "                       "
	

makeThird: $(SOURCE_ROOT_GSLIB)/install
	cd $(SOURCE_ROOT_GSLIB); \
	./install $(CC) $(FC)
	@echo "                       "
	@echo "***********************"
	@echo "*** INSTALLED GSLIB ***"
	@echo "***********************"
	@echo "                       "

clean:
	 rm -r $(SOURCE_ROOT_PPICLF)/libppiclF.a                
	 rm -r $(SOURCE_ROOT_PPICLF)/*.o                        
	 rm -r $(SOURCE_ROOT_PPICLF)/ppiclf_DynamicAllocation.mod
	 rm -r $(SOURCE_ROOT_PPICLF)/ppiclf.f
	 rm -r $(SOURCE_ROOT_PPICLF)/PPICLF
	 rm -r $(SOURCE_ROOT_PPICLF)/PPICLF_tmp/*                
	 rm -r $(SOURCE_ROOT_GSLIB)/gslib                       
	 rm -r $(SOURCE_ROOT_GSLIB)/lib                         
	 rm -r $(SOURCE_ROOT_GSLIB)/include                     
	 rm -r $(SOURCE_ROOT_GSLIB)/*.tar.gz                    
	 rm -r $(INSTALL_LOCATION)/short_tests/*.pyc            
	 rm -r $(INSTALL_LOCATION)/short_tests/lib/*.pyc        
	 rm -r $(INSTALL_LOCATION)/short_tests/test1/test1.log*
	 rm -r $(INSTALL_LOCATION)/short_tests/test1/*.vtu
	 cd $(INSTALL_LOCATION)/short_tests/test1/; \
	make clean

