FC               = mpif77 # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -cpp -w -O3 -std=legacy -g -fmax-errors=1#-fallow-argument-mismatch -fcheck=all
#FFLAGS          += -cpp -w -std=legacy -g #-fcheck=all
CPFLAGS         += -E
#FFLAGS          += -DPPICLC
 	           
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
      $(SOURCE_ROOT_USER)/ppiclf_user_MapProjPart.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_InitZero.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_debug.f     \
      $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.f   \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.f     \
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
      $(SOURCE_ROOT_USER)/ppiclf_user_MapProjPart.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_InitZero.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_debug.o     \
      $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.o   \
      $(SOURCE_ROOT_PPICLF)/ppiclf_comm.o     \
      $(SOURCE_ROOT_PPICLF)/ppiclf_op.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_io.o       \
      $(SOURCE_ROOT_PPICLF)/ppiclf_solve.o    \
      $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.o
PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h      \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h
POS = $(SOURCE_ROOT_PPICLF)/PPICLF
      

# Make commands
default: makeThird preProcess getObjs libObjs 

libObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.o
	@ar crv $(SOURCE_ROOT_PPICLF)/libppiclF.a $(SOURCE_ROOT_PPICLF)/ppiclf.o $(SOURCE_ROOT_GSLIB_OBJ)/*.o
	@echo "                       "
	@echo "***********************"
	@echo "*** LIBRARY SUCCESS ***"
	@echo "***********************"
	@echo "                       "

getObjs: $(SOURCE_ROOT_PPICLF)/ppiclf.f
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
	 rm -r $(SOURCE_ROOT_PPICLF)/*.o                        
	 rm -r $(SOURCE_ROOT_PPICLF)/ppiclf.f
	 rm -r $(SOURCE_ROOT_PPICLF)/PPICLF
	 rm -r $(SOURCE_ROOT_PPICLF)/libppiclF.a                
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

