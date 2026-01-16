FC               = mpifort # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -cpp -w -O3 -std=legacy -g -fmax-errors=1#-fallow-argument-mismatch -fcheck=all
CPFLAGS         += -E
ifeq ($(TEST),1)
FFLAGS += -DTEST
PPICLF_FFLAGS += -DTEST
endif
 	           
INSTALL_LOCATION = .

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=$(INSTALL_LOCATION)/ppiclF/source
SOURCE_ROOT_USER=$(INSTALL_LOCATION)/user_files
BUILD_ROOT_PPICLF=$(INSTALL_LOCATION)/ppiclF/build
BUILD_ROOT_USER=$(BUILD_ROOT_PPICLF)/user
PPICLF_IFLAGS+=-I$(SOURCE_ROOT_PPICLF)/include -I$(SOURCE_ROOT_USER) -I$(SOURCE_ROOT_USER)/include 
PPICLF_FFLAGS+=$(PPICLF_IFLAGS) -cpp -Wall -ffree-line-length-none -Wno-missing-include-dirs -fmax-errors=20
PPICLF_USER_FFLAGS:=$(PPICLF_FFLAGS) -J$(BUILD_ROOT_USER) -I$(BUILD_ROOT_PPICLF)
PPICLF_FFLAGS+=-J$(BUILD_ROOT_PPICLF) 

# GSLIB LIBRARY
SOURCE_ROOT_GSLIB=$(INSTALL_LOCATION)/ppiclF/3rd_party/gslib
SOURCE_ROOT_GSLIB_OBJ=$(SOURCE_ROOT_GSLIB)/gslib/src
GSLIB_IFLAGS+=-I$(SOURCE_ROOT_GSLIB)/include

SRCUSR =    $(SOURCE_ROOT_USER)/ppiclf_user_imp.f90     
UNUSED=            $(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.f       \
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
            $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.f   

SRC =       $(SOURCE_ROOT_PPICLF)/ppiclf_data.f90 \
            $(SOURCE_ROOT_PPICLF)/ppiclf_op.f90       \
            $(SOURCE_ROOT_PPICLF)/ppiclf_commu.f90     \
            $(SOURCE_ROOT_PPICLF)/ppiclf_io.f90       \
			$(SOURCE_ROOT_PPICLF)/ppiclf_initsolve.f90       \
            $(SOURCE_ROOT_PPICLF)/ppiclf_solve.f90    \
			$(SOURCE_ROOT_PPICLF)/ppiclf_user.f90	\
            $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.f90      

OBJ = 	$(patsubst $(SOURCE_ROOT_PPICLF)/%.f90,$(BUILD_ROOT_PPICLF)/%.o,$(SRC))
OBJUSR =  	$(patsubst $(SOURCE_ROOT_USER)/%,$(BUILD_ROOT_USER)/%,$(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SRCUSR))))

PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h      \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h
POS = $(SOURCE_ROOT_PPICLF)/PPICLF
      
default: libppicl
# default: $(BUILD_ROOT_USER)/ppiclf_user_imp.o

libppicl: $(BUILD_ROOT_PPICLF)/libppiclF.a

libppiclcore: $(BUILD_ROOT_PPICLF)/libppiclcore.a

# core library build rule
$(BUILD_ROOT_PPICLF)/libppiclcore.a: $(OBJ) makeThird
	ar crv $(BUILD_ROOT_PPICLF)/libppiclcore.a $(OBJ) $(SOURCE_ROOT_GSLIB_OBJ)/*.o



# module dependencies
$(BUILD_ROOT_PPICLF)/ppiclf_commu.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_commu.o: $(BUILD_ROOT_PPICLF)/ppiclf_io.o
$(BUILD_ROOT_PPICLF)/ppiclf_commu.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o
$(BUILD_ROOT_PPICLF)/ppiclf_commu.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o



$(BUILD_ROOT_PPICLF)/ppiclf_op.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o

$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o
$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o

$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_commu.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o

$(BUILD_ROOT_PPICLF)/%.o : $(SOURCE_ROOT_PPICLF)/%.f90
	$(FC) $< -c -o $@ $(PPICLF_FFLAGS)
# mpifort source/ppiclf_data.f90 -cpp -o ppiclf_data.o -Isource/include -I../user_files -ffree-line-length-none





# full library build rule, adds the user files to the core archive generated above
# this separation allows for the core library to be build without the user files (except the .h file)
$(BUILD_ROOT_PPICLF)/libppiclF.a: $(BUILD_ROOT_PPICLF)/libppiclcore.a $(OBJUSR)
	cp $(BUILD_ROOT_PPICLF)/libppiclcore.a $(BUILD_ROOT_PPICLF)/libppiclF.a
	ar rv $(BUILD_ROOT_PPICLF)/libppiclF.a $(OBJUSR)

$(BUILD_ROOT_USER)/%.o : $(SOURCE_ROOT_USER)/%.f90
	$(FC) $< -c -o $@ $(PPICLF_USER_FFLAGS)







# Old Make commands
defaultold: makeThird preProcess getObjs libObjs 

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
	 -rm -r $(SOURCE_ROOT_PPICLF)/*.o                        
	 -rm -r $(SOURCE_ROOT_PPICLF)/ppiclf.f
	 -rm -r $(SOURCE_ROOT_PPICLF)/PPICLF
	 -rm -r $(SOURCE_ROOT_PPICLF)/libppiclF.a                
	 -rm -r $(SOURCE_ROOT_PPICLF)/PPICLF_tmp/*                
	 -rm -r $(SOURCE_ROOT_GSLIB)/gslib                       
	 -rm -r $(SOURCE_ROOT_GSLIB)/lib                         
	 -rm -r $(SOURCE_ROOT_GSLIB)/include                     
	 -rm -r $(SOURCE_ROOT_GSLIB)/*.tar.gz                    
	 -rm -r $(INSTALL_LOCATION)/short_tests/*.pyc            
	 -rm -r $(INSTALL_LOCATION)/short_tests/lib/*.pyc        
	 -rm -r $(INSTALL_LOCATION)/short_tests/test1/test1.log*
	 -rm -r $(INSTALL_LOCATION)/short_tests/test1/*.vtu
	 -rm -r $(BUILD_ROOT_PPICLF)
	 mkdir $(BUILD_ROOT_PPICLF)
	 mkdir $(BUILD_ROOT_USER)
	 cd $(INSTALL_LOCATION)/ppiclF/short_tests/test1/; \
	make clean

