FC               = mpif77 # Valid MPI Fortran compiler
CC               = mpicc  # Valid MPI C compiler
FFLAGS          += -cpp -w -O3 -std=legacy -g -fmax-errors=1#-fallow-argument-mismatch -fcheck=all
CPFLAGS         += -E
ifeq ($(TEST),1)
    FFLAGS += -DTEST
endif

ifeq ($(PERF),1)
    FFLAGS += -DPERF
endif

PREPROCDEFINES+=-DNDIM=3
 	           
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

SRCUSR =    $(SOURCE_ROOT_USER)/ppiclf_user_MapProjPart.f90    \
			$(SOURCE_ROOT_USER)/ppiclf_m_user_data.f90 \
			$(SOURCE_ROOT_USER)/ppiclf_m_user_RFLUdata.f90 \
			$(SOURCE_ROOT_USER)/ppiclf_user_SetYdot.fpp		\
            $(SOURCE_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_InitZero.fpp     \
			$(SOURCE_ROOT_USER)/ppiclf_user_NearestNeighbor.f90 \
     		$(SOURCE_ROOT_USER)/ppiclf_user_subbinMap.fpp       \
			$(SOURCE_ROOT_USER)/ppiclf_m_user_ForceModels.fpp \
            $(SOURCE_ROOT_USER)/ppiclf_user_BurnRateModel.fpp   \
            $(SOURCE_ROOT_USER)/ppiclf_user_Lift.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_Torque.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_HeatTransfer.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_Fluctuations.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_QuasiSteady.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_ViscousUnsteady.fpp     \
            $(SOURCE_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.fpp   \
            $(SOURCE_ROOT_USER)/ppiclf_user_AddedMass.fpp     \
			
UNUSED=     $(SOURCE_ROOT_USER)/ppiclf_user_unitproblems.fpp     \
            $(SOURCE_ROOT_USER)/ppiclf_user_debug.fpp     \

SRC =       $(SOURCE_ROOT_PPICLF)/ppiclf_data.fpp \
            $(SOURCE_ROOT_PPICLF)/ppiclf_op.fpp       \
            $(SOURCE_ROOT_PPICLF)/ppiclf_m_comm.fpp     \
            $(SOURCE_ROOT_PPICLF)/ppiclf_io.fpp       \
			$(SOURCE_ROOT_PPICLF)/ppiclf_initsolve.fpp       \
            $(SOURCE_ROOT_PPICLF)/ppiclf_solve.fpp    \
			$(SOURCE_ROOT_PPICLF)/ppiclf_user.fpp	\
            $(SOURCE_ROOT_PPICLF)/ppiclf_mxm.fpp	\
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_types.fpp	\
			$(SOURCE_ROOT_PPICLF)/ppiclf_m_particledata.fpp	

OBJ = 	$(patsubst $(SOURCE_ROOT_PPICLF)/%.fpp,$(BUILD_ROOT_PPICLF)/%.o,$(SRC))
OBJUSR =  	$(patsubst $(SOURCE_ROOT_USER)/%,$(BUILD_ROOT_USER)/%,$(patsubst %.f90,%.o,$(patsubst %.fpp,%.o,$(SRCUSR))))

PRE = $(SOURCE_ROOT_PPICLF)/PPICLF_SOLN.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GRID.h     \
      $(SOURCE_ROOT_PPICLF)/PPICLF_OPT.h      \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PARALLEL.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_USER_COMMON.h \
      $(SOURCE_ROOT_PPICLF)/PPICLF_GEOM.h        \
      $(SOURCE_ROOT_PPICLF)/PPICLF_PERFORMANCE.h
POS = $(SOURCE_ROOT_PPICLF)/PPICLF


default: libppicl
test: $(BUILD_ROOT_PPICLF)/ppiclf_m_particledata.f90
# default: $(BUILD_ROOT_USER)/ppiclf_user_imp.o

libppicl: $(BUILD_ROOT_PPICLF)/libppiclF.a

libppiclcore: $(BUILD_ROOT_PPICLF)/libppiclcore.a

# core library build rule
$(BUILD_ROOT_PPICLF)/libppiclcore.a: $(OBJ) makeThird
	ar crv $(BUILD_ROOT_PPICLF)/libppiclcore.a $(OBJ) $(SOURCE_ROOT_GSLIB_OBJ)/*.o

$(BUILD_ROOT_PPICLF)/testfile.f90: $(SOURCE_ROOT_PPICLF)/include/PPICLF_PARTMACROS.fypp

# module dependencies
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_io.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o
$(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o



$(BUILD_ROOT_PPICLF)/ppiclf_op.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o

$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o
$(BUILD_ROOT_PPICLF)/ppiclf_io.o: $(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o

$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_data.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_comm.o
$(BUILD_ROOT_PPICLF)/ppiclf_solve.o: $(BUILD_ROOT_PPICLF)/ppiclf_op.o

$(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o: $(BUILD_ROOT_PPICLF)/ppiclf_user.o
$(BUILD_ROOT_PPICLF)/ppiclf_initsolve.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_particledata.o

$(BUILD_ROOT_PPICLF)/ppiclf_m_particledata.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_types.o

$(BUILD_ROOT_PPICLF)/ppiclf_user.o: $(BUILD_ROOT_PPICLF)/ppiclf_m_types.o

# user file dependencies 
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_m_user_ForceModels.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_user_subbinMap.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_SetYdot.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o: $(BUILD_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.o
$(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_EvalNearestNeighbor.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_Torque.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_Torque.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_Lift.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_Lift.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_HeatTransfer.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_HeatTransfer.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_Fluctuations.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_Fluctuations.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_ViscousUnsteady.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_ViscousUnsteady.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_AddedMass.o: $(BUILD_ROOT_USER)/ppiclf_m_user_RFLUdata.o
$(BUILD_ROOT_USER)/ppiclf_user_AddedMass.o: $(BUILD_ROOT_USER)/ppiclf_m_user_data.o

$(BUILD_ROOT_USER)/ppiclf_user_AddedMass.o: $(BUILD_ROOT_USER)/Added_Mass/ppiclf_user_AddedMass_functions.o


$(BUILD_ROOT_PPICLF)/%.f90 : $(SOURCE_ROOT_PPICLF)/%.fpp
	@mkdir -vp $(BUILD_ROOT_PPICLF)
	@fypp -I $(SOURCE_ROOT_USER) -I $(SOURCE_ROOT_PPICLF)/include $< $@ -m fyppmacros -M $(INSTALL_LOCATION)/ppiclF/tools/ ${PREPROCDEFINES} -n

$(BUILD_ROOT_PPICLF)/%.o : $(BUILD_ROOT_PPICLF)/%.f90
	@mkdir -vp $(BUILD_ROOT_PPICLF)
	@$(FC) $< -c -o $@ $(PPICLF_FFLAGS)
# mpifort source/ppiclf_data.f90 -cpp -o ppiclf_data.o -Isource/include -I../user_files -ffree-line-length-none





# full library build rule, adds the user files to the core archive generated above
# this separation allows for the core library to be build without the user files (except the .h file)
$(BUILD_ROOT_PPICLF)/libppiclF.a: $(BUILD_ROOT_PPICLF)/libppiclcore.a $(OBJUSR)
	cp $(BUILD_ROOT_PPICLF)/libppiclcore.a $(BUILD_ROOT_PPICLF)/libppiclF.a
	ar rv $(BUILD_ROOT_PPICLF)/libppiclF.a $(OBJUSR)

$(BUILD_ROOT_USER)/%.o : $(SOURCE_ROOT_USER)/%.f90
	mkdir -vp $(BUILD_ROOT_USER)
	$(FC) $< -c -o $@ $(PPICLF_USER_FFLAGS)

$(BUILD_ROOT_USER)/%.f90 : $(SOURCE_ROOT_USER)/%.fpp
	@mkdir -vp $(BUILD_ROOT_USER)
	@fypp -I $(SOURCE_ROOT_USER) -I $(SOURCE_ROOT_PPICLF)/include $< $@ -m fyppmacros -M $(INSTALL_LOCATION)/ppiclF/tools/ ${PREPROCDEFINES} -n

$(BUILD_ROOT_PPICLF)/%.o : $(BUILD_ROOT_USER)/%.f90
	@mkdir -vp $(BUILD_ROOT_USER)
	@$(FC) $< -c -o $@ $(PPICLF_FFLAGS)






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
	@rm -rf $(SOURCE_ROOT_GSLIB)/gslib                       
	@rm -rf $(SOURCE_ROOT_GSLIB)/lib                         
	@rm -rf $(SOURCE_ROOT_GSLIB)/include                     
	@rm -rf $(SOURCE_ROOT_GSLIB)/*.tar.gz                    
	@rm -rf $(INSTALL_LOCATION)/short_tests/*.pyc            
	@rm -rf $(INSTALL_LOCATION)/short_tests/lib/*.pyc        
	@rm -rf $(INSTALL_LOCATION)/short_tests/test1/test1.log*
	@rm -rf $(INSTALL_LOCATION)/short_tests/test1/*.vtu
	@rm -rf $(BUILD_ROOT_PPICLF)
	@mkdir $(BUILD_ROOT_PPICLF)
	@mkdir $(BUILD_ROOT_USER)
	@mkdir $(BUILD_ROOT_USER)/Added_Mass
	@cd $(INSTALL_LOCATION)/ppiclF/short_tests/test1/; \
	make clean


.SECONDARY: $(patsubst $(SOURCE_ROOT_PPICLF)/%.fpp,$(BUILD_ROOT_PPICLF)/%.f90,$(SRC)) $(patsubst $(SOURCE_ROOT_USER)/%.fpp,$(BUILD_ROOT_USER)/%.f90,$(SRCUSR))

FORCE: