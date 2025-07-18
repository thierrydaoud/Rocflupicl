ENAME            = PICL_Test
FC               = mpif77
FFLAGS          +=  -z muldefs -w -std=legacy -g -fmax-errors=1

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=../ppiclF/source
NEW_SOURCE_ROOT_PPICLF= ../ppiclF/source
PPICLF_IFLAGS+=-I $(SOURCE_ROOT_PPICLF)
USR_LFLAGS+=-L $(NEW_SOURCE_ROOT_PPICLF) -lppiclF
TEST_SOURCE=source

SRC = $(TEST_SOURCE)/test_main.f    \
      $(TEST_SOURCE)/test_grid.f      \
      $(TEST_SOURCE)/test_particles.f 

 OBJ = $(TEST_SOURCE)/test_main.o   \
      $(TEST_SOURCE)/test_grid.o      \
      $(TEST_SOURCE)/test_particles.o 

# Make commands
default: getObjs linkObjs

linkObjs:  $(OBJ)
	$(FC) $(FFLAGS) -o $(ENAME) $(OBJ) $(USR_LFLAGS)
	@echo "                         "
	@echo "*************************"
	@echo "*** USER LINK SUCCESS ***"
	@echo "*************************"
	@echo "                         "

getObjs: $(SRC)
	$(FC) $(FFLAGS) -c $(SRC) $(PPICLF_IFLAGS)
	mv *.o $(TEST_SOURCE) 
	@echo "                           "
	@echo "***************************"
	@echo "*** USER OBJECT SUCCESS ***"
	@echo "***************************"
	@echo "                           "

# Cleaning stuff
clean:
	rm $(ENAME)
	rm $(OBJ)
	rm run.err
	rm ppiclF_test.txt
	rm core*
