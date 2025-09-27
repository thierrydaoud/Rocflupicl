ENAME            = PICL_Test
FC               = mpif77
FFLAGS          += -cpp           \
 		   -z muldefs     \
 	           -Wall          \
             	   -fbacktrace    \
                   -fallow-argument-mismatch  

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF= ../ppiclF/source
PPICLF_IFLAGS+=-I $(SOURCE_ROOT_PPICLF)
USR_LFLAGS+=-L $(SOURCE_ROOT_PPICLF) -lppiclF

SRC = testCode/picl_test.f    
OBJ = testCode/picl_test.o   

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
	@echo "                           "
	@echo "***************************"
	@echo "*** USER OBJECT SUCCESS ***"
	@echo "***************************"
	@echo "                           "
	@mv *.o testCode
# Cleaning stuff
clean:
	@rm $(ENAME)
	@rm $(OBJ)
