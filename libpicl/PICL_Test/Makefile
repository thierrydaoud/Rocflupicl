ENAME            = PICL_Test
FC               = mpif77
FFLAGS          +=  -z muldefs -w -std=legacy -g -fmax-errors=1

####################
# DO NOT TOUCH BELOW
####################

# PPICLF LIBRARY
SOURCE_ROOT_PPICLF=../ppiclF/source
PPICLF_IFLAGS+=-I $(SOURCE_ROOT_PPICLF)
USR_LFLAGS+=-L $(SOURCE_ROOT_PPICLF) -lppiclF

SRC = main.f
OBJ = main.o

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

# Cleaning stuff
clean:
	@rm $(ENAME) *.o run*
