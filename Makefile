SERIAL_CC= gcc
MPI_CC= mpicc
WOPT_DEFAULT= -Wall -Werror
ADD_OMP= yes
ADD_MPI= yes
DEBUG_VERSION= no
DEBUG_SINGLEPIX= no
LIB_GSL= -L
INC_GSL= -I
LIB_HP= -L/home/bthorne/Codes/Healpix_3.31/lib
INC_HP= -I/home/bthorne/Codes/Healpix_3.31/include/chealpix.h
LIB_FITS= 
INC_FITS= 

WOPT=$(WOPT_DEFAULT) -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF
ifeq ($(ADD_OMP),yes)
	WOPT+= -fopenmp -D_WITH_OMP
endif
ifeq ($(ADD_MPI),yes)
	WOPT+= -D_WITH_MPI
	CC= $(MPI_CC)
else
	CC= $(SERIAL_CC)
endif
ifeq ($(DEBUG_VERSION),yes)
	WOPT+= -D_DEBUG
ifeq ($(DEBUG_SINGLEPIX),yes)
	WOPT+= -D_DEBUG_SINGLEPIX
endif
endif

CFLAGS= $(WOPT) -I./src $(INC_GSL) $(INC_HP) $(INC_FITS)
LIBS= $(LIB_GSL) $(LIB_HP) $(LIB_FITS) $(LIB_SHARP)
LIBS+= -lgsl -lgslcblas -lchealpix -lcfitsio -lm

COMMONO= src/common.o
HEO= src/healpix_extra.o
RNGO= src/rng.o
POWELLO= src/powell.o
BFOREO= src/bfore.o
MAINO= src/main.o
OBJ= $(COMMONO) $(HEO) $(RNGO) $(POWELLO) $(BFOREO) $(MAINO)

EXEC= BFoRe
all: $(EXEC)

BFoRe : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o $@

clean :
	rm -f $(EXEC) $(OBJ)
cleaner :
	rm -f $(EXEC) $(OBJ) *~ src/*~
