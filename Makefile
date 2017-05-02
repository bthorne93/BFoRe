SERIAL_CC= icc
MPI_CC= mpicc
WOPT_DEFAULT= -Wall -O3
ADD_OMP= yes
ADD_MPI= yes
ADD_NILC= yes
DEBUG_VERSION= no
DEBUG_SINGLEPIX= no
LIB_GSL= -L/home/damonge/lib
INC_GSL= -I/home/damonge/include
LIB_HP=
INC_HP=
LIB_FITS=
INC_FITS=
LIB_SHARP=
INC_SHARP=

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

INCS= -I./src $(INC_GSL) $(INC_HP) $(INC_FITS)
LIBDIRS= $(LIB_GSL) $(LIB_HP) $(LIB_FITS) $(LIB_SHARP)
LIBS+= -lgsl -lgslcblas -lchealpix -lcfitsio -lm
ifeq ($(ADD_NILC),yes)
	WOPT+= -D_WITH_SHT -D_WITH_NEEDLET
	INCS+= $(INC_SHARP)
	LIBDIRS+= $(LIB_SHARP)
	LIBS+= -lsharp -lfftpack -lc_utils
endif
CFLAGS= $(WOPT)
CFLAGS+= $(INCS)

COMMONO= src/common.o
HEO= src/healpix_extra.o
RNGO= src/rng.o
POWELLO= src/powell.o
BFOREO= src/bfore.o
MAINO= src/main.o
OBJ= $(COMMONO) $(HEO) $(RNGO) $(POWELLO) $(BFOREO) $(MAINO)

EXEC= BFoRe_db
all: $(EXEC)

BFoRe_db : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBDIRS) $(LIBS) -o $@

clean :
	rm -f $(EXEC) $(OBJ)
cleaner :
	rm -f $(EXEC) $(OBJ) *~ src/*~
