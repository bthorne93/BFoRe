CC= gcc
WOPT= -Wall -fopenmp -O3
#WOPT+= -D_DEBUG
LIB_GSL= -L/home/damonge/lib
INC_GSL= -I/home/damonge/include
LIB_HP=
INC_HP=
LIB_FITS=
INC_FITS=
LIB_SHARP=
INC_SHARP=

CFLAGS= $(WOPT) -I./src $(INC_GSL) $(INC_HP) $(INC_FITS) $(INC_SHARP)
LIBS= $(LIB_GSL) $(LIB_HP) $(LIB_FITS) $(LIB_SHARP)
LIBS+= -lgsl -lgslcblas -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -lm

COMMONO= src/common.o
HEO= src/healpix_extra.o
RNGO= src/rng.o
BFOREO= src/bfore.o
MAINO= src/main.o
OBJ= $(COMMONO) $(HEO) $(RNGO) $(BFOREO) $(MAINO)

EXEC= BFoRe
all: $(EXEC)

BFoRe : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LIBS) -o $@

clean :
	rm -f $(EXEC) $(OBJ)
cleaner :
	rm -f $(EXEC) $(OBJ) *~ src/*~
