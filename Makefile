# MAKEFILE for C-program SCMF

# OPTIONS

NAME     	= TRY

HEADER   	= spline.h string.h

LDFLAGS  	= -lm  -lz -lfftw3 #-static 

INCLUDES 	= 
CC       	=  mpicc    				# standalone (no MPI)
LIBS_MPI 	= 
CADD     	= -O3  -Wall -pedantic

SCMF = spline.o string.o init.o solve.o
scmf:   $(SCMF)
	$(CC) -o $(NAME) $(CFLAGS) $(SCMF) $(LDFLAGS)
spline.o          : spline.c $(HEADER)
string.o          : string.c $(HEADER)
init.o		  : init.c   $(HEADER)
solve.o		  : solve.c  $(HEADER)

#----------------------------------------------------------------------------

new          :
	touch *.c

clean        :
	@rm *.o

#-----------------------------

