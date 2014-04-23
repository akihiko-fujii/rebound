PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)

PREDEF+= -D$(shell basename `readlink ../../../src/gravity.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink ../../../src/boundaries.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink ../../../src/collisions.c` '.c' | tr '[a-z]' '[A-Z]')
# PREDEF+= -D$(shell basename `readlink ../../../src/integrator.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF += -DINTEGRATOR_SEI	# always shearing sheet

CC = gcc
OPT += -std=c99 
OBJECTS_P2D += particle2density.o readfile.o ../../../src/input_params.o
OBJECTS_ACR += autocorrelation.o readfile.o ../../../src/input_params.o
OBJECTS_SPR += showparams.o readfile.o ../../../src/input_params.o

all:
	$(CC) $(PREDEF) readfile.c -c
	$(CC) $(PREDEF) particle2density.c -c
	$(CC) $(PREDEF) autocorrelation.c -c
	$(CC) $(PREDEF) showparams.c -c
	$(CC) $(OBJECTS_P2D) -o particle2density.out
	$(CC) $(OBJECTS_ACR) -o autocorrelation.out
	$(CC) $(OBJECTS_SPR) -o showparams.out
clean:
	rm -vf *.o *.out
