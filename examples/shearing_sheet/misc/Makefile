PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)

PREDEF+= -D$(shell basename `readlink ../../../src/gravity.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink ../../../src/boundaries.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF+= -D$(shell basename `readlink ../../../src/collisions.c` '.c' | tr '[a-z]' '[A-Z]')
# PREDEF+= -D$(shell basename `readlink ../../../src/integrator.c` '.c' | tr '[a-z]' '[A-Z]')
PREDEF += -DINTEGRATOR_SEI	# always shearing sheet

CC = gcc
OPT += -std=c99 
OBJECTS_QVALUE += qvalue.o ../../../src/input_params.o ../p2d/readfile.o

all:
	$(CC) $(PREDEF) -c qvalue.c
	$(CC) $(OBJECTS_QVALUE) -o qvalue.out

clean:
	rm -vf *.o *.out
