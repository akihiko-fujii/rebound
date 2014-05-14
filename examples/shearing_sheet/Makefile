PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export PROFILING=0
export OPT=-O3 
export OPENGL=0
export LIBPNG=0
export OUTPUT=1
export USEGLUT=0
export MPI=0

P2DDIR = p2d
MISCDIR = misc

all:
	# Setup link to different modules
	ln -fs gravity_tree.c ../../src/gravity.c
	ln -fs boundaries_shear.c ../../src/boundaries.c
	ln -fs integrator_sei.c ../../src/integrator.c
	ln -fs collisions_tree.c ../../src/collisions.c
	# Setup link to problem file
	ln -fs ../$(PROBLEMDIR)/problem.c ../../src/problem.c
	# Compile
	$(MAKE) -C ../../src/
	$(MAKE) -C ./$(P2DDIR)/
	$(MAKE) -C ./$(MISCDIR)/
	# Copy result
	cp ../../src/nbody .
	cp ./$(P2DDIR)/*.out .
	rm ./$(P2DDIR)/*.out
	cp ./$(MISCDIR)/*.out .
	rm ./$(MISCDIR)/*.out
doc: all
	$(MAKE) -C ../../src/ doc

clean:
	$(MAKE) -C ../../src/ clean
	$(MAKE) -C ./$(P2DDIR)/ clean
	$(MAKE) -C ./$(MISCDIR)/ clean
	rm -vf nbody
	rm *.out
