PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export PROFILING=0
export OPT=-O3 
export OPENGL=1
export LIBPNG=0
export OUTPUT=1
export USEGLUT=0
export MPI=0

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
	# Copy result
	cp ../../src/nbody .
	$(MAKE) -C ./p2d/
	cp ./p2d/*.out .
	rm ./p2d/*.out
doc: all
	$(MAKE) -C ../../src/ doc

clean:
	$(MAKE) -C ../../src/ clean
	$(MAKE) -C ./p2d/ clean
	rm -vf nbody
	rm *.outg
