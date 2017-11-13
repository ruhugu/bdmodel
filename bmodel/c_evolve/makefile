all: evolution.so dranxor2/dranxor2.o

evolution.so: c_evolve.c c_evolve.h evolve.pyx  dranxor2/dranxor2.o
	python setup.py build_ext --inplace 

dranxor2/dranxor2.o : dranxor2/dranxor2.f
	gfortran dranxor2/dranxor2.f -O3 -fPIC -c -o $@
