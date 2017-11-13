all: evolution.so 

evolution.so: c_evolve.c c_evolve.h evolve.pyx  dranxor2/dranxor2C.o
	python setup.py build_ext --inplace 

dranxor2/dranxor2C.o : dranxor2/dranxor2C.f dranxor2/dranxor2C.h
	gfortran dranxor2/dranxor2C.f -O3 -fPIC -c -o $@
