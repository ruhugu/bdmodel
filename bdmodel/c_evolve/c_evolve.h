#pragma once
#include "dranxor2/dranxor2C.h"
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

void c_evolve(long int* in_heights, long int* out_heights, int length, 
              long int nsteps, long int* pbc, 
              void (*depositfunc)(int, long int*, long int*));
void c_evolveBD(long int* in_heights, long int* out_heights, int length, 
              long int nsteps, long int* pbc);
void c_evolveRD(long int* in_heights, long int* out_heights, int length, 
                long int nsteps, long int* pbc);
void c_evolveRDdiff(double* in_heights, double* out_heights, double ht,
                    int length, long int nMCsteps, long int* pbc);
void c_evolveRelaxRD(long int* in_heights, long int* out_heights, int length, 
                long int nsteps, long int* pbc);
void c_evolveRelaxRDdiff(double* in_heights, double* out_heights, double ht,
                         int length, long int nsteps, long int* pbc);
void c_evolveWalledEW(double* in_heights, double* out_heights, double ht,
                      double a, double p, int length, long int nsteps, 
                      long int* pbc);

void c_depositBD(int j_latt, long int *heights, long int* pbc);
void c_depositRD(int j_latt, long int *heights, long int* pbc);
void c_depositRelaxRD(int j_latt, long int *heights, long int* pbc);

long int max_three(long int a, long int b, long int c);
