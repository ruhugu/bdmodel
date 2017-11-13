#pragma once
#include "dranxor2/dranxor2C.h"
#include <stdlib.h>

void c_evolve(long int* in_heights, long int* out_heights, int length, 
              long int nsteps);
long int max_three(long int a, long int b, long int c);
