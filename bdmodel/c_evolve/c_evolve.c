#include "c_evolve.h"

// Evolve "nsteps" timesteps the system given in in_heights with
// "length" lattice points.
// The evolved state of the system is stored in out_heights.
void c_evolve(long int* in_heights, long int* out_heights, int length, 
              long int nsteps)
{
    int j_latt;
    int* pbc;
    
    // Generate a length+2 size vector to implement periodic 
    // boundary conditions. in_heights[pbc[i]] will return the
    // height of the i-th position, being 1 the first lattice
    // point and length the last. i=0 returns the height of the 
    // last point, while i=lenght+1 gives the first.
    pbc = (int*) malloc((length+2)*sizeof(int*));
    pbc[0] = length - 1;
    pbc[length + 1] = 0;
    for (int i=1; i<length+1; i++) pbc[i] = i - 1;

    //Initialize out_heights with the values in in_heights
    for (int i=0; i<length; i++) out_heights[i] = in_heights[i];

    for (int t=0; t<nsteps; t++)
    {
        // Choose randomly the deposition point. 
        // idran_ generates a random integer between 1 and length.
        j_latt = idran_(&length);
        out_heights[j_latt] = max_three(out_heights[pbc[j_latt - 1]],
                                        out_heights[j_latt] + 1,
                                        out_heights[pbc[j_latt + 1]]);
    }

    free(pbc);
    
    return;
}


long int max_three(long int a, long int b, long int c)
// Return the maximum of the 3 given numbers.
{
    long int ret;    
    
    ret = a > b ? a : b;
    ret = c > ret ? c : ret;
    
    return ret;
}
    
    
