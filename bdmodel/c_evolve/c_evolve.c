#include "c_evolve.h"

// Evolve "nsteps" timesteps the system given in in_heights with
// "length" lattice points.
// The evolved state of the system is stored in out_heights.
void c_evolve(long int* in_heights, long int* out_heights, int length, 
              long int nsteps, long int* pbc)
{
    int j_latt;
    
    //Initialize out_heights with the values in in_heights
    for (int i=0; i<length; i++) out_heights[i] = in_heights[i];

    for (int t=0; t<nsteps; t++)
    {
        // Choose randomly the deposition point. 
        // idran_ generates a random integer between 1 and length.
        j_latt = idran_(&length) - 1;
        
        c_deposit(j_latt, out_heights, pbc);
    }
    
    return;
}

// Throw a particle in the "j_latt" point of the lattice given
// by "heights" and return the resulting lattice heights.
// j_latt must be a number between 0 and length - 1 (the 
// function does not check this, so be careful).
void c_deposit(int j_latt, long int *heights, long int* pbc)
{
        // We use (j_latt + 1) because of the way pbc is
        // defined (heights[pbc[i]] = heights[i-1] with
        // the exceptions heights[pbc[0]] = heights[length-1] 
        // and heights[pbc[length+1]] = heights[0]).
        heights[j_latt] = max_three(heights[pbc[(j_latt + 1) - 1]],
                                    heights[pbc[(j_latt + 1)]] + 1,
                                    heights[pbc[(j_latt + 1) + 1]]);
        return;
}
    


// Return the maximum of the 3 given numbers.
long int max_three(long int a, long int b, long int c)
{
    long int ret;    
    
    ret = a > b ? a : b;
    ret = c > ret ? c : ret;
    
    return ret;
}
    
    
