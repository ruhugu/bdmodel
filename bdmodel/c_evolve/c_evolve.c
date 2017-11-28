#include "c_evolve.h"


// Generic function for system evolution. Deposits "nsteps" particles
// on the lattice with heights "in_heights" with "length" number of
// points. The deposition is done using the function "depositfunc".
// The heights after the deposition are stored in out_heights.
void c_evolve(long int* in_heights, long int* out_heights, int length, 
              long int nsteps, long int* pbc, 
              void (*depositfunc)(int, long int*, long int*))
{
    int j_latt;
    
    //Initialize out_heights with the values in in_heights
    for (int i=0; i<length; i++) out_heights[i] = in_heights[i];

    for (int t=0; t<nsteps; t++)
    {
        // Choose randomly the deposition point. 
        // idran_ generates a random integer between 1 and length.
        j_latt = idran_(&length) - 1;
      
        (*depositfunc)(j_latt, out_heights, pbc);
    }
    
    return;
}

// Evolution function using ballistic deposition model. The model details
// are implemented via the funcion c_depositBD.
void c_evolveBD(long int* in_heights, long int* out_heights, int length, 
                long int nsteps, long int* pbc)
{
    c_evolve(in_heights, out_heights, length, nsteps, pbc, c_depositBD);  
    return;
}


// Evolution function using ballistic deposition model. The model details
// are implemented via the funcion c_depositRD.
void c_evolveRD(long int* in_heights, long int* out_heights, int length, 
                long int nsteps, long int* pbc)
{
    c_evolve(in_heights, out_heights, length, nsteps, pbc, c_depositRD);  
    return;
}

// Evolution solving the stochastic differential equation for the random
// deposition model using Euler method.
// Now the steps are in ht steps 
void c_evolveRDdiff(double* in_heights, double* out_heights, double ht,
                    int length, long int nsteps, long int* pbc)
{
    int j_latt;
    double sqrtht;

    // Define the square root of the timestep outside to improve speed
    sqrtht = sqrt(ht);
    
    //Initialize out_heights with the values in in_heights
    for (int i=0; i<length; i++) out_heights[i] = in_heights[i];

    for (int t=0; t<nsteps; t++)
    {
        for (j_latt=0; j_latt<length; j_latt++)
        {
            out_heights[j_latt] += ht + sqrtht*drang_();
        }
    }
        
    return;
}


// Throw a particle in the "j_latt" point of the lattice given
// by "heights" and return the resulting lattice heights using
// the ballistic deposition model.
// j_latt must be a number between 0 and length - 1 (the 
// function does not check this, so be careful).
void c_depositBD(int j_latt, long int *heights, long int* pbc)
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


// Throw a particle in the "j_latt" point of the lattice given
// by "heights" and return the resulting lattice heights using
// the random deposition model.
// j_latt must be a number between 0 and length - 1 (the 
// function does not check this, so be careful).
void c_depositRD(int j_latt, long int *heights, long int* pbc)
{
        // We use (j_latt + 1) because of the way pbc is
        // defined (heights[pbc[i]] = heights[i-1] with
        // the exceptions heights[pbc[0]] = heights[length-1] 
        // and heights[pbc[length+1]] = heights[0]).
        heights[j_latt] += 1; 
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
    
 
