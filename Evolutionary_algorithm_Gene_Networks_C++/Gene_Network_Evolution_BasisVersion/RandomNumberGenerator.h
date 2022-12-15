//  RandomNumberGenerator.h

#ifndef EmbrionicDevelopment_RandomNumberGenerator_h
#define EmbrionicDevelopment_RandomNumberGenerator_h

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <cfloat>

using namespace std;

inline float UnifRand()           // creates a [0,1)-uniform random number
{
    return rand() / (RAND_MAX+1.0);     // write 1.0 instead of 1 to convert result of addition and division to double
}


inline short ChooseInd(short N)     // randomly & uniformly picks an integer number between 0 and N-1 (using the [0,1) unif random number r) that can be used as an index for an array or vector of dimension N
{
    static double r = rand() / (RAND_MAX+1.0);              // to enhance efficancy use the same random number several times to "extract an index out of it"
    static short count = 0;                    // use static variables r that stores the random number and count that counts how often the random number has already been used here. Only after every 7th time create a new random number (7 doesn't have a special meaning despite it's the virgin number and the number of capital sins and virtues).
    if (count++ > 7)  {
        r = rand() / (RAND_MAX+1.0);
        count = 1;
    }
    
    short s = short(r * N);
    r = (r * N) - s;                      // reset the class's random number r for the next call of ChooseInd thus avoiding to constantly create a new random number with UniRand()
    return s;
}


template<size_t size>
short ChooseInd_respective_weights(const array<size_t,size>& weights)
{
    short sum = 0;
    short i,d;
    
    for(i=0; i<size; ++i)
        sum += weights[i];
    if (sum==0)
        return -1;
    d = ChooseInd(sum);
    
    for (i=0; i<size; ++i) {
        if (d<weights[i])
            return i;
        else
            d -= weights[i];
    }
    exit(3);
}


#endif
