//  Individual.h

#ifndef ____Individual__
#define ____Individual__

#include "Genome.h"

using namespace std;

class Individual : public Genome {
public:
    float Fitness;
    array<float, 2> Scaling_Fitness;
    
    Individual() : Fitness(0), Scaling_Fitness({0,0}) {Genome();};   // trivial Constructor: Initially sets all Fitness values to 0. Fitness values must then be properly initialized in the program that creates the object. A trivial constructor is, however, necessary if e.g. an array of this object is to be created.
    
    template<class T>
    Individual(T FitnessEvaluator) {Genome(); FitnessEvaluator(*this);};                 // non-trivial Constructor (template) that takes a functor 'FitnessEvaluator' that is supposed to calculate Fitness. Constructor applies FitnessEvaluator to the object itself to initialize all Fitness values properly
};


#endif /* defined(____Individual__) */
