//  Individual.h

#ifndef ____Individual__
#define ____Individual__

#include "Genome_Version2.hpp"

using namespace std;

class Individual : public Genome_Version2 {
public:
    float Fitness;
//    array<float, 2> Scaling_Fitness;
    short ancestry_index_1;     // index of the the ancester in the array Individuals (see EvoAlg.h) important for Evolution_Histroy (also see EvoAlg.h) to find the lineage of one individual of the final generation at the end.
    short ancestry_index_2;
    
    Individual() : Fitness(0)   /*, Scaling_Fitness({0,0})*/ {Genome_Version2();};   // trivial Constructor: Initially sets all Fitness values to 0. Fitness values must then be properly initialized in the program that creates the object. A trivial constructor is, however, necessary if e.g. an array of this object is to be created.
    
    template<class T>
    Individual(T FitnessEvaluator) {Genome_Version2(); FitnessEvaluator(*this);};                 // non-trivial Constructor (template) that takes a functor 'FitnessEvaluator' that is supposed to calculate Fitness. Constructor applies FitnessEvaluator to the object itself to initialize all Fitness values properly
};


#endif /* defined(____Individual__) */
