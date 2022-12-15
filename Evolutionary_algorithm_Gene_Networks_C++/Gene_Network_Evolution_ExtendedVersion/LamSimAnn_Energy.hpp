//
//  LamSimAnn_Interface.hpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 10/08/16.
//  Copyright © 2016 Florian Gartner. All rights reserved.
//

#ifndef LamSimAnn_Interface_hpp
#define LamSimAnn_Interface_hpp

#include "LamSimAnn.h"
#include "Dynamics.hpp"
#include "Individual.h"

template<class T, short N_L>
class LamSimAnn_Energy : public Energy_base
{
private:
    vector<double> LB;
    vector<double> UB;
    int N_par;
    
    T& FitnessEvaluator;
    Individual& individual;
    ReactionDiffusionDynamics<N_L>& RDD;
    double lower_bound = FLT_MIN;
    double upper_bound = 20;
    
public:
    LamSimAnn_Energy(T& FitnessEvaluator, Individual& individual, ReactionDiffusionDynamics<N_L>& RDD) : FitnessEvaluator(FitnessEvaluator), individual(individual), RDD(RDD)        // Constructor
    {
        N_par = RDD.get_NumberOfIndependentParameters();
        LB.resize(N_par,lower_bound);
        UB.resize(N_par,upper_bound);
    };
    
    virtual double Energy (double *par)    // return fitness value (being minimized)
    {
        RDD.modify_OdeSyst(par);
        
        return -FitnessEvaluator(individual, RDD);
    };
    
    virtual unsigned int GetN() {return N_par;}; // return number of parameters
    
    virtual double* GetLB() {return &LB[0];}; // return lower boundaries of parameters
    
    virtual double* GetUB() {return &UB[0];}; // return upper boundaries of parameters
    
    ReactionDiffusionDynamics<N_L>& GetDynamicsProvider() {return RDD;};   // return reference to RDD for possibility of accessing data and calling functions therein from EvoAlg.LamSimAnn_Fit
};




#endif /* LamSimAnn_Interface_hpp */
