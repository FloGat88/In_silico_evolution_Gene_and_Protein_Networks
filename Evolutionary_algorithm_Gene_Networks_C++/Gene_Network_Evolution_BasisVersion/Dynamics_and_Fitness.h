//  Dynamics_and_Fitness.h

#ifndef __EmbrionicDevelopment__Dynamics_and_Fitness__
#define __EmbrionicDevelopment__Dynamics_and_Fitness__

#include "Individual.h"
//#include "LamSimAnn_Interface.hpp"
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/odeint.hpp>

// #define state_type_matrix boost::numeric::ublas::matrix< float, boost::numeric::ublas::column_major , boost::numeric::ublas::bounded_array<float,N_tot*N_L> >
// #define state_type_array array<float,N_tot*N_L>

using namespace std;


struct ode1_reduced {      // struct relevant for reaction_diffusion_dynamics_fast; like ode1 just without a address to the respective element. Can indeed be initialized from ode1
    float rate;
    short ind;
    ode1_reduced(): rate(0),ind(-1) {}
    ode1_reduced(ode1& ode1): rate(ode1.rate),ind(ode1.ind) {}
};

struct ode2_reduced {    // same as above, just for ode2
    float rate;
    short ind1;
    short ind2;
    ode2_reduced(): rate(0),ind1(-1),ind2(-1) {}
    ode2_reduced(ode2& ode2): rate(ode2.rate),ind1(ode2.ind1),ind2(ode2.ind2) {}
};



/*
class reaction_diffusion_dynamics     // Functor that is once initialized for a certain individual and can then be used to for the Ode Solver in the standard way
{
private:
//    float D;
//    float dx;
    float D_dxdx;                      // relevant constant D/dx^2 for the dynamics to be evaluated only once during construction for increase in efficency
    array< vector<ode1> , N_tot >& OdeSyst1; // Reference to OdeSyst1 and ~2 of the Individual whos dynamics is to be simulated (being initialized at construction from this individual).
    array< vector<ode2> , N_tot >& OdeSyst2;

public:
    reaction_diffusion_dynamics(Individual&, float, float);  // Constructor
    
    template<short N_L>
    void operator () (const array< array<float, N_L>, N_tot > &x , array< array<float, N_L>, N_tot > &dxdt , const double );    // state type array<array> is not supported by the Ode Solver. Is first implemetation of operator () is therfore no longer relevant
    
    template<short N_L>
    void operator () (const state_type_matrix& , state_type_matrix& , const double);
    
};
*/



template<short N_L>
class reaction_diffusion_dynamics_fast    // Functor similar as the one above that is however assumed to work more efficiently as, when constructed, it copies all relevant entries of OdeSyst1 and OdeSyst2 in a compact one-dimensional array format that is supposed to provide faster access.
{
// template<short N> friend class LamSimAnn_Interface;
    
private:
//    float D;
//    float dx;
    float D_dxdx;
    array<ode1_reduced, OdeSyst1_MaxSize> OdeSyst1;
    array<size_t, N_tot> OdeSyst1_numbers;
    array<ode2_reduced, OdeSyst2_MaxSize> OdeSyst2;
    array<size_t, N_tot> OdeSyst2_numbers;
    
    vector< vector< float* >> Rates_distributor;   // for the interface to LamSimAnn

public:
    reaction_diffusion_dynamics_fast(Individual& Individual, float D, float dx);    // Constructor
    
/*    template<size_t N_L>
    void operator () (const array< array<float, N_L>, N_tot > &x , array< array<float, N_L>, N_tot > &dxdt , const double); // state type array<array> is not supported by the Ode Solver. Is first implemetation of operator () is therfore no longer relevant
    
    template<size_t N_L>
    void operator () (const state_type_matrix& , state_type_matrix& , const double);
*/
    void operator () (const boost::array<float, N_L*N_tot>& , boost::array<float, N_L*N_tot>& , const double /* t */ );    // utrafast: supposedly the fastest provider for the dynamics
};





template<short N_L>                             // template parameter is discretisation number
class FitnessEvaluation_2GaussianPeaks          // Functor that must be initialized and then can be used to evaluate the Fitness of an Individual to create a predifined Pattern of the two Gaussian Peaks a distance of L/3 apart with std deviation (normally sigma = L/6) scaling with the system size
{
private:
    float L;            // length of system
    float dx;           // Discretisation
    float sigma;        // Std deviation for the Gaussian Peaks of the Target Concentration Prifiles
    array<float, N_L> TargetConcentrationProfile1;   // Target Concentration Prifile (Gaussian Peak at L/3 with var = sigma^2) for Essential_Protein[1]
    array<float, N_L> TargetConcentrationProfile2;  // Target Concentration Prifile (Gaussian Peak at 2/3 *L with var = sigma^2) for Essential_Protein[2]
    float TrivialFitnessValue_L1;                     // Fitness Value that a trivial network (with no components) would have in the L1 Norm. Used for Normalisation of Fitness when ProbeProfile has been called with the L1 Norm variant.
    
public:
    FitnessEvaluation_2GaussianPeaks(float length);         // Constructor: input: system length; initialise all member variables, in particular calculate TargetConcentrationProfile1 and -2 and TrivialFitnessValue so that these are only evaluated once for the whole program
    
private:
    float ProbeProfile_2GaussianPeaks_L1Norm (const float* ProbeProfile1, const float* ProbeProfile2); // private member function for internal calculation of the L^1 metric distance between the probe profile and the Target Profile: input: Reference to the arrays containing the Concentration profiles of the two Proteins to be tested (Essential_Protein[1] and ~[2]). Returns: the normed (average) L^1 metric distance between the probe profile and the Target Profile, d_{L^1}(ProbeProfile,TargetProfile) = 1/L * int_0^L(|ProbeProfile-TargetProfile|). Note that the normalization with 1/L is necessary to later directly compare Fitness values for different system sizes.
    
    float ProbeProfile_2GaussianPeaks_L2Norm (const float* ProbeProfile1, const float* ProbeProfile2);  // same as above but with averaged L^2 norm instead of L^1. Can be used instead of the L^1 norm. In that case the function call in the operator () definition must be changed and the calculation of the TrivialFitnessValue must be changed in the constructor  
    
public:
    float operator () (Individual& individual, float D, double T_max); // overloaded function call operator that calculates fitness: input: Reference to an Individual whos Fitness is to be evaluated, Diffusion constant and Time that the simulation of the dynamics shall take. Returns: the Fitness value defined as -d_{L^1}(ProbeProfile,TargetProfile) *100/TrivialFitnessValue + 100, i.e. the L^1 metric distance relative to the TrivialFitnessValue, scaled to %-values and inversed and shifted such that the maximal fitness value is now 100 and that of a trivial network is 0. Fintess of n roughly means "n% close to Target (or to maximal fitness)". Note that negative fitness values are possible (-> worse performance than the trivial network)
};







template<short N_L, short N_profiles>                             // template parameter is discretisation number
class FitnessEvaluation_MatchProfiles          // Functor that must be initialized and then can be used to evaluate the Fitness of an Individual to create a predifined Pattern of the two Gaussian Peaks a distance of L/3 apart with std deviation (normally sigma = L/6) scaling with the system size
{
private:
    float L;            // length of system
    float dx;           // Discretisation
    array<array<float, N_L>, N_profiles> TargetConcentrationProfiles;   // Target Concentration Prifiles
    float TrivialFitnessValue_L1;                     // Fitness Value that a trivial network (with no components) would have in the L1 Norm. Used for Normalisation of Fitness when ProbeProfile has been called with the L1 Norm variant.
    
public:
    FitnessEvaluation_MatchProfiles() {};      // trivial constructor
    FitnessEvaluation_MatchProfiles(array<array<float, N_L>, N_profiles> TargetConcentrationProfiles, float length);         // Constructor: input: Profiles of all proteins that are to be matched, system length; initialise all member variables, in particular calculate TargetConcentrationProfile1 and -2 and TrivialFitnessValue so that these are only evaluated once for the whole program
    
private:
    float ProbeProfiles_L1Norm (array<const float*, N_profiles> ProbeProfiles); // private member function for internal calculation of the L^1 metric distance between the probe profile and the Target Profile: input: Reference to the arrays containing the Concentration profiles of the two Proteins to be tested (Essential_Protein[1] and ~[2]). Returns: the normed (average) L^1 metric distance between the probe profile and the Target Profile, d_{L^1}(ProbeProfile,TargetProfile) = 1/L * int_0^L(|ProbeProfile-TargetProfile|). Note that the normalization with 1/L is necessary to later directly compare Fitness values for different system sizes.
    
    float ProbeProfiles_L2Norm (array<const float*, N_profiles> ProbeProfiles);  // same as above but with averaged L^2 norm instead of L^1. Can be used instead of the L^1 norm. In that case the function call in the operator () definition must be changed and the calculation of the TrivialFitnessValue must be changed in the constructor
    
public:
    float operator () (Individual& individual, float D, double T_max, reaction_diffusion_dynamics_fast<N_L>& RDD); // overloaded function call operator that calculates fitness: input: Reference to an Individual whos Fitness is to be evaluated, Diffusion constant, Time that the simulation of the dynamics shall take as well as a reference to an istance of a dynamics provider. Returns: the Fitness value defined as -d_{L^1}(ProbeProfile,TargetProfile) *100/TrivialFitnessValue + 100, i.e. the L^1 metric distance relative to the TrivialFitnessValue, scaled to %-values and inversed and shifted such that the maximal fitness value is now 100 and that of a trivial network is 0. Fintess of n roughly means "n% close to Target (or to maximal fitness)". Note that negative fitness values are possible (-> worse performance than the trivial network)
};








template<short N_L>
class FitnessEvaluation_2PeakScaling {
    
private:
    FitnessEvaluation_2GaussianPeaks<N_L> FitnessEvaluatorL1;
    FitnessEvaluation_2GaussianPeaks<N_L> FitnessEvaluatorL2;
    float D1;
    float D2;
    double T_max;
    
public:
    FitnessEvaluation_2PeakScaling(float L1, float D1, float L2, float D2, double T_max);
    
    float operator () (Individual& individual);
};






template<short N_L>
class FitnessEvaluation_ComplementaryGradients {
    
private:
    FitnessEvaluation_MatchProfiles<N_L, 2> FitnessEvaluator;
    float D;
    float dx;
    double T_max;
    
public:
    FitnessEvaluation_ComplementaryGradients(float L, float D, double T_max);
    
    float operator () (Individual& individual);     // dynamics provider is instantiated automatically from individual
    float operator () (Individual& individual, reaction_diffusion_dynamics_fast<N_L>& RDD);      // overloaded method which takes reference to an instantiation of the dynamics provider explicitly as argument. For use in LamSimAnn. Only the memeber variables of RDD are varied during optimization by LamSimAnn, members of individual are not touched at the first place for reasons of efficiency. Therefore, also the fitness value is not written in Individual.fitness in contrast to the former ()-method. Only when optimization by LamSimAnn is finished, the member variables of individual are adapted to the optimal fit parameters.  
};



#include "Dynamics_and_Fitness.tpp"

#endif /* defined(__EmbrionicDevelopment__Dynamics_and_Fitness__) */
