//  Dynamics_and_Fitness.h

#ifndef __EmbrionicDevelopment__Dynamics_and_Fitness__
#define __EmbrionicDevelopment__Dynamics_and_Fitness__

#include "Individual.h"
#include "Dynamics.hpp"
#include "Stability_Analysis.h"
#include <boost/numeric/odeint.hpp>
#include "Ode_Integrator.h"
//#include "Implicit_Stepper_Sparse.h"


// Base class for fitness evaluation. Inherit from this class to create easily a derived class that selects for an arbitrary concentration pattern. The derived class only has to specify TargetConcentrationProfile (number of its rows, i.e. number of proteins that form the pattern is the template parameter N_profiles) to select for that pattern in the essential proteins (it is important that the number of essential proteins Number_Essential_Proteins is large enough, i.e. N_profiles + 1 initial protein to make sure that all essential proteins are always available)

template<short N_L, short N_profiles>                             // template parameter is discretisation number
class FitnessEvaluation_MatchProfiles          // Functor that must be initialized and then can be used to evaluate the Fitness of an Individual to create a predifined Pattern of the two Gaussian Peaks a distance of L/3 apart with std deviation (normally sigma = L/6) scaling with the system size
{
private:
    float L;
    float dx;           // Discretisation
    const bool SymmetricMatch;                      // indicates whether also the mirrored target profile should be checked. For example when an initially homogeneous concentration profile is selected to polarise like in 'FitnessEvaluation_Polarisation' and the random initial perturbations decide how the cell polarises, then it makes sense to set this on true
    float TrivialFitnessValue;                     // Fitness Value that a trivial network (with no components) would have in the L1 Norm. Used for Normalisation of Fitness when ProbeProfile has been called with the L1 Norm variant.
protected:
    array<array<float, N_L>, N_profiles> TargetConcentrationProfiles;   // Target Concentration Prifiles
    vector<short> right_boundaries;        // only needed for 'ProbeProfiles_CoarseGrainedL1'. Contains the right boundaries of the intervals in which the local averages are evaluated, e.g. {10,20,30,40} (for N_L = 40); the last element should always be 'N_L';
    vector<float> differences;             // only needed for 'ProbeProfiles_CoarseGrainedL1'. Defines the, either additive or multiplicative minimal differences or ratios, respectively, between the local averages in two successive intervals that are desired and selected for (stronger real differences than described by 'differences' is also ok). The absolute concentration of the first average does not matter, therefore 'differences' contains one element less than 'right_boundaries' (the first entry defines the difference or ratio between the second and first interval and so on, e.g. {3,3,3} for the example of the right_boundaries vector above).
    
    FitnessEvaluation_MatchProfiles(float L, bool SymmetricMatch) : L(L), dx(L/N_L), SymmetricMatch(SymmetricMatch)  {};     // Constructor
    
    void CalculateTrivialFitnessValue();         // Constructor: input: system length, diffusion constant and Time that the simulation of the dynamics shall take;
    
    
private:
    float ProbeProfiles_L1Norm(array<vector<const float*>, N_profiles>& ProbeProfiles); // private member function for internal calculation of the L^1 metric distance between the probe profile and the Target Profile: input: Reference to the arrays containing the Concentration profiles of the two Proteins to be tested (Essential_Protein[1] and ~[2]). Returns: the normed (average) L^1 metric distance between the probe profile and the Target Profile, d_{L^1}(ProbeProfile,TargetProfile) = 1/L * int_0^L(|ProbeProfile-TargetProfile|). Note that the normalization with 1/L is necessary to later directly compare Fitness values for different system sizes. If 'SymmetricMatch' == true also check the mirrored target profile and take the better match of both.
   
    
    template<short p>
    float ProbeProfiles_LPNorm(array<vector<const float*>, N_profiles>& ProbeProfiles);  // same as above but with averaged L^p norm instead of L^1. Can be used instead of the L^1 norm. In that case the function call in the operator () definition must be changed and the calculation of the TrivialFitnessValue must be changed in the constructor
    
    
    float ProbeProfiles_CoarseGrainedL1(array<vector<const float*>, N_profiles>& ProbeProfiles)  {   // calculates first the local averages in the intervals defined by the (right) boundaries in right_boundaries and then evaluates these local averages by comparing the additive or multiplicative differences (or ratios) between local averages in successive intervals with the desired differences (ratios) in 'differences' (the average absolute concentration within the first interval is thereby irrelevant!). It is not punished (by a contribution to 'coarse_grained_match_L1') if the tendency is even stronger than required by 'differences'.
        float coarse_grained_match_L1 = 0;
        float coarse_grained_match_L1_mirrored = 0;
        vector<float> local_averages(right_boundaries.size());         // saves the local averages in the intervals defined by the boundaries in right_boundaries
        
        auto make_localAverages = [&](const float* _ProbeProfile)  {
            short j=0;
            for(short i=0; i<right_boundaries.size(); ++i)  {
                local_averages[i] = 0;               // initialize again with 0
                while(j<right_boundaries[i])
                    local_averages[i] += _ProbeProfile[j++];
                local_averages[i] / ((i==0) ? right_boundaries[0] : right_boundaries[i]-right_boundaries[i-1]);
            }
        };
        auto evaluate_localAverages_additive = [this](vector<float>& averages, float& _coarse_grained_match_L1)  {     // evaluate the local averages stored in 'local_averages' according to 'differences', interpreted as additive differences from the left local average to its right neighbour. The local average at the very left end is thereby arbitrary, only the differences between the local averages count.
            for(short i=0; i<differences.size(); ++i)  {
                float temp = differences[i]-(averages[i+1]-averages[i]);         // temporary variable
                if((differences[i]<0 && temp>0) || (differences[i]>0 && temp<0))     // if the tendency is even stronger, i.e. the real difference between local averages is is more positive (negative) than the desired positive (negative) difference, then do not punish this and set 'temp' which is added to 'coarse_grained_match_L1' to 0
                    temp = 0;
                _coarse_grained_match_L1 += abs(temp);
            }
        };
        auto evaluate_localAverages_multiplicative = [this](vector<float>& averages, float& _coarse_grained_match_L1)  {     // evaluate the local averages stored in 'local_averages' according to 'differences', interpreted now as multiplicative differences or ratios between the right local average and its left neighbour. The local average at the very left end is thereby arbitrary, only the ratios between the local averages count.
            for(short i=0; i<differences.size(); ++i)  {
                float temp = differences[i]*averages[i] - averages[i+1];         // temporary variable
                if((differences[i]>1 && temp<0) || (differences[i]<1 && temp>0))     // if the tendency is even stronger, i.e. the real ratio between local averages is larger (smaller) than the desired ratio >1 (<1), then do not punish this and set 'temp' which is added to 'coarse_grained_match_L1' to 0.
                    temp = 0;
                _coarse_grained_match_L1 += abs(temp);
            }
        };
        
        for (short i=0; i<N_profiles; ++i)  {
            switch(ProbeProfiles[i].size())  {
                case 0:         // in case there is no ProbeProfile for that function, create an array of only 0s that is delivered to the function
                    return TrivialFitnessValue;
                    
                case 1:         // in case there is one ProbeProfile Concentration for that function
                {
                    make_localAverages(ProbeProfiles[i][0]);
                    break;
                }
                default:        // // in case there are several ProbeProfile Concentrations for that function, first calculate the sum of all those in 'Sum_ProbeProfiles' and then let the local averages be calculated for this summarized profile of all the concentrations of the essential elements with this function.
                {
                    array<float,N_L> Sum_ProbeProfiles = {};
                    for(short j=0; j<ProbeProfiles[i].size(); ++j)
                        for(short k=0; k<N_L; ++k)
                            Sum_ProbeProfiles[k] += ProbeProfiles[i][j][k];
                    make_localAverages(&Sum_ProbeProfiles[0]);
                    break;
                }
            }
            
            evaluate_localAverages_multiplicative(local_averages, coarse_grained_match_L1);        // add the contribution of the coarse grained profile for ONE essential function (the function with index i in the loop) to 'coarse_grained_match_L1'.
            if(SymmetricMatch)  {       // if SymmetricMatch == true reverse local_averages and evaluate the mirrord profile, adding the result of the evaluation to the variable coarse_grained_match_L1_mirrored.
                reverse(local_averages.begin(), local_averages.end());
                evaluate_localAverages_multiplicative(local_averages, coarse_grained_match_L1_mirrored);
            }
        }
        if(SymmetricMatch && coarse_grained_match_L1_mirrored < coarse_grained_match_L1)     // return the minimum of either the fit with the profile or the fit with the mirrored profile, each normalised properly (remeber that the first local average is irrelevant, therefore the -1 in the normalisation factor)
            return coarse_grained_match_L1_mirrored/(local_averages.size()-1);
        else
            return coarse_grained_match_L1/(local_averages.size()-1);
    };
    
    
protected:
    float MatchProfiles (Individual& individual, ReactionDiffusionDynamics<N_L>& RDD)  { // function to calculate fitness: input: Reference to an Individual whos Fitness is to be evaluated, as well as a reference to an istance of a dynamics provider (a class that provides the right hand side of the ODE). Returns: the Fitness value defined as -d_{L^1}(ProbeProfile,TargetProfile) *100/TrivialFitnessValue + 100, i.e. the L^1 metric distance relative to the TrivialFitnessValue, scaled to %-values and inversed and shifted such that the maximal fitness value is now 100 and that of a trivial network is 0. Fintess of n roughly means "n% close to Target (or to maximal fitness)". Note that negative fitness values are possible (-> worse performance than the trivial network)

        using namespace boost::numeric::odeint;
        
        short STEPPER_TYPE = 1;
        
        if(individual.Get_NumberElements() == 0)     // in case the network is empty return 0 for the fitness (trivial fitness!). This is necessary since otherwise 'concentrations' is empty and Make_ProbeProfiles would be called with the NULL pointer which later causes BAD_ACCESS
            return 0;
        
        vector<float> concentrations(individual.Get_NumberElements() * N_L, 0);     // define the state type that the Ode Solver deals with as a vector that is initialized with 0
        // set initial protein concentrations
        individual.FillIn_InitialHomogeneousProteinConcentrations(concentrations, N_L);     // maybe make this and Make_ProbeProfiles a member function of RDD as it can only be called after an RDD object has been created !!!!
#ifndef NDEBUG
        individual.Save_InitialConcentrations(concentrations);
#endif
        array<vector<const float*>, N_profiles> ProbeProfiles;    // create the array ProbeProfiles with the addresses of the start positions of the concerning Essential Elements (whos profiles are to be matched by the Algorithm) within concentrations
        individual.Make_ProbeProfiles(ProbeProfiles, &concentrations[0], N_L);
        
        //ReactionDiffusionDynamics<N_L> RDD(individual, dx, D);
        
        typedef runge_kutta_dopri5< vector<float> > stepper_type_explicit;
//        typedef StepperRoss<ReactionDiffusionDynamics<N_L> , N_L*N_tot> stepper_type_implicit;     // !!!!! TODO: adapt implicit stepper
        //typedef runge_kutta_cash_karp54< boost::array<float, N_L*N_tot> > stepper_type;
        double initial_step = 1.0e-4;
        int steps1, steps2;
        if(STEPPER_TYPE==1)         // if STEPPER_TYPE = 1 use the explicit runge kutta routine, for STEPPER_TYPE = 2 use the implicit rosenbrock method for the integration
            steps1 = integrate_adaptive( make_controlled<stepper_type_explicit>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, 0.0, T_max, initial_step ); // matlab default error tolerances are (AbsTol = 1.0e-6, RelTol = 1.0e-3);; matlab also uses runge_kutta_dopri5 as its ode45; alternatively try runge_kutta_cash_karp54
            //size_t steps = integrate_const( runge_kutta4< boost::array<float, N_L*N_tot> >(), boost::ref(RDD) , concentrations , 0.0 , T_max , 0.1 );
        else  {         //if(STEPPER_TYPE==2)
            Output out;         // output object for the method
            Doub t1 = 0.0, t2 = T_max, atol = 1.0e-4, rtol = 1.0e-3, h1 = 1.0e-3, hmin = 0.0;
   //             Odeint<stepper_type_implicit, N_L*N_tot> ode_implicit(concentrations, t1, t2, atol , rtol, h1, hmin, out, RDD);  // !!!!
   //             steps1 = ode_implicit.integrate();        // !!! TODO: adapt implicit solver
            }
#ifndef NDEBUG
        for(float c : concentrations)       // check if concentrations are actually numbers (not nan) and bigger or equal 0
            assert(c>=0);
        individual.Make_SpeciesObject();
        individual.Check_MassConservation(concentrations);
#endif
        float fitness1 = -ProbeProfiles_L1Norm( ProbeProfiles ) * 100./TrivialFitnessValue + 100.; // normalize the metric distance between the Probe Profiles and the Target Profile to a fitness value in [0,100]
        
        double T_check = 0.2 * T_max;
    //    double T_check = 20;           // time span after which a second check is done to ensure that the profile was indeed stationary
        if(STEPPER_TYPE==1)
            steps2 = integrate_adaptive( make_controlled<stepper_type_explicit>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, T_max, T_max + T_check, initial_step );   // integrate the system T_check time longer
        else {         //if(STEPPER_TYPE==2)
            Output out;
            Doub t1 = T_max, t2 = T_max + T_check, atol = 1.0e-4, rtol = 1.0e-3, h1 = 1.0e-3, hmin = 0.0;
        //        Odeint<stepper_type_implicit, N_L*N_tot> ode_implicit(concentrations, t1, t2, atol , rtol, h1, hmin, out, RDD);   // !!!!
        //        steps2 = ode_implicit.integrate();        // !!! TODO: adapt implicit stepper!
            }
        
        cout << steps1 + steps2 << '\n';       // display steps as the sum of the integration steps from both parts.
#ifndef NDEBUG
        for(float c : concentrations)       // check if concentrations are actually numbers (not nan) and bigger or equal 0
            assert(c>=0);
        individual.Check_MassConservation(concentrations);
#endif
        float fitness2 = -ProbeProfiles_L1Norm( ProbeProfiles ) * 100./TrivialFitnessValue + 100.;    // check again for the desired profiles
        cout << "fitness1 = " << fitness1 << "    " << "fitness2 = " << fitness2 << '\n';
        return (fitness1 < fitness2 ? fitness1 : fitness2);     // return the final fitness value as the minimum of the two fitness values
    }

public:
    float operator () (Individual* individual);     // dynamics provider is instantiated automatically from individual

    float operator () (Individual& individual, ReactionDiffusionDynamics<N_L>& RDD);      // overloaded method which takes reference to an instantiation of the dynamics provider explicitly as argument. For use in LamSimAnn. Only the memeber variables of RDD are varied during optimization by LamSimAnn, members of individual are not touched at the first place for reasons of efficiency. Therefore, also the fitness value is not written in Individual.fitness in contrast to the former ()-method. Only when optimization by LamSimAnn is finished, the member variables of individual are adapted to the optimal fit parameters.
};





// Fitness Evaluators
// Functor Classes (overloaded function call operator is inherited from the base class) that select for certain target concentration profiles. Inherit from the class FitnessEvaluation_MatchProfiles and only specify the Target Profile (the array TargetconcentrationProfile in the base class) that is to be selected for.



// #1: Complementary Gradients
template<short N_L>
class FitnessEvaluation_ComplementaryGradients : public FitnessEvaluation_MatchProfiles<N_L, 2>
{
public:
    FitnessEvaluation_ComplementaryGradients(float L);
};


// #2: Two Gaussian Peaks
template<short N_L>
class FitnessEvaluation_2GaussianPeaks : public FitnessEvaluation_MatchProfiles<N_L, 2> {
public:
    FitnessEvaluation_2GaussianPeaks(float L);
};


// #3: DrosophilaGapgenes
template<short N_L>
class FitnessEvaluation_DrosophilaGapgenes : public FitnessEvaluation_MatchProfiles<N_L, 4> {
public:
    FitnessEvaluation_DrosophilaGapgenes(float L);
};


// #4: Polarisation
template<short N_L>          // select for a Theta-Step function that switches from 0 to 1 at L/2 as the concentration profile of a membrane-bound component
class FitnessEvaluation_Polarisation : public FitnessEvaluation_MatchProfiles<N_L, 1>
{
public:
    FitnessEvaluation_Polarisation(float L);
};


// #5: Polarisation Coarse Grained
template<short N_L>          // select for a Theta-Step function that switches from 0 to 1 at L/2 as the concentration profile of a membrane-bound component
class FitnessEvaluation_Polarisation_CoarseGrained : public FitnessEvaluation_MatchProfiles<N_L, 1>
{
public:
    FitnessEvaluation_Polarisation_CoarseGrained(float L);
};


/*
 // #3: Two Gaussian Peaks Scaling
 template<short N_L>
 class FitnessEvaluation_2GaussianPeaks_Scaling {
 
 private:
 FitnessEvaluation_2GaussianPeaks<N_L> FitnessEvaluatorL1;
 FitnessEvaluation_2GaussianPeaks<N_L> FitnessEvaluatorL2;
 float dx;
 
 public:
 FitnessEvaluation_2GaussianPeaks_Scaling(float L1, float L2, float D, double T_max);
 
 float operator () (Individual& individual);     // dynamics provider is instantiated automatically from individual
 float operator () (Individual& individual, ReactionDiffusionDynamics<N_L>& RDD);      // overloaded method which takes reference to an instantiation of the dynamics provider explicitly as argument. For use in LamSimAnn. Only the memeber variables of RDD are varied during optimization by LamSimAnn, members of individual are not touched at the first place for reasons of efficiency. Therefore, also the fitness value is not written in Individual.fitness in contrast to the former ()-method. Only when optimization by LamSimAnn is finished, the member variables of individual are adapted to the optimal fit parameters.
 };
 */





#include "Fitness_Evaluation.tpp"

#endif /* defined(__EmbrionicDevelopment__Dynamics_and_Fitness__) */
