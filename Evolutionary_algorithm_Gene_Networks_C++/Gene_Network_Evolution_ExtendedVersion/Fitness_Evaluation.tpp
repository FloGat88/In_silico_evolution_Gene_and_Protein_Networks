//  Dynamics_and_Fitness.tpp

//#include "Dynamics_and_Fitness.h"
// #include "Individual.h"
//#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
// using namespace boost::numeric::ublas;




template <short N_L, short N_profiles>
void FitnessEvaluation_MatchProfiles<N_L, N_profiles>::CalculateTrivialFitnessValue()
{
    array<float,N_L> ZeroProfile;                          // create an array of zeros. The address of the first element is copied N_profiles times in the array ProbeProfilesZero that is handed over to ProbeProfiles_L1Norm to calculate the TrivialFitnessValue
    fill_n(ZeroProfile.begin(), N_L, 0);
    array<vector<const float*>, N_profiles> ProbeProfilesZero;
    for(short i=0; i<N_profiles; ++i)
        ProbeProfilesZero[i].push_back(&ZeroProfile[0]);
    
    TrivialFitnessValue = ProbeProfiles_L1Norm(ProbeProfilesZero);
    // TrivialFitnessValue = ProbeProfiles_L2Norm(ProbeProfilesZero);         // if L2 Norm is used, calculate TrivialFitnessValue_L2 instead
};


template<short N_L, short N_profiles>
template<short p>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::ProbeProfiles_LPNorm (array<vector<const float*>, N_profiles>& ProbeProfiles)
{
    float profile_match_LPNorm = 0;
    float profile_match_LPNorm_mirrored = 0;       // L1 norm for the mirrored profile if SymmetricMatch == true
    
    for (short i=0; i<N_profiles; ++i)  {
        switch(ProbeProfiles[i].size())  {
                
            case 0:       // in case there is no ProbeProfile for that function, set 0 for the ProbeProfile
                return TrivialFitnessValue;
                
            case 1:     // in case there is one ProbeProfile Concentration for that function
                for (short j=0; j<N_L; ++j)  {
                    profile_match_LPNorm += pow( abs(ProbeProfiles[i][0][j] - TargetConcentrationProfiles[i][j]) , p);
                    if(SymmetricMatch)          // if SymmetricMatch == true also test the mirrored TargetProfile and take the better fit below
                        profile_match_LPNorm_mirrored += pow( abs(ProbeProfiles[i][0][j] - TargetConcentrationProfiles[i][(N_L-1)-j]) , p);
                }
                break;
                
            default:        // in case there are several ProbeProfile Concentrations for that function, first calculate the sum of all those in 'Sum_ProbeProfiles' and then calculate the LP Norm from this sum.
                array<float,N_L> Sum_ProbeProfiles = {};
                for(short j=0; j<ProbeProfiles[i].size(); ++j)
                    for(short k=0; k<N_L; ++k)
                        Sum_ProbeProfiles[k] += ProbeProfiles[i][j][k];
                for (short j=0; j<N_L; ++j)  {
                    profile_match_LPNorm += pow( abs(Sum_ProbeProfiles[j] - TargetConcentrationProfiles[i][j]) , p);
                    if(SymmetricMatch)      // if SymmetricMatch == true also test the mirrored TargetProfile and take the better fit below
                        profile_match_LPNorm_mirrored += pow( abs(Sum_ProbeProfiles[j] - TargetConcentrationProfiles[i][(N_L-1)-j]) , p);
                }
                break;
        }
    }
    if(SymmetricMatch && profile_match_LPNorm_mirrored < profile_match_LPNorm)     //Return the normed (average) L^p metric distance between the probe profile and the Target Profile, d_{L^p}(ProbeProfile,TargetProfile) = (1/L * int_0^L(|ProbeProfile-TargetProfile|^p))^(1/p). return the minimum of the two calculated values (one for the profile and the other one for the mirrored profile)
        return pow( profile_match_LPNorm_mirrored/N_L, 1.0/p );
    else
        return pow( profile_match_LPNorm/N_L, 1.0/p );
}


template<short N_L, short N_profiles>           // in principle this function is obsolete as the same function could also be covered with the previous function. remove this for more compact code
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::ProbeProfiles_L1Norm (array<vector<const float*>, N_profiles>& ProbeProfiles)
{
    float profile_match_L1Norm = 0;
    float profile_match_L1Norm_mirrored = 0;           // L1 norm for the mirrored profile
    for (short i=0; i<N_profiles; ++i)  {
        switch(ProbeProfiles[i].size())  {
                
            case 0:         // in case there is no ProbeProfile for that function, set 0 for the ProbeProfile
                return TrivialFitnessValue;
                
            case 1:         // in case there is one ProbeProfile Concentration for that function
                for (short j=0; j<N_L; ++j)  {
                    profile_match_L1Norm += abs(ProbeProfiles[i][0][j] - TargetConcentrationProfiles[i][j]);
                    if(SymmetricMatch)
                        profile_match_L1Norm_mirrored += abs(ProbeProfiles[i][0][j] - TargetConcentrationProfiles[i][(N_L-1)-j]);
                }
                break;
                
            default:        // // in case there are several ProbeProfile Concentrations for that function, first calculate the sum of all those in 'Sum_ProbeProfiles' and then calculate the LP Norm from this sum.
                array<float,N_L> Sum_ProbeProfiles = {};
                for(short j=0; j<ProbeProfiles[i].size(); ++j)
                    for(short k=0; k<N_L; ++k)
                        Sum_ProbeProfiles[k] += ProbeProfiles[i][j][k];
                for (short j=0; j<N_L; ++j)  {
                    profile_match_L1Norm += abs(Sum_ProbeProfiles[j] - TargetConcentrationProfiles[i][j]);
                    if(SymmetricMatch)
                        profile_match_L1Norm_mirrored += abs(Sum_ProbeProfiles[j] - TargetConcentrationProfiles[i][(N_L-1)-j]);
                }
                break;
        }
    }
    if(SymmetricMatch && profile_match_L1Norm_mirrored < profile_match_L1Norm)     //Return the normed (average) L^1 metric distance between the probe profile and the Target Profile, d_{L^1}(ProbeProfile,TargetProfile) = 1/L * int_0^L(|ProbeProfile-TargetProfile|); return the minimum of the two calculated values (one for the profile and the other one for the mirrored profile)
        return profile_match_L1Norm_mirrored/N_L;
    else
        return profile_match_L1Norm/N_L;
}

/*
template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::MatchProfiles (Individual& individual, ReactionDiffusionDynamics<N_L>& RDD)
{
    using namespace boost::numeric::odeint;
    
    short STEPPER_TYPE = 2;
    
    boost::array<float, N_L*N_tot> concentrations = {};     // define the state type that the Ode Solver deals with as an boost array that is initialized with 0
    
    // set initial conditions
    for (short i=0; i!=N_P; ++i)               // all existing genes have initial concentration of 1
        if (individual.OdeSyst1[N_P+i].size()>0)        // for each Gene check if the corresponding Protein has at least one entry in OdeSyst1, e.g. due to gene expression (thus implying that the Gene exists) and if so...
            for (short j=i*N_L; j!=(i+1)*N_L; ++j)
                concentrations[j] = 1;      // ...set the initial concentration of the Gene to 1
    
    concentrations[(N_P+individual.Essential_Proteins[0])*N_L] = 5;      // create a high concentration of the maternal RNA (EssentialProtein[0]) at the left end of the embryo, so the cell is assumed to be already poloarized/symmetry-broken
    /*  for (short i=0; i!=N_L; ++i)        // alternatively create a linear gradient in essential protein 1 over the entire embryo
     concentrations[(N_P+individual.Essential_Proteins[0])*N_L + i] = 5-i/N_L*5;   */
    /*    for (short i=N_P+individual.Essential_Proteins[1]; i<=N_P+individual.Essential_Proteins[2]; ++i)               // the other two essential proteins have initial concentration of 1
     for (short j=i*N_L; j!=(i+1)*N_L; ++j)
     concentrations[j] = 1;      // ...set the initial concentration of the essentail Proteins to 1
     */
/*
    array<const float*, N_profiles> ProbeProfiles;    // create the array ProbeProfiles with the addresses of the start positions of the concerning Essential Proteins (whos profiles are to be matched by the Algorithm) within concentrations
    for(short i=0; i<N_profiles; ++i)
        ProbeProfiles[i] = &concentrations[(N_P+individual.Essential_Proteins[i+1])*N_L];
    
    //ReactionDiffusionDynamics<N_L> RDD(individual, dx, D);
    
    typedef runge_kutta_dopri5< boost::array<float, N_L*N_tot> > stepper_type_explicit;
    typedef StepperRoss<ReactionDiffusionDynamics<N_L> , N_L*N_tot> stepper_type_implicit;
    //typedef runge_kutta_cash_karp54< boost::array<float, N_L*N_tot> > stepper_type;
    int steps1, steps2;
    if(STEPPER_TYPE==1)         // if STEPPER_TYPE = 1 use the explicit runge kutta routine, for STEPPER_TYPE = 2 use the implicit rosenbrock method for the integration
        steps1 = integrate_adaptive( make_controlled<stepper_type_explicit>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, 0.0, T_max, 0.1 ); // matlab default error tolerances are (AbsTol = 1.0e-6, RelTol = 1.0e-3);; matlab also uses runge_kutta_dopri5 as its ode45; alternatively try runge_kutta_cash_karp54
        //size_t steps = integrate_const( runge_kutta4< boost::array<float, N_L*N_tot> >(), boost::ref(RDD) , concentrations , 0.0 , T_max , 0.1 );
    else    {         //if(STEPPER_TYPE==2)
        Output out;         // output object for the method
        Doub t1 = 0.0, t2 = T_max, atol = 1.0e-6, rtol = 1.0e-4, h1 = 1.0e-3, hmin = 0.0;
        Odeint<stepper_type_implicit, N_L*N_tot> ode_implicit(concentrations, t1, t2, atol , rtol, h1, hmin, out, RDD);
        steps1 = ode_implicit.integrate();
    }

    float fitness1 = -ProbeProfiles_L1Norm( ProbeProfiles ) * 100./TrivialFitnessValue + 100.; // normalize the metric distance between the Probe Profiles and the Target Profile to a fitness value in [0,100]
    
    double T_check = 20;           // time span after which a second check is done to ensure that the profile was indeed stationary
    if(STEPPER_TYPE==2)
        steps2 = integrate_adaptive( make_controlled<stepper_type_explicit>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, T_max, T_max + T_check, 0.1 );   // integrate the system T_check time longer
    else {         //if(STEPPER_TYPE==2)
        Output out;
        Doub t1 = T_max, t2 = T_max + T_check, atol = 1.0e-6, rtol = 1.0e-4, h1 = 1.0e-3, hmin = 0.0;
        Odeint<stepper_type_implicit, N_L*N_tot> ode_implicit(concentrations, t1, t2, atol , rtol, h1, hmin, out, RDD);
        steps2 = ode_implicit.integrate();
    }
    
    cout << steps1 + steps2 << '\n';       // display steps as the sum of the integration steps from both parts.
    
    float fitness2 = -ProbeProfiles_L1Norm( ProbeProfiles ) * 100./TrivialFitnessValue + 100.;    // check again for the desired profiles
    cout << "fitness1 = " << fitness1 << "    " << "fitness2 = " << fitness2 << '\n';
    return (fitness1 < fitness2 ? fitness1 : fitness2);     // return the final fitness value as the minimum of the two fitness values
}
*/
/*
template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::operator () (Individual* individual)      // give individual as pointer in order to be able to parallelise the method in the Evovler.
{
    ReactionDiffusionDynamics<N_L> RDD(*individual, dx);
    individual->Fitness = MatchProfiles(*individual, RDD);
    return individual->Fitness;
};
*/
template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::operator() (Individual* individual)      // give individual as pointer in order to be able to parallelise the method in the Evovler.
{
    ReactionDiffusionDynamics<N_L> RDD(*individual, dx);
    Stability_Analysis<N_L>(*individual, RDD, &(individual->Fitness));
    return individual->Fitness;
};

template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::operator () (Individual& individual, ReactionDiffusionDynamics<N_L>& RDD)
{
    return MatchProfiles(individual, RDD);
};




/*
 template<short N_L>
 FitnessEvaluation_2PeakScaling<N_L>::FitnessEvaluation_2PeakScaling(float L1, float D1, float L2, float D2, double T_max) :
 FitnessEvaluatorL1(FitnessEvaluation_2GaussianPeaks<N_L>(L1)),
 FitnessEvaluatorL2(FitnessEvaluation_2GaussianPeaks<N_L>(L2)),
 D1(D1), D2(D2), T_max(T_max)
 {};
 
 
 template<short N_L>
 float FitnessEvaluation_2PeakScaling<N_L>::operator () (Individual& individual)
 {
 float fit1 = individual.Scaling_Fitness[0] = FitnessEvaluatorL1(individual, D1, T_max);
 //    float fit2 = individual.Scaling_Fitness[1] = FitnessEvaluatorL2(individual, D2, T_max);
 //    individual.Fitness = fit1 < fit2 ? fit1 : fit2;
 individual.Fitness = fit1;   //!!!  remove this line when two system sizes are simulated (see below!)
 
 return individual.Fitness;
 }
 */

/*
 template<short N_L>
 float FitnessEvaluation_2PeakScaling<N_L>::operator () (Individual& individual)
 {
 float fit1 = individual.Scaling_Fitness[0] = FitnessEvaluatorL1(individual, D1, T_max);
 float fit2 = individual.Scaling_Fitness[1] = FitnessEvaluatorL2(individual, D2, T_max);
 individual.Fitness = fit1 < fit2 ? fit1 : fit2;
 return individual.Fitness;
 }
 */




template<short N_L>
FitnessEvaluation_ComplementaryGradients<N_L>::FitnessEvaluation_ComplementaryGradients(float L) : FitnessEvaluation_MatchProfiles<N_L, 2>(L, false)
{
    float dx = L/N_L;
    float xi = L/2;       // decay length of the gradient
    for (short i=0; i<N_L; ++i) {
        this->TargetConcentrationProfiles[0][i] = 2 * exp(-i*dx/xi );
        this->TargetConcentrationProfiles[1][i] = 2 * exp(-(N_L-1-i)*dx/xi );
    }
    this->CalculateTrivialFitnessValue();      // calculate the trivial fitness value form the Target Concentration Profile
};




template<short N_L>
FitnessEvaluation_2GaussianPeaks<N_L>::FitnessEvaluation_2GaussianPeaks(float L) : FitnessEvaluation_MatchProfiles<N_L, 2>(L,false)
{
    float dx = L/N_L;
    float sigma = L/6;
    float stdnormalization = 1/sqrt(2*M_PI*pow(sigma,2));
    for (short i=0; i<N_L; ++i) {
        this->TargetConcentrationProfiles[0][i] = 7 * stdnormalization * exp(-pow(i*dx - L/3.0,2)/(2*pow(sigma,2)) );         // previously there was a facotr of 2 instead of 7. 7 has been chosen to make the profile more similar to ComplementrayGradients so that such a seed can be used efficiently.
        this->TargetConcentrationProfiles[1][i] = 7 * stdnormalization * exp(-pow(i*dx - 2.0/3.0*L,2)/(2*pow(sigma,2)) );
    }
    this->CalculateTrivialFitnessValue();
};



template<short N_L>
FitnessEvaluation_DrosophilaGapgenes<N_L>::FitnessEvaluation_DrosophilaGapgenes(float L) : FitnessEvaluation_MatchProfiles<N_L, 4>(L,false)
{
    extern array<float, 1000> Kni, Kr, Gt, Hb;
    extern float Gt_max;                 // extern variables describing the expression profiles of the four gapgenes (stored in Drosophila_Gapgene_data.cpp), Gt_max is the maximum measured concentration of Gt that is used for the normalization
    
    float Normalization = 3./Gt_max;                       // Normalization factor for the profiles: normalize such that highest maximum is at 3 (corresponds roughly to the concentration levels of ComplementaryGradients and 2GaussianPeaks)
    for (short i=0; i<N_L; ++i) {
        this->TargetConcentrationProfiles[0][i] = Normalization * Kni[int(i*999./(N_L-1))];       // transform the gapgene arrays of length 1000 (max index = 999) to the TargetProfile vectors of length N_L that will be used here
        this->TargetConcentrationProfiles[1][i] = Normalization * Kr[int(i*999./(N_L-1))];
        this->TargetConcentrationProfiles[2][i] = Normalization * Gt[int(i*999./(N_L-1))];
        this->TargetConcentrationProfiles[3][i] = Normalization * Hb[int(i*999./(N_L-1))];
    }
    this->CalculateTrivialFitnessValue();
};


template<short N_L>
FitnessEvaluation_Polarisation<N_L>::FitnessEvaluation_Polarisation(float L) : FitnessEvaluation_MatchProfiles<N_L, 1>(L,true)
{
    for (short i=0; i<N_L; ++i) {       // define a Theta-Step function that switches from 0 to 1 at L/2
        if(i<N_L/2)
            this->TargetConcentrationProfiles[0][i] = 0;
        else
            this->TargetConcentrationProfiles[0][i] = InitialConcentration_Max/2;
    }
    this->CalculateTrivialFitnessValue();      // calculate the trivial fitness value form the Target Concentration Profile
};


template<short N_L>
FitnessEvaluation_Polarisation_CoarseGrained<N_L>::FitnessEvaluation_Polarisation_CoarseGrained(float L) : FitnessEvaluation_MatchProfiles<N_L, 1>(L, true)
{
    this->right_boundaries = {N_L/2, N_L};
    this->differences = {2};

    this->CalculateTrivialFitnessValue();      // calculate the trivial fitness value form the Target Concentration Profile
};

/*
 template<short N_L>
 FitnessEvaluation_2GaussianPeaks_Scaling<N_L>::FitnessEvaluation_2GaussianPeaks_Scaling(float L1, float L2) :
 dx(L/N_L)
 {
 
 };
 
 
 template<short N_L>
 float FitnessEvaluation_2GaussianPeaks<N_L>::operator () (Individual& individual)
 {
 ReactionDiffusionDynamics<N_L> RDD(individual, dx);
 individual.Fitness = FitnessEvaluator(individual, RDD);
 return individual.Fitness;
 };
 
 
 template<short N_L>
 float FitnessEvaluation_2GaussianPeaks<N_L>::operator () (Individual& individual, ReactionDiffusionDynamics<N_L>& RDD)
 {
 return FitnessEvaluator(individual, RDD);
 };
 */




