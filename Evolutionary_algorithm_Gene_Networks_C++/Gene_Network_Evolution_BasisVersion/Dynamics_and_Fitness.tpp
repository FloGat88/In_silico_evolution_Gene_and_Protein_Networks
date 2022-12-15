//  Dynamics_and_Fitness.tpp

//#include "Dynamics_and_Fitness.h"
// #include "Individual.h"
//#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;
// using namespace boost::numeric::ublas;

/* short N_L1 = (short)(L1/delta_x);
 short N_L2 = (short)(L2/delta_x);   */




/*
reaction_diffusion_dynamics::reaction_diffusion_dynamics(Individual& Individual, float D, float dx): OdeSyst1(Individual.OdeSyst1), OdeSyst2(Individual.OdeSyst2), D_dxdx(D/pow(dx,2)) {};  // Constructor
    
template<short N_L>
void reaction_diffusion_dynamics::operator () (const array< array<float, N_L>, N_tot > &x , array< array<float, N_L>, N_tot > &dxdt , const double)
{
    short i,j,k;
    short size_OdeSyst1, size_OdeSyst2;
        
    for (i=0;i<N_tot;++i)  {
        size_OdeSyst1 = OdeSyst1[i].size();
        if (size_OdeSyst1 > 0)  {                               // note that if OdeSyst1[i].size() = 0 then also OdeSyst2[i].size() = 0
            size_OdeSyst2 = OdeSyst2[i].size();
                
            // dynamics due to diffusion
            if ( i<N_PPC && (i<N_P || i>=N_GPC) )                                    // for Genes and GPComp there is no diffusion, only initialize the corresponding entries in dxdt with 0
                for (j=0;j<N_L;++j)
                    dxdt[i][j] = 0;
            else {
                dxdt[i][0] = D_dxdx * (- x[i][0] + x[i][1]);        // left boundary term
                for (j=1;j<N_L-1;++j)                                           // bulk terms
                    dxdt[i][j] = D_dxdx * (x[i][j-1] - 2*x[i][j] + x[i][j+1]);
                dxdt[i][N_L-1] = D_dxdx * (x[i][N_L-2] - x[i][N_L-1]);  // right boundary term
            }
            // reaction dynamics
            for (k=0;k<size_OdeSyst1;++k)                                   // OdeSyst1
                for (j=0;j<N_L;++j)
                    dxdt[i][j] += OdeSyst1[i][k].rate * x[OdeSyst1[i][k].ind][j];
                
            for (k=0;k<size_OdeSyst2;++k)                           // OdeSyst2
                for (j=0;j<N_L;++j)
                    dxdt[i][j] += OdeSyst2[i][k].rate * x[OdeSyst2[i][k].ind1][j] * x[OdeSyst2[i][k].ind2][j];
        }
    }
}


template<short N_L>
void reaction_diffusion_dynamics::operator () (const state_type_matrix &x , state_type_matrix &dxdt , const double)
{
    short i,j,k;
    short size_OdeSyst1, size_OdeSyst2;
    
    for (i=0;i<N_tot;++i)  {
        size_OdeSyst1 = OdeSyst1[i].size();
        if (size_OdeSyst1 > 0)  {                               // note that if OdeSyst1[i].size() = 0 then also OdeSyst2[i].size() = 0
            size_OdeSyst2 = OdeSyst2[i].size();
            
            // dynamics due to diffusion
            if ( i<N_PPC && (i<N_P || i>=N_GPC) )                                    // for Genes and GPComp there is no diffusion, only initialize the corresponding entries in dxdt with 0
                for (j=0;j<N_L;++j)
                    dxdt(i,j) = 0;
            else {
                dxdt(i,0) = D_dxdx * (- x(i,0) + x(i,1));        // left boundary term
                for (j=1;j<N_L-1;++j)                                           // bulk terms
                    dxdt(i,j) = D_dxdx * (x(i,j-1) - 2*x(i,j) + x(i,j+1));
                dxdt(i,N_L-1) = D_dxdx * (x(i,N_L-2) - x(i,N_L-1));  // right boundary term
            }
            // reaction dynamics
            for (k=0;k<size_OdeSyst1;++k)                                   // OdeSyst1
                for (j=0;j<N_L;++j)
                    dxdt(i,j) += OdeSyst1[i][k].rate * x(OdeSyst1[i][k].ind,j);
            
            for (k=0;k<size_OdeSyst2;++k)                           // OdeSyst2
                for (j=0;j<N_L;++j)
                    dxdt(i,j) += OdeSyst2[i][k].rate * x(OdeSyst2[i][k].ind1,j) * x(OdeSyst2[i][k].ind2,j);
        }
    }
}
*/


/*
template<short N_L>
void reaction_diffusion_dynamics_fast::operator () (const array< array<float, N_L>, N_tot > &x , array< array<float, N_L>, N_tot > &dxdt , const double)
{
    short i,j,k;
    short m = 0,n = 0;
        
    for (i=0;i<N_tot;++i)  {
        if (OdeSyst1_numbers[i] > 0)  {                               // note that if OdeSyst1[i].size() = 0 then also OdeSyst2[i].size() = 0
                
            // dynamics due to diffusion
            if ( i<N_PPC && (i<N_P || i>=N_GPC) )                                    // for Genes and GPComp there is no diffusion, only initialize the corresponding entries in dxdt with 0
                for (j=0;j<N_L;++j)
                    dxdt[i][j] = 0;
            else {
                dxdt[i][0] = D_dxdx * (- x[i][0] + x[i][1]);        // left boundary term
                for (j=1;j<N_L-1;++j)                                           // bulk terms
                    dxdt[i][j] = D_dxdx * (x[i][j-1] - 2*x[i][j] + x[i][j+1]);
                dxdt[i][N_L-1] = D_dxdx * (x[i][N_L-2] - x[i][N_L-1]);  // right boundary term
            }
    
            // reaction dynamics
            for (k=0;k<OdeSyst1_numbers[i];++k,++m)                                     // OdeSyst1
                for (j=0;j<N_L;++j)
                    dxdt[i][j] += OdeSyst1[m].rate * x[OdeSyst1[m].ind][j];
                
            for (k=0;k<OdeSyst2_numbers[i];++k,++n)                                   // OdeSyst2
                for (j=0;j<N_L;++j)
                    dxdt[i][j] += OdeSyst2[n].rate * x[OdeSyst2[n].ind1][j] * x[OdeSyst2[n].ind2][j];
        }
    }
};


template<short N_L>
void reaction_diffusion_dynamics_fast::operator () (const state_type_matrix &x , state_type_matrix &dxdt , const double)
{
    short i,j,k;
    short m = 0,n = 0;
    
    for (i=0;i<N_tot;++i)  {
        if (OdeSyst1_numbers[i] > 0)  {                               // note that if OdeSyst1[i].size() = 0 then also OdeSyst2[i].size() = 0
            
            // dynamics due to diffusion
            if ( i<N_PPC && (i<N_P || i>=N_GPC) )                                    // for Genes and GPComp there is no diffusion, only initialize the corresponding entries in dxdt with 0
                for (j=0;j<N_L;++j)
                    dxdt(i,j) = 0;
            else {
                dxdt(i,0) = D_dxdx * (- x(i,0) + x(i,1));        // left boundary term
                for (j=1;j<N_L-1;++j)                                           // bulk terms
                    dxdt(i,j) = D_dxdx * (x(i,j-1) - 2*x(i,j) + x(i,j+1));
                dxdt(i,N_L-1) = D_dxdx * (x(i,N_L-2) - x(i,N_L-1));  // right boundary term
            }
            
            // reaction dynamics
            for (k=0;k<OdeSyst1_numbers[i];++k,++m)                                     // OdeSyst1
                for (j=0;j<N_L;++j)
                    dxdt(i,j) += OdeSyst1[m].rate * x(OdeSyst1[m].ind,j);
            
            for (k=0;k<OdeSyst2_numbers[i];++k,++n)                                   // OdeSyst2
                for (j=0;j<N_L;++j)
                    dxdt(i,j) += OdeSyst2[n].rate * x(OdeSyst2[n].ind1,j) * x(OdeSyst2[n].ind2,j);
        }
    }
};
*/




















template<short N_L>
reaction_diffusion_dynamics_fast<N_L>::reaction_diffusion_dynamics_fast(Individual& Individual, float D, float dx): D_dxdx(D/pow(dx,2))  // Constructor
{
    short i,j;
    short k = 0;
    for (i=0; i<N_tot; ++i) {
        OdeSyst1_numbers[i] = Individual.OdeSyst1[i].size();
        for (j=0; j<OdeSyst1_numbers[i]; ++j)
            OdeSyst1[k++] = Individual.OdeSyst1[i][j];
    }
    if (k > OdeSyst1_MaxSize) {
        cout << "ERROR: OdeSyst1 exceeds maximum size OdeSyst1_MaxSize" << '\n';
        exit(4);
    }

    k = 0;
    for (i=0; i<N_tot; ++i) {
        OdeSyst2_numbers[i] = Individual.OdeSyst2[i].size();
        for (j=0; j<OdeSyst2_numbers[i]; ++j)
            OdeSyst2[k++] = Individual.OdeSyst2[i][j];
    }
    if (k > OdeSyst2_MaxSize) {
        cout << "ERROR: OdeSyst2 exceeds maximum size OdeSyst2_MaxSize" << '\n';
        exit(4);
    }

};


template<short N_L>
void reaction_diffusion_dynamics_fast<N_L>::operator () (const boost::array<float, N_L * N_tot> &x , boost::array<float, N_L * N_tot> &dxdt , const double )
{
    size_t elem_index,i,j,j1,j2,j3,k;
    float rate;
    size_t m = 0,n = 0;
    
    for (i = 0, elem_index = 0; i!=N_tot; ++i, elem_index += N_L)  {
        
        if (OdeSyst1_numbers[i] == 0)          // if the size of the corresponding OdeSyst1 vector is 0 then only initialize the corresponding entries in dxdt with 0; this is however necessary otherwise dxdt if not initialized properly makes weird things
            for (j=elem_index; j!=elem_index+N_L; ++j)
                dxdt[j] = 0;
        else {                               // note that if OdeSyst1[i].size() = 0 then also OdeSyst2[i].size() = 0
            
            // dynamics due to diffusion
            if ( i<2*N_P+N_GPC && (i<N_P || i>=2*N_P) )           // for Genes and GPComp there is no diffusion, only initialize the corresponding entries in dxdt with 0
                for (j=elem_index; j!=elem_index+N_L; ++j)
                    dxdt[j] = 0;
            else {
                dxdt[elem_index] = D_dxdx * (- x[elem_index] + x[elem_index+1]);        // left boundary term
                for (j=elem_index+1; j!=elem_index+N_L-1; ++j)                                           // bulk terms
                    dxdt[j] = D_dxdx * (x[j-1] - 2*x[j] + x[j+1]);
                dxdt[j] = D_dxdx * (x[j-1] - x[j]);  // right boundary term; j is now equal to elem_index+N_L-1 after the for loop has finished
            }
            
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[i];++k,++m)                                     // OdeSyst1
                for (j1=elem_index, j2=OdeSyst1[m].ind*N_L, rate=OdeSyst1[m].rate; j1!=elem_index+N_L; ++j1, ++j2)
                    dxdt[j1] += rate * x[j2];
            
            for (k=0;k!=OdeSyst2_numbers[i];++k,++n)                                   // OdeSyst2
                for (j1=elem_index, j2=OdeSyst2[n].ind1*N_L, j3=OdeSyst2[n].ind2*N_L, rate=OdeSyst2[n].rate; j1!=elem_index+N_L; ++j1, ++j2, ++j3)
                    dxdt[j1] += rate * x[j2] * x[j3];
        }
    }
};





template <short N_L>
FitnessEvaluation_2GaussianPeaks<N_L>::FitnessEvaluation_2GaussianPeaks(float length): L(length), dx(L/N_L), sigma(L/6)                // Constructor
{
    float stdnormalization = 1/sqrt(2*M_PI*pow(sigma,2));
    for (short i=0; i<N_L; ++i) {
        TargetConcentrationProfile1[i] = 2 * stdnormalization * exp(-pow(i*dx - L/3.0,2)/(2*pow(sigma,2)) );
        TargetConcentrationProfile2[i] = 2 * stdnormalization * exp(-pow(i*dx - 2.0/3.0*L,2)/(2*pow(sigma,2)) );
    }
        TrivialFitnessValue_L1 = 0;
        for (short i=0; i<N_L; ++i)     // calculate the trivial fitness value as the L^2 norm of the Target Profiles
            TrivialFitnessValue_L1 += abs(TargetConcentrationProfile1[i]) + abs(TargetConcentrationProfile2[i]);
        TrivialFitnessValue_L1 = TrivialFitnessValue_L1 / N_L;
    
};
    

template<short N_L>
float FitnessEvaluation_2GaussianPeaks<N_L>::ProbeProfile_2GaussianPeaks_L2Norm (const float* ProbeProfile1, const float* ProbeProfile2)
{
    float profile_match_L2Norm = 0;
    for (short i=0; i<N_L; ++i)
        profile_match_L2Norm += pow( ProbeProfile1[i] - TargetConcentrationProfile1[i] , 2) + pow( ProbeProfile2[i] - TargetConcentrationProfile2[i] , 2);
    return sqrt(profile_match_L2Norm/N_L);    //Return the normed (average) L^2 metric distance between the probe profile and the Target Profile, d_{L^2}(ProbeProfile,TargetProfile) = sqrt(1/L * int_0^L(|ProbeProfile-TargetProfile|^2))
}


template<short N_L>
float FitnessEvaluation_2GaussianPeaks<N_L>::ProbeProfile_2GaussianPeaks_L1Norm (const float* ProbeProfile1, const float* ProbeProfile2)
{
    float profile_match_L1Norm = 0;
    for (short i=0; i<N_L; ++i)
        profile_match_L1Norm += abs(ProbeProfile1[i] - TargetConcentrationProfile1[i]) + abs(ProbeProfile2[i] - TargetConcentrationProfile2[i]);
    return profile_match_L1Norm/N_L;    //Return the normed (average) L^1 metric distance between the probe profile and the Target Profile, d_{L^1}(ProbeProfile,TargetProfile) = 1/L * int_0^L(|ProbeProfile-TargetProfile|)
}


template<short N_L>
float FitnessEvaluation_2GaussianPeaks<N_L>::operator () (Individual& individual, float D, double T_max)
{
    using namespace boost::numeric::odeint;
    
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
    reaction_diffusion_dynamics_fast<N_L> RDD(individual, D, dx);
    
    typedef runge_kutta_dopri5< boost::array<float, N_L*N_tot> > stepper_type;
    //typedef runge_kutta_cash_karp54< boost::array<float, N_L*N_tot> > stepper_type;
    size_t steps = integrate_adaptive( make_controlled<stepper_type>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, 0.0, T_max, 0.5 ); // matlab default error tolerances are (AbsTol = 1.0e-6, RelTol = 1.0e-3);; matlab also uses runge_kutta_dopri5 as its ode45; alternatively try runge_kutta_cash_karp54
    //size_t steps = integrate_const( runge_kutta4< boost::array<float, N_L*N_tot> >(), boost::ref(RDD) , concentrations , 0.0 , T_max , 0.1 );
    
    cout << steps << '\n';
    
    return -ProbeProfile_2GaussianPeaks_L1Norm( &concentrations[(N_P+individual.Essential_Proteins[1])*N_L], &concentrations[(N_P+individual.Essential_Proteins[2])*N_L] ) * 100./TrivialFitnessValue_L1 + 100.; // normalize the metric distance between the Probe Profile and the Target Profile to a fitness value in [0,100]
}







template <short N_L, short N_profiles>
FitnessEvaluation_MatchProfiles<N_L, N_profiles>::FitnessEvaluation_MatchProfiles(array<array<float, N_L>, N_profiles> TargetConcentrationProfiles, float length): L(length), dx(L/N_L), TargetConcentrationProfiles(TargetConcentrationProfiles)          // Constructor
{
/*    float sigma = L/6;
    float stdnormalization = 1/sqrt(2*M_PI*pow(sigma,2));
    for (short i=0; i<N_L; ++i) {
        TargetConcentrationProfiles[1][i] = 2 * stdnormalization * exp(-pow(i*dx - L/3.0,2)/(2*pow(sigma,2)) );
        TargetConcentrationProfiles[2][i] = 2 * stdnormalization * exp(-pow(i*dx - 2.0/3.0*L,2)/(2*pow(sigma,2)) );
    }
*/
        TrivialFitnessValue_L1 = 0;
        for (short i=0; i<N_profiles; ++i)
            for (short j=0; j<N_L; ++j)     // calculate the trivial fitness value as the L^2 norm of the Target Profiles
                TrivialFitnessValue_L1 += abs(TargetConcentrationProfiles[i][j]);
        TrivialFitnessValue_L1 = TrivialFitnessValue_L1 / N_L;
};


template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::ProbeProfiles_L2Norm (array<const float*, N_profiles> ProbeProfiles)
{
    float profile_match_L2Norm = 0;
    for (short i=0; i<N_profiles; ++i)
        for (short j=0; j<N_L; ++j)
            profile_match_L2Norm += pow( ProbeProfiles[i][j] - TargetConcentrationProfiles[i][j] , 2);
    return sqrt(profile_match_L2Norm/N_L);    //Return the normed (average) L^2 metric distance between the probe profile and the Target Profile, d_{L^2}(ProbeProfile,TargetProfile) = sqrt(1/L * int_0^L(|ProbeProfile-TargetProfile|^2))
}


template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::ProbeProfiles_L1Norm (array<const float*, N_profiles> ProbeProfiles)
{
    float profile_match_L1Norm = 0;
    for (short i=0; i<N_profiles; ++i)
        for (short j=0; j<N_L; ++j)
            profile_match_L1Norm += abs(ProbeProfiles[i][j] - TargetConcentrationProfiles[i][j]);
    return profile_match_L1Norm/N_L;    //Return the normed (average) L^1 metric distance between the probe profile and the Target Profile, d_{L^1}(ProbeProfile,TargetProfile) = 1/L * int_0^L(|ProbeProfile-TargetProfile|)
}


template<short N_L, short N_profiles>
float FitnessEvaluation_MatchProfiles<N_L, N_profiles>::operator () (Individual& individual, float D, double T_max, reaction_diffusion_dynamics_fast<N_L>& RDD)
{
    using namespace boost::numeric::odeint;
    
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
    
    array<const float*, N_profiles> ProbeProfiles;    // create the array ProbeProfiles with the addresses of the start positions of the concerning Essential Proteins (whos profiles are to be matched by the Algorithm) within concentrations
    for(short i=0; i<N_profiles; ++i)
        ProbeProfiles[i] = &concentrations[(N_P+individual.Essential_Proteins[i+1])*N_L];

    //reaction_diffusion_dynamics_fast<N_L> RDD(individual, D, dx);
    
    typedef runge_kutta_dopri5< boost::array<float, N_L*N_tot> > stepper_type;
    //typedef runge_kutta_cash_karp54< boost::array<float, N_L*N_tot> > stepper_type;
    size_t steps1 = integrate_adaptive( make_controlled<stepper_type>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, 0.0, T_max, 0.5 ); // matlab default error tolerances are (AbsTol = 1.0e-6, RelTol = 1.0e-3);; matlab also uses runge_kutta_dopri5 as its ode45; alternatively try runge_kutta_cash_karp54
    //size_t steps = integrate_const( runge_kutta4< boost::array<float, N_L*N_tot> >(), boost::ref(RDD) , concentrations , 0.0 , T_max , 0.1 );
    
    float fitness1 = -ProbeProfiles_L1Norm( ProbeProfiles ) * 100./TrivialFitnessValue_L1 + 100.; // normalize the metric distance between the Probe Profiles and the Target Profile to a fitness value in [0,100]
    
    double T_check = 20;           // time span after which a second check is done to ensure that the profile was indeed stationary
    size_t steps2 = integrate_adaptive( make_controlled<stepper_type>( 1.0e-6 , 1.0e-4 ), boost::ref(RDD), concentrations, T_max, T_max + T_check, 0.5 );   // integrate the system T_check time longer
    cout << steps1 + steps2 << '\n';       // display steps as the sum of the integration steps from both parts.
    
    float fitness2 = -ProbeProfiles_L1Norm( ProbeProfiles ) * 100./TrivialFitnessValue_L1 + 100.;    // check again for the desired profiles
    cout << "fitness1 = " << fitness1 << "    " << "fitness2 = " << fitness2 << '\n';
    return (fitness1 < fitness2 ? fitness1 : fitness2);     // return the final fitness value as the minimum of the two fitness values
}









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
FitnessEvaluation_ComplementaryGradients<N_L>::FitnessEvaluation_ComplementaryGradients(float L, float D, double T_max) :
D(D), T_max(T_max), dx(L/N_L)
{
    array<array<float, N_L>, 2> TargetConcentrationProfiles;
    float xi = L/2;       // decay length of the gradient
    for (short i=0; i<N_L; ++i) {
        TargetConcentrationProfiles[0][i] = 2 * exp(-i*dx/xi );
        TargetConcentrationProfiles[1][i] = 2 * exp(-(N_L-1-i)*dx/xi );
    }
    FitnessEvaluator = FitnessEvaluation_MatchProfiles<N_L, 2>(TargetConcentrationProfiles, L);
};


template<short N_L>
float FitnessEvaluation_ComplementaryGradients<N_L>::operator () (Individual& individual)
{
    reaction_diffusion_dynamics_fast<N_L> RDD(individual, D, dx);
    individual.Fitness = FitnessEvaluator(individual, D, T_max, RDD);
    return individual.Fitness;
};


template<short N_L>
float FitnessEvaluation_ComplementaryGradients<N_L>::operator () (Individual& individual, reaction_diffusion_dynamics_fast<N_L>& RDD)
{
    return FitnessEvaluator(individual, D, T_max, RDD);
};


