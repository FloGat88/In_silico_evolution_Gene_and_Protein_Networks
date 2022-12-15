//
//  Stability_Analysis.h
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 20/04/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#ifndef Stability_Analysis_h
#define Stability_Analysis_h

#include "Individual.h"
#include "Dynamics.hpp"
#include <boost/numeric/odeint.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Eigenvalues>
//#include <levmar.h>
// #include <libraries/Levenberg_Marquardt_Algorithm_levmar-2.6>

using namespace Eigen;


/*
template<class Doub>
void Dynamics_Fwd(Doub* c, Doub* y, int N1, int N2, void* adata) {
    static_cast< ReactionDiffusionDynamics<N_discretisation>* >(adata)->ReactionDynamics_WellMixed<Doub>(c, y, N1, N2);
}

template<class Doub>
void Jac_Fwd(Doub* c, Doub* y, int N1, int N2, void* adata) {
    static_cast< ReactionDiffusionDynamics<N_discretisation>* >(adata)->Jacobian_WellMixed<Doub>(c, y, N1, N2);
}
*/


struct ReactionDynamics_Jacobian_WellMixed
{
    vector<ode1>& OdeSyst1;
    vector<short>& OdeSyst1_numbers;
    vector<ode2>& OdeSyst2;
    vector<short>& OdeSyst2_numbers;
    int N_elems;
    int N_constraints;
    Matrix<double,Dynamic,Dynamic> A;         // constraints due to mass conservation
    Map<VectorXf> b;         // total masses
    
    ReactionDynamics_Jacobian_WellMixed(ReactionDiffusionDynamics<N_discretisation>& RDD,vector<vector<short>>& SpeciesObject, vector<float>& TotalMasses) : OdeSyst1(RDD.OdeSyst1),OdeSyst1_numbers(RDD.OdeSyst1_numbers),OdeSyst2(RDD.OdeSyst2),OdeSyst2_numbers(RDD.OdeSyst2_numbers),N_elems(RDD.N_elems),N_constraints(SpeciesObject.size()),b(TotalMasses.data(),TotalMasses.size())
    {
        A.setZero(SpeciesObject.size(),N_elems);
        for(short i=0; i<SpeciesObject.size(); ++i)
            for(short ind : SpeciesObject[i])
                A(i,ind) += 1.;
    };
    
    int inputs() {return N_elems;}
    int values() {return N_elems+N_constraints;}
        
    int operator()(const VectorXd &x, VectorXd &F) const
    {
        
        short elem_index,k;
        short m = 0, n = 0;
        
        F.setZero(N_elems+N_constraints);
        for (elem_index = 0; elem_index!=N_elems; ++elem_index)  {
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[elem_index];++k,++m)                                     // OdeSyst1
                F[elem_index] += OdeSyst1[m].rate * x[OdeSyst1[m].ind];
            
            for (k=0;k!=OdeSyst2_numbers[elem_index];++k,++n)                                   // OdeSyst2
                F[elem_index] += OdeSyst2[n].rate * x[OdeSyst2[n].ind1] * x[OdeSyst2[n].ind2];
        }
        F.bottomRows(N_constraints) = A*x - b.cast<double>();
        return 0;
    }
    
    int df(const VectorXd &x, MatrixXd &J) const
    {
        short elem_index,k;
        short m = 0, n = 0;
        
        J.setZero(N_elems+N_constraints,N_elems);          // initialize J with 0
        for (elem_index = 0; elem_index!=N_elems; ++elem_index)  {
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[elem_index];++k,++m)                                     // OdeSyst1
                J(elem_index,OdeSyst1[m].ind) += OdeSyst1[m].rate;
            
            for (k=0;k!=OdeSyst2_numbers[elem_index];++k,++n)   {                                  // OdeSyst2
                J(elem_index,OdeSyst2[n].ind1) += OdeSyst2[n].rate * x[OdeSyst2[n].ind2];
                J(elem_index,OdeSyst2[n].ind2) += OdeSyst2[n].rate * x[OdeSyst2[n].ind1];
            }
        }
        J.bottomRows(N_constraints) = A;
        return 0;
    }
};

template<short N_L>
short Stability_Analysis(Individual& individual, ReactionDiffusionDynamics<N_L>& RDD, float* Fitness=NULL)  {
    using namespace boost::numeric::odeint;
    using namespace std::placeholders;
    typedef double Doub;
    constexpr double PI = 3.14159265359;
    
    constexpr double T_precond = T_max/3;   ///6;    // !!!!!!!!!!!!!
    short N_elems = individual.Get_NumberElements();
    short N_constraints = individual.Get_NumberElements(Prot);
    int info;           // output of eige::Levenberg-Marquardt
    
    if(N_elems==0)  {
        cout << "number of elements is zero\n";
        if(Fitness!=NULL)
            *Fitness=-1000.;
        return -1;
    }
// create initial concentration vector
    typedef boost::numeric::ublas::vector<Doub> vector_type;
    vector_type concentrations_wellmixed(N_elems,0.);
//    vector<Doub> concentrations_wellmixed(N_elems, 0);
    individual.FillIn_InitialHomogeneousProteinConcentrations(concentrations_wellmixed, 1, false);
    
// precondition by simulating well mixed system for time T_precond
    /*
// (explicit solver)
    typedef runge_kutta_dopri5< vector_type > stepper_type_explicit;
//    function<void(const vector<Doub>&, vector<Doub>&, const double)>   // boost::ref(RDD.WellMixed(const vector<Doub>& x, vector<Doub>& dxdt, const double t)),
    auto fp=bind(&ReactionDiffusionDynamics<N_L>::template WellMixed<vector_type>,RDD,_1,_2,_3);
    size_t steps = integrate_adaptive( make_controlled<stepper_type_explicit>( 1.0e-6 , 1.0e-4 ), fp,  concentrations_wellmixed, 0.0, T_precond, 1.0e-4);
    cout << "Steps needed: " << steps << endl;
*/
    // (implicit solver)
    auto fp=bind(&ReactionDiffusionDynamics<N_L>::template WellMixed<vector_type>,RDD,_1,_2,_3);
    auto Jp=bind(&ReactionDiffusionDynamics<N_L>::template Jacobian_WellMixed<vector_type,boost::numeric::ublas::matrix<Doub>>,RDD,_1,_2,_3,_4);
    size_t steps = integrate_adaptive( make_controlled< rosenbrock4<Doub> >( 1.0e-6 , 1.0e-4 ), make_pair(fp, Jp), concentrations_wellmixed, 0.0, T_precond, 1.0e-2);
    cout << "Steps needed: " << steps << endl;
    
// find fixed point by least-square-minimization via Levenberg-Marquardt
    VectorXd x(N_elems);             // transform concentrations_well mixed to type VectorXd
    for(short i=0; i<N_elems; ++i)
        x[i] = double(concentrations_wellmixed[i]);
    individual.Make_SpeciesObject();
    ReactionDynamics_Jacobian_WellMixed RD_Jac_WellMixed(RDD, individual.SpeciesObject, individual.TotalMasses);
    LevenbergMarquardt<ReactionDynamics_Jacobian_WellMixed> lm(RD_Jac_WellMixed);
    info = lm.minimize(x);
    cout << "Info: " << info << "  Iterations: " << lm.iter << "  function evaluations: " << lm.nfev << "  jacobi evaluations: " << lm.njev << endl;
//    if(lm.fvec.norm() >= 1e-8)  {
    if(lm.fnorm >= 1e-8)  {
        cout << "Fixed point of the reaction system could not be found\n";
        if(Fitness!=NULL)
            *Fitness = -1000.;
        return -1;
    }
/*
    cout << x << endl << endl << endl;   // !!!!!
    VectorXd g;   //!!!
    RD_Jac_WellMixed(x,g);  //!!!
    cout << g << endl;  //!!!
    MatrixXd popo;  //!!!
    RD_Jac_WellMixed.df(x,popo);        //!!!
    cout << popo - lm.fjac << endl;     //!!!
*/
    
// determine eigenvalues of J(x)-k^2*Diff
    VectorXd Diffusion(N_elems);
    for(short i=0; i<N_elems; ++i)
        Diffusion[i] = (RDD.states[i]==cytosolic)  ?  D_cytosolic : D_membrane;
    
    MatrixXd jacobian;          // attention: lm.fjac does not contain the correct jacobian at the minimum x. Probably the algorithm stops if it realizes that the function has norm very close to 0 (which then must be a minimum) without calculating the jacobian any more. Therefore we have to recalculate the jacobian at minimum x here again using RD_Jac_WellMixed
    RD_Jac_WellMixed.df(x,jacobian);
//    jacobian = jacobian.topRows(N_elems);
    
    auto MaxEigenvalueRealPart = [&jacobian, &Diffusion, N_elems](double k)  {
        EigenSolver<MatrixXd> es;
        MatrixXd Jac_Dk2 = jacobian.topRows(N_elems);
        Jac_Dk2.diagonal() -= Diffusion*pow(k,2);
        es.compute(Jac_Dk2, false);
        cout << es.info() << endl;
//        cout << es.eigenvalues() << endl << endl << endl;
        return es.eigenvalues().real().maxCoeff();
    };
    
    double k = PI/L_small;
    double maxEV_RealPart = MaxEigenvalueRealPart(k);
    if(Fitness!=NULL)
        *Fitness = maxEV_RealPart;
    cout << "maximum real part: " << maxEV_RealPart << endl;
    return 0;
}




























#endif /* Stability_Analysis_h */
