//
//  Dynamics.hpp
//  EmbrionicDevelopment
//
//  Created by Florian Gartner on 21/04/17.
//  Copyright © 2017 Florian Gartner. All rights reserved.
//

#ifndef Dynamics_hpp
#define Dynamics_hpp


#include "Individual.h"
//#include "LamSimAnn_Interface.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>           // for input and output of ublas vectors and matrices
//#include <boost/numeric/ublas/storage.hpp>


// #define state_type_matrix boost::numeric::ublas::matrix< float, boost::numeric::ublas::column_major , boost::numeric::ublas::bounded_array<float,N_tot*N_L> >
// #define state_type_array array<float,N_tot*N_L>

using namespace std;


template<short N_L>
class ReactionDiffusionDynamics   {   // Functor similar as previous imlementations that is however assumed to work more efficiently as, when constructed, it copies all relevant entries of OdeSyst1 and -2 in a compact one dimensional array format that is supposed to provide faster access.
    friend struct ReactionDynamics_Jacobian_WellMixed;
private:
    short N_elems;              // saves the number of elements in individual
    float D_cytosolic_dxdx;     // cytosolic diffusion constant divided by dx^2
    float D_membrane_dxdx;      // membranic diffusion constant divided by dx^2
public:
    vector<short> states;     // saves the diffusion constants of only the existing elements in a row (leaves out all the -1 in individual.DiffConsts); only relevnat if MutateDiffusionConstants == true, otherwise D_dxdx which is fixed in the constructor is the universal diffusion constnat for all elements
    vector<ode1> OdeSyst1;
    vector<short> OdeSyst1_numbers;
    vector<ode2> OdeSyst2;
    vector<short> OdeSyst2_numbers;
private:
    vector< vector< float* >> RatesDistributor;   // for the interface to LamSimAnn: vector of vectors that in each row composes the addresses of all rate entries in OdeSyst1 and OdeSyst2 that belong to the same reaction/arrow so that RatesDistributor[i][:] is a vector of addresses to equal rates that belong to reaction i. A paramter vector of size RatesDistributor.size() can thereby be used to modify all the entries in OdeSyst1 and -2 by distributing it via RatesDistributor.
    Individual& individual;                       // for convey_parameters to hand optimized parameters over to individual
    
    //    bool jacobian_initialized;                      // bool that states whether the Jacobian has already been called for this system (true) and has thus been initialized once or if not (false). In the former case it is only necessary to update the sparse Jacobian matrix for OdeSyst2 because the entries of the Jacobian deriving from OdeSyst1 are independent of the concentrations and therefore remain unchanged throughout the simulation via the implicit solver. Jacobian_initialized is initialized in the constructor as 'false' and set to 'true' after the first call of 'Jacobian'.

///////////////////////////////////////////////// Constructor ////////////////////////////////////////////////////////////////////////////////
public:
    ReactionDiffusionDynamics(Individual& individual, float dx = 0, bool TransferDirectly = true) : individual(individual)
    //, jacobian_initialized(false)    // Constructor: leave out the arguments dx and D in the well-mixed case, leave out only D if mutation of diffusion constants is switched on (then the diffusion constants are taken from DiffConsts)
    {
        auto Transfer_OdeSyst = [&]()  {
            // resize the vectors to their appropriate size
            OdeSyst1_numbers.resize(N_elems);
            OdeSyst2_numbers.resize(N_elems);
            OdeSyst1.clear();
            OdeSyst2.clear();
            OdeSyst1.reserve(individual.Get_NumberOde1());
            OdeSyst2.reserve(individual.Get_NumberOde2());
            
            for (short i=0; i<N_elems; ++i) {
                OdeSyst1_numbers[i] = individual.OdeSyst1[i].size();
                for (short j=0; j<OdeSyst1_numbers[i]; ++j)
                    OdeSyst1.push_back(individual.OdeSyst1[i][j]);
            }
            assert(OdeSyst1.size()==individual.Get_NumberOde1());
            
            for (short i=0; i<N_elems; ++i) {
                OdeSyst2_numbers[i] = individual.OdeSyst2[i].size();
                for (short j=0; j<OdeSyst2_numbers[i]; ++j)
                    OdeSyst2.push_back(individual.OdeSyst2[i][j]);
            }
            assert(OdeSyst2.size()==individual.Get_NumberOde2());
        };
        
        N_elems = individual.Get_NumberElements();
        D_cytosolic_dxdx = D_cytosolic/pow(dx,2);
        D_membrane_dxdx = D_membrane/pow(dx,2);
        
        if(TransferDirectly)
            individual.Transfer_directly_to_OdeSystCompressed(states, OdeSyst1, OdeSyst1_numbers, OdeSyst2, OdeSyst2_numbers);
        else  {
            individual.Make_OdeSyst(states);
            Transfer_OdeSyst();
        }
    };
    
////////////////////////////////// Dynamics Provider for Ode Solver  ///////////////////////////////////////////////////////////
///////////////////////////////////////////////// Dynamics (spatial) //////////////////////////////////////////////
    template<class vector_type>             // template for the vector type of x and dxdt. The explicit Runge Kutta solver from Odeint needs boost::array<float, N_L * N_tot> whereas the implicit solver from Numerical Recipes needs a mtl4 vector type
    void operator () (const vector_type& x, vector_type& dxdt, const double )  // !!!!!! put this into .tpp file again and outcomment t again   // utrafast: supposedly the fastest provider for the dynamics
    {
        short elem_index,i,j,j1,j2,j3,k;
        float rate;
        short m = 0, n = 0;
        
        for (i = 0, elem_index = 0; i!=N_elems; ++i, elem_index += N_L)  {
            
            // diffusion dynamics
            float D_dxdx_loc = (states[i]==cytosolic)  ?  D_cytosolic_dxdx : D_membrane_dxdx;
            dxdt[elem_index] = D_dxdx_loc * (- x[elem_index] + x[elem_index+1]);        // left boundary term
            for (j=elem_index+1; j!=elem_index+N_L-1; ++j)                                           // bulk terms
                dxdt[j] = D_dxdx_loc * (x[j-1] - 2*x[j] + x[j+1]);
            dxdt[j] = D_dxdx_loc * (x[j-1] - x[j]);  // right boundary term; j is now equal to elem_index+N_L-1 after the for loop has finished
            
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[i];++k,++m)                                     // OdeSyst1
                for (j1=elem_index, j2=OdeSyst1[m].ind*N_L, rate=OdeSyst1[m].rate; j1!=elem_index+N_L; ++j1, ++j2)
                    dxdt[j1] += rate * x[j2];
            
            for (k=0;k!=OdeSyst2_numbers[i];++k,++n)                                   // OdeSyst2
                for (j1=elem_index, j2=OdeSyst2[n].ind1*N_L, j3=OdeSyst2[n].ind2*N_L, rate=OdeSyst2[n].rate; j1!=elem_index+N_L; ++j1, ++j2, ++j3)
                    dxdt[j1] += rate * x[j2] * x[j3];
        }
    };
    
/////////////////////////////////////////////////// Dynamics (well-mixed) ///////////////////////////////////////////////
    template<class vector_type>
    void WellMixed (const vector_type &x , vector_type &dxdt , const double )     // dynamics provider for well-mixed systems, i.e. N_L = 1
    {
        short elem_index,k;
        short m = 0, n = 0;
        
        for (elem_index = 0; elem_index!=N_elems; ++elem_index)  {
            dxdt[elem_index] = 0;
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[elem_index];++k,++m)                                     // OdeSyst1
                dxdt[elem_index] += OdeSyst1[m].rate * x[OdeSyst1[m].ind];
            
            for (k=0;k!=OdeSyst2_numbers[elem_index];++k,++n)                                   // OdeSyst2
                dxdt[elem_index] += OdeSyst2[n].rate * x[OdeSyst2[n].ind1] * x[OdeSyst2[n].ind2];
        }
    };
    
/////////////////////////////////////////////////// Jacobian (well-mixed) ///////////////////////////////////////////////
    template<class vector_type, class matrix_type>        // provider of the Jacobian for odeint (for implicit steppers)
    void Jacobian_WellMixed(const vector_type &x, matrix_type &J , const double /* t */ , vector_type &dfdt )
    {
        short elem_index,k;
        short m = 0, n = 0;
        
        fill(J.begin1(), J.end1(), 0.0);              // initialize J with 0
        for (elem_index = 0; elem_index!=N_elems; ++elem_index)  {
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[elem_index];++k,++m)                                     // OdeSyst1
                J(elem_index,OdeSyst1[m].ind) += OdeSyst1[m].rate;
            
            for (k=0;k!=OdeSyst2_numbers[elem_index];++k,++n)   {                                  // OdeSyst2
                J(elem_index,OdeSyst2[n].ind1) += OdeSyst2[n].rate * x[OdeSyst2[n].ind2];
                J(elem_index,OdeSyst2[n].ind2) += OdeSyst2[n].rate * x[OdeSyst2[n].ind1];
            }
        }
        fill_n(dfdt.begin(),N_elems, 0.0);
    };
    
/////////////////////////////////////////////////// Jacobian (spatial) ///////////////////////////////////////////////
    /*
     // functions important for the implicit solver
     void Jacobian(const typename StepperBase<N_L*N_tot>::vector_type& x, typename StepperBase<N_L*N_tot>::matrix_type& jacobian)
     {
     static typename StepperBase<N_L*N_tot>::matrix_type basis_jacobian(N_L*N_tot,N_L*N_tot);
     
     if(jacobian_initialized==false) {
     basis_jacobian = 0;
     
     mtl::mat::inserter<typename StepperBase<N_L*N_tot>::matrix_type> ins(basis_jacobian,15);
     
     short elem_index,i,j,j1,j2,k;
     float rate;
     short d = 0, m = 0;
     for (i = 0, elem_index = 0; i!=N_tot; ++i, elem_index += N_L)  {
     
     if (OdeSyst1_numbers[i] != 0)  {                               // note that if OdeSyst1[i].size() = 0 then also OdeSyst2[i].size() = 0
     
     // dynamics due to diffusion
     if ( i<2*N_P+N_GPC && (i>=N_P && i<2*N_P) ) {
     float D_dxdx_loc;                   // use a local variable for the diffusion constant since with compiler optimization (that will optimize it away) this is supposedly faster (if DiffConsts_dxdx is used) than assigning the constants each time to the classe's D_dxdx where the assignment really has to be executed.
     if(MutateDiffusionConstants)                       //if mutation of diffusion constants is switched on, use one entry of DiffConsts_dxdx after the other for the diffusion constant (save this on D_dxdx) DiffConsts_dxdx contains only the diffusin constants of the existing elements as these lines are only executed for existing elements (and only then is the index d incremented) (if mutation of diffusion constants is off, leave D_dxdx unchanged so that its universal constructor-defined value is used instead)
     D_dxdx_loc = DiffConsts_dxdx[d++];
     else
     D_dxdx_loc = D_dxdx;
     
     ins[elem_index][elem_index] = -D_dxdx_loc;              // left boundary term
     ins[elem_index][elem_index+1] = D_dxdx_loc;
     for (j=elem_index+1; j!=elem_index+N_L-1; ++j)   {      // bulk terms
     ins[j][j-1] = D_dxdx_loc;
     ins[j][j] = -2 * D_dxdx_loc;
     ins[j][j+1] = D_dxdx_loc;
     }
     ins[j][j-1] = D_dxdx_loc;                               // right boundary term; j is now equal to elem_index+N_L-1 after the for loop has finished
     ins[j][j] = -D_dxdx_loc;
     }
     
     // reaction dynamics of OdeSyst1
     for (k=0;k!=OdeSyst1_numbers[i];++k,++m)
     for (j1=elem_index, j2=OdeSyst1[m].ind*N_L, rate=OdeSyst1[m].rate; j1!=elem_index+N_L; ++j1, ++j2)
     ins[j1][j2] += rate;
     }
     }
     jacobian_initialized = true;
     }
     
     short elem_index,i,j1,j2,j3,k, n=0;
     float rate;
     
     jacobian = basis_jacobian;
     mtl::mat::inserter<typename StepperBase<N_L*N_tot>::matrix_type> ins(jacobian,10);
     for (i = 0, elem_index = 0; i!=N_tot; ++i, elem_index += N_L)
     for (k=0;k!=OdeSyst2_numbers[i];++k,++n)                                   // OdeSyst2
     for (j1=elem_index, j2=OdeSyst2[n].ind1*N_L, j3=OdeSyst2[n].ind2*N_L, rate=OdeSyst2[n].rate; j1!=elem_index+N_L; ++j1, ++j2, ++j3)     {
     ins[j1][j2] += rate * x[j3];
     ins[j1][j3] += rate * x[j2];
     }
     };
     */
    
/////////////////////////////////////////////////// Interface to Levenberg-Marquardt //////////////////////////////////////////
/////////////////////////////////////////////////// Dynamics and Jacobian (well-mixed) /////////////////////////////////////////
/*    template<class Doub>
    void ReactionDynamics_WellMixed (Doub* c, Doub* y, int N1, int N2)
    {
        short elem_index,k;
        short m = 0, n = 0;
        
        for (elem_index = 0; elem_index!=N2; ++elem_index)  {
            y[elem_index] = 0;
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[elem_index];++k,++m)                                     // OdeSyst1
                y[elem_index] += OdeSyst1[m].rate * c[OdeSyst1[m].ind];
            
            for (k=0;k!=OdeSyst2_numbers[elem_index];++k,++n)                                   // OdeSyst2
                y[elem_index] += OdeSyst2[n].rate * c[OdeSyst2[n].ind1] * c[OdeSyst2[n].ind2];
        }
    };
    
    template<class Doub>
    void Jacobian_WellMixed (Doub* c, Doub* J, int N1, int N2)
    {
        short elem_index,k;
        short m = 0, n = 0;
        
        fill_n(J, N1*N2, (Doub)0.0);          // initialize J with 0
        for (elem_index = 0; elem_index!=N2; ++elem_index)  {
            short row_index = N1*elem_index;
            // reaction dynamics
            for (k=0;k!=OdeSyst1_numbers[elem_index];++k,++m)                                     // OdeSyst1
                J[row_index + OdeSyst1[m].ind] += OdeSyst1[m].rate;
            
            for (k=0;k!=OdeSyst2_numbers[elem_index];++k,++n)   {                                  // OdeSyst2
                J[row_index + OdeSyst2[n].ind1] += OdeSyst2[n].rate * c[OdeSyst2[n].ind2];
                J[row_index + OdeSyst2[n].ind2] += OdeSyst2[n].rate * c[OdeSyst2[n].ind1];
            }
        }
    };
*/
/////////////////////////////////////////////////// Interface to LamSimAnn /////////////////////////////////////////////////////
     // functinos enabling the interface to LamSimAnn
/*
     void make_RatesDistributor(vector<double>& start_par);     // compose the RatesDistributor object from OdeSyst1 and 2 (searching for equal rates therein) hand over a vector reference start_par in which the initial parameter set that is extracted from OdeSyst1 and 2 while composing RatesDistributor is saved.
     int get_NumberOfIndependentParameters()  {return RatesDistributor.size();};   // returns RatesDistributor.size(), i.e the number of paramters that completely describe the system that is needed e.g in LamSimAnn_Energy. Only returns a reasonable value if the RatesDistributor object has been created before by calling make_RatesDistributor().
     
     void modify_OdeSyst(double* par);       // hand over a vector of parameters (resp. pointer to initial element of a block of size RatesDistributor.size(), the way in which LamSimAnn deals with the parameters) that fully describes the Ode_system of the Individual under consideration and that with the help of RatesDistributor is distributed to OdeSyst1 and -2.
     void modify_OdeSyst(vector<double> par);   // same as above just with vector input for the final optimized parameterset that is given by LamSimAnn as a vector of doubles. The optimized parameters are again used to modify OdeSyst1 and -2 that will then be conveyed to OdeSyst1 and -2 in Individual via convey_parameters();
     
     void convey_parameters(); // conveys the parameters in OdeSyst1 and -2 in ReactionDiffusionDynamics to OdeSyst1 and -2 in Individual.
*/
};





/*
 template<short N_L>
 void ReactionDiffusionDynamics<N_L>::make_RatesDistributor(vector<double>& start_par)  {
     if(start_par.size()!=0)
        cerr << "ERROR in make Rates_Distributor: delivered start_par is not empty\n";
 
     short i,j;
     for(i=0; i<OdeSyst1.size(); ++i)   {
        if(OdeSyst1[i].ind==-1)             // if ind ==-1 i.e. all non-empty entries of OdeSyst1 are already through, break up
            break;
         for(j=0; j<RatesDistributor.size(); ++j)
             if( abs(OdeSyst1[i].rate) == abs(*(RatesDistributor[j][0]))  )  {
                 RatesDistributor[j].push_back( &OdeSyst1[i].rate );
                 break;
             }
         if(j==RatesDistributor.size())  {
             start_par.push_back( abs(OdeSyst1[i].rate) );
             RatesDistributor.push_back({&OdeSyst1[i].rate});
         }
     }
 
     short RatesDistributor_size1 = RatesDistributor.size();
     for(i=0; i<OdeSyst2.size(); ++i)   {
         if(OdeSyst2[i].ind1==-1)             // if ind1 ==-1 i.e. all non-empty entries of OdeSyst2 are already through, break up
             break;
         for(j=RatesDistributor_size1; j<RatesDistributor.size(); ++j)
             if( abs(OdeSyst2[i].rate) == abs(*(RatesDistributor[j][0])) )  {
                 RatesDistributor[j].push_back( &OdeSyst2[i].rate );
                 break;
             }
         if(j==RatesDistributor.size())  {
             start_par.push_back( abs(OdeSyst2[i].rate) );
             RatesDistributor.push_back({&OdeSyst2[i].rate});
         }
     }
 };


 template<short N_L>
 void ReactionDiffusionDynamics<N_L>::modify_OdeSyst(double* par)  {
     for(short i=0; i<RatesDistributor.size(); ++i)  {
         for(short j=0; j<RatesDistributor[i].size(); ++j)
             if(*RatesDistributor[i][j]<0)
                 *RatesDistributor[i][j] = -par[i];
             else
                 *RatesDistributor[i][j] = par[i];
     }
 }


 template<short N_L>
 void ReactionDiffusionDynamics<N_L>::modify_OdeSyst(vector<double> par)  {
     for(short i=0; i<RatesDistributor.size(); ++i)  {
         for(short j=0; j<RatesDistributor[i].size(); ++j)
             if(*RatesDistributor[i][j]<0)
                 *RatesDistributor[i][j] = -par[i];
             else
                 *RatesDistributor[i][j] = par[i];
     }
 }

template<short N_L>
void ReactionDiffusionDynamics<N_L>::convey_parameters()  {
    for(short i=0, k=0; i<N_tot; ++i)           // convey the optimized rates from RDD.OdeSyst1 to individual.OdeSyst1
        for(short j=0; j<OdeSyst1_numbers[i]; ++j, ++k)
            individual.OdeSyst1[i][j].rate = OdeSyst1[k].rate;
 
    for(short i=0, k=0; i<N_tot; ++i)           // convey the optimized rates from RDD.OdeSyst2 to individual.OdeSyst2
        for(short j=0; j<OdeSyst2_numbers[i]; ++j, ++k)
            individual.OdeSyst2[i][j].rate = OdeSyst2[k].rate;
}
*/




























#endif /* Dynamics_hpp */
