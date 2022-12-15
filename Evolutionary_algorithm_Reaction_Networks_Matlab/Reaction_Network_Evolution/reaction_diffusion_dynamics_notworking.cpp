#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;
    const TypedArray<double> OdeTransitions_rates;
    const TypedArray<double> OdeReactions_rates;
    const TypedArray<double> OdeDecays_rates;
    const TypedArray<short> OdeTransitions;
    const TypedArray<short> OdeReactions;
    const TypedArray<short> OdeDecays;
    const TypedArray<double> D_eff;
    const short N_elems;
    const short N_L;
   
public:
    MexFunction()  {
        mexLock();
        OdeTransitions_rates = inputs[2];
        OdeReactions_rates = inputs[3];
        OdeDecays_rates = inputs[4];
        OdeTransitions = inputs[5];
        OdeReactions = inputs[6];
        OdeDecays = inputs[7];
        D_eff = inputs[8];
        param_N = inputs[9];
        N_elems = param_N[0];
        N_L = param_N[1];
    };
    ~MexFunction() {
        mexUnlock();
    }
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        const TypedArray<double> t = inputs[0];
        const TypedArray<double> x = inputs[1];
        TypedArray<double>  dxdt = factory.createArray<double>( {N_L*N_elems,1} );
// diffusion dynamics
        int i,j; 
        for(j=0,i=0; i<N_elems; ++i, ++j)
            dxdt_[j] = D_eff[i]*(-x_[j]+x_[j+N_elems]);
        for(i=0; j<(N_L-1)*N_elems; ++j,++i)  {
             if (i==N_elems) i=0;
             dxdt[j] = D_eff[i]*(x[j-N_elems]-2*x[j]+x[j+N_elems]);
        }
        for(j=(N_L-1)*N_elems,i=0; i<N_elems; ++j,++i)
             dxdt[j] = D_eff[i]*(x[j-N_elems]-x[j]);      
        
// reaction dynamics
        double* x_, dxdt_;
        for(x_=&x[0], dxdt_=&dxdt[0]; x_!=&x[0]+N_L*N_elems; x_+=N_elems,dxdt_+=N_elems)  {
            for(i=0,j=0; i<sizeof(OdeTranstions_rates); ++i,j+=2)   {      // contribution of Transitions
                double s = OdeTransitions_rates[i] * x_[OdeTransitions[j]];
                dxdt_[OdeTransitions[j]] -= s;
                dxdt_[OdeTransitions[j+1]] += s;
            }
            for(i=0,j=0; i<sizeof(OdeReactions_rates); ++i,j+=3)   {      // contribution of Reactions
                double s = OdeReactions_rates[i] * x_[OdeReactions[j]] * x_[OdeReactions[j+1]];
                dxdt_[OdeReactions[j]] -= s;
                dxdt_[OdeReactions[j+1]] -= s;
                dxdt_[OdeReactions[j+2]] += s;
            }
            for(i=0,j=0; i<sizeof(OdeDecays_rates); ++i,j+=3)   {       // contribution of Decays
                double s = OdeDecays_rates[i] * x_[OdeDecays[j]];
                dxdt_[OdeDecays[j]] -= s;
                dxdt_[OdeDecays[j+1]] += s;
                dxdt_[OdeDecays[j+2]] += s;
            }
        }
        outputs[0] = dxdt;
    };

};