#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using namespace matlab::data;

class MexFunction : public matlab::mex::Function {
    ArrayFactory factory;

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        const TypedArray<double> t = inputs[0];
        const TypedArray<double> x = inputs[1];
        const TypedArray<double> OdeTransitions_rates = inputs[2];
        const TypedArray<double> OdeReactions_rates = inputs[3];
        const TypedArray<double> OdeDecays_rates = inputs[4];
        const TypedArray<short> OdeTransitions = inputs[5];
        const TypedArray<short> OdeReactions = inputs[6];
        const TypedArray<short> OdeDecays = inputs[7];
        const TypedArray<double> D_eff = inputs[8];
        const TypedArray<short> paraN = inputs[9];
        short N_elems = paraN[0];
        short N_L = paraN[1];
        TypedArray<double>  dxdt = factory.createArray<double>( {std::size_t(N_L*N_elems),1} );
            
// diffusion dynamics
        int i,j; 
        for(j=0,i=0; i<N_elems; ++i, ++j)
            dxdt[j] = D_eff[i]*(-x[j]+x[j+N_elems]);
        for(i=0; j<(N_L-1)*N_elems; ++j,++i)  {
             if (i==N_elems) i=0;
             dxdt[j] = D_eff[i]*(x[j-N_elems]-2*x[j]+x[j+N_elems]);
        }
        for(j=(N_L-1)*N_elems,i=0; i<N_elems; ++j,++i)
             dxdt[j] = D_eff[i]*(x[j-N_elems]-x[j]);  
                
// reaction dynamics
        auto x_=x.begin();
        auto dxdt_=dxdt.begin();
        for(; x_!=x.begin()+N_L*N_elems; x_+=N_elems,dxdt_+=N_elems)  {
            for(i=0; i<OdeTransitions_rates.getNumberOfElements(); ++i)   {      // contribution of Transitions
                double s = OdeTransitions_rates[i] * x_[OdeTransitions[0][i]];
                dxdt_[OdeTransitions[0][i]] -= s;
                dxdt_[OdeTransitions[1][i]] += s;
            }
            for(i=0; i<OdeReactions_rates.getNumberOfElements(); ++i)   {      // contribution of Reactions
                double s = OdeReactions_rates[i] * x_[OdeReactions[0][i]] * x_[OdeReactions[1][i]];
                dxdt_[OdeReactions[0][i]] -= s;
                dxdt_[OdeReactions[1][i]] -= s;
                dxdt_[OdeReactions[2][i]] += s;
            }
            for(i=0; i<OdeDecays_rates.getNumberOfElements(); ++i)   {       // contribution of Decays
                double s = OdeDecays_rates[i] * x_[OdeDecays[0][i]];
                dxdt_[OdeDecays[0][i]] -= s;
                dxdt_[OdeDecays[1][i]] += s;
                dxdt_[OdeDecays[2][i]] += s;
            }      
        }          
        outputs[0] = dxdt;
    };
};