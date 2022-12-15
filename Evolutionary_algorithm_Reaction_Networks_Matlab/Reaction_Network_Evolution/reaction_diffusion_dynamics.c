
#include "mex.h"
#include "math.h"
#include "mex.h"

// The computational routine 

// Leschi model
void dynamics(double t, double *x, double *Transtions_rates, double *Reactions_rates, double *Decays_rates, short *Transitions, short *Reactions, short *Decays, mwSize N, double *dxdt)
        // parameters = [mu, a, alpha, gamma0, delta, D]
{
    
    
    
    // dynamics of rho 
    double prefactor0 = parameters[5]/pow(parameters[1],2);
    dxdt[0] = parameters[3] - parameters[4]*x[0] + prefactor0*(x[N-1]-2*x[0]+x[1]) + noise_para[0];
    for (mwSize i=1; i<N-1; i++) 
        dxdt[i] = parameters[3] - parameters[4]*x[i] + prefactor0*(x[i-1]-2*x[i]+x[i+1]) + noise_para[i];
    dxdt[N-1] = parameters[3] - parameters[4]*x[N-1] + prefactor0*(x[N-2]-2*x[N-1]+x[0]) + noise_para[N-1];
  //  mexPrintf("\n%f  %f  %f %f \n%f  %f  %f  %f \n%f  %f  %f  %f\n", x[N-1], x[0], x[1], x[250], dxdt[N-1],dxdt[0],dxdt[1],dxdt[250],noise_para[N-1], noise_para[0],noise_para[1],noise_para[250]);
  //  mexEvalString("drawnow;");
// dynamics of h 
    double prefactor2 = parameters[2]/pow(parameters[1],2);  
   
    dxdt[N] = prefactor2*(x[2*N-1]-2*x[N]+x[N+1]) + x[0];
    for(mwSize i=N+1; i<2*N-1; ++i)
        dxdt[i] = prefactor2*(x[i-1]-2*x[i]+x[i+1]) + x[i-N];
    dxdt[2*N-1] = prefactor2*(x[2*N-2]-2*x[2*N-1]+x[N]) + x[N-1];
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    size_t N;                   /* size of matrix */
    double *output;              /* output matrix */
    
    N = mxGetM(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    output = mxGetPr(plhs[0]);
    dynamics(mxGetScalar(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]), mxGetPr(prhs[3]), mxGetM(prhs[1]), output);
}
