/*
 *  LamSimAnn.h
 *  SimAnn
 *
 *  Created by Patrick Hillenbrand on 5/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _LAMSIMANN_
#define _LAMSIMANN_

#include <vector>
#include "boost/random.hpp"

/**********************************************************************************************************************
*********** base class for energy function ****************************************************************************
**********************************************************************************************************************/

/*
 DEFINE YOUR TARGET FUNCTION (to be minimized by LamSimAnn)
 
 - The target function is defined in class that inherits publicly from Energy_base [declare: class YourClass : public Energy_base {...};]
 - Inheriting from virtual class (like that one) means that you must implement all functions defined in the virtual class (here four)
 - function 'Energy' should implement your target function and return its value for parameters *par (explained below)
 - LamSimAnn always minimizes Energy; if you want to maximize, let Energy return the negative of the target function
 - The target function class must know the number of parameters N that are optimized; that number should be returned by GetN()
 - The Energy function takes a pointer as argument. LamSimAnn makes sure that *par points to the first element of a memory block with N doubles, containing the proposed parameter values. So while implementing the target function, either read out the parameters sequentially (double var1 = *par; par++; double var2 = *par; p++ ....) or access the parameter values directly (double value = par[0]*par[1]-par[2]; ...)
 - LamSimAnn requires the parameter values to have a upper and lower boundary. These boundaries should be defined in the target function class and returned by GetLB() (lower boundary) and GetUB() (upper boundary). The number each boundaries must (of course) be N.
 - Remark: you can simply store the boundaries in vectors and let GetLB()/GetUB() return a pointer to the first element (return &LB[0];). Not super clean but works.
 
 */

class Energy_base // derive your fitness class from this template
{
public:
	virtual double Energy (double *par) = 0; // return fitness value (being minimized)
	virtual unsigned int GetN () = 0; // return number of parameters
	virtual double* GetLB () = 0; // return lower boundaries of parameters
	virtual double* GetUB () = 0; // return upper boundaries of parameters
};

/**********************************************************************************************************************
*********** class simulated annealing *********************************************************************************
**********************************************************************************************************************/

class LamSimAnn
{

public:
	// constructor / destructor
	explicit LamSimAnn (double lambda=0.0001, unsigned int tau=100, double theta_m_0=0.1, double frozen_N=10, double frozen_eps=10e-10);
    /*
     Contructor arguments:
     
     ** lambda: the only parameter that determines the speed of the cool-down. The smaller, the more thourogh the parameter space is searched but the longer the run-time. The 'right' value is a trade-off between speed and accuracy, dependent on how long your target function takes to evaluate. Should not be bigger than 0.1
     ** tau/theta_m_0: leave those at their standard values
     ** frozen_N/frozen_eps: determine the algorithm termination criterion. The algorithm stops if the state of the target function does not change more than frozen_eps in frozen_N consecutive steps. Decrease frozen_N to speed up the algorithm a little. Decrease frozen_eps if the target function return high values or is stochastic.
     */
	virtual ~LamSimAnn ();
	
	// interface
	void Run (Energy_base *loss_fct, double *start_par=NULL);
    /*
     Main function:
     Systematically varies parameters within the given boundaries and evaluates the energy function in loss_fct until the termination criterion is fulfilled. Per default a random set of parameters is used as a start point, but if a good guess is available it is also possible to define the initial parameters start_par
     */
	
	std::vector<double> GetBestPar() {return best_par;} // Returns the set of the best parameter found during the optimization
	double GetBestEnergy () {return best_energy;} // Returns the smallest value of the target function found found during the optimization (corresponding to the best parameters)
    
    // functions that return timeseries of quantities pertaining the optimization algorithm. Only good for optimizing the optimization parameters
	std::vector<double> GetAccRatioHist() {return acc_ratio_hist;} // not really infomative
	std::vector<double> GetMeanEnergyHist() {return mean_energy_hist;} // returns (locally averaged) values of the target function the algorithm went through. Good to check if the algorithm takes too long to terminate
	std::vector<double> GetTHist() {return temperature_hist;} // returns the history of temperature values used by the algorithm. (Always starts at high temperature and cools down to low temperatures)
	
	// data
private:
	// metaparameters
	unsigned int N;
	double m_lambda, m_K, m_theta_m_0;
	unsigned int m_tau;
    int m_frozen_N;
    double m_frozen_eps;
	double *m_LB, *m_UB;
	
	double w_a, w_b;
	
	// random number generator
	boost::mt19937 m_mt;
	boost::variate_generator< boost::mt19937& , boost::uniform_real<> > *m_rng;
	
	// auxiliary variables needed in member functions
	int i,l,k;
	int frozen_cnt;
	double rnd, rnd1, rnd2;
	Energy_base *m_anchor;
	
	// algorithm variables
	double current_energy;
	double m_mean, m_vari;
	std::vector<double> current_par;
	double proposed_par;
	std::vector<double> tmp_par;
	std::vector<double> m_theta;
	std::vector<double> acc_ratio;
	
	// Temperature and associated variables
	double S,dS;
	double estimate_mean, estimate_std;
	double A,B,D,E;
	
	double usxx,usyy,usxy,usx,usy,usum;
	double vsxx,vsyy,vsxy,vsx,vsy,vsum;
	double alpha;
	
	// results
	double best_energy;
	std::vector<double> best_par;
	std::vector<double> mean_energy_hist;
	std::vector<double> acc_ratio_hist;
	std::vector<double> temperature_hist;
	
	// auxiliary functions
private:
	double make_move(int par_n);
	void update_S();
	void update_parameter();
	void initialize_parameter();
	void update_control();
	bool is_frozen();

};

#endif