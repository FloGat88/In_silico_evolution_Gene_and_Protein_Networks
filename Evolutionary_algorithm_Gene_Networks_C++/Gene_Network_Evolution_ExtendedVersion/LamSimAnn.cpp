/*
 *  LamSimAnn.cpp
 *  SimAnn
 *
 *  Created by Patrick Hillenbrand on 5/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "LamSimAnn.h"
#include <fstream>
#include <iostream>
#include <stdexcept>


/**********************************************************************************************************************
*********** class simulated annealing implementation ******************************************************************
**********************************************************************************************************************/


LamSimAnn::LamSimAnn (double lambda, unsigned int tau, double theta_m_0, double frozen_N, double frozen_eps)
: m_lambda(lambda), m_tau(tau), m_theta_m_0(theta_m_0), m_frozen_N(frozen_N), m_frozen_eps(frozen_eps)
{
	m_K = 3.;
	
	// initialize random number generator
	m_rng = new boost::variate_generator< boost::mt19937& , boost::uniform_real<> >(m_mt, boost::uniform_real<>(0,1));
	
	unsigned int seed;
	std::ifstream dev_random("/dev/random", std::ios::binary|std::ios::in);
	dev_random.read(reinterpret_cast<char*>(&seed), sizeof(seed));
	dev_random.close();
	m_rng->engine().seed(seed);
}

/*********************************************************************************************************************/

LamSimAnn::~LamSimAnn ()
{
	delete m_rng;
}

/*********************************************************************************************************************/

void LamSimAnn::Run (Energy_base *loss_fct, double *start_par)
{
	// initialize variables
	m_anchor = loss_fct;
	N = loss_fct->GetN();
	m_LB = loss_fct->GetLB();
	m_UB = loss_fct->GetUB();
	acc_ratio_hist.clear();
	temperature_hist.clear();
	mean_energy_hist.clear();
	frozen_cnt=0;
	
	best_energy = 1.0e300;
	
	double proposed_energy;
	acc_ratio.resize(N);
	
	m_theta.resize(N);
	for (i = 0; i < N; i++) m_theta[i]=m_theta_m_0*(m_UB[i]-m_LB[i]);
	
	// choose either random initial parameter set or set to assigned value
	current_par.resize(N);
	if (start_par == NULL) {
        for (i = 0; i < N; i++) current_par[i]=(*m_rng)()*(m_UB[i]-m_LB[i])+m_LB[i];
    } else {
        for (i = 0; i < N; i++) {
            if (start_par[i] < m_LB[i]) {
                current_par[i] = m_LB[i];
            } else if (start_par[i] > m_UB[i]) {
                current_par[i] = m_UB[i];
            } else {
                current_par[i] = start_par[i];
            }
        }
    }
	current_energy = m_anchor->Energy(&current_par[0]);
	
	// length of burn-in runs
	int M = 1000;
	
	// decorrelation run
	m_mean = 0.;
	for (i = 0; i < M; i++)
	{
		for (k = 0; k < N; k++) {
			proposed_energy = make_move(k);
			if (proposed_energy <= current_energy) {
				current_energy = proposed_energy;
				current_par[k] = proposed_par;
			}
			m_mean += current_energy;
		}
	}
	m_mean /= (double) (M*N);
	S = 1./m_mean;
	
	// gather initial statistics
	m_mean = m_vari = 0.0;
	
	for (i = 0; i < M; i++)
	{
		for (k = 0; k < N; k++) {
			proposed_energy = make_move(k);
			rnd = (*m_rng)();
			if ((proposed_energy - current_energy <= 0.0) || exp(-S*(proposed_energy-current_energy))>rnd) {
				current_energy = proposed_energy;
				current_par[k] = proposed_par;
				acc_ratio[k] += 1.;
				if (current_energy < best_energy) { 
					best_energy=current_energy;
					best_par = current_par;
				}
			}
			m_mean += current_energy;
			m_vari += pow(current_energy - estimate_mean,2);
		}
	}
	m_mean /= (double) (M*N);
	m_vari = m_vari / ((double) (M*N)) - m_mean * m_mean;
	for (i = 0; i < N; i++) {acc_ratio[i] /= ((double) M);}
	
	initialize_parameter();
	
	// main loop
	
		// continue algorithm until stopping criterion (frozen state) is fulfilled
	while (true) {
		m_mean = m_vari = 0.0;
		for (i = 0; i < N; i++) {
			acc_ratio[i] = 0.;
		}
		
			// do m_tau parameter sweeps until parameter update
		for (i = 0; i < m_tau; i++) {
			
				// toggle through all parameters seperately
			for (k = 0; k < N; k++) {
				proposed_energy = make_move(k);
				rnd = (*m_rng)();
				if ((proposed_energy - current_energy <= 0.0) || exp(-S*(proposed_energy-current_energy))>rnd) 
				{
					current_energy = proposed_energy;
					current_par[k] = proposed_par;
					
					acc_ratio[k] += 1.;
					
					if (current_energy < best_energy) {
						best_energy = current_energy;
						best_par = current_par;
					}
				}
				m_mean += current_energy;
				m_vari += pow(current_energy - estimate_mean,2);
			}
			
			update_S(); //update the temperature
			
		}
		
		m_mean /= (double) (m_tau*N);
		m_vari /= (double) (m_tau*N);
		for (i = 0; i < N; i++) {acc_ratio[i] /= ((double) m_tau);}
		
		// check termination condition
		if (is_frozen()) break;
		
		// update parameters
		update_parameter();
		update_control();
	}
	
	m_LB = NULL;
	m_UB = NULL;
	m_anchor = NULL;
}

/*********************************************************************************************************************/

double LamSimAnn::make_move (int par_n)
{	
	tmp_par = current_par;
    int n = 0;
	while (1) {
		rnd1 = (*m_rng)();
		rnd2 = (*m_rng)();
		if (rnd1<=0.5) {
			proposed_par=current_par[par_n]+m_theta[par_n]*log(rnd2);
		} else {
			proposed_par=current_par[par_n]-m_theta[par_n]*log(rnd2);
		}
		tmp_par[par_n] = proposed_par;
		if (proposed_par >= m_LB[par_n] && proposed_par <= m_UB[par_n]) {
			return m_anchor->Energy(&tmp_par[0]);
		} else n++;
        if (n == 100) {
            std::cerr << "stuck!\n";
            throw std::runtime_error("..");
        }
	}
}

/*********************************************************************************************************************/

void LamSimAnn::update_S()
{
	S += dS;
	estimate_mean = 1.0 / (A * S + B);
	estimate_std = 1.0 / (D * S + E);
	dS = m_lambda * alpha / (pow(S * estimate_std,2) * estimate_std);
}

/*********************************************************************************************************************/

void LamSimAnn::update_parameter()
{
	double d;
	
	d = 1.0 / m_mean;
	usyy *= w_a;
	usxy *= w_a;
	usy *= w_a;
	usx *= w_a;
	usxx *= w_a;
	usum *= w_a;
	usyy += d * d;
	usxy += S*d;
	usy += d;
	usx += S;
	usxx += S * S;
	usum += 1.0;
	A = (usum * usxy - usx * usy) / (usum * usxx - usx *usx);
	B = ( usy - A * usx ) /usum;
	estimate_mean = 1.0 / (A * S + B);
	if (m_vari > 0.0) 
	{
		d = 1.0 / sqrt(m_vari);
		vsyy *= w_b;
		vsxy *= w_b;
		vsy *= w_b;
		vsx *= w_b;
		vsxx *= w_b;
		vsum *= w_b;
		vsyy += d * d;
		vsxy += S*d;
		vsy += d;
		vsx += S;
		vsxx += S * S;
		vsum += 1.0;
		D=(vsum*vsxy-vsx*vsy)/(vsum*vsxx-vsx*vsx);
		E=(vsy-D*vsx)/vsum;
	}
	estimate_std = 1.0 / (D * S + E);
	double tt=0.;
	for (i = 0; i < N; i++) {
		tt += acc_ratio[i];
	}
	tt /= (double) N;
	acc_ratio_hist.push_back(tt);
	temperature_hist.push_back(1./S);
	d = (1.0 - tt) / (2.0 - tt);
	alpha = 4.0 * tt * d * d;
}

/*********************************************************************************************************************/

void LamSimAnn::initialize_parameter ()
{
	double d;
	
	estimate_std = sqrt(m_vari);
	estimate_mean = m_mean;
	A = estimate_std * estimate_std / (estimate_mean * estimate_mean);
	B = 1.0 / estimate_mean-A*S;
	D = estimate_std / estimate_mean;
	E = 1.0 / estimate_std-D*S;
	usum = vsum = 1.0;
	usxy = usxx = usx = 0.0;
	usy = 1.0 / estimate_mean;
	usyy = usy * usy;
	vsxy = vsxx = vsx = 0.0;
	vsy = 1.0 / estimate_std;
	vsyy = vsy * vsy;
	
	double tt=0.;
	for (i = 0; i < N; i++) {
		tt += acc_ratio[i];
	}
	tt /= (double) N;
	acc_ratio_hist.push_back(tt);
	temperature_hist.push_back(1./S);
	d = (1.0 - tt) / (2.0 - tt);
	alpha = 4.0 * tt * d * d;
	dS=m_lambda*alpha/(estimate_std*estimate_std*estimate_std*S*S);
	
	w_a = 2./m_lambda;
	w_a = 1.0 - m_tau / w_a;
	if (w_a < 0.0) w_a = 0.0;
	w_b = 100./m_lambda;
	w_b = 1.0 - m_tau / w_b;
	if (w_b < 0.0) w_b = 0.0;
}

/*********************************************************************************************************************/

void LamSimAnn::update_control()
{
	for (i=0; i < N; i++) {
		m_theta[i]=m_theta[i]*exp(m_K*(acc_ratio[i]-0.44));
		if (m_theta[i] > (m_UB[i]-m_LB[i])) m_theta[i] = (m_UB[i]-m_LB[i]);
	}
}

/*********************************************************************************************************************/

bool LamSimAnn::is_frozen()
{
	if (!mean_energy_hist.empty()) {
		if ((m_mean - mean_energy_hist.back())/m_mean < m_frozen_eps) {
			frozen_cnt++;
		} else {
			frozen_cnt = 0;
		}
	}
	mean_energy_hist.push_back(m_mean);
	return (frozen_cnt >= m_frozen_N);
}

