#ifndef __EmbrionicDevelopment__Ode_Integrator__
#define __EmbrionicDevelopment__Ode_Integrator__

// the following code is adapted from 'Numerical Recipes in C++'


#include "Implicit_Stepper_Sparse.h"
// #include <boost/numeric/mtl/mtl.hpp>
// #include <boost/numeric/itl/itl.hpp>
// #include <boost/array.hpp>



struct Output {
    typedef mtl::dense_vector<Doub> VecDoub;
    typedef mtl::dense2D<Doub> MatDoub;             // matrix type for ysave to save output data in y
    
	Int kmax;
	Int nvar;
	Int nsave;
	bool dense;
	Int count;
	Doub x1,x2,xout,dxout;
	VecDoub xsave;                              // one could easily change this and the next to containers from std library e.g. vector and vactor of vector
	MatDoub ysave;
	Output() : kmax(-1),dense(false),count(0) {}
	Output(const Int nsavee) : kmax(500),nsave(nsavee),count(0),xsave(kmax) {
		dense = nsave > 0 ? true : false;
	}
	void init(const Int neqn, const Doub xlo, const Doub xhi) {
		nvar=neqn;
		if (kmax == -1) return;
		ysave.change_dim(nvar,kmax);           // change to .resize() when using a library different from mtl4
		if (dense) {
			x1=xlo;
			x2=xhi;
			xout=x1;
			dxout=(x2-x1)/nsave;
		}
	}
	void resize() {
		Int kold=kmax;
		kmax *= 2;
		VecDoub tempvec(xsave);
		xsave.change_dim(kmax);             // change to .resize() when using a library different from mtl4
		for (Int k=0; k<kold; k++)          // this would be obsolete when using containers from the std library
			xsave[k]=tempvec[k];
		MatDoub tempmat(ysave);
		ysave.change_dim(nvar,kmax);       // change to .resize() when using a library different from mtl4
		for (Int i=0; i<nvar; i++)          // this would be obsolete when using containers from the std library
			for (Int k=0; k<kold; k++)
				ysave[i][k]=tempmat[i][k];
	}
	template <class Stepper>
	void save_dense(Stepper &s, const Doub xout, const Doub h) {
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=s.dense_out(i,xout,h);
		xsave[count++]=xout;
	}
    template<class VecDoub_IO>      // template here also for the input vector type because the Odeint class usually uses a different Vector type than Output class
	void save(const Doub x, const VecDoub_IO &y) {
		if (kmax <= 0) return;
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=y[i];
		xsave[count++]=x;
	}
	template <class Stepper, class VecDoub_IO>     // similar as for save, use a template here for the input vector type
	void out(const Int nstp,const Doub x,const VecDoub_IO &y,Stepper &s,const Doub h) {
		if (!dense)
			throw("dense output not set in Output!");
		if (nstp == -1) {
			save(x,y);
			xout += dxout;
		} else {
			while ((x-xout)*(x2-x1) > 0.0) {
				save_dense(s,xout,h);
				xout += dxout;
			}
		}
	}
};





template<class Stepper, short syst_size>
struct Odeint {
    typedef boost::array<Doub_IO, syst_size> VecDoub_IO;        // call the type of the initial vector VecDoub_I which can be different from the actual type VecDoub with which the routine works. E.G other solvers in odeint want boost::array, so this solver can be used interchangeably with the solvers of odeint on the same initial vector type.
    
    typedef typename StepperBase<syst_size>::vector_type VecDoub;       // get the internal vector type from StepperBase class. Should be a vector type of the mtl4 library, so it can be used in the bicgstab and preconditioer method.
    
	static const Int MAXSTP=500000;
	Doub EPS;
	Int nok;
	Int nbad;
	Int nvar;
	Doub x1,x2,hmin;
	bool dense;
	VecDoub y,dydx;
	VecDoub_IO &ystart;
	Output &out;
	typename Stepper::Dtype &derivs;
	Stepper s;
	Int nstp;
	Doub x,h;
	Odeint(VecDoub_IO &ystartt,const Doub xx1,const Doub xx2,
		const Doub atol,const Doub rtol,const Doub h1,
		const Doub hminn,Output &outt,typename Stepper::Dtype &derivss);
	Int integrate();
};

template<class Stepper, short syst_size>
Odeint<Stepper, syst_size>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
	const Doub atol, const Doub rtol, const Doub h1, const Doub hminn,
	Output &outt,typename Stepper::Dtype &derivss) : nvar(syst_size),
	ystart(ystartt),x(xx1),nok(0),nbad(0),                                             // no need to set size of y and dydx when sizes are fixed by type
	x1(xx1),x2(xx2),hmin(hminn),dense(outt.dense),out(outt),derivs(derivss),
	s(y,dydx,x,atol,rtol,dense) {
        EPS=std::numeric_limits<Doub>::epsilon();
	h=SIGN(h1,x2-x1);
	for (Int i=0;i<nvar;i++) y[i]=ystart[i];
	out.init(s.neqn,x1,x2);
}

template<class Stepper, short syst_size>
Int Odeint<Stepper, syst_size>::integrate() {
	derivs(y,dydx,x);
	if (dense)
		out.out(-1,x,y,s,h);
	else
		out.save(x,y);
	for (nstp=0;nstp<MAXSTP;nstp++) {
		if ((x+h*1.0001-x2)*(x2-x1) > 0.0)
			h=x2-x;
		s.step(h,derivs);
		if (s.hdid == h) ++nok; else ++nbad;
		if (dense)
			out.out(nstp,x,y,s,s.hdid);
		else
			out.save(x,y);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (Int i=0;i<nvar;i++) ystart[i]=y[i];
			if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS)
				out.save(x,y);
			return nok + nbad;
		}
   //     std::cout << std::abs(s.hnext) << std::endl;
		if (abs(s.hnext) <= hmin) throw("Step size too small in Odeint");
		h=s.hnext;
	}
	throw("Too many steps in routine Odeint");
}

#endif
