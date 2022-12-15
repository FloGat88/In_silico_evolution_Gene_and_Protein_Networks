#ifndef __EmbrionicDevelopment__Implicit_Stepper_Sparse__
#define __EmbrionicDevelopment__Implicit_Stepper_Sparse__

// the following code is adapted from 'Numerical Recipes in C++'

#include <limits>
#include <algorithm>
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <boost/array.hpp>


typedef int Int;
typedef float Doub;
typedef float Doub_IO;          // floating point number type in which the rhs of the system and the jacobian are delivered. Therefore, the floating point number type of vector_type and matrix_type in StepperBase and Odeint. This may be different from the Doub, the data type in which internal calculations are performed and the constants of StepperBase are stored. It is however reasonable to use the same data type for Doub and Doub_IO


// four utility functions from nr3.h (macro-like inline functions)
template<class T>                                       // with help of these one can replace the std::min and std::max again by MIN and MAX as in the original form
inline const T &MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

template<class T>
inline const T &MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

template<class T>
inline T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline T SQR(const T a) {return a*a;}
//-------------------------------------------------------------------------------------------------------------


// linear solver from mtl4 (templated for all kinds of matrices and vectors)
template < class LinearOperator, class HilbertSpaceX, class HilbertSpaceB, class Preconditioner>
void solve(const LinearOperator& A, HilbertSpaceX& x, const HilbertSpaceB& b,
           const Preconditioner& L, int max_iteractions =500)
{
    // Termination criterion: r < 1e-6 * b or N iterations
    itl::basic_iteration< Doub > iter( b , max_iteractions , 1e-6 );      // considier changing floating point type to double here as precision is essential for stability of iterative solvers.
    itl::bicgstab( A , x , b , L , iter );
    
}
//-------------------------------------------------------------------------------------------------------------


template<short syst_size>
struct StepperBase {
    
//    typedef boost::array<Doub_IO, syst_size> vector_type;     // does not work with boost::array because bicgstab and preconditioner need mtl4 vector types
    typedef mtl::dense_vector<Doub_IO, mtl::vec::parameters<mtl::tag::col_major, mtl::vec::fixed::dimension<syst_size>>> vector_type;
    //    typedef mtl::vector<Doub_IO, dim<syst_size>> vector_type;    // equivalently created with vector type generator of mtl4 (col_major is default) if the macros MTL_WITH_TEMPLATE_ALIAS, MTL_WITH_VARIADIC_TEMPLATE and MTL_WITH_STATICASSERT are enabled in the header
    typedef mtl::compressed2D<Doub_IO, mtl::mat::parameters<mtl::tag::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, unsigned short>> matrix_type;
    // matrix<float, sparse, as_size_type<unsigned short> >           // equivalent creation by the matrix type generator of mtl4 (again if three macros are enabled)
    
    Doub &x;
    Doub xold;
    vector_type &y,&dydx;
    Doub atol,rtol;
    bool dense;
    Doub hdid;
    Doub hnext;
    Doub EPS;
    Int n,neqn;
    vector_type yout,yerr;
    StepperBase(vector_type &yy, vector_type &dydxx, Doub &xx, const Doub atoll,
                const Doub rtoll, bool dens) : x(xx),y(yy),dydx(dydxx),atol(atoll),
    rtol(rtoll),dense(dens),n(syst_size),/*n(y.size())*/neqn(n)/*,yout(n),yerr(n)*/ {}
};


struct Ross_constants {
	static const Doub c2,c3,c4,bet2p,bet3p,bet4p,d1,d2,d3,d4,a21,a31,a32,
		a41,a42,a43,a51,a52,a53,a54,c21,c31,c32,c41,c42,c43,c51,c52,
		c53,c54,c61,c62,c63,c64,c65,gam,d21,d22,d23,d24,d25,d31,d32,
		d33,d34,d35;
};




template <class D, short syst_size>
struct StepperRoss : StepperBase<syst_size>, Ross_constants {
//-------------------------------------------------------------------------------------
    using typename StepperBase<syst_size>::vector_type;                      // make all typedefs and data types of the template base class StepperBase accessible in the derived class (this is neccessary because base class is template base class; alternatives would be to call members of the base class as StepperBase::a or this->a)
    using typename StepperBase<syst_size>::matrix_type;
    
    using StepperBase<syst_size>::x;
    using StepperBase<syst_size>::xold;
    using StepperBase<syst_size>::y;
    using StepperBase<syst_size>::dydx;
    using StepperBase<syst_size>::atol;
    using StepperBase<syst_size>::rtol;
    using StepperBase<syst_size>::dense;
    using StepperBase<syst_size>::hdid;
    using StepperBase<syst_size>::hnext;
    using StepperBase<syst_size>::EPS;
    using StepperBase<syst_size>::n;
    using StepperBase<syst_size>::neqn;
    using StepperBase<syst_size>::yout;
    using StepperBase<syst_size>::yerr;
    using StepperBase<syst_size>::StepperBase;
//-------------------------------------------------------------------------------------
    
    typedef D Dtype;
    
    matrix_type dfdy = matrix_type(syst_size, syst_size);          // as the size of the matix is fixed, set the size directly here in the declaration
//	vector_type dfdx;
	vector_type k1,k2,k3,k4,k5,k6;
	vector_type cont1,cont2,cont3,cont4;
    matrix_type a = matrix_type(syst_size, syst_size);
	StepperRoss(vector_type &yy, vector_type &dydxx, Doub &xx, const Doub atoll,
		const Doub rtoll, bool dens);
	void step(const Doub htry,D &derivs);
	void dy(const Doub h,D &derivs);
	void prepare_dense(const Doub h,const vector_type &dydxnew);
	Doub dense_out(const Int i, const Doub x, const Doub h);
	Doub error();
	struct Controller {
		Doub hnext;
		bool reject;
		bool first_step;
		Doub errold;
		Doub hold;
		Controller();
		bool success(Doub err, Doub &h);
	};
	Controller con;
};

template <class D, short syst_size>
StepperRoss<D,syst_size>::StepperRoss(vector_type &yy, vector_type &dydxx, Doub &xx,
	const Doub atoll,const Doub rtoll, bool dens) :
	StepperBase<syst_size>(yy,dydxx,xx,atoll,rtoll,dens)/*,dfdy(n,n)/*,dfdx(n)*//*,k1(n),k2(n),             // no need to set the sizes and dimensions of all the vectors and matrices, as these are already predifined in the data type or in the declaration.
	k3(n),k4(n),k5(n),k6(n),cont1(n),cont2(n),cont3(n),cont4(n),a(n,n)*/ {
        EPS=std::numeric_limits<Doub>::epsilon();
}

template <class D, short syst_size>
void StepperRoss<D,syst_size>::step(const Doub htry,D &derivs) {
	vector_type dydxnew;
	Doub h=htry;
	derivs.Jacobian(/*x,*/ y /*,dfdx*/,dfdy);
	for (;;) {
		dy(h,derivs);
		Doub err=error();
		if (con.success(err,h)) break;
		if (abs(h) <= abs(x)*EPS)
			throw("stepsize underflow in StepperRoss");
	}
	derivs(yout,dydxnew,double(x+h));       // convert to double as ReactionDiffusionDynamics takes time as double type (as is required from the Odeint library) and we want to avoid using a template on ReactionDiffusionDynamics
	if (dense)
		prepare_dense(h,dydxnew);
	dydx=dydxnew;
	y=yout;
	xold=x;
	x += (hdid=h);
	hnext=con.hnext;
}

template<class D, short syst_size>
void StepperRoss<D,syst_size>::dy(const Doub h,D &derivs) {
	vector_type ytemp,dydxnew;          // remove initialization of vecotrs by sizes as the size of vector_type is already fixed by its definition
	Int i;
//	for (i=0;i<n;i++) {
//		for (Int j=0;j<n;j++) a[i][j] = -dfdy[i][j];
    
    a = 1.0/(gam*h);
    a = a - dfdy;
    
//    for (i=0;i<n;i++)
//		a[i][i] += 1.0/(gam*h);
	
	itl::pc::ilu_0<matrix_type, float> L( a );      // use float precision here independent of the floating number type elsewhere (Doub) as for preconditioning precision is not essential for stability
//    LUdcmp alu(a);
	for (i=0;i<n;i++)
        ytemp[i]=dydx[i];      //+h*d1*dfdx[i];
 	solve(a,k1,ytemp,L);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+a21*k1[i];
	derivs(ytemp,dydxnew,x+c2*h);
	for (i=0;i<n;i++)
		ytemp[i]=dydxnew[i]+c21*k1[i]/h;   // +h*d2*dfdx[i]
	solve(a,k2,ytemp,L);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+a31*k1[i]+a32*k2[i];
	derivs(ytemp,dydxnew,x+c3*h);
	for (i=0;i<n;i++)
		ytemp[i]=dydxnew[i]+(c31*k1[i]+c32*k2[i])/h;    // +h*d3*dfdx[i]
	solve(a,k3,ytemp,L);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+a41*k1[i]+a42*k2[i]+a43*k3[i];
	derivs(ytemp,dydxnew,x+c4*h);
	for (i=0;i<n;i++)
		ytemp[i]=dydxnew[i]+(c41*k1[i]+c42*k2[i]+c43*k3[i])/h;      // +h*d4*dfdx[i]
	solve(a,k4,ytemp,L);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i];
	Doub xph=x+h;
	derivs(ytemp,dydxnew,xph);
	for (i=0;i<n;i++)
		k6[i]=dydxnew[i]+(c51*k1[i]+c52*k2[i]+c53*k3[i]+c54*k4[i])/h;
	solve(a,k5,k6,L);
	for (i=0;i<n;i++)
		ytemp[i] += k5[i];
	derivs(ytemp,dydxnew,xph);
	for (i=0;i<n;i++)
		k6[i]=dydxnew[i]+(c61*k1[i]+c62*k2[i]+c63*k3[i]+c64*k4[i]+c65*k5[i])/h;
	solve(a,yerr,k6,L);
	for (i=0;i<n;i++)
		yout[i]=ytemp[i]+yerr[i];
}


template <class D, short syst_size>
void StepperRoss<D,syst_size>::prepare_dense(const Doub h,const vector_type &dydxnew) {
	for (Int i=0;i<n;i++) {
		cont1[i]=y[i];
		cont2[i]=yout[i];
		cont3[i]=d21*k1[i]+d22*k2[i]+d23*k3[i]+d24*k4[i]+d25*k5[i];
		cont4[i]=d31*k1[i]+d32*k2[i]+d33*k3[i]+d34*k4[i]+d35*k5[i];
	}
}

template <class D, short syst_size>
Doub StepperRoss<D,syst_size>::dense_out(const Int i,const Doub x,const Doub h) {
	Doub s=(x-xold)/h;
	Doub s1=1.0-s;
	return cont1[i]*s1+s*(cont2[i]+s1*(cont3[i]+s*cont4[i]));
}

template <class D, short syst_size>
Doub StepperRoss<D,syst_size>::error() {
	Doub err=0.0,sk;
	for (Int i=0;i<n;i++) {
        sk=atol+rtol*std::max(abs(y[i]),abs(yout[i]));
		err += SQR(yerr[i]/sk);
	}
	return sqrt(err/n);
}

template <class D, short syst_size>
StepperRoss<D,syst_size>::Controller::Controller() : reject(false), first_step(true) {}

template <class D, short syst_size>
bool StepperRoss<D, syst_size>::Controller::success(Doub err, Doub &h) {
	static const Doub safe=0.9,fac1=5.0,fac2=1.0/6.0;
    Doub fac=std::max(fac2,std::min(fac1,Doub(pow(err,0.25)/safe)));         // change MAX and MIN from numerical recipes (nr3.h) to max and min, repectively, from stdlib.
	Doub hnew=h/fac;
	if (err <= 1.0) {
		if (!first_step) {
			Doub facpred=(hold/h)*pow(err*err/errold,0.25)/safe;
			facpred=std::max(fac2,std::min(fac1,facpred));               // change MAX and MIN from numerical recipes to max and min, repectively, from stdlib.
			fac=std::max(fac,facpred);                               // change MAX and MIN from numerical recipes to max and min, repectively, from stdlib.
			hnew=h/fac;
		}
		first_step=false;
		hold=h;
		errold=std::max(Doub(0.01),err);                                      // change MAX and MIN from numerical recipes to max and min, repectively, from stdlib.
		if (reject)
            hnew=(h >= 0.0 ? std::min(hnew,h) : std::max(hnew,h));        // change MAX and MIN from numerical recipes to max and min, repectively, from stdlib.
		hnext=hnew;
		reject=false;
		return true;
	} else {
		h=hnew;
		reject=true;
		return false;
	}
}



#endif
