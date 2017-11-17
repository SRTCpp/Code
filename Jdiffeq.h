#include <valarray>
#include "Jcrap.h"

template <class Datatype>
class Jcrap::diffeq
{
	public:
		diffeq();
		virtual ~diffeq() = 0;
		virtual void odeint(Datatype, Datatype);	// called with starting and ending x values
		void graph(string="");
		void setup();

	protected:
// NR odeint variables
		long int maxsteps;
		Datatype tiny;
		int kmax;						// desired max number of intermediate steps to save
		int nvar;			// number of equations to be solved simultaneously
		
		Datatype dxsav;				// approximate interval between intermediate saves (in x)
		Jcrap::unit xunit;			// units x
		Jcrap::unit yunit;			// units y
		vector<Datatype> ystart;	// starting y values		
		Datatype eps;					// accuracy, defaults to 1e-3
		Datatype h1;					// guessed first stepsize
		Datatype hmin;					// minimum stepsize
		vector<Datatype> yscal,y,dydx;

		// dude this shit is weird
// new NR routines
		virtual void derivs(Datatype, vector<Datatype>&, vector<Datatype>&)=0; 
		virtual Datatype func(Datatype d) {return d;}
		virtual Datatype auxfunc(Datatype d) {return d;}
		
		void rkqs(vector<Datatype>&, vector<Datatype>&, Datatype&, Datatype, 
				Datatype&, Datatype&);
		void rkck(Datatype, Datatype, vector<Datatype>&, vector<Datatype>&);	
		
		Datatype qromb(Datatype, Datatype);
		Datatype trapzd(Datatype, Datatype, int);
		Datatype qtrap(Datatype, Datatype);

	public:	
		Datatype zriddr(Datatype, Datatype, Datatype);
		pair<Datatype, Datatype> polint(Datatype *, Datatype *, int, Datatype);
		unsigned long locate(vector<Datatype> xx, Datatype x);
//	protected:
		int kount;						//   should be deprecated in the final writing I believe


// variables for the NR routines		
	// for odeint	
		vector<Datatype> xp;			// storage of intermediate x values
		vector<vector<Datatype> > yp;			// storage of intermediate y values	
	protected:	
		int nok;							// # of good steps taken
		int nbad;						// # of bad steps taken
	// for qromb
		Datatype EPS;
		int JMAX;
		int JMAXP;
		int K;
};

template <class Datatype>
Jcrap::diffeq<Datatype>::diffeq()
{
//	cout << "Diffeq constructor.\n";
	kmax=10000;
	maxsteps=10000;
	setup();
//	cout << "End Diffeq constructor\n";
}

template <class Datatype>
void Jcrap::diffeq<Datatype>::setup()
{
	tiny=1e-50;
	nvar=1;
	eps=1e-3;
	EPS=eps;
	JMAX=20;
	JMAXP=JMAX+1;
	K=5;
	yscal.resize(nvar);
	y.resize(nvar);
	dydx.resize(nvar);
}

template <class Datatype>
Jcrap::diffeq<Datatype>::~diffeq() {}
#include "Jdiffeq.c++"
