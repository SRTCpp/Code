#include <vector>

#include "Jcrap.h"

using namespace std;

namespace Jcrap {

template <class Datatype>
class nDfunction
{
	public:
		nDfunction(int);
			
		Datatype the_function(Datatype x); 
		virtual Datatype the_function(vector<Datatype> x) {}
		Datatype neg_function(vector<Datatype> x) { return (-the_function(x));  }
//		pair<Datatype,Datatype> derivitive(Datatype, Datatype=-1.);
	
// methods
//		Datatype findzero(Datatype,Datatype,Datatype);
//		pair<Datatype,Datatype> findminimum(Datatype,Datatype,Datatype,Datatype);
//		pair<Datatype,Datatype> findmaximum(Datatype,Datatype,Datatype,Datatype);
//		Datatype integrate_smooth(Datatype,Datatype,Datatype=1e-6, int=5);
//		Datatype integrate_rough (Datatype,Datatype,Datatype=1e-6, int=0);
//		cube plot(Datatype, Datatype, int=100);
		cube plot(const cube&);
		
		Datatype operator()(vector<Datatype> x) { return the_function(x); }

//  Tunable parameters
		int jMAX;
		
	public:
		int dimensions;
//		Datatype trapzd(Datatype, Datatype, int);
//		Datatype s; // for trapzd.
		

};

template <class Datatype>
class nDmultifunction : public nDfunction<Datatype>
{		
	public:
		nDmultifunction(int d, int na=0);
		nDmultifunction(int d, vector<Datatype>);
		
// function related	
		virtual Datatype the_function(const vector<Datatype>& x, 
				vector<Datatype>& v)=0;
		Datatype the_function(vector<Datatype> x) {return the_function(x,a);}		
		Datatype operator()(vector<Datatype> x) {return the_function(x,a);}
		virtual vector<Datatype> the_derivitives(vector<Datatype> x, vector<Datatype>& aa);
		
		Datatype funcs(vector<Datatype>, vector<Datatype>&, vector<Datatype>&);
		Datatype aderiv(int, vector<Datatype>, vector<Datatype>, Datatype);
		virtual Datatype aconvert(vector<Datatype>&, vector<int>&, int);
		
// additional methods
		pair<vector<Datatype>, Datatype> fit(const cube&);
		pair<vector<Datatype>, Datatype> fit(const cube&, const cube&);
		pair<vector<Datatype>, Datatype> fit(const cube&, const cube&, vector<Datatype>);
		pair<vector<Datatype>, Datatype> fit(const cube&, const cube&, vector<Datatype>, vector<int>);

// some NR methods for fit
		Datatype Jmrqmin(const cube&, const cube&,
				vector<Datatype>&, vector<int>&);
		void Jcovsrt(vector<vector<Datatype> > &covar, int ma, vector<int> &ia, int mfit);
		void Jgaussj(vector<vector<Datatype> > &a, int n, vector<vector<Datatype> > &b, int m);
		Datatype Jmrqcof(const cube&, const cube&, vector<Datatype>&, 
				vector<int>&, int, vector<vector<Datatype> > &, vector<Datatype>&);
			
// additional data for fitting routines
		vector<Datatype> a;
		vector<string> param_names;
		vector<Datatype> dyda;
		vector<int> ia;
		vector<vector<Datatype> > covar;	
		vector<vector<Datatype> > alpha;	
		Datatype alamda;
		Datatype lamda_up;
		Datatype lamda_down;
		Datatype upperstop;
		Datatype initialalamda;
		Datatype chidiffstop;

};

}

#include "JnDfunction.c++"
#include "JnDmultifunction.c++"
#include "JnDmultifunction.implimentations.c++"
