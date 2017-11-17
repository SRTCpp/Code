#include <vector>
#include <algorithm>

#define _USESTDVECTOR_ 1
#include "../../nr/nr3/Jnr3.h"

#include "Jcrap.h"


#ifndef JFUNC
#define JFUNC 1


namespace Jcrap {

template <class Datatype>
class function
{
	public:
		function();
			
		virtual Datatype the_function(Datatype x) {
			Dimensionality = 1;
			vector<Datatype> d(1);
			d[0] = x;
			cout << "in creepy virtual the_function in Jcrap::function<Datatype>\n";
			return the_function(d);}
		virtual Datatype the_function(vector<Datatype> x) {
			cout << "not using the multifunction for the_function(vector<>)\n";
			return the_function(x[0]);   // will cause infinite loop if neither is overridden
		}
		Datatype neg_function(Datatype x) { return (-the_function(x));  }
		Datatype neg_function(vector<Datatype> x) { return (-the_function(x));  }
		pair<Datatype,Datatype> derivative(Datatype, Datatype=-1.);
		pair<Datatype,Datatype> derivative_dumb(Datatype, Datatype=-1.);
		pair<Datatype,Datatype> derivative2nd_dumb(Datatype, Datatype=-1.);
	
// methods
		Datatype findzero(Datatype,Datatype,Datatype,Datatype=0.,
				string=string(""));
		Datatype findzero(Datatype,Datatype,Datatype,string);
		pair<Datatype,Datatype> findminimum(Datatype,Datatype,Datatype,Datatype);
		pair<Datatype,Datatype> findmaximum(Datatype,Datatype,Datatype,Datatype);
		vector<pair<Datatype, Datatype> > bracketminimum(Datatype,Datatype);
		Datatype integrate_smooth(Datatype,Datatype,Datatype=1e-6, int=5, string=string(""));
		Datatype integrate_rough (Datatype,Datatype,Datatype=1e-6, int=0);
		Datatype integrate(Datatype a, Datatype b, Datatype c=1.e-6, int d=5, string s=string("")) {
			return integrate_smooth(a,b,c,d,s);}
		cube plot(Datatype, Datatype, int=100);
		cube plot(const cube&);
		
		Datatype operator()(Datatype x) { return the_function(x); }
		Datatype operator()(vector<Datatype> x) { return the_function(x); }

		void goldenratio(Datatype d) { GOLD=d;}
		Datatype goldenratio() {return GOLD;}
		void goldenlimit(Datatype d) { GLIMIT=d;}
		Datatype goldenlimit() {return GLIMIT;}
		void tiny(Datatype d) { TINY=d;}
		Datatype tiny() {return TINY;}

//  Tunable parameters
		int jMAX;
		
	protected:
		int Dimensionality;
			
	private:
		Datatype trapzd(Datatype, Datatype, int);
		Datatype s; // for trapzd.
		Datatype GOLD,GLIMIT,TINY;
		

};

template <class Datatype>
class multifunction : public function<Datatype>
{		
	public:
		multifunction(int na=0);
		multifunction(vector<Datatype>);
		
// function related	
		virtual Datatype the_function(Datatype, vector<Datatype>&, unsigned int=0)=0;
		Datatype the_function(Datatype x) {if ((unsigned int)(thisdim)<a.size()&&thisdim>=0) a.at(thisdim)=x; return the_function(x,a);}
		Datatype the_function(vector<Datatype> x) {return the_function(0.,x);}
		Datatype operator()(Datatype x) {return the_function(x);}	
		Datatype operator()(vector<Datatype> x) {return the_function(x);}
		virtual vector<Datatype> the_derivitives(Datatype, vector<Datatype>&, unsigned int=0);
		
		Datatype funcs(Datatype, vector<Datatype>&, vector<Datatype>&, unsigned int=0);
		Datatype aderiv(int, Datatype, vector<Datatype>, Datatype, unsigned int=0);
		Datatype dumbderiv(int, Datatype, vector<Datatype>, Datatype, unsigned int=0);
		virtual Datatype aconvert(vector<Datatype>&, vector<int>&, int);
		pair<vector<Datatype>,Datatype> findmaximum();  // all-automatic version
		pair<vector<Datatype>,Datatype> findminimum();  // all-automatic version
		pair<vector<Datatype>,Datatype> findminimum(vector<Datatype>, Datatype, Datatype=1e-5);
		pair<vector<Datatype>,Datatype> findminimum(vector<Datatype>, vector<Datatype>, Datatype=1e-5);
		pair<vector<Datatype>,Datatype> findminimum(MatDoub_I &pp, Datatype=1e-5);
		pair<vector<Datatype>,Datatype> findmaximum(MatDoub_I &pp, Datatype=1e-5);
		pair<vector<Datatype>,Datatype> findmaximum(vector<Datatype>, Datatype);
		pair<vector<Datatype>,Datatype> findminimum_dumb(vector<Datatype>);
		pair<vector<Datatype>,Datatype> findminimum_dumb();
		pair<Datatype,Datatype> findminimum(int,Datatype,Datatype,Datatype,Datatype);
		
// additional methods
		pair<vector<Datatype>, Datatype> fit(cube);
		pair<vector<Datatype>, Datatype> fit(cube, vector<Datatype>);
		pair<vector<Datatype>, Datatype> fit(cube, vector<Datatype>, vector<int>);
		
		Datatype chisquared(cube&);
		Datatype reducedchisquared(cube&);
		vector<Datatype> finderrors() const;
		
// MCMC methods
		Datatype Pstate(cube&, vector<Datatype>);
		Datatype Pstate(cube&);
		vector<Datatype> Proposal();
		Datatype MCMCstep(cube&);

// some NR methods for fit
		Datatype Jmrqmin(vector<Datatype>&, vector<Datatype>&, vector<Datatype>&, 
				vector<Datatype>&, vector<int>&, vector<unsigned int>&);
		void Jcovsrt(vector<vector<Datatype> > &covar, int ma, vector<int> &ia, int mfit);
		void Jgaussj(vector<vector<Datatype> > &a, int n, vector<vector<Datatype> > &b, int m);
		Datatype Jmrqcof(vector<Datatype>&, vector<Datatype>&, vector<Datatype>&, 
				vector<Datatype>&, vector<int>&, int,vector<vector<Datatype> > &, vector<Datatype>&, vector<unsigned int>&);
			
// additional data for fitting routines
		vector<Datatype> a;
		vector<string> param_names;
		vector<Datatype> dyda;
		vector<int> ia;
		vector<double> h;
		vector<int> quantized;
		vector<vector<Datatype> > covar;	
		vector<vector<Datatype> > alpha;	
		Datatype alamda;
		Datatype lamda_up;
		Datatype lamda_down;
		Datatype upperstop;
		Datatype initialalamda;
		Datatype chidiffstop;		
		
		struct Jamoeba;
	
		int thisdim;
		bool skiperroranalysis;
};


}  // (end namespace)

#include "Jfunction.c++"
#include "Jmultifunction.implimentations.c++"
#include "Jmultifunction.c++"

#endif
