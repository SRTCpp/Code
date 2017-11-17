#include "Jnr.h"

#define MAXIT 60

template <class Datatype>
Datatype Jcrap::nDfunction<Datatype>::the_function(Datatype x)
{
	vector<Datatype> xv; 
	xv[0]=x; 
	return the_function(xv);
}
		
		
template <class Datatype>
Datatype JMAX(Datatype x, Datatype y)
{
	if (x > y) return x;
	else return y;
}

template <class Datatype>
Jcrap::nDfunction<Datatype>::nDfunction(int d) : dimensions(d) {jMAX=200;}

template <class Datatype>
cube Jcrap::nDfunction<Datatype>::plot(const cube& c)
{
	cube answer(c);
	cout << "Plotting --  00%";

	if (dimensions > 3) {
		cout << "Attempt to plot too many dimensions\n";
		return answer;
	}
	
	for (int x=0;x<answer.N(X);x++) {
		printpercent(x, c.N(X)-1);
		for (int y=0;y<answer.N(Y);y++) {
			for (int z=0;z<answer.N(Z);z++) {
				vector<Datatype> xvec(dimensions);
				if (dimensions >= 1) xvec.at(0) = c.Axis(X,x);
				if (dimensions >= 2) xvec.at(1) = c.Axis(Y,y);
				if (dimensions >= 3) xvec.at(2) = c.Axis(Z,z);
				
				answer(x,y,z) = the_function(xvec);
			}
		}
	}
	cout << "\n";
	return answer;
}

