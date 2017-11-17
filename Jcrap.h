#ifndef Jcrap_h
#define Jcrap_h 1

#include <string>
#include <strstream>
#include <vector>
#include <list>
#include <map>


using namespace std;


namespace Jcrap {
	const double pi(3.141592653589793238462643383279502884197);

	class bigcube;
	
	class unitelement;
	class unit;
	class diffeq;
//	template <class Datatype> class function;
//	template <class Datatype> class multifunction;
	class date;
	void printfitserror(int);

// for new Jcubeview, 8/2004
	class imagespace;
	class slidedial;
	class cubeview;
	class cubeviewer;
	class controlpanel;
	
// binary tree for classification 2005/09
	class pixelnode;

	double polygonarea(vector<pair<double, double> >);
}

using namespace Jcrap;


#include "Jvalue.h"
#include "Jcube.h"
#include "Jaddons.h"
#include "Jassociated.h"

#endif
