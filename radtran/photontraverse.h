#include "atmospheres.h"
#include "../Jfunction.h"
#include "photons.h"

class photontraverse;

class segment : public Jcrap::function<double>
{
	public:
		segment(double, atmozone*, photontraverse*);
	
		double Start();
		atmozone* Zone();
		photontraverse* Traverse();
		
				
		double the_function(double);  // returns extinction as a function of d along
												// the photon's path.
			
	private:
		double _start;
		atmozone* _tozone;
		photontraverse* _totraverse;
};		

class photontraverse : public Jcrap::function<double>
{
	public:
		photontraverse(photon&, atmosphere&);
	
		list<segment>& Segments() {return _segments;}
		
		double the_function(double);  // this function returns the integral
		                              // of extinction (i.e., Tau_total(d))
												// along the photon's flight path
		
		photon Photon();
		
	private:
			
		list<segment> _segments;
		photon _thisphoton;
	
};

class photontraverse_extinction : public photontraverse
//  This class is just a photontraverse that has a different the_function,
//  this time the_function returns the extinction, instead of the integral of extinction
{
	public:
		photontraverse_extinction(photon&, atmosphere&);	
		
		double the_function(double);  // now returning local extinction from the segments
};
