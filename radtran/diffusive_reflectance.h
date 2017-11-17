#include "../Jcrap.h"
#include "../Jfunction.h"

class diffusive_reflectance : public Jcrap::function<double>
{	// taken from Hapke section 10.6, this approximation gives the reflectance at the top of a two-layer system
	// with a finite optically thick top layer and a semi-infinite, optically thick bottom layer
	
	public:
		diffusive_reflectance(double, double); // cumulative tau,surface albedo
		double cumulativetau(){return _cumulativetau;}
		double surfacealbedo(){return _surfacealbedo;}
		double gammasurf(){return _gammasurf;}
		
		double the_function(double);
		
		
	private:
		double _cumulativetau, _surfacealbedo, _gammasurf;
};

class diffusive_reflectance1D : public Jcrap::function<double>
{	// taken from Hapke section 8.7, this approximation gives the reflectance at the top of a 1-layer system
	
	public:
		diffusive_reflectance1D(double); 
		
		double the_function(double);
			
};
