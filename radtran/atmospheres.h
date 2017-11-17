#include "atmolayers.h"

class atmosphere : public list<atmolayer>
{
	
	private:		
		
	public:
		/*Constructors*/
		atmosphere();
	
		// here are some ways to generate canned atmospheres
		static atmosphere orangerind(cube, value, cube);
		static atmosphere orangerind(double, value=value(300., "km"), double=1.);
		static atmosphere orangerind(double, value, cube);
		static atmosphere SebastienValidateDISR08();
		static atmosphere SebastienValidateDISR08_layered();
		static atmosphere DISR08();
		static atmosphere DISR(); 

		/*Getters*/
		double totaltau(value wavelength, geomvector latlonloc=geomvector(2575.,
			geomvector::geolatitude(angle(0., angle::DEG)), 
			geomvector::geolongitude(angle(0., angle::DEG))));
		atmolayer* layerataltitude(value);
		
		/*Setters*/
		
		/*Member Functions*/
		value outeredge(); //return the other edge of the atmosphere
		double outeredge_km(); //return the other edge of the atmosphere
		value surfacealtitude();
		double surfacealtitude_km(); // not implimented yet
		void addlayer(atmolayer);
		
		
		/*Atmosphere test functions*/
		void normal_extinction();
		cube extinction_profile(photon&);	
};

ostream& operator<<(ostream&, atmosphere&);

class DISR08_tau_function : public Jcrap::function<double>
{
	public:
		DISR08_tau_function(double, double, value, value);
	
		double the_function(double);		
	
	private:
			
		double _constant;
		double _exponent;
		double _depth_km;
		double _scaleheight_km;
};



class DISR_tau_function : public Jcrap::function<double>
{
	public:
		DISR_tau_function(value);
	
		double the_function(atmozone*,double);	
			
		
		void calculate_scaleheight_fromtau(atmozone*,double);
		
		
	private:
		double _depth_km;
		double _scaleheight_km;	
		
	
		//Tau function parameters: see Table 2 of Doose et al. 2015	
		static constexpr double _a=17.5;
		static constexpr double _b=6.;
		static constexpr double _c=45.;
		static constexpr double _d=55.;
		static constexpr double _h=1.35;
		static constexpr double _f=2.;
		static constexpr double _g=0.1; 
};
