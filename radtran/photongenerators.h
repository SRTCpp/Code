#include "photons.h"

class photongenerator
{
	public:
		photongenerator(value, long long, value, vector<value>);
		photongenerator(value=value(2575., "km"), long long=1000,
				value=value(10., "AU"), value=value(0.5, "um"));
		
		vector<value> wavelengths();
		
	// VIRTUAL FUNCTIONS THAT *MUST* BE OVERRIDDEN!!
		virtual photon generatephoton(long long) = 0;
		virtual long long size() = 0;
		virtual double photons_per_square_km() = 0;
	
		
	protected:
		value _lengthscale;
		long long _number;
		value _distance;
		vector<value> _wavelengths;
				
};

class photongenerator_square : public photongenerator
{
	public:
		photongenerator_square(value, long long, value, vector<value>);
			
		photon generatephoton(long long);
		long long size();
		double photons_per_square_km();
		void photondirection(geomvector);// (2 Sept 2016 by smack)
		void fieldcenter(geomvector); // (2 Sept 2016 by smack)
		void inc_offaxis(angle); // INC (12 April 2017 by smack)
		
		void makehole(value); // (24 Aug 2016 by smack)
		bool isinsidehole(geomvector); // (24 Aug 2016 by smack)
				
	private:
		geomvector _fieldcenter;
		geomvector _zvector, _yvector;
		geomvector _photondirection;
		long long _oneside;
		bool _makehole; // (24 Aug 2016, smack)
		geomvector _holevector; // (24 Aug 2016, smack)
		
};

class photongenerator_diffuselamp : public photongenerator
{ 	// the "searchlight" case of a diffuse source on the surface 
	// (smack 2016.11.16)
	
	public:
		photongenerator_diffuselamp(long long, vector<value>, geomvector,geomvector); // n photons, wavelengths, source position, source direction
	
		photon generatephoton(long long);
		geomvector sourcelocation(){return _sourcelocation;}
		long long size();
		double photons_per_square_km(){return size()*Jcrap::pi;} 
		
	private:
		geomvector _sourcelocation;
		geomvector _sourcepointing;
	
};

class photongenerator_laser : public photongenerator 
{  // photons released in a particular direction in a single beam from a location
	// (smack 2016.12.2)
	
	public:
		photongenerator_laser(long long, vector<value>, geomvector,geomvector); //number of phots, vector of wavelengths, position, diretion
		
		photon generatephoton(long long);
		long long size();
		double photons_per_square_km(){return 1;} //this metric doesn't make sense for our ground lamp
																//but needs to be 1 in order to avoid division by zero in the detector normalization
		
		
	private:
		geomvector _sourcelocation;
		geomvector _photondirection;
		
};

class photongenerator_multitarget : public photongenerator
{
	public:
		photongenerator_multitarget(value, long long, value, vector<value>);
			
		photon generatephoton(long long);
		long long size();
		double photons_per_square_km();
				
		void addfieldcenter(geomvector); 
	private:
		vector<photongenerator_square> _target_areas;
		
};
