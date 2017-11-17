#ifndef PHOTONS_H
#define PHOTONS_H 1

#include"../Jcrap.h"
#include "geomvector.h"

using namespace std;


class atmozone;

class photon
{
	
	private:
		geomvector _photposition;
		geomvector _photdirection;
		double _photaltitude_km; //smack 2017.june.8.;
		double _amplitude;
		double _lambda_um;
		unsigned int atmoscatterers,surfscatterers; //number of scattering events off atmo,surf for each photon
		atmozone* _scattering_zone;
		
		pair<double, double> _ranthetaphi; //for testing groundlamp
		geomvector _randirection; //for testing groundlamp
			
	public:
		/*Constructors and Annihilators*/
		photon(geomvector pos, geomvector direc, 
			double A, value inlambda=value(0.5, "um"));
		
		/*Getters*/
		double amplitude() const {return _amplitude;}
		value lambda(){return value(_lambda_um, "um");}
		double lambda_um() const {return _lambda_um;}
		geomvector position() const {return _photposition;}
		geomvector direction() const {return _photdirection;}
		static double generate_random_tau();
		atmozone* scattering_zone() {
			if (_scattering_zone==0) 
				cout << "oh crap pointer to a scattering zone not set in photon.  Segfault imminent\n";
			return _scattering_zone;}
		pair<double, double> ranthetaphi(){return _ranthetaphi;}
		geomvector randirection(){return _randirection;}	
			
		/*Setters*/
		void amplitude(double A){_amplitude= A ;} //function to get amplitude (private)
		void lambda(value l){ _lambda_um = l.convert("km");}
		void position(geomvector g){_photposition=g;}
		void direction(geomvector g){_photdirection=g;}
		void set_scattering_zone(atmozone* z){_scattering_zone=z;}
		void ranthetaphi(double, double );
		void randirection(geomvector g){_randirection=g;}
				
		/*Member Functions*/
		void setrandtau();
		geomvector project_km(double d);  // project a photon d km along its flight path
		photon scatter(atmozone*);
		
		/*Friend Functions*/
		friend std::ostream& operator<< (std::ostream&, photon&);
				
	
};



#endif
