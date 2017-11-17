#ifndef GEOMVECTOR_H
#define GEOMVECTOR_H 1

#include"../Jangle.h"

using namespace std;

enum coordinatesystem{CARTESIAN=0,SPHERICALPOLAR=1,XYZ=0,CART=0,RECT=0,POLAR=1,SPHERICAL=1,
		SPHERE=1,RTHETAPHI=1};



// Note that member functions are declared in the class construc,
// but are defined in the .c++ file
class geomvector
{
	
	private:
		double _x,_y,_z;   //centered on center of planet
									// X=to the right
									// Y=away from observer
									// Z=up
		angle _lat, _lon;  	// lat measured northward from equator
									// lon measured CCW from X axis as viewed from NP (East Longitude)

		double _r;
		coordinatesystem _system;
		 
		
	public:	
		class geolatitude : public angle {
			public:
				geolatitude(angle a) : angle(a) {};
		};
		class geolongitude : public angle {
			public:
				geolongitude(angle a): angle(a) {};
		};
		
		static geolatitude geocolatitude(angle);
		
		/*Constructors*/
		geomvector(double=0., double=0., double=0.);
		geomvector(geolatitude, geolongitude, double);
		geomvector(geolongitude, geolatitude, double);
		geomvector(double, geolatitude, geolongitude);
		geomvector(double, geolongitude, geolatitude);
			
		/*Getters*/
		double x() const;
		double y() const;
		double z() const;
		angle lat() const;
		angle latitude() const;
		angle colatitude() const;
		angle longitude() const;
		angle lon() const;
		double r() const;
		
		/*Get simple modifications*/
		geomvector unitvector() const;
		
		/*Setters*/
		/* NONE -- on PURPOSE!  These would give undefined behavior depending
		on the _system, so if you want to do something funny update your own
		damned components by creating a new geomvector*/
		
		/*Functions*/
		geomvector as(coordinatesystem) const;
		geomvector to(coordinatesystem);
		
		geomvector transform_to(const geomvector) const;  // transform to and from surface coordinate systems
		geomvector transform_from(const geomvector) const; // i.e. rotate *this as at the input lat/lon so that the local zenith points up the Z axis

		angle angulardistance(const geomvector) const;
		angle angulardistancel(const geomvector) const;  // uses long doubles for precision
				
		/*Operators*/
		geomvector operator-() const;
		geomvector operator*(double) const;
		geomvector operator/(double) const;
		geomvector operator+(geomvector) const;
		geomvector operator-(geomvector) const;
		
		/*Would be an Operator but it's ambiguous*/
		double dot(const geomvector) const;
		
		/*Manipulators*/
		geomvector euler_ZXZ(angle, angle, angle) const;
		geomvector euler(angle, angle, angle) const;
		
		friend std::ostream& operator<< (std::ostream&, const geomvector&);
			
		
};

#endif
