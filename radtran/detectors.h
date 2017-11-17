#include <bitset>
#include "photons.h"

class detector
{
	public:
		detector(geomvector);
			
		virtual void detect(photon)=0;
		virtual void detect(list<photon>);
		virtual void write(string)=0;
		
		void boresight(geomvector);
		
		virtual detector& operator/=(double);
		virtual void keyword(string,string);
		virtual double pixelarea_in_square_km()=0;
		virtual angle iFOV()=0;
		
		geomvector position();
		geomvector boresight();
		
		string name() {return _name;}
		void name(string s) {_name=s;}
		
		cube& CCD();
		cube& geo();
				
	protected:
		geomvector _position;
		geomvector _boresight_direction;
		string _name;
		cube _CCD;
		cube _geo;
};

class colorCCD : public detector
{
	public:
		colorCCD(geomvector, int, angle, vector<value>);
		colorCCD(geomvector, value, value, vector<value>);
		
		void write(string);
		
			
		double km_per_pixel();
		double fieldofview_km();
		double pixelarea_in_square_km();
		angle iFOV();
		void detect(photon);
		

	protected:
		angle _iFOV;  // per pixel
};

class elephant : public colorCCD 
// because it never forgets the entire photon history! 
// elephants have vectors of vectors of ccd cubes 
{   
	public:
		elephant(colorCCD, int=5); // last int=scatter depth
		
		void write(string);
		elephant& operator/=(double);
		void keyword(string,string);
		
		cube& CCD(unsigned int, bitset<64>);
		
		void detect(list<photon>);
		
		
	private:
			
		vector<vector<colorCCD> > _CCD_vv;  // the elephant(s) in the room	
};
		
