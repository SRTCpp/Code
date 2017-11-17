#ifndef ATMOZONES_H
#define ATMOZONES_H 1

#include "../Jcrap.h"
#include "../Jfunction.h"
#include "phasefunctions.h"
#include "photons.h"

class atmolayer;

class atmozone
{

	protected:
			

		vector<pair<double, cube> > _ss_albedos;  // pair is double alt_km, albedocube
								// in the cube x=lon, y=lat, z=wave, by default
	
		angle _left_lon, _right_lon, _top_lat, _bottom_lat;
		atmolayer* _toparentlayer;
		
		phasefunction* _zonephasefunction;
		
		static atmozone* _tospace;
		
		void initialsetup();
		
			
	public:					
		/*Constructors*/
		atmozone(phasefunction*, vector<pair<double, cube> >);
		atmozone(phasefunction*, cube);
		atmozone(phasefunction*, double=1.);
		atmozone(vector<pair<double, cube> >);
		atmozone(cube);
		atmozone(double=1.);
	
			
		/*Getters*/ 
		angle left_lon(){return _left_lon;}
		angle right_lon(){return _right_lon;}
		angle top_lat(){return _top_lat;}
		angle bottom_lat(){return _bottom_lat;}
		virtual double ss_albedo_latdeg_londeg_rkm_waveum(double, 
			double, double, double);
		virtual double totaltau(value);
		vector<pair<double, cube> > ss_albedos();
		phasefunction* zonephasefunction();
		
		atmolayer* toparentlayer(){return _toparentlayer;}
		
		/*Setters*/
		void left_lon(angle a){_left_lon=a;}
		void right_lon(angle a){_right_lon=a;}
		void top_lat(angle a){_top_lat=a;}
		void bottom_lat(angle a){_bottom_lat=a;}
		void toparentlayer(atmolayer* thispointer){_toparentlayer=thispointer;}
		void ss_albedo(cube);
		void ss_albedo(double,cube);
		void zonephasefunction(phasefunction*);
		
		//default to DISR08 versions
		virtual double extinction_km_um(double, double);
		virtual double extinction_gas_km_um(double, double);
		virtual double extinction_haze_km_um(double, double)=0;
		
		bool issurface(); 
		string layername();
		
		
		geomvector scatter(photon);
		
		static atmozone& Space();
		
		static cube bitmap_ew0_to_cylmap(cube);
		static cube bitmap_west180_to_cylmap(cube);
		static cube VIMScube_to_albedomap(cube, vector<value>);
		static cube bitmap_to_cylmap(cube, latangle, latangle, lonangle, lonangle);
		static cube bitmap_to_cylmap(cube, latangle, lonangle, angle);

		
};

ostream& operator<<(ostream&, atmozone&);
 
class space_atmozone : public atmozone
{
	public:
		space_atmozone(double);
	
		double extinction_haze_km_um(double, double);
};

class surface_atmozone : public atmozone
{
	public:
		surface_atmozone(double);
		surface_atmozone(cube);
			
		double extinction_haze_km_um(double, double);
};

class DISR_Tomasko_atmozone : public atmozone
{
	private:
		cube _haze_extinction_0; // just haze dTau/dz at bottom of this layer
//		cube _gas_extinction_0; // gas dTau/dz at bottom of this layer
			
	public:
			
		DISR_Tomasko_atmozone(vector<pair<double, cube> >);
		DISR_Tomasko_atmozone(cube);
			
		double extinction_haze_km_um(double, double);
		
	
		void zone_haze_extinction_0(cube c){_haze_extinction_0 = c;}
//		void zone_gas_extinction_0(cube c){_gas_extinction_0 = c;}
	
};

class DISR_Doose_atmozone : public atmozone
{
	private:
			
	public:
		DISR_Doose_atmozone(cube);
			
		double extinction_haze_km_um(double, double);
	
		
		static constexpr  double _a=17.5;
		static constexpr  double _b=6.;
		static constexpr  double _c=45.;
		static constexpr  double _d=55.;
		static constexpr  double _h=1.35;
		static constexpr  double _f=2.;
		static constexpr  double _g=0.1;
};

class exponential_atmozone : public atmozone
{
	cube _extinction_at_bottom;
	
	public:
		exponential_atmozone(cube, cube);
			
		cube extinction_at_bottom(){return _extinction_at_bottom;}
		void extinction_at_bottom(cube e){_extinction_at_bottom=e;}
		
		double extinction_haze_km_um(double, double);
};

class orangerind_atmozone : public atmozone
{
	cube _extinctions;  // z = wave
	public:
			
		orangerind_atmozone(cube);
		orangerind_atmozone(vector<pair<double, cube> >);
			
		double extinction_haze_km_um(double, double);
			
		cube extinctions(){return _extinctions;}
		void extinctions(cube e){_extinctions=e;}
};


#endif
