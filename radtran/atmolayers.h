#include "atmozones.h"

class atmolayer : public vector<atmozone*>
{	
	//
	private:
		double _depth_km;
		double _scaleheight_km;
		phasefunction* _phase_func;
		atmolayer* _toabove;
		atmolayer* _tobelow;
		string _layername;
		
	public:
		/*Constructors*/
		atmolayer(value indepth=value(1., "R_Earth"), phasefunction* infunc=0);
	
		
		/*Getters*/ 
		value depth(){return value(_depth_km,"km");};
		double depth_km(){return _depth_km;};
		value scaleheight(){return value(_scaleheight_km,"km");};
		double scaleheight_km() {return _scaleheight_km;}
		string layername();
		phasefunction* tophasefunction(){return _phase_func;}
		double totaltau(value, geomvector);
		atmolayer* toabove(){return _toabove;};
		atmolayer* tobelow(){return _tobelow;};
		
		/* Setters*/	
		void depth(value indepth){_depth_km=indepth.convert("km");}
		void layername(string);
		void scaleheight(value inscaleheight){_scaleheight_km=inscaleheight.convert("km");};
		void toabove(atmolayer* inabove);
		void tobelow(atmolayer* inbelow);
		void phasefunc(phasefunction *p){_phase_func=p;}
		
		/*Member functions*/
		value bottom(); // calculate bottoms
		double bottom_km();
		value top(); // calculate tops
		double top_km();
		bool issurface();
		static atmolayer default_surface(value=value(2575., "km"));
		void addzone(atmozone *thiszone);
		
};


ostream& operator<<(ostream&, atmolayer&);
  
 
class cubeinterpfunction : public Jcrap::function<double>
{

	
	public:	
		cubeinterpfunction(cube incube, double inlambda): _thiscube(incube),
			_lambda(inlambda) {;}
		
		double the_function(double);
	
		double lambda(){return _lambda;}
		void lambda(double inlambda){_lambda=inlambda;}
		
	
	private:
		cube _thiscube;
		double _lambda;
};
