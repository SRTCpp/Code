#include "atmolayers.h"


atmolayer::atmolayer(value indepth, phasefunction* infunc)
{
 	
	_depth_km=indepth.convert("km");
 	_phase_func=infunc;
	
}


double cubeinterpfunction::the_function (double x){

	return(_thiscube(_lambda,x,0,mLINEAR));

}



value atmolayer::bottom()
{
// each layer needs to calculate its own bottom
	return value(bottom_km(), "km"); 	 
}

double atmolayer::bottom_km()
{
	double answer(0.);
	
	if (tobelow() == 0) answer = 0.;
	else answer = tobelow()->depth_km()+tobelow()->bottom_km();	
	
	return answer;	
}


atmolayer atmolayer::default_surface(value radius)
// default surface
{
	static phasefunction_lambertian phasefunc_l;
	surface_atmozone* surfacezone=new surface_atmozone(1.);
	
	atmolayer answer(radius, &phasefunc_l);
	string name("Default surface of radius ");
	name += radius.tostring();
	answer.layername(name);
		
	answer.addzone(surfacezone);
	
	return answer;
}

value atmolayer::top()
{	
	return value(top_km(), "km");	 
}

double atmolayer::top_km()
{	
	if (tobelow() == 0) return depth();
	else return depth_km()+bottom_km();
}

void atmolayer::addzone(atmozone *thiszone)
{
	push_back(thiszone);
	back()->toparentlayer(this);
}

bool atmolayer::issurface()
{
	bool answer;
	
	if (bottom_km()==0.) answer = true;
	else answer = false;
	
	return answer;
}

string atmolayer::layername()
{
	return _layername;
}

void atmolayer::layername(string inname)
{
	_layername = inname;
}

void atmolayer::toabove(atmolayer* inabove){_toabove = inabove;}

void atmolayer::tobelow(atmolayer* inbelow){_tobelow = inbelow;}

ostream& operator<<(ostream& out, atmolayer& a)
{
	out << "\tAtmolayer number " << &a << ":  " << a.size() << " zones, name " << a.layername() << ":\n";
	out << "\tissurface? " << a.issurface() << "  depth = " << a.depth() << "\n";
	out << "\tscaleheight " << a.scaleheight() << ", bottom is " << a.bottom() << "\n";
	out << "\tphasefunction name is " << a.tophasefunction()->name() << ", address " << a.tophasefunction() << "\n";
	out << "\tZones:\n";
	int n(0);
	for (atmolayer::iterator i(a.begin());i!=a.end();i++) {
		out << "\tZone #" << n << ":\n";
		out << **i;
		n++;
	}
	return out << "\n";
}

double atmolayer::totaltau(value w, geomvector g)
{
	double answer(0.);
	
	answer = at(0)->totaltau(w);
	
	return answer;
}
