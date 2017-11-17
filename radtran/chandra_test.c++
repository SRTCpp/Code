#include "chandra_test.h"

double Hfunction::the_function(double mu) 
{
	double gamma(sqrt(1.-ssalbedo));
	if (!(gamma>0. || gamma<=0.)) gamma=0.;
	double answer(1.+2.*mu);
	answer /= 1+2*mu*gamma;
	return answer;
}

double Hfunction_8_56::the_function(double mu) // Hapke eq. 8.56
{
	double gamma(sqrt(1.-ssalbedo));
	if (!(gamma>0. || gamma<=0.)) gamma=0.;
	
	double r0((1.-gamma)/(1.+gamma));
	
	
	double answer(1.-ssalbedo*mu* 
		(r0 + (1.-2.*r0*mu) /2. *log((1+mu)/mu) )  );
	answer = 1./answer;
	return answer;
}

double chandraslab::the_function(double ssalbedo) 
{
	double answer(ssalbedo/4.);
	double mu_0(incidence.cos());
	double mu(emission.cos());
	answer *= mu_0/(mu+mu_0);
	Hfunction_8_56 H;
	H.ssalbedo=ssalbedo;
	answer *= H(mu)*H(mu_0);
			
	return answer;
}

azimuthintegrator::azimuthintegrator(double _a, angle _p, double _r, testslab* _slab) : 
		ssalbedo(_a), phase(_p), radius(_r), slab(_slab) {}

double azimuthintegrator::the_function(double az_radians) 
{
	angle azimuth(az_radians, angle::RADIANS);
	double answer(0.);
			
	geomvector subsolarpoint(geomvector::geolatitude(angle(0., angle::DEG)), 
			geomvector::geolongitude(angle(0., angle::DEG)), 1.);
	geomvector subspacecraftpoint(geomvector::geolatitude(angle(0., angle::DEG)),
			geomvector::geolongitude(phase), 1.); 
	
	slab->emission = angle(asin(radius), angle::RADIANS);
	slab->incidence = angle_acos(slab->emission.cos()*phase.cos() +
									slab->emission.sin()*phase.sin()*azimuth.cos());
	
	if (slab->incidence.degrees()>90.)
		answer = 0.;
	else
		answer = slab->the_function(ssalbedo);
	return answer;
}


diskintegrator::diskintegrator(angle _p, testslab* _slab) : 
		phase(_p), slab(_slab) {}

double diskintegrator::the_function(double radius) 
{
	double answer(0.);
	
	azimuthintegrator aint(ssalbedo, phase, radius, slab);
	
	if (phase.degrees()==0.) {
		answer = (aint(0.)*2.*Jcrap::pi*radius);
	} else {
		answer = aint.integrate(0., 2.*Jcrap::pi)*radius;
	}
		
	return answer;
}

disktest::disktest(testslab* inslab, angle inphase) : 
		phase(inphase), slab(inslab) {}

double disktest::the_function(double ssalbedo) 
{
	diskintegrator d(phase,slab);
	d.ssalbedo = ssalbedo;
	return d.integrate(0., 1.)/Jcrap::pi;
}
