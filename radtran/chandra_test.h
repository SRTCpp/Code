#include "../Jcrap.h"
#include "../Jfunction.h"
#include "../Jangle.h"
#include "geomvector.h"

class Hfunction : public Jcrap::function<double>
{
	public:
		double ssalbedo;
		double the_function(double);
};

class Hfunction_8_56 : public Jcrap::function<double>  // Hapke eq. 8.56
{
	public:
		double ssalbedo;
		double the_function(double);
};

class testslab : public Jcrap::function<double>
{
	public:
		angle incidence, emission;
		virtual double the_function(double)=0;
};	

class chandraslab : public testslab
{
	public:
		double the_function(double);
		
};

class azimuthintegrator : public Jcrap::function<double>
{
	public:
		azimuthintegrator(double, angle, double, testslab*);
	
		double ssalbedo;
		angle phase;
		double radius;
		testslab* slab;
		
		double the_function(double);
};

class diskintegrator : public Jcrap::function<double>
{
	public:
		diskintegrator(angle, testslab*);
	
		double ssalbedo;
		angle phase;
		testslab* slab;
		
		double the_function(double);
};

class disktest : public Jcrap::function<double>
{
	public:
		disktest(testslab*, angle = angle(0., angle::DEG));
	
		angle phase;
		testslab* slab;	
		double the_function(double);
};
