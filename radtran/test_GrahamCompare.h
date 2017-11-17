#include "SRTC++_testsuite.h"

static double accuracy(1.e-4); 
static int min_number(6);

class Graham10 : public Jcrap::function<double>
{
	public:
		atmosphere* a;
		angle incidence, emission, phase;
		double wavelength_um;
		
		double the_function (double theta_prime_rad);
};


class Graham10_azimuthal : public Jcrap::function<double>
{
	public:
		Graham10 *parent;
		angle theta_prime;
	
		double the_function(double azimuthal_angle_rad);
};


class Graham01 : public Jcrap::function<double>
{
	public:
		atmosphere* a;
		angle incidence, emission, phase;
		double wavelength_um;
		bool makegraphs;
		
		double the_function (double theta_prime_rad);
};


class Graham01_azimuthal : public Jcrap::function<double>
{
	public:
		Graham01 *parent;
		angle theta_prime;
	
		double the_function(double azimuthal_angle_rad);
};

