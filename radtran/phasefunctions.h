#include "../Jcrap.h"
#include "../Jfunction.h"
#include "photons.h"

class phasefunction 
{
	protected:
		cube _Cumulative_Distribution_Function;
		vector<cube> _inverse_CDFs;
		
		cube _test_distribution_function;
		bool _keeptrack;
		
		cube _integrates_to;
		
		virtual double the_phasefunction(double, angle phase, angle incidence, angle emission) = 0;

	public:
		double normalized_phasefunction(double, angle, angle, angle);
		
		phasefunction();
		void initialize();

		cube& integrates_to();
		double integrates_to(double);
		void precalculate_CDF(double);

		void keeptrack(bool);
		bool keeptrack();
		cube test_distribution_function();
		cube test_distribution_normalized();
		void add_to_test_distribution(geomvector);
		
		virtual string name()=0;
		virtual bool is_surface()=0;
		geomvector randomscatter(const photon);
		virtual cube CDF_calculate(cube)=0;
		
		virtual void check_that_phase_function_integrates(vector<value>)=0;

		vector<cube> invert_CDF(cube&);
		string CDFname();


};

class surface_phasefunction : public phasefunction
{
	protected:
				
		class phase_sinemission : public Jcrap::function<double>
		{
			public:
				phase_sinemission(surface_phasefunction* f, double wavelength_um) : _thephasefunc(f), _wavelength_um(wavelength_um) {}
					
				double the_function(double emission) {
					angle emissionangle(emission, angle::DEGREES);
					angle phase(0., angle::DEGREES);  // not dealing with non-azimuthally symmetric cases yet
					angle incidence(0., angle::DEGREES);  // assuming this doesn't matter so far
					return _thephasefunc->the_phasefunction(_wavelength_um,phase,incidence,emissionangle)*emissionangle.sin()*(2.*Jcrap::pi);
				}
				
				surface_phasefunction *_thephasefunc;
				double _wavelength_um;
		};
	
	public:
		void check_that_phase_function_integrates(vector<value>); 
		cube CDF_calculate(cube);

		bool is_surface() {return true;}
};

class atmospheric_phasefunction : public phasefunction
{
	protected:
		class phase_sintheta : public Jcrap::function<double>
		{
			public:
				phase_sintheta(atmospheric_phasefunction* f, double wavelength_um) : _thephasefunc(f), _wavelength_um(wavelength_um) {}
					
				double the_function(double phase) {
					angle phaseangle(phase, angle::DEGREES);
					angle emission(0., angle::DEGREES);  // presumably this doesn't matter
					angle incidence(0., angle::DEGREES);  // this too
					return _thephasefunc->the_phasefunction(_wavelength_um,phaseangle,incidence,emission)*phaseangle.sin()*(2.*Jcrap::pi);
				}
				
				atmospheric_phasefunction *_thephasefunc;
				double _wavelength_um;
		};

	
	public:			
		void check_that_phase_function_integrates(vector<value>);
		cube CDF_calculate(cube);
		
		bool is_surface() {return false;}
};

class phasefunction_lambertian : public surface_phasefunction
{
		
	public:
		
		phasefunction_lambertian();
			
		double the_phasefunction(double, angle,angle,angle);
		
		string name() {return string("Lambertian");}
};

class phasefunction_surfaceisotropic : public surface_phasefunction
{
	public:
		phasefunction_surfaceisotropic();
	
		double the_phasefunction(double, angle, angle, angle);
		
		string name() {return string("surface_isotropic");}
};

class phasefunction_cube : public atmospheric_phasefunction
{
	private:
		cube _thiscube;
		angle::Jangleunit _radiansordegrees;
	
	public:
		phasefunction_cube(cube, angle::Jangleunit);
		double the_phasefunction(double, angle, angle=angle(0., angle::DEG), angle=angle(0., angle::DEG));
		
		string name();
};

class phasefunction_isotropic : public atmospheric_phasefunction
{

	public:
		phasefunction_isotropic() : atmospheric_phasefunction() {initialize();}
		double the_phasefunction(double, angle phase, angle=angle(0., angle::DEG), angle=angle(0., angle::DEG));
	
		
		string name() {return string("Isotropic");}
};

class phasefunction_forward : public atmospheric_phasefunction
{

	public:
		phasefunction_forward() : atmospheric_phasefunction(){initialize();}
		double the_phasefunction(double, angle phase, angle=angle(0., angle::DEG), angle=angle(0., angle::DEG));
	
		
		string name() {return string("forward");}
};
