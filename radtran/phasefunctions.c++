#include "phasefunctions.h"

phasefunction::phasefunction() : _keeptrack(false), _test_distribution_function(1,1,0)
{
	_integrates_to=0.;
	_keeptrack = false;
	_test_distribution_function.keyword("title", "Actual scatter angle distribution(deg)");


}

void phasefunction::initialize()
{
	if (FILE *file = fopen(CDFname().c_str(), "r")) {
		fclose(file);
		_Cumulative_Distribution_Function = cube(CDFname().c_str());
		_inverse_CDFs = invert_CDF(_Cumulative_Distribution_Function);
	}
}


void phasefunction::keeptrack(bool _k)
{
	_keeptrack = _k;
}

bool phasefunction::keeptrack()
{
	return _keeptrack;
}

string phasefunction::CDFname()
{
	string answer(name());
	answer += string(".CDF.Jcube");
	return answer;
}

cube phasefunction::test_distribution_function()
{
	cube answer(_test_distribution_function.histogram(80, -200., 200.));
	answer.keyword("title", "Actual scatter angle distribution(deg)");
	return answer;
}

cube phasefunction::test_distribution_normalized()
{
	cube answer(test_distribution_function());
	answer /= answer.sum()/2.;
	answer /= (answer.Axis(Z,1)-answer.Axis(Z,0));
	answer *= 180./Jcrap::pi/Jcrap::pi;
	for (int i(0);i<answer.N(Z);i++) {
		answer(0,0,i) /= angle(answer.Axis(Z,i), angle::DEGREES).sin();
		if (angle(answer.Axis(Z,i), angle::DEGREES).sin()<1.e-8) answer(0,0,i) = 0.;
	}
	answer.keyword("title", "Per-steradian scatter angle distribution(deg)");
	return answer;
}

void phasefunction::add_to_test_distribution(geomvector scattered_angle)
{	
#pragma omp critical
{
	_test_distribution_function = _test_distribution_function.blocks(Z,cube(1,1,1,scattered_angle.colatitude().degrees()));
}
}

double phasefunction::normalized_phasefunction(double thiswavelength_um, angle phase, angle incidence, angle emission)
// This is the public interface to a phasefunction.  It automatically enforces that each
// phasefunction integrates to 1 by dividing by whatever crazy phase function that a user
// comes up with actually numerically integrates to.  The advantage to this is that any
// (or just about any) phase function input becomes a valid phasefunction that won't lead
// to completely incomprehensible behavior later on.
{	
	if (integrates_to().N()==0) {
		cout << "You haven't run check_that_phase_function_integrates() so I can't normalize\n";
		exit(1);
	}

	double answer(the_phasefunction(thiswavelength_um, phase, incidence, emission));
	
	answer /= integrates_to(thiswavelength_um);

	return answer;
}

cube& phasefunction::integrates_to() {return _integrates_to;}

double phasefunction::integrates_to(double inwav) {return _integrates_to(inwav,0,0);}
 
void phasefunction::precalculate_CDF(double lambda_um)
{
	bool havewedonethisonebefore(false);
	for (int i(0);i<_Cumulative_Distribution_Function.N(X);i++)
		if (lambda_um == _Cumulative_Distribution_Function.Axis(X,i)) 
			havewedonethisonebefore=true;
	if (!havewedonethisonebefore) {
		cube new_CDF(1,1,101);
		new_CDF.Axis(X,0) = lambda_um;
		new_CDF.keyword("scattering_name", name());
				
		cout << "Precalculating CDF at " << lambda_um << "um:  00%";
		new_CDF = CDF_calculate(new_CDF);  // calls different CDF_calculate methods depending if surface_ or atmospheric_phasefunction
		printpercent(10,10);  cout << "\n";
		
		new_CDF.keyword("title", double2str(lambda_um)+name()+string(" CDF"));
		double normalizeby(new_CDF(0,0,new_CDF.N(Z)-1));
		if (osuppress<-3) cout << double2str(lambda_um) << " normalized by " << normalizeby << "\n";
		new_CDF /= normalizeby;
		if (osuppress < -1) new_CDF.graph();
		
		if (_Cumulative_Distribution_Function.N(Z)==1)
			_Cumulative_Distribution_Function = new_CDF;
		else
			_Cumulative_Distribution_Function = _Cumulative_Distribution_Function.blocks(X,new_CDF).dirinc(X);
		
		_Cumulative_Distribution_Function.write(CDFname().c_str());
		
		_inverse_CDFs = invert_CDF(_Cumulative_Distribution_Function);
	}
}
		
cube atmospheric_phasefunction::CDF_calculate(cube templateCDF)
{
	cube answer(templateCDF);

#pragma omp parallel for schedule(dynamic)
	for (int i=0;i<answer.N(Z);i++) {			
		phase_sintheta inphasefunc(this, templateCDF.Axis(X,0));
		double thisangle_deg(180.*double(i)/double(answer.N(Z)-1));
		answer.Axis(Z,i)=thisangle_deg;
		answer(0,0,i) = inphasefunc.integrate(0.,thisangle_deg, 1.e-8, 16);
		printpercent(i,answer.N(Z));	
	}
	
	return answer;
}
		
cube surface_phasefunction::CDF_calculate(cube templateCDF)
{
	cube answer(templateCDF);
	
#pragma omp parallel for schedule(dynamic)
	for (int i=0;i<answer.N(Z);i++) {
		phase_sinemission inphasefunc(this, templateCDF.Axis(X,0));	
		double thisangle_deg(90.*double(i)/double(answer.N(Z)-1));
		answer.Axis(Z,i)=thisangle_deg;
		answer(0,0,i) = inphasefunc.integrate(0.,thisangle_deg, 1.e-8, 16);	
		printpercent(i,answer.N(Z));	
	}
	
	return answer;
}
		
phasefunction_lambertian::phasefunction_lambertian(){initialize();}

double phasefunction_lambertian::the_phasefunction(double, angle,angle,angle in_emission){
	
	return in_emission.cos()/Jcrap::pi;

}

phasefunction_surfaceisotropic::phasefunction_surfaceisotropic(){initialize();}

double phasefunction_surfaceisotropic::the_phasefunction(double,angle,angle,angle){
	
	return 1./(2.*Jcrap::pi);

}

phasefunction_cube::phasefunction_cube(cube thiscube, angle::Jangleunit thisunit)
{	
	_thiscube= thiscube;
	_radiansordegrees= thisunit;
	initialize();
}
	
string phasefunction_cube::name()
{
	return _thiscube.keyword("scattering_name");
}

void atmospheric_phasefunction::check_that_phase_function_integrates(vector<value> these_wavelengths)
// All atmospheric phasefunctions need necessarily to integrate to 1 when integrated over 4pi steradians
{
	
	cube integrates_to(these_wavelengths.size(), 1, 1);
	for (int i(0);i<integrates_to.N(X);i++) {
		integrates_to.Axis(X,i) = these_wavelengths.at(i);
		
		phase_sintheta inphasefunc(this, double(these_wavelengths.at(i).convert("um")));
		
		double answer(inphasefunc.integrate(0., 180., 1.e-6,12));
		if (osuppress < -1) inphasefunc.plot(0., 180., 100).graph();
		answer *= Jcrap::pi/180.;
		if (osuppress < -1) cout << "This phase function, " << inphasefunc._thephasefunc->name() << " at " << integrates_to.Axis(X,i);
		if (osuppress < -1) cout << " integrates to " << answer << "\n";
		integrates_to(i,0,0) = answer;
	}
	if (osuppress < -1) cout << "cube _integrates_to: " << integrates_to << "\n";
	_integrates_to = integrates_to;
}


void surface_phasefunction::check_that_phase_function_integrates(vector<value> these_wavelengths)
{
	cout << "PHASEFUNCTION CHECKS FOR SURFACE_PHASEFUNCTION STILL ASSUME AZIMUTHAL SYMMETRY\n";
	
	
	cube integrates_to(these_wavelengths.size(), 1, 1);
	for (int i(0);i<integrates_to.N(X);i++) {
		integrates_to.Axis(X,i) = these_wavelengths.at(i);
		
		phase_sinemission inphasefunc(this, double(these_wavelengths.at(i).convert("um")));
		
		double answer(inphasefunc.integrate(0., 90., 1.e-6,12));
		if (osuppress < -1) inphasefunc.plot(0., 180., 100).graph();
		answer *= Jcrap::pi/180.;
		if (osuppress < -1) cout << "This phase function, " << inphasefunc._thephasefunc->name() << " at " << integrates_to.Axis(X,i);
		if (osuppress < -1) cout << " integrates to " << answer << "\n";
		integrates_to(i,0,0) = answer;
	}
	if (osuppress < -1) cout << "cube _integrates_to: " << integrates_to << "\n";
	_integrates_to = integrates_to;
	
	
	_integrates_to = cube(1,1,1,1.0);
}

vector<cube> phasefunction::invert_CDF(cube& CDF)
{
	vector<cube> answer;
	
	for (int w(0);w<CDF.N(X);w++) {
		cube inverse_CDF = CDF.plane(X,0) * 0.;
		for (int i(0);i<inverse_CDF.N(Z);i++) {
			inverse_CDF.Axis(Z,i) = CDF(w,0,i);
			inverse_CDF(0,0,i) = _Cumulative_Distribution_Function.Axis(Z,i);
		}
		answer.push_back(inverse_CDF);
	}
	
	return answer;
}

double phasefunction_cube::the_phasefunction(double thiswavelength_um, angle a1, angle, angle){
	
	//understood that incoming cube Y = phase, X = wavelength
	double answer;	
	double convertedangle(a1.as(_radiansordegrees).number());
	
	int w(_thiscube.Jclosest_dumb(X, thiswavelength_um));
	
	answer=_thiscube(w,convertedangle,0,mSLOW);
	
	return answer;	
}

geomvector phasefunction::randomscatter(const photon inphoton)
{
	double lambda_um(inphoton.lambda_um());
	
	double ran_1=cube::Jrandom(), ran_2=cube::Jrandom();
	
	if (osuppress<-2)cout << "_CDF info:\n" << _Cumulative_Distribution_Function.info() << "\n"; cout.flush();
	int w(-1);
	for (int i(0);i<_Cumulative_Distribution_Function.N(X);i++)
		if (_Cumulative_Distribution_Function.Axis(X,i) == lambda_um)
			w = i;
	if (w==-1) {
		// okay, have to go generate a new CDF cube for this precise wavelength
		precalculate_CDF(lambda_um);
		w = _Cumulative_Distribution_Function.N(X)-1;
	}
	
	if (osuppress<-4)cout << "Jclosest_dumb found to be " << w << "\n";
	
	
	angle colatitude(_inverse_CDFs.at(w)(0,0,ran_1,mSLOW), angle::DEGREES);
	if (osuppress<-4)cout << "about to go from colatitude(" << colatitude << "; from ran1=" << ran_1 << ") to latitude\n";
	angle latitude(angle(90., angle::DEGREES)-colatitude);
	angle phi(ran_2*360., angle::DEGREES);
	
	geomvector answer(geomvector::geolatitude(latitude), geomvector::geolongitude(phi), 1.);
	if (osuppress<-2) cout << "This photon experiences a scatter directional deviation of " << answer << "\n";
	return answer;
}
		
double phasefunction_isotropic::the_phasefunction(double, angle phase, angle, angle)
{
	double answer(1.);
	
	return answer;
}

double phasefunction_forward::the_phasefunction(double, angle phase, angle, angle)
{
	double answer(0.);
	if (phase < angle(60., angle::DEGREES)) answer = 4.;
	return answer;  // should integrate to 2
}
