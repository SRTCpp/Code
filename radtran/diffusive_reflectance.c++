#include "diffusive_reflectance.h"

diffusive_reflectance::diffusive_reflectance(double totaltau, double surfacealbedo) :
	_cumulativetau(totaltau), _surfacealbedo(surfacealbedo), _gammasurf(sqrt(1-surfacealbedo))
{
}


double diffusive_reflectance::the_function(double ssa)
{
	double answer(0.);
	double gammaatmo(sqrt(1-ssa));
	double r_U((1-gammaatmo)/(1+gammaatmo)),r_L((1-gammasurf())/(1+gammasurf()));
	if (surfacealbedo() <=0.5) r_L=0; //low albedo case
	answer=1 + (1/r_U)*(r_L-r_U)/(1-r_L*r_U)*exp(-4*gammaatmo*cumulativetau());
	answer/=1+r_U*(r_L-r_U)/(1-r_L*r_U)*exp(-4*gammaatmo*cumulativetau());
	answer*=r_U;
	return answer;//translating to I/F
}

diffusive_reflectance1D::diffusive_reflectance1D(double tau)
{
	if (tau<1.) cout <<"Note that 1D diffusive reflectance ASSUMES medium of infinite thickness\n";
}

double diffusive_reflectance1D::the_function(double ssa)
{
	double answer(0.);
	double gammaatmo=sqrt(1-ssa);
	answer=(1-gammaatmo)/(1+gammaatmo);
	return answer;
}
