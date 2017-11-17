#include "SRTC++_testsuite.h"

cube diskintegrate(cube incube, double planetradius_in_km)
{
	cube answer(incube.skewer(Z,0,0)*0.);
	
	for (int z(0);z<answer.N(Z);z++) {
		int n(0);
		for (int x(0);x<incube.N(X);x++) {
			for (int y(0);y<incube.N(Y);y++) {
				double d(incube.Axis(X,x)*incube.Axis(X,x) + 
					incube.Axis(Y,y)*incube.Axis(Y,y));
				d = sqrt(d);
				if (d < planetradius_in_km) {
					n++;
					answer(0,0,z) += incube(x,y,z);
				}
			}
		}
		answer(0,0,z)/= double(n);
	}
	
	return answer;
}

double test_LambertSurface()
{
// micro-thin black atmosphere so as not to obstruct surface (much)
	atmosphere a(atmosphere::orangerind(1.e-6, value(30., "km"), 0.));

// have to put in a wavelength, tho doesn't really matter in this case		
	vector<value> wavelengths;
	wavelengths.push_back(value(2.0,"um"));
	
// setting up the surface
	cube albedoone(1,2,1);
	albedoone(0,0,0)=1.0; 
	albedoone.Axis(Y,0)=-1.;
	albedoone(0,1,0)=1.0; 
	albedoone.Axis(Y,0)=1.;  // surface albedo is all one.
	atmozone *uniformsurface=new surface_atmozone(albedoone);
 	a.front().addzone(uniformsurface);

// setting up the Sun	
	photongenerator_square hv(value(4200., "km"), 42*qualityfactor, value(4200., "km"), wavelengths);	
	geomvector incomingdirection(-1.,0.,0.);  // set the photons to come in initially toward -x
	hv.photondirection(incomingdirection);
	hv.fieldcenter(geomvector(4200.,0.,0.));
	
// Detector		
	vector<detector*> Detectors;
	double x(1.e9); 
	double y(0.);
	colorCCD *CCD_detector;
	CCD_detector = new colorCCD(geomvector(x,y,0.),value(3200., "km"),value(100., "km"),wavelengths);
	Detectors.push_back(CCD_detector);
	
// RUN SRTC++!
	SRTC S(a, &hv, Detectors);
	S.run();	
	CCD_detector->CCD().write("LambertSurface.Jcube");
	
// Analyze result
	return diskintegrate(CCD_detector->CCD())(0,0,0);
}

