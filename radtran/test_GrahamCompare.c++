#include "test_GrahamCompare.h"
 
extern int graphno;
const double acceptancecone(5.); // degrees from specified geom that we accept

double Graham10_azimuthal::the_function(double azimuthal_angle_rad) 
{
	double answer;
	
	angle azimuthal_angle(azimuthal_angle_rad, angle::RADIANS);
	double emitted_prime_atmos(1./theta_prime.cos());
	double incident_atmos(1./parent->incidence.cos());
	
	geomvector::geolongitude observation_azimuth(angle_acos(parent->phase.cos()-(parent->incidence.cos()*parent->emission.cos())  /  
			parent->incidence.sin()*parent->emission.sin()));	
	geomvector emission_vec(geomvector::geocolatitude(parent->emission), geomvector::geolongitude(observation_azimuth), 1.);
	
	geomvector::geolongitude Sun_azimuth(angle(0., angle::DEGREES));  // obs_az above defined relative to Sun az, so this is zero by convention
	geomvector incidence_vec(geomvector::geocolatitude(parent->incidence), Sun_azimuth, 1.);
	
	geomvector Omega_prime_vec(geomvector::geocolatitude(theta_prime), geomvector::geolongitude(azimuthal_angle), 1.);
	
	angle phase_incident_to_prime(incidence_vec.angulardistancel(Omega_prime_vec));
	angle phase_prime_to_emitted(emission_vec.angulardistancel(Omega_prime_vec));
	
	angle doesntmatter(0., angle::RADIANS);
	
	answer = parent->a->front().back()->zonephasefunction()->normalized_phasefunction(value(parent->wavelength_um, "um"), 
			phase_prime_to_emitted, theta_prime, parent->emission);  // surface scatter
	answer*= incident_atmos / (incident_atmos-emitted_prime_atmos);
	answer*= parent->a->back().front()->zonephasefunction()->normalized_phasefunction(value(parent->wavelength_um, "um"), 
			phase_incident_to_prime, doesntmatter, doesntmatter);   // atmospheric scatter
	answer*= exp(-(parent->a->totaltau(value(parent->wavelength_um, "um"))*emitted_prime_atmos)) -
			   exp(-(parent->a->totaltau(value(parent->wavelength_um, "um"))*incident_atmos));

	return answer;
}

double Graham10::the_function (double theta_prime_rad)
{  // theta_prime is zenith angle as viewed from surface
	double answer;
	
	Graham10_azimuthal az_func;
	az_func.parent = this;
	az_func.theta_prime = angle(theta_prime_rad, angle::RADIANS);

/*	cube graph10_az(az_func.plot(0., 2.*Jcrap::pi, 100));
	graph10_az.keyword("title", (string("10 az ") +
			double2fstr(az_func.theta_prime.degrees(), 5)).c_str());
	graph10_az.graph();*/
	answer = az_func.theta_prime.sin()*az_func.integrate(0., 2.*Jcrap::pi, accuracy, min_number, "Graham10 internal azimuthal");

	return answer;
}

double Graham01_azimuthal::the_function(double azimuthal_angle_rad) 
{
	double answer;
	
	angle azimuthal_angle(azimuthal_angle_rad, angle::RADIANS);
	double emitted_prime_atmos(1./theta_prime.cos());
	double emitted_atmos(1./parent->emission.cos());
	
	angle phase_incident_to_prime(parent->phase); // phase is wrong but won't matter for Lambertian -- update later to full spherical geometry if non-Lambertian
	
	double azimuth_calculation_argument((parent->phase.cos()-parent->incidence.cos()*parent->emission.cos())  /  
			(parent->incidence.sin()*parent->emission.sin()));
	geomvector::geolongitude observation_azimuth(angle_acos(azimuth_calculation_argument));
	if (isnan(observation_azimuth.degrees())) {
//		cout << "Observation_azimuth NaN -- setting to zero\n";
//		cout << "incidence " << parent->incidence << "   emission " << parent->emission << "   phase " << parent->phase << "  argument " << 
//				(parent->phase.cos()-parent->incidence.cos()*parent->emission.cos())  /  (parent->incidence.sin()*parent->emission.sin()) 
//				<< "\n";
		if (azimuth_calculation_argument > 1.)
			observation_azimuth = angle(0., angle::DEGREES);
		else if (azimuth_calculation_argument < -1.)
			observation_azimuth = angle(180., angle::DEGREES);
		else observation_azimuth = angle(0., angle::DEGREES);
	}
	if(parent->makegraphs) cout << "observation_azimuth is " << observation_azimuth.as(angle::DEGREES) << "\n";
	geomvector emission_vec(geomvector::geocolatitude(parent->emission), geomvector::geolongitude(observation_azimuth), 1.);	
	if(parent->makegraphs) cout << "emission_vec is " << emission_vec << "\n";
	geomvector Omega_prime_vec(geomvector::geocolatitude(theta_prime), geomvector::geolongitude(azimuthal_angle), 1.);	
	if(parent->makegraphs) cout << "Omega_prime_vec is " << Omega_prime_vec << "\n";
	angle phase_prime_to_emitted(emission_vec.angulardistance(Omega_prime_vec));
	
//	cout << "In Graham01, phase_prime_to_emitted is " << phase_prime_to_emitted.as(angle::DEGREES) << "\n";
//	cout << "when emission equals " << emission_vec << " and Omega_prime equals " << Omega_prime_vec << "\n";
	
	angle doesntmatter(0., angle::RADIANS);
	
	answer = parent->a->front().back()->zonephasefunction()->normalized_phasefunction(value(parent->wavelength_um, "um"), 
			phase_incident_to_prime, parent->incidence, theta_prime);  // surface scatter
	answer*= emitted_prime_atmos / (emitted_atmos-emitted_prime_atmos);
	answer*= parent->a->back().front()->zonephasefunction()->normalized_phasefunction(value(parent->wavelength_um, "um"), 
			phase_prime_to_emitted, doesntmatter, doesntmatter);  // atmospheric scatter
//	if(parent->makegraphs) {
//		cout << "Well, for theta_prime " << theta_prime.as(angle::DEGREES) << " the atmo phase function is " << parent->a->back().front()->zonephasefunction()->normalized_phasefunction(value(parent->wavelength_um, "um"), 
//			phase_prime_to_emitted, doesntmatter, doesntmatter) << " for phase_prime_to_emitted " << phase_prime_to_emitted.as(angle::DEGREES) << "\n";
//	}
	answer*= exp(-(parent->a->totaltau(value(parent->wavelength_um, "um"))*emitted_prime_atmos)) -
			   exp(-(parent->a->totaltau(value(parent->wavelength_um, "um"))*emitted_atmos));

	return answer;
}

double Graham01::the_function (double theta_prime_rad)
{  // theta_prime is zenith angle as viewed from surface
	double answer;
	
	Graham01_azimuthal az_func;
	az_func.parent = this;
	az_func.theta_prime = angle(theta_prime_rad, angle::RADIANS);
	
	answer = az_func.theta_prime.sin()*az_func.integrate(0., 2.*Jcrap::pi, accuracy, min_number, "Graham01 internal azimuthal");

	if (makegraphs) { 
		cube azgraph01(az_func.plot(0., 2.*Jcrap::pi, 100));
		azgraph01.keyword("title", (string("Azplot01 theta ")+int2str(int(az_func.theta_prime.degrees()))));
		azgraph01.graph();
	}
	
	return answer;
}

double test_GrahamCompare()
{
// phase angles to test
	vector<angle> emissions;
	emissions.push_back(angle( 0., angle::DEGREES));
	emissions.push_back(angle(10., angle::DEGREES));
	emissions.push_back(angle(20., angle::DEGREES));
	emissions.push_back(angle(30., angle::DEGREES));
	emissions.push_back(angle(40.01, angle::DEGREES));
	emissions.push_back(angle(50., angle::DEGREES)); 
	emissions.push_back(angle(60., angle::DEGREES));
	emissions.push_back(angle(70., angle::DEGREES));
	emissions.push_back(angle(80., angle::DEGREES));
	emissions.push_back(angle(90., angle::DEGREES));
	emissions.push_back(angle(100., angle::DEGREES));
	emissions.push_back(angle(110., angle::DEGREES));
	emissions.push_back(angle(120., angle::DEGREES));
	emissions.push_back(angle(130., angle::DEGREES));
	emissions.push_back(angle(140., angle::DEGREES));
	emissions.push_back(angle(150., angle::DEGREES));
	emissions.push_back(angle(160., angle::DEGREES));
	emissions.push_back(angle(170., angle::DEGREES));

	bool backwards(false);
	
// 
	vector<value> wavelengths;
	wavelengths.push_back(value(0.933078,"um"));
	wavelengths.push_back(value(1.08183,"um"));
	wavelengths.push_back(value(1.27813,"um"));
	wavelengths.push_back(value(1.59018, "um"));
	wavelengths.push_back(value(2.01781,"um"));
	wavelengths.push_back(value(2.69620,"um"));
	wavelengths.push_back(value(2.79889,"um"));
	wavelengths.push_back(value(5.00576, "um"));
	
	vector<value> backwards_wavelengths;
	cout << "b_w size is " << backwards_wavelengths.size() << "\n";
	for (int i(wavelengths.size()-1);i>=0;i--) {
		backwards_wavelengths.push_back(wavelengths.at(i));
		cout << "b_w size is " << backwards_wavelengths.size() << ", i = " << i << "\n";
	}
	if (backwards) wavelengths = backwards_wavelengths;
	cout << "b_w size is " << backwards_wavelengths.size() << "\n";
	
	
// Comparing to a generic, Titan-ish orangerind just to test
	cube taus(1,1,wavelengths.size());
	taus(0,0,0) = 3.2;
	taus(0,0,1) = 2.5;
	taus(0,0,2) = 2.1;
	taus(0,0,3) = 1.45;
	taus(0,0,4) = 1.02;
	taus(0,0,5) = 0.8;
	taus(0,0,6) = 0.8;
	taus(0,0,7) = 0.3;
	for (int i(0);i<taus.N(Z);i++) taus.Axis(Z,i) = double(wavelengths.at(i).convert("um"));
	

	cube backtau(1,1,0);
	for (int i(taus.N(Z)-1);i>=0;i--) backtau = backtau.blocks(Z,taus.plane(Z,i));
	if (backwards) taus = backtau;
	
//	taus*=0.; taus+=1.0;
	
	cube ssalbedos(taus);
	ssalbedos(0,0,0) = 1.00;
	ssalbedos(0,0,1) = 0.99;
	ssalbedos(0,0,2) = 0.98;
	ssalbedos(0,0,3) = 0.96;
	ssalbedos(0,0,4) = 0.77;
	ssalbedos(0,0,5) = 0.507;
	ssalbedos(0,0,6) = 0.463;
	ssalbedos(0,0,7) = 0.4998;
	
//	ssalbedos*=0.;  ssalbedos+=1.0;
	
	cout << "About to create atmosphere\n";
	atmosphere a1(atmosphere::orangerind(taus, value(10., "km"), ssalbedos));
	cout << "After atmo create\n";
	
	cube hpf_Above(cube("/vims/DISRhazephasedata_above80km.Jcube"));
	hpf_Above.keyword("scattering_name", "hpf_Above");
	phasefunction_cube hpfcube_Above(hpf_Above, angle::DEGREES);
	
	cube hpf_Below(cube("/vims/DISRhazephasedata_below80km.Jcube"));
	hpf_Below.keyword("scattering_name", "hpf_Below");		
	phasefunction_cube hpfcube_Below(hpf_Below, angle::DEGREES);
	
	cube Below_but_not_wavedependent(hpf_Below.plane(X,0.93));
	Below_but_not_wavedependent.Axis(X,0) = wavelengths.at(0);
	for (int i(0);i<wavelengths.size();i++) {
		Below_but_not_wavedependent=Below_but_not_wavedependent.blocks(X,hpf_Below.plane(X,0.93));
		Below_but_not_wavedependent.Axis(X,i) = wavelengths.at(i);
	}
	Below_but_not_wavedependent.keyword("scattering_name", "Below_but_not_wavedependent_0.93um");
	phasefunction_cube hpfcube_Below_but_not_wavedependent(Below_but_not_wavedependent, angle::DEGREES);	
			
	cube hpf_allforward("/vims/forwardscatteringhazephasedata.txt.Jcube");
	hpf_allforward.keyword("scattering_name", "allforward_cube");	
	hpf_allforward += 2.*0.;
	phasefunction_cube forwardscattering(hpf_allforward, angle::DEGREES);
	
	cube hpf_allbackward(hpf_allforward);
	for (int x(0);x<hpf_allforward.N(X);x++) {
		for (int i(0);i<hpf_allforward.N(Y);i++) {
			hpf_allbackward(x,i,0) = hpf_allforward(x,hpf_allforward.N(Y)-1-i,0);
		}
	}
	hpf_allbackward.keyword("scattering_name", "allbackward_cube");	
	hpf_allbackward.write("/vims/backwardscatteringhazephasedata.txt.Jcube");
	phasefunction_cube backscattering(hpf_allbackward, angle::DEGREES); 

	cube hpf_isotropic("/vims/isotropichazephasedata.txt.Jcube");
	hpf_isotropic.keyword("scattering_name", "isotropic_cube");
	phasefunction_cube isotropicphasefunction_cube(hpf_isotropic, angle::DEGREES);
	
	phasefunction_isotropic isotropicphasefunction;		// usual analytic isotropic		
	phasefunction_surfaceisotropic surfaceisotropicphasefunction;
	phasefunction_lambertian surfaceLambertian;
	
	a1.back().at(0)->zonephasefunction(&hpfcube_Below);
	
	a1.back().at(0)->zonephasefunction()->check_that_phase_function_integrates(wavelengths);
	cube phasefunctiongraph(wavelengths.size(),1,2001,0.);
	for (int w(0);w<wavelengths.size();w++) {
		for (int i(0);i<phasefunctiongraph.N(Z);i++) {
			phasefunctiongraph.Axis(Z,i) = double(i)/phasefunctiongraph.N(Z)*180.;
			angle phase(phasefunctiongraph.Axis(Z,i), angle::DEG);
			angle incidence(phase);
			angle emission(phase);
			phasefunctiongraph(w,0,i) = a1.back().at(0)->zonephasefunction()->normalized_phasefunction(double(wavelengths.at(w)), phase, incidence, emission);
//			cout << "phase function is " << a1.back().at(0)->zonephasefunction()->name() << " " << a1.back().at(0)->zonephasefunction() << "\n";
		}
	}
	phasefunctiongraph *= 4.*Jcrap::pi;
	phasefunctiongraph.keyword("title", "Phasefunction of 1st atmospheric layer's zone");
	phasefunctiongraph.graph();
	phasefunctiongraph.write("phasefunctiongraph.Jcube");
	
// setting up the surface
	cube one(1,1,1); 
	one(0,0,0) =1.0;
	atmozone *uniformsurface1=new surface_atmozone(one);
 	a1.front().addzone(uniformsurface1);
	a1.front().at(0)->zonephasefunction(&surfaceLambertian);
	a1.front().at(1)->zonephasefunction(&surfaceLambertian);
	
	
	cout << "Atmosphere for GrahamCompare1: " << a1 << "\n";

// setting up the Sun	
	photongenerator_square hv(value(4200., "km"), 420.*qualityfactor, value(4200., "km"), wavelengths);	
	geomvector incomingdirection(-1.,0.,0.);  // set the photons to come in initially toward -x
	hv.photondirection(incomingdirection);
	hv.fieldcenter(geomvector(4200.,0.,0.));
	
// Detector		
	vector<detector*> Detectors1;
	for (int i(0);i<emissions.size();i++) {
		angle phase(emissions.at(i));
		double x(1.e9*phase.cos()); 
		double y(1.e9*phase.sin());

		colorCCD *CCD_detector1;
		CCD_detector1 = new elephant(colorCCD(geomvector(x,y,0.),value(3200., "km"), value(100., "km"),wavelengths), 3);
		CCD_detector1->name(string("Graham_e")+int2str(i));
		Detectors1.push_back(CCD_detector1);
		
	}
	
// RUN SRTC++!  for zero albedo case.
	SRTC S1(a1, &hv, Detectors1);
	cout << "Running GrahamTest unity-surface-albedo case:\n";
	S1.run();
	
// write out the results (for now)
	cout << "Entering writing loop\n"; cout.flush();
	for (int i(0);i<Detectors1.size();i++) {
		Detectors1.at(i)->write(Detectors1.at(i)->name());
	}
	
// Generate Graham results
	vector<cube> Graham__0(Detectors1.size());  // surface only
	vector<cube> Graham__1(Detectors1.size());  // atmosphere first-scatter
	vector<cube> Graham_01(Detectors1.size());  // surface, then atmo
	vector<cube> Graham_10(Detectors1.size());  // atmo, then surface
	
	for (int i(0);i<Detectors1.size();i++) {
		Graham__0.at(i) = Detectors1.at(0)->CCD()*0.;
		Graham__1.at(i) = Detectors1.at(0)->CCD()*0.;
		Graham_01.at(i) = Detectors1.at(0)->CCD()*0.;
		Graham_10.at(i) = Detectors1.at(0)->CCD()*0.;
	}
	

	for (int d(0);d<Detectors1.size();d++) {
		string prefix("Graham_analytical_phase");
		prefix += int2str(int(emissions.at(d).degrees()));
		if (FILE *file = fopen((prefix+string("__0.Jcube")).c_str(), "r")) {
			fclose(file);
			Graham__0.at(d) = cube((prefix+string("__0.Jcube")).c_str());
			Graham__1.at(d) = cube((prefix+string("__1.Jcube")).c_str());
			Graham_10.at(d) = cube((prefix+string("_10.Jcube")).c_str());
			Graham_01.at(d) = cube((prefix+string("_01.Jcube")).c_str());
		} else {	
			cout << "Calculating Vixie et al. (2015) based analytical results for phase " << emissions.at(d).as(angle::DEGREES) << " --  00%";
			for (int x(0);x<Detectors1.at(0)->CCD().N(X);x++) {
				printpercent(x,Detectors1.at(0)->CCD().N(X)-1);
//#pragma omp parallel for schedule(dynamic)
				for (int y=0;y<Detectors1.at(0)->CCD().N(Y);y++) {
				
					if (Detectors1.at(d)->geo()(x,y,0) >=-90. &&
							 	Detectors1.at(d)->geo()(x,y,0) <= 90.  &&  Detectors1.at(d)->geo()(x,y,5) < 90.) {// i.e. not NaN's
						angle incidence(Detectors1.at(d)->geo()(x,y,5), angle::DEG);
						angle emission(Detectors1.at(d)->geo()(x,y,6), angle::DEG);
						angle phase(Detectors1.at(d)->geo()(x,y,4), angle::DEG);  
						angle atmo_phase(angle(180.,angle::DEG)-phase);
						
						double incident_atmospheres(1./incidence.cos());
						double emitted_atmospheres(1./emission.cos());
						
						
						
						for (int w(0);w<Detectors1.at(d)->CCD().N(Z);w++) {
							bool calculate(true);
							for (int yy(0);yy<y;yy++) {// check to see if this is a previously solved geometry
								if (Detectors1.at(d)->geo()(x,y,5)  ==  Detectors1.at(d)->geo()(x,yy,5)    &&
									 Detectors1.at(d)->geo()(x,y,6)  ==  Detectors1.at(d)->geo()(x,yy,6)    &&
									 Detectors1.at(d)->geo()(x,y,4)  ==  Detectors1.at(d)->geo()(x,yy,4)) {
									calculate = false;
									Graham__0.at(d)(x,y,w) = Graham__0.at(d)(x,yy,w);
									Graham__1.at(d)(x,y,w) = Graham__1.at(d)(x,yy,w);
									Graham_10.at(d)(x,y,w) = Graham_10.at(d)(x,yy,w);
									Graham_01.at(d)(x,y,w) = Graham_01.at(d)(x,yy,w);
								}
							}
							if (calculate) {
							
								Graham__0.at(d)(x,y,w) = exp(-(a1.totaltau(wavelengths.at(w)))*(incident_atmospheres+emitted_atmospheres));
								Graham__0.at(d)(x,y,w)*= a1.front().back()->zonephasefunction()->normalized_phasefunction(Detectors1.at(d)->CCD().Axis(Z,w), phase, incidence, emission);
								Graham__0.at(d)(x,y,w)*= Jcrap::pi;
								
								if (emission.degrees()==0.) cout << "emission is zero, phasefunction is " <<
									a1.front().back()->zonephasefunction()->normalized_phasefunction(Detectors1.at(d)->CCD().Axis(Z,w), phase, incidence, emission) << "\n";
								
								Graham__1.at(d)(x,y,w) = ssalbedos(0,0,w)*Jcrap::pi;
								Graham__1.at(d)(x,y,w)*= a1.back().at(0)->zonephasefunction()->normalized_phasefunction(Detectors1.at(d)->CCD().Axis(Z,w),atmo_phase, incidence, emission);
								Graham__1.at(d)(x,y,w)*= emitted_atmospheres / (incident_atmospheres+emitted_atmospheres);
								Graham__1.at(d)(x,y,w)*= 1.-exp(-(a1.totaltau(wavelengths.at(w))*(incident_atmospheres+emitted_atmospheres)));
								
								Graham10 Zd_integrator_10;
								Zd_integrator_10.a = &a1;
								Zd_integrator_10.incidence = incidence;
								Zd_integrator_10.emission = emission;
								Zd_integrator_10.phase = phase;
								Zd_integrator_10.wavelength_um = double(wavelengths.at(w));
								 
								Graham_10.at(d)(x,y,w) = 1.*Jcrap::pi*ssalbedos(0,0,w); 
								Graham_10.at(d)(x,y,w)*= exp(-(a1.totaltau(wavelengths.at(w))*emitted_atmospheres));
//								if () {
//									cube graph10(Zd_integrator_10.plot(0., Jcrap::pi/2., 90));
//									graph10.keyword("title", string("10 integration").c_str());
//									graph10.graph();
//								}
								Graham_10.at(d)(x,y,w)*= Zd_integrator_10.integrate(0., Jcrap::pi/2., accuracy, min_number, "Graham10 ZD integrator");
								
								
								Graham01 Zd_integrator_01;
								Zd_integrator_01.a = &a1;
								Zd_integrator_01.incidence = incidence;
								Zd_integrator_01.emission = emission;
								Zd_integrator_01.phase = phase;
								Zd_integrator_01.makegraphs = false;
								Zd_integrator_01.wavelength_um = double(wavelengths.at(w));
							
								Graham_01.at(d)(x,y,w) = 1.*Jcrap::pi*ssalbedos(0,0,w); 
								Graham_01.at(d)(x,y,w)*= exp(-(a1.totaltau(wavelengths.at(w))*incident_atmospheres));						
								Graham_01.at(d)(x,y,w)*= Zd_integrator_01.integrate(0., Jcrap::pi/2., accuracy/10., min_number, "Graham01 ZD integrator");			
								if (Graham_01.at(d)(x,y,w) >= 1. /*(x==10 && y==28 && w==-1) || (x==11 && y==28 && w==-1)*/ ) {
									cout << "Plotting an unusual point, sir!\n";
								//	Zd_integrator_01.makegraphs=true;
									cube graph01(Zd_integrator_01.plot(0.001, Jcrap::pi/2.-0.001, 1000));
									cout << "point plotted.\n";
									string title(string("01 integration i=")+int2str(int(incidence.degrees()))+string("  e=")+int2str(int(emission.degrees())));
									title+= string("  int=");
									title+= double2fstr(Zd_integrator_01.integrate(0., Jcrap::pi/2., accuracy/10., min_number, "Graham01 ZD integrator"),6);
									title+= string("  x=");
									title+= int2str(x);
									graph01.keyword("title", title.c_str());
									graph01.graph("", oX);
			//						graph01.graph("", oGIF);
									
									cout << "point graphed.\n";
									if (w==Detectors1.at(d)->CCD().N(Z)-1) {
										graph01.graph("8   -C ");
//										char c; cin>>c;
									}
								}	
							}		
						}
					}
				}
			}		
			printpercent(100,100);
			cout << "\n";	
			Graham__0.at(d).write((prefix+string("__0.Jcube")).c_str());
			Graham__1.at(d).write((prefix+string("__1.Jcube")).c_str());
			Graham_01.at(d).write((prefix+string("_01.Jcube")).c_str());
			Graham_10.at(d).write((prefix+string("_10.Jcube")).c_str());
		}
	}


// Analyze results
	cube result_template(emissions.size()/2,emissions.size()/2,ssalbedos.N(Z), 0.);
	result_template.copyaxis(ssalbedos,Z);
	for (int i(0);i<result_template.N(X);i++) result_template.Axis(X,i) = emissions.at(i).degrees();
	for (int i(0);i<result_template.N(Y);i++) result_template.Axis(Y,i) = emissions.at(i).degrees();
	cube rt__0(result_template);
	cube rt__1(result_template);
	cube rt_01(result_template);
	cube rt_10(result_template);
	cube Gr__0(result_template);
	cube Gr__1(result_template);
	cube Gr_01(result_template);
	cube Gr_10(result_template);
	
	cube weights(result_template);
	cout << "Analyzing results:   00%";
	for (int x(0);x<Detectors1.at(0)->CCD().N(X);x++) {
		printpercent(x,Detectors1.at(0)->CCD().N(X)-1);
		for (int y(0);y<Detectors1.at(0)->CCD().N(Y);y++) {
			for (int d(0);d<Detectors1.size();d++) {
				if (Detectors1.at(d)->geo()(x,y,0) >=-90. &&
					 	Detectors1.at(d)->geo()(x,y,0) <= 90.) {// i.e. not NaN's
					double thisincidence(Detectors1.at(d)->geo()(x,y,5));
					double thisemission(Detectors1.at(d)->geo()(x,y,6));
					for (int i(0);i<rt__0.N(X);i++) {  // incidence
						double incidence(rt__0.Axis(X,i));
						if (fabs(thisincidence-incidence)<acceptancecone) {
							for (int e(0);e<rt__0.N(Y);e++) {	// emission
								double emission(rt__0.Axis(Y,e));
								if (fabs(thisemission-emission)<acceptancecone) {
									for (int w(0);w<rt__0.N(Z);w++) {	// wavelength
										if (isnormal(Detectors1.at(d)->CCD()(x,y,w)) || Detectors1.at(d)->CCD()(x,y,w)==0.) {
											weights(i,e,w) += 1.;
											rt__0(i,e,w) += ((elephant*)(Detectors1.at(d)))->CCD(1,bitset<64>(string("0")))(x,y,w);
											rt__1(i,e,w) += ((elephant*)(Detectors1.at(d)))->CCD(1,bitset<64>(string("1")))(x,y,w);
											rt_01(i,e,w) += ((elephant*)(Detectors1.at(d)))->CCD(2,bitset<64>(string("01")))(x,y,w);
											rt_10(i,e,w) += ((elephant*)(Detectors1.at(d)))->CCD(2,bitset<64>(string("10")))(x,y,w);
											Gr__0(i,e,w) += Graham__0.at(d)(x,y,w);
											Gr__1(i,e,w) += Graham__1.at(d)(x,y,w);
											Gr_01(i,e,w) += Graham_01.at(d)(x,y,w);
											Gr_10(i,e,w) += Graham_10.at(d)(x,y,w);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	printpercent(12,12);
	cout << "\n";
	rt__0/=weights;
	rt__1/=weights;
	rt_01/=weights;
	rt_10/=weights;
	Gr__0/=weights;
	Gr__1/=weights;
	Gr_01/=weights;
	Gr_10/=weights;
	weights.write("weights.Jcube");
	
	cube rterror__0(result_template);
	cube rterror__1(result_template);
	cube rterror_01(result_template);
	cube rterror_10(result_template);
	cout << "Analyzing errors:   00%";
	for (int x(0);x<Detectors1.at(0)->CCD().N(X);x++) {
		printpercent(x,Detectors1.at(0)->CCD().N(X)-1);
		for (int y(0);y<Detectors1.at(0)->CCD().N(Y);y++) {
			for (int d(0);d<Detectors1.size();d++) {
				if (Detectors1.at(d)->geo()(x,y,0) >=-90. &&
					 	Detectors1.at(d)->geo()(x,y,0) <= 90.) {// i.e. not NaN's
					double thisincidence(Detectors1.at(d)->geo()(x,y,5));
					double thisemission(Detectors1.at(d)->geo()(x,y,6));
					for (int i(0);i<rt__0.N(X);i++) {  // incidence
						double incidence(rt__0.Axis(X,i));
						if (fabs(thisincidence-incidence)<acceptancecone) {
							for (int e(0);e<rt__0.N(Y);e++) {	// emission
								double emission(rt__0.Axis(Y,e));
								if (fabs(thisemission-emission)<acceptancecone) {
									for (int w(0);w<rt__0.N(Z);w++) {	// wavelength
										if (isnormal(Detectors1.at(d)->CCD()(x,y,w)) || Detectors1.at(d)->CCD()(x,y,w)==0.) {
											rterror__0(i,e,w) += pow(((elephant*)(Detectors1.at(d)))->CCD(1,bitset<64>(string("0")))(x,y,w)-rt__0(i,e,w), 2.);
											rterror__1(i,e,w) += pow(((elephant*)(Detectors1.at(d)))->CCD(1,bitset<64>(string("1")))(x,y,w)-rt__1(i,e,w), 2.);
											rterror_01(i,e,w) += pow(((elephant*)(Detectors1.at(d)))->CCD(2,bitset<64>(string("01")))(x,y,w)-rt_01(i,e,w), 2.);
											rterror_10(i,e,w) += pow(((elephant*)(Detectors1.at(d)))->CCD(2,bitset<64>(string("10")))(x,y,w)-rt_10(i,e,w), 2.);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	printpercent(12,12);
	cout << "\n";
	weights-=1.;
	rterror__0 /= weights;
	rterror__1 /= weights;
	rterror_10 /= weights;
	rterror_01 /= weights;
	rterror__0 = rterror__0.pow(0.5);
	rterror__1 = rterror__1.pow(0.5);
	rterror_10 = rterror_10.pow(0.5);
	rterror_01 = rterror_01.pow(0.5);
	weights+=1.;
	weights=weights.pow(0.5);
	rterror__0 /= weights;  // to get to std dev of the mean
	rterror__1 /= weights;
	rterror_10 /= weights;
	rterror_01 /= weights;
	
	
	
	
	
	rt__0.write("rt__0.Jcube");
	rt__1.write("rt__1.Jcube");
	rt_01.write("rt_01.Jcube");
	rt_10.write("rt_10.Jcube");
	Gr__0.write("Gr__0.Jcube");
	Gr__1.write("Gr__1.Jcube");
	Gr_01.write("Gr_01.Jcube");
	Gr_10.write("Gr_10.Jcube");
	rterror__0.write("rterror__0.Jcube");
	rterror__1.write("rterror__1.Jcube");
	rterror_10.write("rterror_10.Jcube");
	rterror_01.write("rterror_01.Jcube");	
	
	string command;
// produce output plots
	
	command = string("graph -T fig -X \"Wavelength\" -Y \"I/F\" -L \"__0 Graham Result ");
	command+= a1.back().at(0)->zonephasefunction()->name();
	command+= string("\" -C ");
	command+= string(" -x 0.8 5.2 -y 0 0.6 ");
	command+= string("--pen-colors \"1=red:2=purple:3=blue:4=green:5=black\" ");
	for (int i(0);i<Gr__0.N(X)/2;i++) {
		Gr__0.skewer(Z,i,i).graph("", oGIF);
		rt__0.keyword("WITHERROR", "yes");
		rt__0.skewer(Z,i,i).blocks(X,rterror__0.skewer(Z,i,i)).graph("", oGIF);
		command+= string(" -S 0 -m ");
		command+= int2str(i%5+1);
		command+= string(" -I a graph");
		command+= int2str(graphno-2);
		command+= string(".out -S 3 -m -");
		command+= int2str(i%5+1);
		command+= string(" -I e graph");
		command+= int2str(graphno-1);
		command+= string(".out ");
	}
	command += string(" > __0.fig");
	{
		ofstream gc("graph.commands", ios::app);
		gc << command << "\n";
	}
	system(command.c_str());
	
	
	command = string("graph -T fig -X \"Wavelength\" -Y \"I/F\" -L \"__1 Graham Result ");
	command+= a1.back().at(0)->zonephasefunction()->name();
	command+= string("\" -C ");
	command+= string(" -x 0.8 5.2 -y 0 0.12 ");
	command+= string("--pen-colors \"1=red:2=purple:3=blue:4=green:5=black\" ");
	for (int i(0);i<Gr__0.N(X)/2;i++) {
		Gr__1.skewer(Z,i,i).graph("", oGIF);
		rt__1.keyword("WITHERROR", "yes");
		rt__1.skewer(Z,i,i).blocks(X,rterror__1.skewer(Z,i,i)).graph("", oGIF);
		command+= string(" -S 0 -m ");
		command+= int2str(i%5+1);
		command+= string(" -I a graph");
		command+= int2str(graphno-2);
		command+= string(".out -S 3 -m -");
		command+= int2str(i%5+1);
		command+= string(" -I e graph");
		command+= int2str(graphno-1);
		command+= string(".out ");
	}
	command+= string(" > __1.fig");
	{
		ofstream gc("graph.commands", ios::app);
		gc << command << "\n";
	}
	system(command.c_str());
	
	command = string("graph -T fig -X \"Wavelength\" -Y \"I/F\" -L \"_10 Graham Result ");
	command+= a1.back().at(0)->zonephasefunction()->name();
	command+= string("\" -C ");
	command+= string(" -x 0.8 5.2 -y 0 0.1 ");
	command+= string("--pen-colors \"1=red:2=purple:3=blue:4=green:5=black\" ");
	for (int i(0);i<Gr__0.N(X)/2;i++) {
		Gr_10.skewer(Z,i,i).graph("", oGIF);
		rt_10.keyword("WITHERROR", "yes");
		rt_10.skewer(Z,i,i).blocks(X,rterror_10.skewer(Z,i,i)).graph("", oGIF);
		command+= string(" -S 0 -m ");
		command+= int2str(i%5+1);
		command+= string(" -I a graph");
		command+= int2str(graphno-2);
		command+= string(".out -S 3 -m -");
		command+= int2str(i%5+1);
		command+= string(" -I e graph");
		command+= int2str(graphno-1);
		command+= string(".out ");
	}
	command+= string(" > _10.fig");	
	{
		ofstream gc("graph.commands", ios::app);
		gc << command << "\n";
	} 
	system(command.c_str());	
		
	command = string("graph -T fig -X \"Wavelength\" -Y \"I/F\" -L \"_01 Graham Result ");
	command+= a1.back().at(0)->zonephasefunction()->name();
	command+= string("\" -C ");
	command+= string(" -x 0.8 5.2 -y 0 0.1 ");
	command+= string("--pen-colors \"1=red:2=purple:3=blue:4=green:5=black\" ");
	for (int i(0);i<Gr__0.N(X)/2;i++) {
		Gr_01.skewer(Z,i,i).graph("", oGIF);
		rt_01.keyword("WITHERROR", "yes");
		rt_01.skewer(Z,i,i).blocks(X,rterror_01.skewer(Z,i,i)).graph("", oGIF);
		command+= string(" -S 0 -m ");
		command+= int2str(i%5+1);
		command+= string(" -I a graph");
		command+= int2str(graphno-2);
		command+= string(".out -S 3 -m -");
		command+= int2str(i%5+1);
		command+= string(" -I e graph");
		command+= int2str(graphno-1);
		command+= string(".out ");
	}
	command+= string(" > _01.fig");
	{
		ofstream gc("graph.commands", ios::app);
		gc << command << "\n";
	}
	system(command.c_str());	

	cube diff__0((Gr__0-rt__0)/rterror__0);
	diff__0.write("diff__0.Jcube");
	cube diff__1((Gr__1-rt__1)/rterror__1);
	diff__1.write("diff__1.Jcube");
	cube diff_10((Gr_10-rt_10)/rterror_10);
	diff_10.write("diff_10.Jcube");
	cube diff_01((Gr_01-rt_01)/rterror_01);
	diff_01.write("diff_01.Jcube");
	
	double diff__0mean(0.), diff__1mean(0.), diff_01mean(0.), diff_10mean(0.);
	unsigned int d__0n(0), d__1n(0), d_01n(0), d_10n(0);
	for (int x(0);x<diff__0.N(X);x++)
		for (int y(0);y<diff__0.N(Y);y++)
			for (int z(0);z<diff__0.N(Z);z++) {
				if (isnormal(diff__0(x,y,z))) {
					d__0n++;
					if (diff__0(x,y,z) < 1.e3)  diff__0mean+=diff__0(x,y,z);
				}
				if (isnormal(diff__1(x,y,z))) {
					d__1n++;
					if (diff__1(x,y,z) < 1.e3) diff__1mean+=diff__1(x,y,z);
				}
				if (isnormal(diff_01(x,y,z))) {
					d_01n++;
					if (diff_01(x,y,z) < 1.e3) diff_01mean+=diff_01(x,y,z);
				}
				if (isnormal(diff_10(x,y,z))) {
					d_10n++;
					if (diff_10(x,y,z) < 1.e3) diff_10mean+=diff_10(x,y,z);
				}
			}
	cout << "diff__0=" << diff__0 << "\n diff__0mean is " << diff__0mean << " with d__0n=" << d__0n << "\n";
	cout << "diff__1=" << diff__1 << "\n diff__1mean is " << diff__1mean << " with d__1n=" << d__1n << "\n";
	cout << "diff_01=" << diff_01 << "\n diff_01mean is " << diff_01mean << " with d_01n=" << d_01n << "\n";
	cout << "diff_10=" << diff_10 << "\n diff_10mean is " << diff_10mean << " with d_10n=" << d_10n << "\n";
			
	return (diff__0mean+diff__1mean+diff_01mean+diff_10mean)/(d__0n+d__1n+d_01n+d_10n);
}

