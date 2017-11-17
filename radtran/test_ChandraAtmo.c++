#include "SRTC++_testsuite.h" 
#include "chandra_test.h"

extern int graphno;

double test_ChandraAtmo()
{
// phase angles to test
	vector<angle> phases;
	phases.push_back(angle(0., angle::DEGREES));
	phases.push_back(angle(45., angle::DEGREES));
	phases.push_back(angle(90., angle::DEGREES));
	phases.push_back(angle(135., angle::DEGREES));
	
// Calculate chandra analytical result
	vector<cube> analyticalresults;
	chandraslab c;
	for (int i(0);i<phases.size();i++) {
		disktest d(&c,phases.at(i));
		cube analyticalresult(d.plot(0.01, 0.99, 101));
		analyticalresult.graph("", outputtype);
		analyticalresults.push_back(analyticalresult);
	}

// setup single-scattering albedos such that a useful graph results
	cube ssalbedos(1,1,13);
	for (int i(0);i<10;i++) {
		ssalbedos(0,0,i) = double(i)/10.;
		ssalbedos.Axis(Z,i) = ssalbedos(0,0,i);
	}
	ssalbedos(0,0,10) = 0.95;
	ssalbedos.Axis(Z,10) = 0.95;
	ssalbedos(0,0,11) = 0.98;
	ssalbedos.Axis(Z,11) = 0.98;
	ssalbedos(0,0,12) = 1.00;
	ssalbedos.Axis(Z,12) = 1.00;
	
// abusing the wavelengths to do those different ssalbedos -- need a tau cube to correspond.
	cube taus(ssalbedos*0.+(3.+sqrt(qualityfactor)));
	
// create the atmosphere
	atmosphere a(atmosphere::orangerind(taus, value(30., "km"), ssalbedos));

// have to put in a wavelength, tho doesn't really matter in this case		
	vector<value> wavelengths;	
	for (int i(0);i<ssalbedos.N(Z);i++) {
		wavelengths.push_back(value(ssalbedos(0,0,i), "um")); // yes now have a wavelength of 0. hopefully doesn't adversely affect things.
	}	
	 
// setting up the surface
	cube albedoone(1,2,1);
	albedoone(0,0,0)=1.0; 
	albedoone.Axis(Y,0)=-1.;
	albedoone(0,1,0)=1.0; 
	albedoone.Axis(Y,0)=1.;  // surface albedo is all one.
	albedoone.Axis(Z,0)=1.5;  
	atmozone *uniformsurface=new surface_atmozone(albedoone);
 	a.front().addzone(uniformsurface);
	
// setting up the Sun	
	photongenerator_square hv(value(4200., "km"), 42*qualityfactor, value(4200., "km"), wavelengths);	
	geomvector incomingdirection(-1.,0.,0.);  // set the photons to come in initially toward -x
	hv.photondirection(incomingdirection);
	hv.fieldcenter(geomvector(4200.,0.,0.));
	
// Detector		
	vector<detector*> Detectors;
	for (int i(0);i<phases.size();i++) {
		angle phase(phases.at(i));
		double x(1.e9*phase.cos()); 
		double y(1.e9*phase.sin());
		colorCCD *CCD_detector;
		CCD_detector = new colorCCD(geomvector(x,y,0.),value(3200., "km"),value(100., "km"),wavelengths);
		Detectors.push_back(CCD_detector);
	}
	
// RUN SRTC++!
	SRTC S(a, &hv, Detectors);
	S.run();

	(++a.begin())->at(0)->zonephasefunction()->test_distribution_function().write("isotropic_dist_function.Jcube");
	(++a.begin())->at(0)->zonephasefunction()->test_distribution_normalized().write("isotropic_dist_normalized.Jcube");
	
// Analyze results
	vector<cube> calculatedresults;
	for (int i(0);i<Detectors.size();i++) {
		string outname("ChandraAtmo_");
		outname += int2str(int(phases.at(i).degrees()));
		outname += ".Jcube";
		Detectors.at(i)->CCD().write(outname);
		
		cube calculatedresult(diskintegrate(Detectors.at(i)->CCD()));
		calculatedresult.graph("", outputtype);
		calculatedresults.push_back(calculatedresult);
	}	

// Generate nice plot
	string command(string("graph -T "));
	command+=asstring(outputtype);
	command+=string(" -X \"Single-Scattering Albedo\" -Y \"Disk-Integrated I/F\" -L \"SRTC++ Comparison to Chandra (");
	for (int i(0);i<phases.size();i++) {
		command+=int2str(int(phases.at(i).degrees()));
		if (i<phases.size()-1) command+=string(",");
	}
	command+=string(" Phase)\" ");
	command+=string("--pen-colors \"1=red:2=purple:3=blue:4=green:5=brown\" "); 
	for (int i(0);i<phases.size();i++) {
		command+=string("-m 4 graph");	
		command+=int2str(graphno-2*int(phases.size())+i);
		command+=string(".out ");
	}	
	command+=string(" -S 3 -C -m 0 ");
	for (int i(0);i<phases.size();i++) {
		command+=string(" graph");	
		command+=int2str(graphno-int(phases.size())+i);
		command+=string(".out ");
	}
	command+=string(" > Chandracompare.");
	command+=asstring(outputtype);
	system(command.c_str());
	cout << command << "\n";
	
	{
		ofstream graphcommands("graph.commands", ios_base::app);
		graphcommands << command << "\n";
	}


// compare analytical and numerical to evaluate goodness of fit
	cube diff(calculatedresults.at(0));
	for (int p(0);p<diff.N(X);p++){
		if (p!=0) diff=diff.blocks(X,calculatedresults.at(p));
		for (int i(0);i<diff.N(Z);i++){
			cout << "dividing " << diff(p,0,i) << " by " << analyticalresults.at(p)(0,0,diff.Axis(Z,i)) << "\n";
			if (analyticalresults.at(p)(0,0,diff.Axis(Z,i) != 0.))
				diff(p,0,i) /= analyticalresults.at(p)(0,0,diff.Axis(Z,i));
			else
				diff(p,0,i) = 1.; 		
		}
	}
	diff -= 1.;
	return diff.mean();
}

