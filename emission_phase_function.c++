#include "SRTC++.h"
 
const double acceptancecone(2.); // degrees from specified geom that we accept
double km_per_pixel(25.);
const double qualityfactor(10.);

const otype outputtype(oFIG);

main(int argn, char **argv)
{
	vector<angle> incidences;  // name lies:  there should be only one
	incidences.push_back(angle(60., angle::DEGREES));
	
// phase angles to test
	vector<angle> phases;
	phases.push_back(incidences.at(0)+angle(-90., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-80., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-70., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-60., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-50., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-40., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-30., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-20., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle(-10., angle::DEGREES));
	phases.push_back(incidences.at(0));
	phases.push_back(incidences.at(0)+angle( 10., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 20., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 30., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 40., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 50., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 60., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 70., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 80., angle::DEGREES));
	phases.push_back(incidences.at(0)+angle( 90., angle::DEGREES));
	
// emission angles to test in azimuth=90 geometry
	vector<angle> emissions;
	emissions.push_back(angle( 0., angle::DEGREES));
	emissions.push_back(angle(10., angle::DEGREES));
	emissions.push_back(angle(20., angle::DEGREES));
	emissions.push_back(angle(30., angle::DEGREES));
	emissions.push_back(angle(40., angle::DEGREES));
	emissions.push_back(angle(50., angle::DEGREES));
	emissions.push_back(angle(60., angle::DEGREES));
	emissions.push_back(angle(70., angle::DEGREES));
	emissions.push_back(angle(80., angle::DEGREES));
	emissions.push_back(angle(90., angle::DEGREES));
	

// Comparing to the Sebastien model just to test
	atmosphere a(atmosphere::SebastienValidateDISR08_layered());

	
	
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

	 
// setting up the surface
	cube surfalbedo(1,1,1);
	surfalbedo(0,0,0)=0.2;
	atmozone *uniformsurface=new surface_atmozone(surfalbedo);
 	a.front().addzone(uniformsurface);
	
	
	cout << "Atmosphere for emission phase function test: " << a << "\n";

// setting up the Sun	
	photongenerator_multitarget hv(value(4200./42., "km"), 42.*qualityfactor, value(4200., "km"), wavelengths);
	
	for (int i(0);i<incidences.size();i++) {
		hv.addfieldcenter(geomvector(4200.,2575.*incidences.at(i).sin(),0.));
	}
	cout << "Sun set up -- now setting up detectors\n";
	
// Detector		
	vector<detector*> Detectors;
	for (int i(0);i<phases.size();i++) {
		angle phase(phases.at(i));
		double x(1.e9*phase.cos()); 
		double y(1.e9*phase.sin());
		
		colorCCD *CCD_detector_az0;
		CCD_detector_az0 = new colorCCD(geomvector(x,y,0.),value(3200., "km"), value(km_per_pixel, "km"),wavelengths);
		CCD_detector_az0->name(string("epf_az0phaseangle")+int2str(int(phases.at(i).degrees())));
		Detectors.push_back(CCD_detector_az0);
	}
	cout << "First set of detectors set up.\n";
	
	for (int i(0);i<emissions.size();i++) {
		angle emission(emissions.at(i));
		double x(1.e9*incidences.at(0).cos()*emission.cos()); 
		double y(1.e9*incidences.at(0).sin()*emission.cos());
		double z(1.e9*emission.sin());
		
		colorCCD *CCD_detector_az90;
		CCD_detector_az90 = new colorCCD(geomvector(x,y,z),value(3200., "km"), value(km_per_pixel, "km"),wavelengths);
		CCD_detector_az90->name(string("epf_az90emissionangle")+int2str(int(emissions.at(i).degrees())));
		Detectors.push_back(CCD_detector_az90);
	}
	
	cout << "Detectors set up.  Now constructing SRTC\n";
	
// RUN SRTC++!
	SRTC S0(a, &hv, Detectors);
	cout << "Running emission_phase_function :\n";
	S0.run();
	
// write out the results (for now)
	cout << "Entering writing loop\n"; cout.flush();
	for (int i(0);i<Detectors.size();i++) {
		Detectors.at(i)->write(Detectors.at(i)->name());
	}


// Analyze results
	cube epf(10,3,wavelengths.size());
	for (int i(0);i<10;i++)
		epf.Axis(X,i) = double(i)*10.;
	for (int i(0);i<epf.N(Z);i++)
		epf.Axis(Z,i) = wavelengths.at(i);
	epf.Axis(Y,0) = 0.;
	epf.Axis(Y,1) = 90.;
	epf.Axis(Y,2) = 180.;
	cube weights(epf.plane(Z,0)*0.);	
	
	vector<cube> whichused;	
	for (int i(0);i<Detectors.size();i++) {
		whichused.push_back(Detectors.at(i)->CCD().plane(Z,0)*0.-1.);
		for (int j(0);j<incidences.size();j++) {
			if (j!=0) whichused.at(i) = whichused.at(i).blocks(Z,Detectors.at(i)->CCD().plane(Z,0)*0.-1.);
			whichused.at(i).Axis(Z,j)=incidences.at(j).degrees();
		}
	}
	cout << "Analyzing results:   00%";
	for (int x(0);x<Detectors.at(0)->CCD().N(X);x++) {
		printpercent(x,Detectors.at(0)->CCD().N(X)-1);
		for (int y(0);y<Detectors.at(0)->CCD().N(Y);y++) {
			for (int d(0);d<Detectors.size();d++) {
				if (Detectors.at(d)->geo()(x,y,0) >=-90. &&
					 	Detectors.at(d)->geo()(x,y,0) <= 90.) {// i.e. not NaN's
					double thisphase(Detectors.at(d)->geo()(x,y,4));
					double thisincidence(Detectors.at(d)->geo()(x,y,5));
					double thisemission(Detectors.at(d)->geo()(x,y,6));
					double incidence(incidences.at(0).degrees());
					double thislatitude(Detectors.at(d)->geo()(x,y,0));
					
					int xpt(0);
					if (d<=9) xpt=0;
					else if (d>9 && d<=18) xpt=1;
					else if (d>=19) xpt=2;
					
					if (xpt<=1 && fabs(thislatitude)>acceptancecone) continue;
					if (xpt==0 && Detectors.at(0)->CCD().Axis(X,x)<0. && thisphase<3.) continue;
					
					int e(0);
					if (xpt==0) e=9-d;
					else if (xpt==1) e=d-9;
					else if (xpt==2) e=d-19;
					double emission(epf.Axis(X,e));
					
					double desiredphase(0.);
					if (xpt<=1) desiredphase = fabs(phases.at(d).degrees());
					else {
						angle inc_ang(incidence, angle::DEGREES);
						angle emi_ang(emission, angle::DEGREES);
						desiredphase = angle_acos(inc_ang.cos()*emi_ang.cos()).degrees();
					}
//					cout << "xpt=" << xpt << ", e=" << e << ", d="<<d<<"\n";
//					cout << "desiredphase = " << desiredphase << ", would have been " << phases.at(e).degrees() << "\n";
//					cout << "for incidence = " << incidence << ", emission = " << emission << "\n";
									
					if (abs(thisincidence-incidence)<acceptancecone) {
						if (abs(thisemission-emission)<acceptancecone) {
							if (abs(thisphase-desiredphase)<acceptancecone) {
								if (xpt<=1 || (xpt==2 && Detectors.at(d)->CCD().Axis(Y,y)>0.)) {
//									cout << "x="<<x<<", y="<<y<<"\n";
//									cout << "e="<<e<<", d="<<d<<", xpt="<<xpt<<"\n";
//									cout << "Emission=" << emission << ", phase = " << desiredphase << ", incidence=" << incidence<<"\n";
//									cout << "this_e=" << thisemission << ", this_p="<< thisphase<<", this_i=" << thisincidence<<"\n"; 
									weights(e,xpt,0) += 1.;
									whichused.at(d)(x,y,0)=incidences.at(0).degrees();
									for (int w(0);w<epf.N(Z);w++) {	// wavelength
										epf(e,xpt,w) += Detectors.at(d)->CCD()(x,y,epf.Axis(Z,w));
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
	weights.write("weights.Jcube");
	epf/=weights;
	for (int w(0);w<epf.N(Z);w++) epf(0,1,w) = epf(0,0,w);  // emission 0 is used in both 0 and 180 azimuth epf's
	epf.write("epf.Jcube");
	for (int i(0);i<whichused.size();i++) {
		if (i<=18)
			whichused.at(i).write(string("whichused_p")+int2str(int(phases.at(i).degrees()))+string(".Jcube"));
		else
			whichused.at(i).write(string("whichused_e")+int2str(int(emissions.at(i-19).degrees()))+string(".Jcube"));
	}
	
// estimate errors
	cube epferror(epf*0.);
	for (int x(0);x<Detectors.at(0)->CCD().N(X);x++) {
		for (int y(0);y<Detectors.at(0)->CCD().N(Y);y++) {
			for (int d(0);d<Detectors.size();d++) {				
				if (whichused.at(d)(x,y,0)>=0.) {
					int xpt(0);
					if (d<=9) xpt=0;
					else if (d>9 && d<=18) xpt=1;
					else if (d>=19) xpt=2;
				
					int e(0);
					if (xpt==0) e=9-d;
					else if (xpt==1) e=d-9;
					else if (xpt==2) e=d-19;
				
					for (int w(0);w<epf.N(Z);w++) {	// wavelength
						epferror(e,xpt,w) += pow(Detectors.at(d)->CCD()(x,y,epf.Axis(Z,w))-epf(e,xpt,w) ,2.);
					}						
				}
			}
		}
	}
	epferror /= weights;
	epferror  = epferror.pow(0.5);
	weights   = weights.pow(0.5);
	epferror /= weights;  // to get to std dev of the mean	
	for (int w(0);w<epferror.N(Z);w++) epferror(0,1,w) = epferror(0,0,w);  // emission 0 is used in both 0 and 180 azimuth epf's

	epferror.write("epferror.Jcube");
	
	
	
// Generate nice plots


	epf = epf(0,8,-1,-1,-1,-1);
	epferror = epferror(0,8,-1,-1,-1,-1);
	
	for (int w(0);w<epf.N(Z);w++) {
		for (int i(0);i<epf.N(Y);i++) {
			cube thisresult(epf.skewer(X,i,w));
			thisresult = thisresult.blocks(Z,epferror.skewer(X,i,w));
			thisresult.keyword("WITHERROR", "yes");
			thisresult = thisresult.yt2xy();
			thisresult.graph("", oGIF);
		}
				
		string command(string("graph -T "));
		command+=asstring(outputtype);
		command+=string(" -X \"Emission Angle (degrees)\" -Y \"I/F\" -L \"EPF i=60 ");
		command+=string(double2fstr(epf.Axis(Z,w),4));
		command+=string("\" --pen-colors \"1=red:2=darkgreen:3=blue:4=darkgreen:5=black\" -C ");	
		command+=string(" -x -3 93  -y 0 0.26 "); 				

		for (int e(0);e<epf.N(Y);e++) {
			command+=string(" -S 0 -m ");
			command+=int2str(e+1);
			command+=string(" -I e graph");	
			command+=int2str(w*epf.N(Y)+e);
			command+=string(".out ");
		}
			
		command+=string(" > emission_phase_function_i=60_");
		command+=string(double2fstr(epf.Axis(Z,w),4));
		command+=string(".");
		command+=asstring(outputtype);
		system(command.c_str());
		cout << command << "\n";
		
		{
			ofstream graphcommands("graph.commands", ios_base::app);
			graphcommands << command << "\n";
		}
	}
/*		
	command = string("graph -T ");
	command+=asstring(outputtype);
	command+=string(" -X \"Wavelength (um)\" -Y \"I/F\" -L \"SebastienCompare A=1\" ");
	command+=string(" -x 0.8 5.2  -y 0.0 1.1 ");
	
	command+=string("--pen-colors \"1=red:2=purple:3=blue:4=darkgreen:5=black\" -C "); 
	for (int e(0);e<incidences.size();e++) {
		command+=string(" -S 0 -m ");
		command+=int2str((e%5)+1);
		command+=string(" -I a graph");	
		command+=int2str(int(graphno-4*incidences.size()+3+4*e));
		command+=string(".out ");
		
		command+=string(" -S 3 -m -");
		command+=int2str((e%5)+1);
		command+=string(" -I e graph");	
		command+=int2str(int(graphno-4*incidences.size()+2+4*e));
		command+=string(".out ");
	}
	
	command+=string(" > Sebastiencompare_a1");
	command+=string(".");
	command+=asstring(outputtype);
	system(command.c_str());
	cout << command << "\n";
	
	{
		ofstream graphcommands("graph.commands", ios_base::app);
		graphcommands << command << "\n";
	}



// compare analytical and numerical to evaluate goodness of fit

	cube diff0((result0-pp0)/rterror0);
	diff0.write("diff0.Jcube");
	cube diff1((result1-pp1)/rterror1);
	diff1.write("diff1.Jcube");

*/
}

