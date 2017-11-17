#include "SRTC++_testsuite.h"
 
extern int graphno;
const double acceptancecone(2.); // degrees from specified geom that we accept
double km_per_pixel(50.);

double test_SebastienCompare()
{
	
//	ingest Sebastien's results
	cube Seb(SRTC::datapath()+string("/data/inputfiles/Test1_Shannon.txt"));
	cube pp0(6,6,16);  // planeparallel results with 0 albedo
	for (int i(0);i<pp0.N(X);i++) pp0.Axis(X,i)=double(i)*10.;
	for (int e(0);e<pp0.N(X);e++) pp0.Axis(Y,e)=double(e)*10.;
	pp0.Axis(Z,0) = 0.933078;  // these are all just what Sebastien sent us
	pp0.Axis(Z,1) = 1.08183;
	pp0.Axis(Z,2) = 1.27813;
	pp0.Axis(Z,3) = 1.59018;
	pp0.Axis(Z,4) = 2.01781;
	pp0.Axis(Z,5) = 2.69620;
	pp0.Axis(Z,6) = 2.79889;
	pp0.Axis(Z,7) = 4.98896;
	pp0.Axis(Z,8) = 5.00576;
	pp0.Axis(Z,9) = 5.02240;
	pp0.Axis(Z,10)= 5.04078;
	pp0.Axis(Z,11)= 5.05734;
	pp0.Axis(Z,12)= 5.07402;
	pp0.Axis(Z,13)= 5.09106;
	pp0.Axis(Z,14)= 5.10680;
	pp0.Axis(Z,15)= 5.12250;	
	cube pp1(pp0);  // and with albedo 1
	for (int i(0);i<pp0.N(X);i++) {  // incidence 
		for (int e(0);e<pp0.N(Y);e++) {	// emission
			for (int w(0);w<pp0.N(Z);w++) {	// wavelength
				pp0(i,e,w) = Seb(w+3,6*i+e,0);
				pp1(i,e,w) = Seb(w+3,6*i+e,1);
			}
		}
	}
	pp0.write("planeparallel0.Jcube");
	pp1.write("planeparallel1.Jcube");
	
	
// phase angles to test
	vector<angle> phases;
	phases.push_back(angle( 0., angle::DEGREES));
//	phases.push_back(angle(10., angle::DEGREES));
	phases.push_back(angle(20., angle::DEGREES));
//	phases.push_back(angle(30., angle::DEGREES));
	phases.push_back(angle(40., angle::DEGREES));
//	phases.push_back(angle(50., angle::DEGREES));
	phases.push_back(angle(60., angle::DEGREES));
//	phases.push_back(angle(70., angle::DEGREES));
	phases.push_back(angle(80., angle::DEGREES));
//	phases.push_back(angle(90., angle::DEGREES));
	phases.push_back(angle(100., angle::DEGREES));
	
	vector<angle> incidences;
	incidences.push_back(angle( 0., angle::DEGREES));
	incidences.push_back(angle(10., angle::DEGREES));
	incidences.push_back(angle(20., angle::DEGREES));
	incidences.push_back(angle(30., angle::DEGREES));
	incidences.push_back(angle(40., angle::DEGREES));
	incidences.push_back(angle(50., angle::DEGREES));

// Comparing to the Sebastien model just to test
	atmosphere a0(atmosphere::SebastienValidateDISR08_layered());
	atmosphere a1(atmosphere::SebastienValidateDISR08_layered());

	
	
// 
	vector<value> wavelengths;
	wavelengths.push_back(value(0.933078,"um"));
	wavelengths.push_back(value(1.08183,"um"));
	wavelengths.push_back(value(1.27813,"um"));
	wavelengths.push_back(value(1.59018, "um"));
	wavelengths.push_back(value(2.01781,"um"));
	wavelengths.push_back(value(2.69620,"um"));
	wavelengths.push_back(value(2.79889,"um"));
//	wavelengths.push_back(value(4.98896, "um"));
	wavelengths.push_back(value(5.00576, "um"));
//	wavelengths.push_back(value(5.02240, "um"));
//	wavelengths.push_back(value(5.04078, "um"));
//	wavelengths.push_back(value(5.05734, "um"));
//	wavelengths.push_back(value(5.07402, "um"));
//	wavelengths.push_back(value(5.09106, "um"));
//	wavelengths.push_back(value(5.10680, "um"));
//	wavelengths.push_back(value(5.12250, "um"));

	 
// setting up the surface
	cube zero(1,1,1);
	zero(0,0,0)=0.0;
	cube one(1,1,1); 
	one(0,0,0) =1.0;
	atmozone *uniformsurface0=new surface_atmozone(zero);
 	a0.front().addzone(uniformsurface0);
	atmozone *uniformsurface1=new surface_atmozone(one);
 	a1.front().addzone(uniformsurface1);
	
	
	cout << "Atmosphere for SebastienCompare0: " << a0 << "\n";
	cout << "Atmosphere for SebastienCompare1: " << a1 << "\n";

// setting up the Sun	
	photongenerator_multitarget hv(value(4200./42., "km"), 42.*qualityfactor, value(4200., "km"), wavelengths);
	
	for (int i(0);i<incidences.size();i++) {
		hv.addfieldcenter(geomvector(4200.,2575.*incidences.at(i).sin(),0.));
	}
	
// Detector		
	vector<detector*> Detectors0;
	vector<detector*> Detectors1;
	for (int i(0);i<phases.size();i++) {
		angle phase(phases.at(i));
		double x(1.e9*phase.cos()); 
		double y(1.e9*phase.sin());
		
		colorCCD *CCD_detector0;
		CCD_detector0 = new colorCCD(geomvector(x,y,0.),value(3200., "km"), value(km_per_pixel, "km"),wavelengths);
		CCD_detector0->name(string("zero_e")+int2str(int(incidences.at(i).degrees())));
		Detectors0.push_back(CCD_detector0);
		
		colorCCD *CCD_detector1;
		CCD_detector1 = new colorCCD(geomvector(x,y,0.),value(3200., "km"), value(km_per_pixel, "km"),wavelengths);
		CCD_detector1->name(string("one_e")+int2str(int(incidences.at(i).degrees())));
		Detectors1.push_back(CCD_detector1);
		
	}
	
// RUN SRTC++!  for zero albedo case.
	SRTC S0(a0, &hv, Detectors0);
	cout << "Running SebastienValidateDISR08_Layered zero-surface-albedo case:\n";
	S0.run();	
	SRTC S1(a1, &hv, Detectors1);
	cout << "Running SebastienValidateDISR08_Layered unity-surface-albedo case:\n";
	S1.run();
	
// write out the results (for now)
	cout << "Entering writing loop\n"; cout.flush();
	for (int i(0);i<Detectors0.size();i++) {
		Detectors0.at(i)->write(Detectors0.at(i)->name());
		Detectors1.at(i)->write(Detectors1.at(i)->name());
	}

// Analyze results
	cube result0(pp0*0.), result1(pp1*0.);
	cube weights(pp0.plane(Z,0)*0.);	
	vector<cube> whichused;	
	for (int i(0);i<Detectors0.size();i++) {
		whichused.push_back(Detectors0.at(i)->CCD().plane(Z,0)*0.-1.);
		for (int j(0);j<incidences.size();j++) {
			if (j!=0) whichused.at(i) = whichused.at(i).blocks(Z,Detectors0.at(i)->CCD().plane(Z,0)*0.-1.);
			whichused.at(i).Axis(Z,j)=incidences.at(j).degrees();
		}
	}
	cout << "Analyzing results:   00%";
	for (int x(0);x<Detectors0.at(0)->CCD().N(X);x++) {
		printpercent(x,Detectors0.at(0)->CCD().N(X)-1);
		for (int y(Detectors0.at(0)->CCD().N(Y)/2-1);y<=Detectors0.at(0)->CCD().N(Y)/2+1;y++) {
			for (int d(0);d<Detectors0.size();d++) {
				if (Detectors0.at(d)->geo()(x,y,0) >=-90. &&
					 	Detectors0.at(d)->geo()(x,y,0) <= 90.) {// i.e. not NaN's
					double thisphase(Detectors0.at(d)->geo()(x,y,4));
					double thisincidence(Detectors0.at(d)->geo()(x,y,5));
					double thisemission(Detectors0.at(d)->geo()(x,y,6));
					for (int i(0);i<pp0.N(X);i++) {  // incidence
						double incidence(result0.Axis(X,i));
						if (abs(thisincidence-incidence)<acceptancecone) {
							for (int e(0);e<pp0.N(Y);e++) {	// emission
								double emission(result0.Axis(Y,e));
								if (abs(thisemission-emission)<acceptancecone) {
									double desiredphase(incidence + emission);
									if (abs(thisphase-desiredphase)<acceptancecone  &&
										 abs(thisincidence-incidence)<acceptancecone) {
										weights(i,e,0) += 1.;
										for (int w(0);w<result0.N(Z);w++) {	// wavelength
											if (i==e) whichused.at(d)(x,y,i)=incidences.at(i).degrees();
											result0(i,e,w) += Detectors0.at(d)->CCD()(x,y,result0.Axis(Z,w));
											result1(i,e,w) += Detectors1.at(d)->CCD()(x,y,result0.Axis(Z,w));
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
	weights.write("weights.Jcube");
	result0/=weights;
	result1/=weights;
	result0.write("result0.Jcube");
	result1.write("result1.Jcube");
	for (int i(0);i<whichused.size();i++) 
		whichused.at(i).write(string("whichused_p")+int2str(int(phases.at(i).degrees()))+string(".Jcube"));
	
	
// estimate errors
	cube rterror0(result0*0.);
	cube rterror1(result1*0.);
	for (int x(0);x<Detectors0.at(0)->CCD().N(X);x++) {
		for (int y(Detectors0.at(0)->CCD().N(Y)/2-1);y<=Detectors0.at(0)->CCD().N(Y)/2+1;y++) {
			for (int d(0);d<Detectors0.size();d++) {
				if (Detectors0.at(d)->geo()(x,y,0) >=-90. &&
					 	Detectors0.at(d)->geo()(x,y,0) <= 90.) {// i.e. not NaN's
					double thisphase(Detectors0.at(d)->geo()(x,y,4));
					double thisincidence(Detectors0.at(d)->geo()(x,y,5));
					double thisemission(Detectors0.at(d)->geo()(x,y,6));
					for (int i(0);i<pp0.N(X);i++) {  // incidence
						double incidence(result0.Axis(X,i));
						if (abs(thisincidence-incidence)<acceptancecone) {
							for (int e(0);e<pp0.N(Y);e++) {	// emission
								double emission(result0.Axis(Y,e));
								if (abs(thisemission-emission)<acceptancecone) {
									double desiredphase(incidence + emission);
									if (abs(thisphase-desiredphase)<acceptancecone  &&
										 abs(thisincidence-incidence)<acceptancecone) {
										for (int w(0);w<result0.N(Z);w++) {	// wavelength
											rterror0(i,e,w) += pow(Detectors0.at(d)->CCD()(x,y,result0.Axis(Z,w))-result0(i,e,w) ,2.);
											rterror1(i,e,w) += pow(Detectors1.at(d)->CCD()(x,y,result0.Axis(Z,w))-result1(i,e,w) ,2.);
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
	rterror0 /= weights;
	rterror1 /= weights;
	rterror0 =  rterror0.pow(0.5);
	rterror1 =  rterror1.pow(0.5);
	weights = weights.pow(0.5);
	rterror0 /= weights;  // to get to std dev of the mean
	rterror1 /= weights; 
	rterror0.write("rterror0.Jcube");
	rterror1.write("rterror1.Jcube");
	
	
	
// Generate nice plots



	for (int i(0);i<incidences.size();i++) {
		cube a0result(result0.skewer(Z,i,i));
		a0result = a0result.blocks(X,rterror0.skewer(Z,i,i));
		a0result.keyword("WITHERROR", "yes");
		a0result = a0result.chunk(Z,0,6).blocks(Z,a0result.plane(Z,11));
		a0result.graph("", oGIF);
		pp0.skewer(Z,i,i).graph("", oGIF);
		
		cube a1result(result1.skewer(Z,i,i));
		a1result = a1result.blocks(X,rterror1.skewer(Z,i,i));
		a1result.keyword("WITHERROR", "yes");
		a1result = a1result.chunk(Z,0,6).blocks(Z,a1result.plane(Z,11));
		a1result.graph("", oGIF);
		pp1.skewer(Z,i,i).graph("", oGIF);
	}
			
	string command(string("graph -T "));
	command+=asstring(outputtype);
	command+=string(" -X \"Wavelength (um)\" -Y \"I/F\" -L \"SebastienCompare A=0\" ");
	command+=string(" --pen-colors \"1=red:2=purple:3=blue:4=darkgreen:5=black\" -C ");	
	command+=string(" -x 0.8 5.2  -y 0 0.26 "); 		

	for (int e(0);e<incidences.size();e++) {
		command+=string(" -S 0 -m ");
		command+=int2str((e%5)+1);
		command+=string(" -I a graph");	
		command+=int2str(int(graphno-4*incidences.size()+1+4*e));
		command+=string(".out ");
		
		command+=string(" -S 3 -m -");
		command+=int2str((e%5)+1);
		command+=string(" -I e graph");	
		command+=int2str(int(graphno-4*incidences.size()+0+4*e));
		command+=string(".out ");
	}
		
	command+=string(" > Sebastiencompare_a0");
	command+=string(".");
	command+=asstring(outputtype);
	system(command.c_str());
	cout << command << "\n";
	
	{
		ofstream graphcommands("graph.commands", ios_base::app);
		graphcommands << command << "\n";
	}
	
	cout << "outputtype = " << outputtype << "\n";
		
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

	
	double diff0mean(0.), diff1mean(0.);
	unsigned int d0n(0), d1n(0);
	for (int x(0);x<diff0.N(X);x++)
		for (int y(0);y<diff0.N(Y);y++)
			for (int z(0);z<diff0.N(Z);z++) {
				if (isnormal(diff0(x,y,z))) {
					d0n++;
					diff0mean+=diff0(x,y,z);
				}
				if (isnormal(diff1(x,y,z))) {
					d1n++;
					diff1mean+=diff1(x,y,z);
				}
			}
	return (diff0mean+diff1mean)/double(d0n+d1n);
}

