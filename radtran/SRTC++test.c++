#include "SRTC++.h"

main(int argn, char **argv)
{	
	osuppress = 0;
	SRTC::debug = 0;
	SRTC::ompdebug = 0; 
	
	cout << "about to create atmosphere()\n";
//	atmosphere a(atmosphere::timetrialorangerind());
//	atmosphere a(atmosphere::mateorangerind());
//	atmosphere a(atmosphere::orangerind(1., value(100., "km"), 1.));
	atmosphere a(atmosphere::DISR08());
	a.normal_extinction();
	exit(0);
	
//	cout << "after atmosphere creation\n";
	
	vector<value> wavelengths;

	wavelengths.push_back(value(0.93,"um"));
   //wavelengths.push_back(value(1.06,"um"));
  // wavelengths.push_back(value(1.28,"um"));
  // wavelengths.push_back(value(1.60, "um"));
 //  wavelengths.push_back(value(2.,"um"));
  // wavelengths.push_back(value(2.75,"um"));
//   wavelengths.push_back(value(5., "um"));
		

	
	cout <<" Size of wavelength vector: "<<wavelengths.size()<<endl;	

	
	//// 	North hemisphere albedo 1, south hemisphere albedo 0 
/* 	cube albedoone(1,2,1);
	albedoone(0,0,0) = 1;
	albedoone(0,1,0) = 0;
	albedoone.Axis(Y,0) = -1.;
	albedoone.Axis(Y,1) = 1.;
	 */
//	cube hirespattern("testpattern.jpg")
//	cube backmap();
//	backmap = atmozone::bitmap_west180_to_cylmap(backmap);
	
	//// uniform surface albedo 
	//cube albedoone(1,1,1);
	//albedoone(0,0,0)=0.2; //titan grey  
	//albedoone.Axis(Y,0)=1.;
	
	//// hibob surface
	
	latangle top(angle(-90.,angle::DEG),NORTH_NEGATIVE), bottom(angle(90.,angle::DEG),NORTH_NEGATIVE);
	lonangle left(angle(-181.,angle::DEG),CENTER_0), right(angle(180.,angle::DEG),CENTER_0);
	
	cube albedoone(1,2,1);
	albedoone(0,0,0)=1.0; 
	albedoone.Axis(Y,0)=-1.;
	albedoone(0,1,0)=0.0; 
	albedoone.Axis(Y,0)=1.;
	
	cube hibob("hibob_map.jpg");
 	////hibob = atmozone::bitmap_west180_to_cylmap(hibob);
 	hibob*= 0.17;
 	hibob+= 0.03;
 	hibob = hibob.plane(Z,0);
 	hibob.Axis(Z,0) = 2.;
 	//hibob = atmozone::bitmap_ew0_to_cylmap(hibob);
 	hibob = atmozone::bitmap_to_cylmap(hibob,top,bottom,left,right); 	
 	//hibob.write("hibob_map.Jcube"); 
	//hibob=hibob.cylindereuler(0., 0., 0., mLINEAR);
	//hibob.write("hibob_map_rotated.Jcube"); 
	
 	//hibob = atmozone::bitmap_to_cylmap(hibob,top,bottom,left,right); 	
 	surface_atmozone *hibobzone=new surface_atmozone(albedoone);
 	a.front().addzone(hibobzone);
			
	//// VIMS map
	/*cube VIMSmap("/vims/maps/Titan.T8T9T34.cyl.Jcube.albedo.jpg");
	VIMSmap+= 0.03;
	VIMSmap = VIMSmap.plane(Z,0);
	VIMSmap.Axis(Z,0) = 2.;
	VIMSmap.write("TitanT8T9T34_bitmapped.Jcube");*/	

	//// VIMS map coadded
	/*cube VIMSmap("/vims/maps/Titan.simplecoadddeg2bestAdiri.cyl.Jcube");
	VIMSmap = atmozone::VIMScube_to_albedomap(VIMSmap,wavelengths);
	VIMSmap.write("VIMScube_to_albedomap.Jcube");	*/
			
	//// VIMS hand drawn		
	/* cube handdrawn("VIMSalbedomap_byhand.jpg");
	handdrawn = atmozone::bitmap_west180_to_cylmap(handdrawn);
	handdrawn*= 0.2;
	handdrawn+= 0.03;
	handdrawn = handdrawn.plane(Z,0);
	handdrawn.Axis(Z,0) = 2.;
	handdrawn.write("Titanglobalcoadd_bitmapped.Jcube"); 	 
	
	latangle top(angle(-90,angle::DEG),NORTH_NEGATIVE), bottom(angle(90,angle::DEG),NORTH_NEGATIVE);
	lonangle left(angle(-360,angle::DEG),CENTER_NEGATIVE180), right(angle(0,angle::DEG),CENTER_NEGATIVE180);
	cube test("testpattern.jpg");
	latangle centerlat(angle(0.,angle::DEG),NORTH_NEGATIVE);
	lonangle centerlon(angle(-180.,angle::DEG),CENTER_NEGATIVE180);
	angle pixelres(2.5e-11,angle::RAD);
	test = atmozone::bitmap_to_cylmap(test,centerlat,centerlon,pixelres, CENTER_NEGATIVE180);
	test*= 0.2;
	test+= 0.03;
	test = test.plane(Z,0);
	test.Axis(Z,0) = 2.;
	test.write("TEST_bitmapped.Jcube"); 	  */
	
	/* Setting up the surface layer*/
  	/* 
	cube black("testbar_black.jpg");
	latangle topb(angle(-90.,angle::DEG),NORTH_NEGATIVE), bottomb(angle(90.,angle::DEG),NORTH_NEGATIVE);
	lonangle leftb(angle(-361.,angle::DEG),CENTER_NEGATIVE180), rightb(angle(0.,angle::DEG),CENTER_NEGATIVE180);
	cout <<"about to read in black bitmap\n";
	black = atmozone::bitmap_to_cylmap(black,topb,bottomb,leftb,rightb);
	cout <<"setting atmozone to black cube\n";
	atmozone blackzone(black);
	cout <<"successfuly set atmozone black!\n";
	
	cube white("testboxwhite.jpg");
	latangle top(angle(10.,angle::DEG),NORTH_NEGATIVE), bottom(angle(-10.,angle::DEG),NORTH_NEGATIVE);
	lonangle left(angle(180.,angle::DEG),CENTER_NEGATIVE180), right(angle(-170.,angle::DEG),CENTER_NEGATIVE180);
	
	latangle centerlat(angle(0.,angle::DEG),NORTH_NEGATIVE);
	lonangle centerlon(angle(-160.,angle::DEG),CENTER_0);
	angle pixelres(0.27,angle::DEG);
	
	cout <<"about to read in white bitmap\n";
	white = atmozone::bitmap_to_cylmap(white,centerlat,centerlon,pixelres);//top,bottom,left,right);
	white.write("whitecube.Jcube");
	cout <<"setting atmozone to white cube\n";
	atmozone whitezone(white);
	cout <<"successfuly set atmozone white!\n";
	
	cube test("testpattern.jpg");	 		 
	centerlon=lonangle(angle(-180.,angle::DEG),CENTER_0);
	
	
	top = latangle(angle(-10.,angle::DEG),NORTH_NEGATIVE);
	bottom = latangle(angle(10.,angle::DEG),NORTH_NEGATIVE);
	left = lonangle(angle(-190.,angle::DEG),CENTER_NEGATIVE180);
	right = lonangle(angle(-170.,angle::DEG),CENTER_NEGATIVE180);
	
	cout <<"about to read in test bitmap\n";
//	test = atmozone::bitmap_to_cylmap(test,centerlat,centerlon,pixelres);
	test = atmozone::bitmap_to_cylmap(test,top, bottom, left, right);
	cout <<"test bitmap created\n";	
	test.write("testpattern_bitmapped.jpg");	
	
	atmozone testzone(test);
			
	a.front().addzone(blackzone);
	a.front().addzone(whitezone);
	a.front().addzone(testzone); */
				 	
	cout << "This atmosphere size = " << a.size()<<" "<<a.front().size() << "\n";
		
//	a.front().at(0).ss_albedo(albedoone);
//	a.front().at(0).ss_albedo(hirespattern);
//	a.front().at(1).ss_albedo(backmap);
//	a.front().at(0).ss_albedo(hibob);
// a.front().at(0).ss_albedo(VIMSmap);
 //a.front().at(0).ss_albedo(test);
// cout<<"created surface layer with ss_albedo from cube TEST\n";
 //atmozone coarsemap(albedoone);
 //cout <<"created atmozone from albedoone cube\n";
 //a.front().push_back(coarsemap);
 
	
	cout << "set surface albedo map \n";
	

	
/*	geomvector h(12., 15., 17.);
	cout << h << "\n";
	
	geomvector g(15., angle(45., angle::DEGREES), angle(90., angle::DEGREES));	
	cout << g << "\n";
	g.to(CARTESIAN);
	cout << " or " << g << "\n";*/
			
			
	bool hires(false);		
	bool onedetector(false);	
	

	//photongenerator_square hv(value(4200., "km"), 1200.,value(4200., "km"), wavelengths);	
	//	photongenerator_square hv(value(12002.,"km"),1000.,value(4200.,"km"),wavelengths); //image total SNR = 10
	//photongenerator_square hv(value(0.2,"km"),1270.,value(4200.,"km"),wavelengths); //50 SNR per 25x25 m pixel
	
	photongenerator_square hv(value(.8,"km"),430.,value(4200.,"km"),wavelengths);
	//hv.makehole(value(0.8,"km"));
	geomvector incomingdirection(-1.,0.,0.);  // set the photons to come in initially toward -x
	
	
	if (!hires)	
		hv =  photongenerator_square(value(3250., "km"), 260L,value(4200., "km"), wavelengths);
//		hv =  photongenerator_square(value(4200., "km"), 140,value(4200., "km"), wavelengths);
	
	
	
	hv.photondirection(incomingdirection);
	hv.fieldcenter(geomvector(4200.,0.,0.));
	
	cout <<" Size of wavelength vector: "<<wavelengths.size()<<endl;	
	vector<detector*> Detectors;
	
	bool detecthistory(true); //set to true to use elephant detectors
	bool detectnormaltoo(true); // if on while detecthistory is on, do both
	int ndetectors(18); 
	if (onedetector) ndetectors =0;
	for (int i(0);i<=ndetectors;i++) {
		angle a(double(i)*10., angle::DEG);
		double x(1.e9*a.cos()); 
		double y(1.e9*a.sin());
		
		colorCCD *CCD_detector;
		
		CCD_detector = new colorCCD(geomvector(x,y,0.),value(3200., "km"),value(100., "km"),wavelengths);
				
		string s("phase_");
		s+=int2str(i,2);
		CCD_detector->name(s);
		
		if (detecthistory) {
			elephant *even_bigger_detector;
			even_bigger_detector = new elephant(*CCD_detector);
			Detectors.push_back(even_bigger_detector);
		}		
		if (!detecthistory || detectnormaltoo) 
			Detectors.push_back(CCD_detector);
	}
	
	
	
	/* 
	colorCCD zerophase_detector(geomvector(-1.e9,0.,0.),value(3100., "km"), value(10., "km"), wavelengths);
	//colorCCD zerophase_detector(geomvector(-1.e9,0.,0.),1, angle(2e-11,angle::RAD), wavelengths);
	
	if (!hires)
		zerophase_detector = colorCCD(geomvector(-1.e9,0.,0.),value(3000., "km"), value(600., "km"), wavelengths);
	zerophase_detector.name("phase_000");
	Detectors.push_back(&zerophase_detector);
	
	
	colorCCD fortyfivephase_detector(geomvector(-1.e9,-1.e9,0.),value(3100., "km"), value(10., "km"), wavelengths);
 	if (!hires)
 		fortyfivephase_detector= colorCCD(geomvector(-1.e9,-1.e9,0.),value(3000., "km"), value(600., "km"), wavelengths);
 	fortyfivephase_detector.name("phase_045");
 	if (!onecamera) Detectors.push_back(&fortyfivephase_detector);
 			
 	colorCCD ninetyphase_detector(geomvector(0.,-1.e9,0.),value(3100., "km"), value(10., "km"), wavelengths);
 	if (!hires)
 		ninetyphase_detector= colorCCD(geomvector(0.,-1.e9,0.),value(3000., "km"), value(600., "km"), wavelengths);
	ninetyphase_detector.name("phase_090");
 	if (!onecamera) Detectors.push_back(&ninetyphase_detector);
 	
 	colorCCD onethirtyfivephase_detector(geomvector(1.e9,-1.e9,0.),value(3100., "km"), value(10., "km"), wavelengths);
 	if (!hires)
 		onethirtyfivephase_detector= colorCCD(geomvector(1.e9,-1.e9,0.),value(3000., "km"), value(600., "km"), wavelengths);
 	onethirtyfivephase_detector.name("phase_135");
 	if (!onecamera) Detectors.push_back(&onethirtyfivephase_detector);	
 	
 	colorCCD oneeightyphase_detector(geomvector(1.e9,0.,0.),value(3100., "km"), value(10., "km"), wavelengths);
	if (!hires)
 		oneeightyphase_detector= colorCCD(geomvector(1.e9,0.,0.),value(3000., "km"), value(600., "km"), wavelengths);
 	oneeightyphase_detector.name("phase_180");
 	if (!onecamera) Detectors.push_back(&oneeightyphase_detector); */
	
	
	SRTC S(a, &hv, Detectors);
	
	S.run();
	
	cout << "Entering writing loop\n"; cout.flush();
	for (int i(0);i<Detectors.size();i++)
		Detectors.at(i)->write(Detectors.at(i)->name());
}
