#include "SRTC++.h"

atmosphere::atmosphere() 
	: list<atmolayer>()
{ 
	cout << "size = " << size() << "\n";
}

void atmosphere::addlayer(atmolayer thislayer)
{
	push_back(thislayer);
	for (int i(0);i<thislayer.size();i++)
		back().at(i)->toparentlayer(&back());
	
	back().toabove(0);
	back().tobelow(0);
	if (size()>1) {
		list<atmolayer>::iterator i(end());
		i--; // last one
		list<atmolayer>::iterator j(end());
		j--; j--; // second-to-last one

		i->tobelow(&(*j));
		j->toabove(&(*i));
	}
}

double atmosphere::totaltau(value w, geomvector g)
{
	double answer(0.);
	
	for (iterator i(begin());i!=end();i++) {
		if (!(i->issurface())) {
			answer += i->totaltau(w,g);
		}
	}
	
	return answer;
}

atmosphere atmosphere::orangerind(cube taus, value vertical_extent, cube ssalbedos)
{
	static phasefunction_isotropic isotropicphasefunction;		// assumed analytic isotropic	-- override with zone as nec
	
	atmolayer rind(vertical_extent, &isotropicphasefunction);
	rind.layername(string("Orange Rind Tau=")+double2str(taus(0,0,0))+string(" extent=")+vertical_extent.tostring());
	
	cube extinctions(ssalbedos);
	for (int z(0);z<extinctions.N(Z);z++)
		extinctions(0,0,z) = taus(0,0,z)/rind.depth_km();
 
	cout <<"Extinction is\t"<<extinctions(0,0,0)<<endl;
	orangerind_atmozone *whitestuff=new orangerind_atmozone(ssalbedos);  // single-scattering albedo
	whitestuff->extinctions(extinctions);
	rind.addzone(whitestuff);
	
	atmosphere answer;
	answer.addlayer(atmolayer::default_surface(value(2575., "km")));
	answer.addlayer(rind);
	
	return answer;	
	
	
}

atmosphere atmosphere::orangerind(double tau, value vertical_extent, cube ssalbedos)
{	
	cube taus(1,1,1,tau);
	
	return orangerind(taus, vertical_extent, ssalbedos);
}


atmosphere atmosphere::orangerind(double tau, value vertical_extent, double ssalbedo)
{	
	cube ssalbedo_cube(1,1,1);
	ssalbedo_cube.Axis(Z,0) = 1.0;
	ssalbedo_cube(0,0,0) = ssalbedo;
	
	
	return orangerind(tau, vertical_extent, ssalbedo_cube);
}

atmosphere atmosphere::SebastienValidateDISR08_layered()
{
	/* Atmosphere created with properties from the aerosol model 
	   derived from DISR obseravations by Tomasko et al. 2008    
		With SSA albedoes from Hirtzig 2013 extrapolation */
	
	/* used in test_SebastienCompare */
		
	/* Constructing zones from different altitude bins of single scatter albedo*/
	// the X value is the BOTTOMMOST extent of the layer: N(X) = 79
	
	cube ssa_cube(SRTC::datapath()+string("/data/inputfiles/ssa_2008Tomasko_HirtzigExtrapolated.txt.Jcube"));
	double collapsefactor(100.);
	for (int x(0);x<ssa_cube.N(X);x++) {
		ssa_cube.Axis(X,x)/=collapsefactor;
		ssa_cube.Axis(X,x)+=2575.;
	}
	ssa_cube = ssa_cube.xt2xy();
	ssa_cube = ssa_cube.dirinc();
	ssa_cube.write("manhandled_ssacube.Jcube");
	
	vector<pair<double, cube> > ssa_cubes_with_altitude;
	for (int x(0);x<ssa_cube.N(X);x++)
		ssa_cubes_with_altitude.push_back(pair<double,cube>(ssa_cube.Axis(X,x), ssa_cube.plane(X,x)));
	
	
	cube wave_cube(SRTC::datapath()+string("/data/inputfiles/wavelength.cal.Jcube")); // all the wavelengths for calc tau later in the loop

	cube hpf_Above(cube(SRTC::datapath()+string("/data/inputfiles/DISRhazephasedata_above80km.Jcube")));
	hpf_Above.keyword("scattering_name", "hpf_Above");			
	static phasefunction_cube hpfcube_Above(hpf_Above, angle::DEGREES);
	
	cube hpf_Below(cube(SRTC::datapath()+string("/data/inputfiles/DISRhazephasedata_below80km.Jcube")));
	hpf_Below.keyword("scattering_name", "hpf_Below");			
	static phasefunction_cube hpfcube_Below(hpf_Below, angle::DEGREES);				
			
	
	atmolayer bottomlayer(value(30./collapsefactor,"km"), &hpfcube_Below);
	value bottomlayerscaleheight(50./collapsefactor, "km");  // this numer is irrelevant as DISR_Tomasko is uniform	
	bottomlayer.scaleheight(bottomlayerscaleheight);	
	bottomlayer.layername("troposphere");
	
	atmolayer middlelayer(value(50./collapsefactor, "km"), &hpfcube_Below);
	value middlelayerscaleheight(50./collapsefactor, "km");  // this numer is irrelevant as DISR_Tomasko is uniform	
	middlelayer.scaleheight(middlelayerscaleheight);
	middlelayer.layername("lower stratosphere");
	
	atmolayer toplayer(value(260./collapsefactor, "km"), &hpfcube_Above);
	value toplayerscaleheight(65./collapsefactor, "km");  // From Tomasko
	toplayer.scaleheight(toplayerscaleheight);
	toplayer.layername("outer atmosphere");
	
	if(osuppress<-4) cout << "About to add zones.\n";
	
	
// setting up the bottom zone.
	cout <<"\n----About to initialize bottomzone\n";
	orangerind_atmozone *bottomzone=new orangerind_atmozone(ssa_cubes_with_altitude);
	
	DISR08_tau_function bottom_tau_function(6.27e2, -0.9706, bottomlayer.depth(), bottomlayer.scaleheight());  // values from Tomasko et al. 2008
	cube thiszone_extinctions_0(bottom_tau_function.plot(300., 6000., 1000));
	
	for (int z(0);z<thiszone_extinctions_0.N(Z);z++)
		thiszone_extinctions_0.Axis(Z,z)/=1000.;
	thiszone_extinctions_0 /= bottomlayer.depth_km();
	thiszone_extinctions_0.write("bottomzone_extinctions.Jcube");
	if(osuppress<-4) cout << "about to zonetau()\n";
	bottomzone->extinctions(thiszone_extinctions_0);
	
	if(osuppress<-4) cout << "zonetau(setter) complete.  Now about to addzone to bottomlayer\n";
	bottomlayer.addzone(bottomzone);
	if(osuppress<-4) cout << "bottom zone addition complete!\n";
	if(osuppress<-4) cout << " bottom zone depth = " << bottomlayer.depth() << ", scaleheight = " << bottomlayer.scaleheight() << "\n";
	
	
// setting up the middle zone.
	cout <<"\n----About to initialize middlezone\n";
	orangerind_atmozone *middlezone=new orangerind_atmozone(ssa_cubes_with_altitude);
	
	DISR08_tau_function middle_tau_function(2.029e4, -1.409, middlelayer.depth(), middlelayer.scaleheight());  // values from Tomasko et al. 2008
	thiszone_extinctions_0 = middle_tau_function.plot(300., 6000., 1000);
	for (int z(0);z<thiszone_extinctions_0.N(Z);z++)
		thiszone_extinctions_0.Axis(Z,z)/=1000.;
	thiszone_extinctions_0 /= middlelayer.depth_km();
	thiszone_extinctions_0.write("middlezone_extinctions.Jcube");
	middlezone->extinctions(thiszone_extinctions_0);
	middlelayer.addzone(middlezone);
	if(osuppress<-4) cout << "Middle zone achieved\n";
	
	
// setting up the top zone.
	cout <<"\n----About to initialize topzone\n";
	DISR_Tomasko_atmozone *topzone=new DISR_Tomasko_atmozone(ssa_cubes_with_altitude);
	DISR08_tau_function top_tau_function(1.012e7, -2.339, toplayer.depth(), toplayer.scaleheight());  // values from Tomasko et al. 2008
	thiszone_extinctions_0 = top_tau_function.plot(300., 6000., 1000);
	for (int z(0);z<thiszone_extinctions_0.N(Z);z++)
		thiszone_extinctions_0.Axis(Z,z)/=1000.;
	
		
	thiszone_extinctions_0 /= (toplayer.scaleheight_km()*(1.-exp(-(toplayer.depth_km())/toplayer.scaleheight_km())));

	
	thiszone_extinctions_0.write("topzone_dTau.Jcube");
	topzone->zone_haze_extinction_0(thiszone_extinctions_0);
	toplayer.addzone(topzone);

	
	if(osuppress<-4) cout << "Zones added in to layers\n";
				
	atmosphere answer;
	// ALWAYS START FROM BOTTOM
	
	if(osuppress<-4) cout <<"atmosphere created\n";
	answer.addlayer(atmolayer::default_surface(value(2575.,"km")));	
	answer.addlayer(bottomlayer);	
	answer.addlayer(middlelayer);	
	answer.addlayer(toplayer);
		
	return answer;
}

atmosphere atmosphere::DISR08()
{
	/* Atmosphere created with properties from the aerosol model 
	   derived from DISR obseravations by Tomasko et al. 2008    */
	
	/* Constructing zones from different altitude bins of single scatter albedo*/
	cube ssa_cube(SRTC::datapath()+string("/data/inputfiles/SingleScatteringAlbedo_DISR_withHirtzigExtrapolation.txt.Jcube"));
	
	double collapsefactor(100.);
	
	ssa_cube.Axis(X,0)/=collapsefactor;
	ssa_cube.Axis(X,1)/=collapsefactor;
	ssa_cube.Axis(X,2)/=collapsefactor;
	
	ssa_cube = ssa_cube.xt2xy();
	
	cube wave_cube(SRTC::datapath()+string("/data/inputfiles/wavelength.cal.Jcube")); // all the wavelengths for calc tau later in the loop

	cube hpf_Above(SRTC::datapath()+string("/data/inputfiles/DISRhazephasedata_above80km.Jcube"));
	hpf_Above.keyword("scattering_name", "hpf_Above");			
	static phasefunction_cube hpfcube_Above(hpf_Above, angle::DEGREES);
	
	cube hpf_Below(SRTC::datapath()+string("/data/inputfiles/DISRhazephasedata_below80km.Jcube"));
	hpf_Below.keyword("scattering_name", "hpf_Below");			
	static phasefunction_cube hpfcube_Below(hpf_Below, angle::DEGREES);				
	
	static phasefunction_isotropic isotrop_phasefunction;
	
	cube hpf_allforward(SRTC::datapath()+string("/data/inputfiles/forwardscatteringhazephasedata.txt.Jcube"));
	hpf_allforward.keyword("scattering_name", "allforward_cube");
	static phasefunction_cube forwardscattering(hpf_allforward, angle::DEGREES);
	
	cube hpf_allbackward(hpf_allforward);
	for (int i(0);i<hpf_allforward.N(Z);i++) {
		hpf_allbackward(0,0,i) = hpf_allforward(0,0,hpf_allforward.N(Z)-1-i);
	}
	static phasefunction_cube backscattering(hpf_allbackward, angle::DEGREES);
		
	
	atmolayer pseudosurface(value(1., "km"), &isotrop_phasefunction);
	value pseudosurfacescaleheight(50., "km");
	pseudosurface.scaleheight(pseudosurfacescaleheight);
	pseudosurface.layername("pseudosurface");
				
	
	
	atmolayer bottomlayer(value(30./collapsefactor,"km"), &hpfcube_Below);
	value bottomlayerscaleheight(50., "km");  // this numer is irrelevant as DISR_Tomasko is uniform	
	bottomlayer.scaleheight(bottomlayerscaleheight);	
	bottomlayer.layername("troposphere");
	
	atmolayer middlelayer(value(50./collapsefactor, "km"), &hpfcube_Below);
	value middlelayerscaleheight(50., "km");  // this numer is irrelevant as DISR_Tomasko is uniform	
	middlelayer.scaleheight(middlelayerscaleheight);
	middlelayer.layername("lower stratosphere");
	
	atmolayer toplayer(value(260./collapsefactor, "km"), &hpfcube_Above);
	value toplayerscaleheight(65./collapsefactor, "km");  // From Tomasko
	toplayer.scaleheight(toplayerscaleheight);
	toplayer.layername("outer atmosphere");
	
	if(osuppress<-4) cout << "About to add zones.\n";
	
	
	orangerind_atmozone *pseudosurfacezone=new orangerind_atmozone(cube(1,1,1,1.0));
	pseudosurfacezone->extinctions(cube(1,1,1,10.));
	pseudosurface.addzone(pseudosurfacezone);
	
	
	
	
// setting up the bottom zone.
	cout <<"\n----About to initialize bottomzone\n";
	orangerind_atmozone *bottomzone=new orangerind_atmozone(ssa_cube.chunk(X,2,2));
	DISR08_tau_function bottom_tau_function(6.27e2, -0.9706, bottomlayer.depth(), bottomlayer.scaleheight());  // values from Tomasko et al. 2008
	cube thiszone_extinctions_0(bottom_tau_function.plot(300., 6000., 1000));
	
	for (int z(0);z<thiszone_extinctions_0.N(Z);z++)
		thiszone_extinctions_0.Axis(Z,z)/=1000.;
	thiszone_extinctions_0 /= bottomlayer.depth_km();
	thiszone_extinctions_0.write("bottomzone_extinctions.Jcube");
	if(osuppress<-4) cout << "about to zonetau()\n";
	bottomzone->extinctions(thiszone_extinctions_0);
	
	
	cout <<" thiszone_haze_extinctions_0 size:\t"<<thiszone_extinctions_0.info()<<endl;
	cout <<" thiszone_haze_extinctions_0 cube:\t"<<thiszone_extinctions_0<<endl;
	if(osuppress<-4) cout << "zonetau(setter) complete.  Now about to addzone to bottomlayer\n";
	bottomlayer.addzone(bottomzone);
	if(osuppress<-4) cout << "bottom zone addition complete!\n";
	if(osuppress<-4) cout << " bottom zone depth = " << bottomlayer.depth() << ", scaleheight = " << bottomlayer.scaleheight() << "\n";
	
	
// setting up the middle zone.
	cout <<"\n----About to initialize middlezone\n";
	orangerind_atmozone *middlezone=new orangerind_atmozone(ssa_cube.chunk(X,1,1));
	DISR08_tau_function middle_tau_function(2.029e4, -1.409, middlelayer.depth(), middlelayer.scaleheight());  // values from Tomasko et al. 2008
	thiszone_extinctions_0 = middle_tau_function.plot(300., 6000., 1000);
	for (int z(0);z<thiszone_extinctions_0.N(Z);z++)
		thiszone_extinctions_0.Axis(Z,z)/=1000.;
	thiszone_extinctions_0 /= middlelayer.depth_km();
	thiszone_extinctions_0.write("middlezone_extinctions.Jcube");
	middlezone->extinctions(thiszone_extinctions_0);
	middlelayer.addzone(middlezone);
	if(osuppress<-4) cout << "Middle zone achieved\n";
	
	
// setting up the top zone.
	cout <<"\n----About to initialize topzone\n";
	DISR_Tomasko_atmozone *topzone=new DISR_Tomasko_atmozone(ssa_cube.chunk(X,0,0));
	DISR08_tau_function top_tau_function(1.012e7, -2.339, toplayer.depth(), toplayer.scaleheight());  // values from Tomasko et al. 2008
	thiszone_extinctions_0 = top_tau_function.plot(300., 6000., 1000);
	for (int z(0);z<thiszone_extinctions_0.N(Z);z++)
		thiszone_extinctions_0.Axis(Z,z)/=1000.;
	
		
	thiszone_extinctions_0 /= (toplayer.scaleheight_km()*(1.-exp(-(toplayer.depth_km())/toplayer.scaleheight_km())));

	
	thiszone_extinctions_0.write("topzone_dTau.Jcube");
	topzone->zone_haze_extinction_0(thiszone_extinctions_0);
	toplayer.addzone(topzone);

	
	if(osuppress<-4) cout << "Zones added in to layers\n";
	
	
	
			
	atmosphere answer;
	// ALWAYS START FROM BOTTOM
	
	if(osuppress<-4) cout <<"atmosphere created\n";
	answer.addlayer(atmolayer::default_surface(value(2575.,"km")));
	answer.addlayer(bottomlayer);
	answer.addlayer(middlelayer);
	answer.addlayer(toplayer);
		
	return answer;
		
}

value atmosphere::outeredge()
{
	return value(outeredge_km(), "km");
}

double atmosphere::outeredge_km()
{
	return back().top_km();
}

value atmosphere::surfacealtitude()
{
	return value(surfacealtitude_km(), "km");
}

double atmosphere::surfacealtitude_km()
{
	return front().depth_km();
}


DISR08_tau_function::DISR08_tau_function(double inconstant, double inexponent, value indepth, value inscaleheight) :
		_constant(inconstant), _exponent(inexponent)
{	
	
	_depth_km = double(indepth.convert("km"));
	_scaleheight_km = double(inscaleheight.convert("km"));
}

double DISR08_tau_function::the_function(double lambda_nm)
{
	double answer(0.);
	
	double Tau_total = _constant * pow(lambda_nm, _exponent);
	
	return Tau_total;
}

void atmosphere::normal_extinction()
{
	geomvector start_at_surface(surfacealtitude_km(),0.,0.);
	geomvector head_outward(1., 0., 0.);
	
	vector<double> wavelengths_to_test;			
	wavelengths_to_test.push_back(0.400);
	wavelengths_to_test.push_back(0.500);
	wavelengths_to_test.push_back(0.700);
	wavelengths_to_test.push_back(1.000);
	wavelengths_to_test.push_back(1.600);
	
	cube e;
	for (int i(0);i<wavelengths_to_test.size();i++) {
		photon outwardbound(start_at_surface, head_outward, 1., value(wavelengths_to_test.at(i), "um"));
			
		cout << "Running atmosphere tests for wavelength " << outwardbound.lambda() << ".\n";
		
		e = extinction_profile(outwardbound);
		e.graph();
	}
}

cube atmosphere::extinction_profile(photon& _inphoton)
{
	cout << "creating extinction profile.  photon wavelength is " << _inphoton.lambda() << "\n";
	photontraverse_extinction p_e(_inphoton, *this);
	
	return p_e.plot(0., 201., 201);
}


atmolayer* atmosphere::layerataltitude(value alt)
{
	atmolayer* answer(&front());
	
	double alt_km = double(alt.convert("km"));
	double dist_km = front().top_km()+alt_km;
	for (atmosphere::iterator i(begin()++);i!=end();i++) {
		if (dist_km>=i->bottom_km()  &&  i->top_km()>=dist_km)
			answer = &(*i);
	}
	if (dist_km>back().top_km()) answer = &back();
	
	return answer;
}


DISR_tau_function::DISR_tau_function(value inde)
{
	_depth_km = double(inde.convert("km"));
	
}


double DISR_tau_function::the_function(atmozone *thiszone, double lambda_um)
{
	double answer(0.);
	double Tau_total;
		if ( thiszone->toparentlayer()->bottom_km()> 55) { 
			
			Tau_total= _a / ( (_g + lambda_um *_f)*(_b+ exp(_depth_km/_c)) );
		
		} else {
		
			double E_55= (_a/_c *exp(_d/_c)) / ((_g +lambda_um*_f)*pow(_b + exp(_d/_c),2));
			Tau_total= _a / ( (_g + lambda_um *_f)*(_b+ exp(_d/_c)) ) + _h*E_55*(_d-_depth_km);
		
		}
	
}

void DISR_tau_function::calculate_scaleheight_fromtau(atmozone *thiszone, double wavelength_um){
	
	double answer;
	
	if ( thiszone->toparentlayer()->bottom_km() > 55) {
		
		answer= -exp(1)*pow(_b + exp(thiszone->toparentlayer()->bottom_km()/_c),2)*(_g+ pow(wavelength_um,_f))*_c;
		answer/=_a*exp(thiszone->toparentlayer()->bottom_km()/_c);
				
	} else {
		
		double E_55= (_a/_c *exp(_d/_c)) / ((_g +wavelength_um*_f)*pow(_b + exp(_d/_c),2));
		answer= - exp(1)/(_h*E_55*thiszone->toparentlayer()->bottom_km());
				
	}
	cout <<"Computed scale height at "<< thiszone->toparentlayer()->bottom_km() << "km to be "<< answer<<endl;
	_scaleheight_km=answer;
}

ostream& operator<<(ostream& out, atmosphere& a)
{
	out << "Atmosphere number " << &a << ":  " << a.size() << " layers:\n";
	int n(0);
	for (atmosphere::iterator i(a.begin());i!=a.end();i++) {
		out << "\tHere's layer #" << n << "\n";
		out << *i;
		n++;
	}
	return out << "------------------------------\n";
}
