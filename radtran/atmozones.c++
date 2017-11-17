#include "atmolayers.h"

atmozone* atmozone::_tospace=new space_atmozone(0.);

atmozone::atmozone(phasefunction* inp, double inalbedo) :
		_zonephasefunction(inp),
		_toparentlayer(0) 
{
	cube ssalbedocube(1,1,1);
	ssalbedocube(0,0,0) = inalbedo;
	ssalbedocube.Axis(Z,0) = 0.5;
	ss_albedo(inalbedo);
	initialsetup();
}

atmozone::atmozone(phasefunction* inp, cube incube) :
		_zonephasefunction(inp),
		_toparentlayer(0) 
{
	ss_albedo(incube);
	initialsetup();
}
	
atmozone::atmozone(phasefunction* inp, vector<pair<double, cube> > inSSA) :
		_zonephasefunction(inp),
		_toparentlayer(0) 
{
	for (int i(0);i<inSSA.size();i++)
		ss_albedo(inSSA.at(i).first,inSSA.at(i).second);
	initialsetup();
}
	
atmozone::atmozone(vector<pair<double, cube> > inSSA) :
		_zonephasefunction(0),
		_toparentlayer(0) 
{
	for (int i(0);i<inSSA.size();i++)
		ss_albedo(inSSA.at(i).first,inSSA.at(i).second);
	initialsetup();
}
		
atmozone::atmozone(double inalbedo) :
		_zonephasefunction(0),
		_toparentlayer(0) 
{
	cube ssalbedocube(1,1,1);
	ssalbedocube(0,0,0) = inalbedo;
	ssalbedocube.Axis(Z,0) = 0.5;
	ss_albedo(ssalbedocube);
	initialsetup();
}

atmozone::atmozone(cube incube) :
		_zonephasefunction(0),
		_toparentlayer(0) 
{  //updated 2016 August 29 by smack; now checks for if incube has lat/lon values stored along axes
	// the complicated if statement is to avoid selecting the DISR ssalbedo cube axes values (which correspond to 
	// altitude, not lat/lon)
	
	ss_albedo(incube);
	
	if (incube.Axis(X,incube.N(X)-1)!=0.  &&  incube.Axis(X,incube.N(X)-1)!=incube.Axis(X,0)){
		_left_lon=angle(incube.Axis(X,0), angle::DEG);
		_right_lon=angle(incube.Axis(X,incube.N(X)-1),angle::DEG);
		_top_lat=angle(incube.Axis(Y,0),angle::DEG);
		_bottom_lat=angle(incube.Axis(Y,incube.N(Y)-1),angle::DEG);
	}
	else initialsetup();
}

void atmozone::initialsetup()
{
	cout <<"ssalbedo cube does not have axes labeled in lat/lon. Assuming global image.\n";
	_left_lon=angle(-180., angle::DEG);
	_right_lon=angle(180., angle::DEG);
	_top_lat=angle(-90., angle::DEG);
	_bottom_lat=angle(90., angle::DEG);
}

bool atmozone::issurface()
{
	bool answer;
	if (toparentlayer())
		answer = toparentlayer()->issurface();
	else
		return false;
	return answer;
}

string atmozone::layername()
{
	string answer;
	
	if (toparentlayer()) 
		answer = toparentlayer()->layername();
	else
		answer = string("Space");
	
	return answer;
}

double atmozone::extinction_km_um(double radialdistance_km, double wavelength_um)
// This function computes the optical depth per km at a given radial
// distance from the center of the planet.  Presumably you have sent in
// a radialdistance for which this zone is valid.
{
	double answer(0.);
	
	
	if (&Space()!=this) {
		if ((radialdistance_km<toparentlayer()->bottom_km()-1  ||   
				radialdistance_km>toparentlayer()->top_km()+1)  &&
				!issurface()) {
			cout << "FAIL:  you've asked this zone for a height outside valid range\n";
			
			cout << "this radialdistance = " << radialdistance_km << " km; ";
			cout << "this layer's bottom is " << toparentlayer()->bottom_km() << "km; ";
			cout << "top is " << toparentlayer()->top_km() <<" km\n";
		}
		answer = extinction_gas_km_um(radialdistance_km, wavelength_um) + extinction_haze_km_um(radialdistance_km, wavelength_um);
	}
	
	return answer;
}






atmozone& atmozone::Space()
{
	return *_tospace;
}

phasefunction* atmozone::zonephasefunction()
{
	phasefunction *answer(_zonephasefunction);
	
	if (answer==0) answer = toparentlayer()->tophasefunction();
		
	return answer;
}

void atmozone::zonephasefunction(phasefunction *p)
{
	_zonephasefunction = p;
}




cube atmozone::bitmap_ew0_to_cylmap(cube incube)
{
	cube answer(incube.N(X),incube.N(Y),incube.N(Z));
	for (int x(0);x<incube.N(X);x++) answer.Axis(X,x) = -180.+double(x)*360./double(incube.N(X)-1);
	for (int y(0);y<incube.N(Y);y++) answer.Axis(Y,y) = -90.+double(y)*180./double(incube.N(Y)-1);
	
	for (int x(0);x<answer.N(X);x++)
		for (int y(0);y<answer.N(Y);y++)
			for (int z(0);z<answer.N(Z);z++)
				answer(x,y,z) = incube(x,y,z);
	
	return answer;
	
}

cube atmozone::bitmap_west180_to_cylmap(cube incube)
{
	cube answer(incube.N(X),incube.N(Y),incube.N(Z));
	for (int x(0);x<incube.N(X);x++) answer.Axis(X,x) = -180.+double(x)*360./double(incube.N(X)-1);
	for (int y(0);y<incube.N(Y);y++) answer.Axis(Y,y) = -90.+double(y)*180./double(incube.N(Y)-1);
	
	for (int x(0);x<answer.N(X)/2;x++)
		for (int y(0);y<answer.N(Y);y++)
			for (int z(0);z<answer.N(Z);z++) {
				answer(x+answer.N(X)/2,y,z) = incube(x,y,z);
				answer(x,y,z) = incube(x+answer.N(X)/2,y,z);
			}
	cout <<"bitmap_west180_to_cylmap incube info:\n";
	cout<<incube.info()<<endl;
	cout<<answer.info()<<endl;
	return answer;
}

cube atmozone::bitmap_to_cylmap(cube incube, latangle lat1, latangle lat2, lonangle lon1, lonangle lon2)
{ 
	// Takes cube created from jpg and maps it out as a proper Jcube, 
	// based on user input for the lat limits, lon limits (see Jangle.h) 
	// output cube will be in CENTER_0 longitude convention to match ssalbedo of atmospheric layers
	// added 26.Aug.2016 by smack
	cube answer(incube.N(X),incube.N(Y),incube.N(Z));
	
	cout <<"bitmapping cube to map:\n"<<endl;
	
	latangle top(lat1), bottom(lat2);
	lonangle left(lon1) , right(lon2);	
	lonconvention readinas(left.isstoredas());	
	cout <<"\t\ttop= "<<top.north_negative().degrees()<< endl;
	cout <<"left="<<left.center_negative180().degrees()<<"\t\t\t"<<" right "<<right.center_negative180().degrees()<<endl;
	cout <<"\t\t bottom= "<<bottom.north_negative().degrees()<<endl;
	
	if (lat1.north_positive().degrees()<lat2.north_positive().degrees()){
		top=lat2;
		bottom=lat1;	
	} 
	if ((left.center_0().degrees() <= -180. && right.center_0().degrees() >= 180.) ||
			(left.center_negative180().degrees() <=-360. && right.center_negative180().degrees() >= 0.) ){
		// checking for degenerate case
		left=lonangle(angle(-181.,angle::DEG),CENTER_0);
		right=lonangle(angle(180.,angle::DEG),CENTER_0);	
	}		
		
	double ydegreesperpixel(abs(top.north_positive().degrees()-bottom.north_positive().degrees())/(incube.N(Y)-1));
	double xdegreesperpixel(abs(right.center_0().degrees()-left.center_0().degrees())/(incube.N(X)-1));
	
	cout <<"\nx/y pixels per degree:\t"<<	xdegreesperpixel<<"\t"<<ydegreesperpixel<<endl;
		
	
	double leftmostlimit(left.center_0().degrees());
	int xlimit(answer.N(X)/2);
	if (readinas == CENTER_0) xlimit=answer.N(X);

	
	for (int x(0); x<incube.N(X);x++)
		answer.Axis(X,x)= leftmostlimit + double(x)*xdegreesperpixel  ;
	
	for (int y(0); y<incube.N(Y);y++)
		answer.Axis(Y,y)=top.north_negative().degrees()+double(y)*ydegreesperpixel;
	
	for (int x(0);x<xlimit;x++)
		for (int y(0);y<answer.N(Y);y++)
			for (int z(0);z<answer.N(Z);z++) {
				if (readinas == CENTER_0) answer(x,y,z)=incube(x,y,z);
				else {
					answer(x+xlimit,y,z) = incube(x,y,z);
					answer(x,y,z) = incube(x+xlimit,y,z);
				}
			}
		
	cout <<"\n after cube axes labeling, answer limits are:\n";
	cout <<"\t\t "<<answer.Axis(Y,0)<<endl;
	cout <<"\t"<<answer.Axis(X,0)<<"\t\t\t"<<answer.Axis(X,answer.N(X)-1)<<endl;
	cout <<"\t\t"<<answer.Axis(Y,answer.N(Y)-1)<<endl;

	return answer;
}

cube atmozone::bitmap_to_cylmap(cube incube, latangle centerlat, lonangle centerlon, angle pixelangularresolution)
{ 	// Takes cube created from jpg and maps it out as a proper Jcube, 
	// based on user input for center lat/lon, resolution, and longitude convention.
	// Lat/lon extents are determined from the resolution and given center lat/lon, then fed into other version of bitmap_to_cylmap
	// resolution assumes SQAURE pixels
	// added 26.Aug.2016 by smack
	cout <<"\n\n............. bitmap with centerlatlon ...............\n";
	cube answer(incube.N(X),incube.N(Y),incube.N(Z));
	cout <<"N(X) N(Y)="<<answer.N(X)<<" "<<answer.N(Y)<<endl;
	cout <<"Angluar resolution in degrees per pixel is: "<<pixelangularresolution.degrees()<<endl;
	cout <<"centerlon "<<centerlon.center_0().degrees()<<centerlon.center_negative180().degrees()<<endl;
	angle del_lat(incube.N(Y)*pixelangularresolution.degrees(),angle::DEG);
	angle del_lon(incube.N(X)*pixelangularresolution.degrees(),angle::DEG);
			
	latangle top(centerlat.north_negative()-del_lat*(0.5),NORTH_NEGATIVE);
	latangle bottom(centerlat.north_negative()+del_lat*(0.5),NORTH_NEGATIVE);
	lonangle left(centerlon.center_0()-del_lon*(0.5),CENTER_0);
	lonangle right(centerlon.center_0()+del_lon*(0.5),CENTER_0);
	
	
	cout<<"\nlimits after calculation\n"<<endl;
	cout <<"\t\t "<<setprecision(11)<<top.north_negative().degrees()<<endl;
	cout <<"\t"<<setprecision(11)<<left.center_0().degrees()<<"\t\t\t"<<setprecision(11)<<right.center_0().degrees()<<endl;
	cout <<"\t\t"<<setprecision(11)<<bottom.north_negative().degrees()<<endl;
 	
	answer=atmozone::bitmap_to_cylmap(incube,top,bottom,left,right);
	
	
	cout <<"\t\t "<<setprecision(11)<<answer.Axis(Y,0)<<endl;
	cout <<"\t"<<setprecision(11)<<answer.Axis(X,0)<<"\t\t\t"<<setprecision(11)<<answer.Axis(X,answer.N(X)-1)<<endl;
	cout <<"\t\t"<<setprecision(11)<<answer.Axis(Y,answer.N(Y)-1)<<endl;
			
	return answer;
}



cube atmozone::VIMScube_to_albedomap(cube incube, vector<value> inwaves)
{
	cube answer(incube.N(X),incube.N(Y),inwaves.size());
	
	cout <<"answer has "<<answer.N(X)<<" in X, "<< answer.N(Y) <<" in Y, and "<<answer.N(Z)<<" in Z."<<endl;

	// if cube is not in 0 to -360, change to appropriate coordinates	
	for (int x(0);x<incube.N(X);x++) answer.Axis(X,x) = -180.+double(x)*360./double(incube.N(X)-1);
	for (int y(0);y<incube.N(Y);y++) answer.Axis(Y,y) = -90.+double(y)*180./double(incube.N(Y)-1);
	
	
	
	
	// slice out the appropriate windows
	for (int w(0); w<answer.N(Z);w++){
		
		cube newplane(incube), sumcube(incube.N(X),incube.N(Y),1);
			;
		cout <<"Newplane has "<<newplane.N(X)<<" in X, "<< newplane.N(Y) <<" in Y, and "<<newplane.N(Z)<<" in Z."<<endl;
		int counter(0);
		cout <<" ---------------------- "<<endl;
		
		if (inwaves.at(w) >= value(0.9,"um") | inwaves.at(w) <= value(0.94,"um")){		
			for (int a(2); a<=3; a++){
				sumcube+=newplane.plane(Z,a)*3;
				++counter;
			}
			sumcube/=counter;
			sumcube.Axis(Z,0)=0.9;
		}
		if (inwaves.at(w) >= value(1.0,"um") | inwaves.at(w) <= value(1.07,"um")){		
			for (int a(11); a<=12; a++){
				sumcube+=newplane.plane(Z,a)*2.3;
				++counter;
			}
			sumcube/=counter;
			sumcube.Axis(Z,0)=1.06;
		}	
		if (inwaves.at(w) >= value(1.22,"um") | inwaves.at(w) <= value(1.3,"um")){		
			for (int a(38); a<=43; a++){
				sumcube+=newplane.plane(Z,a)*2;
				++counter;
			}
			sumcube/=counter;
			sumcube.Axis(Z,0)=1.3;
		}if (inwaves.at(w) > value(1.58,"um") & inwaves.at(w) <value(1.7,"um")){		
			for (int a(23); a<=24; a++){
				sumcube+=newplane.plane(Z,a)*2;
				++counter;
			}
			sumcube/=counter;
			sumcube.Axis(Z,0)=1.6;
		}			
		if(inwaves.at(w) == value(2,"um") | inwaves.at(w) == value(2.015,"um")){
					
			for (int a(55); a<=72; a++){
				sumcube+=newplane.plane(Z,a)*1.5;
				++counter;
			}
			
			sumcube/=counter;
			sumcube.Axis(Z,0)=2;
			
		}	 
		if(inwaves.at(w) == value(2.8,"um") | inwaves.at(w) == value(2.78,"um")){
			answer.plane(Z,w)=newplane.plane(Z,114,116)/2;
						
			for (int a(114); a<=116; a++){
				sumcube+=newplane.plane(Z,a)*1.3;
				++counter;
			}
			
			sumcube/=counter;
			sumcube.Axis(Z,0)=2.8;
		}
		if(inwaves.at(w) == value(5.,"um")){
						
			for (int a(240); a<=255; a++){
				sumcube+=newplane.plane(Z,a);
				++counter;
			}
			
			sumcube/=counter;
			sumcube.Axis(Z,0)=5.;
		}	 	 
		
		for (int x(0); x<incube.N(X); x++){
			for (int y(0); y<incube.N(Y); y++){
				
				answer(x,y,w)=sumcube(x,y,0);	
				
			}
		}
	
	}
	
	
	return answer;
}

		
double atmozone::ss_albedo_latdeg_londeg_rkm_waveum(double latdeg, double londeg, double radialdistancekm,
		double wavelength)
{
	double answer(0.);
	
	int radialdistance_channel(0);
	double closest_distance(1.e99);
	for (int i(0);i<_ss_albedos.size();i++) {
		if (fabs(_ss_albedos.at(i).first-radialdistancekm) < closest_distance) {
			closest_distance = fabs(_ss_albedos.at(i).first-radialdistancekm);
			radialdistance_channel = i;
		}
	}
	if (osuppress<-6) cout << "Closest altitude to " << radialdistancekm << " is " << _ss_albedos.at(radialdistance_channel).first;
	if (osuppress<-6) cout << " at channel " << radialdistance_channel << "\n";
	
	
	long int wavechannel(_ss_albedos.at(radialdistance_channel).second.Jclosest_dumb(Z,wavelength));
	
	if (osuppress<-6) cout << "ssalbedo cube at time of Jclosest_dumb looks like " << _ss_albedos.at(radialdistance_channel).second;
	answer = _ss_albedos.at(radialdistance_channel).second(londeg, -latdeg, wavechannel, mNN);  // no interp for now
	if (osuppress<-4) cout << "Single scattering albedo here (" << wavelength << ", " << wavechannel << ") determined to be "
		 << answer << "\n";
	return answer;
}


double atmozone::extinction_gas_km_um(double, double wavelength_um)
// not used for now
{
	double answer(0.);
	
	return answer;
}

space_atmozone::space_atmozone(double _d) : atmozone(_d) {}

double space_atmozone::extinction_haze_km_um(double, double)
{
	return 0.;
}

surface_atmozone::surface_atmozone(double _d) : atmozone(_d) {}

surface_atmozone::surface_atmozone(cube _c) : atmozone(_c) {}

double surface_atmozone::extinction_haze_km_um(double, double)
{
	return 0.;
}

DISR_Tomasko_atmozone::DISR_Tomasko_atmozone(vector<pair<double, cube> > inSSA) : atmozone::atmozone(inSSA)
{
}

DISR_Tomasko_atmozone::DISR_Tomasko_atmozone(cube inSSA) : atmozone::atmozone(inSSA)
{
}

double DISR_Tomasko_atmozone::extinction_haze_km_um(double radialdistance_km, double wavelength_um)
{
	double answer(0.);
	
	if (&Space()!=this) {
		
		if ((radialdistance_km<toparentlayer()->bottom_km()-1  ||   
				radialdistance_km>toparentlayer()->top_km()+1)  &&
				!issurface()) {
			cout << "FAIL:  you've asked this zone for a height outside valid range\n";
			
			cout << "this radialdistance = " << radialdistance_km << " km; ";
			cout << "this layer's bottom is " << toparentlayer()->bottom_km() << "km; ";
			cout << "top is " << toparentlayer()->top_km() <<" km\n";
		}
		
		double H(toparentlayer()->scaleheight_km());
		double h0(toparentlayer()->bottom_km());
		double h(radialdistance_km);
		double Tau_0(_haze_extinction_0(0,0,wavelength_um));
		answer = Tau_0*exp(-(h-h0)/H);
		
	}
		
	return answer;
}

DISR_Doose_atmozone::DISR_Doose_atmozone(cube incube) : atmozone::atmozone(incube)
{
}


double DISR_Doose_atmozone::extinction_haze_km_um(double altitude_km, double wavelength_um)
{
	double answer(0.);
	
	if (&Space()!=this) {
		
		if ((altitude_km<toparentlayer()->bottom_km()-1  ||   
				altitude_km>toparentlayer()->top_km()+1)  &&
				!issurface()) {
			cout << "FAIL:  you've asked this zone for a height outside valid range\n";
			
			cout << "this radialdistance = " << altitude_km << " km; ";
			cout << "this layer's bottom is " << toparentlayer()->bottom_km() << "km; ";
			cout << "top is " << toparentlayer()->top_km() <<" km\n";
		}
		
		
		
		if(altitude_km > 55.){
			answer = -_a/(_g + pow(wavelength_um,_f)) * exp(altitude_km/_c) / _c * pow(_b + exp(altitude_km/_c),-2);
		}
		else {
			
			double E_55 = (_a/_c *exp(_d/_c)) / ((_g +wavelength_um*_f)*pow(_b + exp(_d/_c),2));
			answer = -_h*E_55*altitude_km;
		}
	}
	
	cout << "DISR_Doose_atmozone extinction is " << answer << " at " << altitude_km << "km and " << wavelength_um << "um\n";
	
	return answer;
}


exponential_atmozone::exponential_atmozone(cube in_ssalbedo, cube in_extinction_at_bottom) :
		atmozone(in_ssalbedo)
{
	extinction_at_bottom(in_extinction_at_bottom);
}

double exponential_atmozone::extinction_haze_km_um(double radialdistance_km, double wavelength_um)
// This function computes the optical depth per km at a given radial
// distance from the center of the planet.  Presumably you have sent in
// a radialdistance for which this zone is valid.
{
	double answer(0.);
		
	if (&Space()!=this) {
		if ((radialdistance_km<toparentlayer()->bottom_km()-1.  ||   
				radialdistance_km>toparentlayer()->top_km()+1.)  &&
				!issurface())
			cout << "FAIL:  you've asked this zone for a height outside valid range\n";
		double H(toparentlayer()->scaleheight_km());
		double h0(toparentlayer()->bottom_km());
		double h(radialdistance_km);
		
		double extinction_at_bottom(_extinction_at_bottom(0,0,wavelength_um));
		answer = extinction_at_bottom*exp(-(h-h0)/H);
	} //else {cout << "dTau/dh in space is zero\n";}
	
	return answer;
}

orangerind_atmozone::orangerind_atmozone(cube incube) : atmozone(incube)
{
	cout << "In orangerind_atmozone constructor.  At this time, the ssalbedos"
			" look like this " << incube << "\n";
}

orangerind_atmozone::orangerind_atmozone(vector<pair<double, cube> > inSSA) : atmozone(inSSA)
{
}


double orangerind_atmozone::extinction_haze_km_um(double alt, double wavelength_um)
// This function computes the optical depth per km at a given radial
// distance from the center of the planet.  Presumably you have sent in
// a radialdistance for which this zone is valid.
{
	double answer(0.);
		
	long int wavechannel(_extinctions.Jclosest_dumb(Z,wavelength_um));
	answer = _extinctions(0,0, wavechannel);  // no interp for now

	if (osuppress<-2) {
		cout << "Extinction here (" << alt << ") determined to be " << answer;
		cout << " at " << wavelength_um << "um (" << wavechannel << ") microns. \n";
	}

	
	return answer;
}

double atmozone::totaltau(value wavelength)
{
	double answer(0.);
	
	answer = toparentlayer()->depth_km()*extinction_haze_km_um(
			(toparentlayer()->top_km()+toparentlayer()->bottom_km())/2.,
			double(wavelength.convert("um")));
	
	return answer;
}

void atmozone::ss_albedo(cube ssa)
{
	ss_albedo(2575., ssa);
}

void atmozone::ss_albedo(double altitude_km, cube ssa)
{
	_ss_albedos.push_back(pair<double, cube>(altitude_km, ssa));
}

vector<pair<double, cube> > atmozone::ss_albedos()
{
	return _ss_albedos;  // return a copy
}

ostream& operator<<(ostream& out, atmozone& a)
{
	vector<pair<double, cube> > ssa_vec(a.ss_albedos());
	out << "\t \tAtmozone address " << &a << ":\n";
	out << "\t \tAltitudes for which we have SSA cubes:  ";
	for (int i(0);i<ssa_vec.size();i++) out << ssa_vec.at(i).first << "km ";
	out << "\t \tFirst single-Scattering albedo cube characteristics (" << ssa_vec.front().second.N(X) << ", ";
	out << ssa_vec.front().second.N(Y) << ", " << ssa_vec.front().second.N(Z) << ") min=" << ssa_vec.front().second.min();
	out << "; max=" << ssa_vec.front().second.max() << "; mean=" << ssa_vec.front().second.mean();
	out << ", median=" << ssa_vec.front().second.median() << "\n";
	return out << "\t \tPhasefunction address " << a.zonephasefunction() << ", named " << a.zonephasefunction()->name() << "\n";

}

