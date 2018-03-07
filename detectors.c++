#include "omp.h"

#include "SRTC++.h"
#include "atmozones.h"

detector::detector(geomvector inpos) : _position(inpos) 
{	
	_boresight_direction = -inpos;  // look nadir by default
}

void detector::detect(list<photon> scatter_history)
// While any detector may wish to know about the entire scatter history,
// it is only at any particular calling of detect() registering the result
// of the LAST event in the history, or the most recent scattering.
// this structure in the parent class detector ensures this:  daughter classes
// are only REQUIRED to  implement detect(photon), NOT detect(list<photon>), and 
// hence they don't really need to understand this subtlety
{
	detect(scatter_history.back());
}

void detector::boresight(geomvector newboresight)
{
	_boresight_direction = newboresight;
}

geomvector detector::position()
{
	return _position;
}

geomvector detector::boresight()
{
	return _boresight_direction;
}

colorCCD::colorCCD(geomvector inloc, int npix, angle FOVpix, vector<value> colorvector) : detector(inloc)
{ //npix is the number of pixels along each axis, NOT the total number of pixels
	_iFOV=FOVpix;

	if (osuppress<-4) cout <<"\ncreating a colorCCD with ";
	
	if (npix%2==0) npix++;
	if (osuppress<-4) {
		cout <<npix<<" pixels of IFOV equivalent to square radius on the surface of ";
		cout <<_iFOV*inloc.r() <<endl;
	}
	_CCD = cube(npix, npix, colorvector.size());
	int mid((npix-1)/2);
	if (osuppress<-4) cout << "km_per_pixel() \t"<<km_per_pixel()<<endl;
	_CCD.keyword("km_per_pixel",double2str(km_per_pixel()));
	
	for (int i(0);i<=mid;i++) {
		_CCD.Axis(X,mid+i) = km_per_pixel()*double(i);
		_CCD.Axis(X,mid-i) = -km_per_pixel()*double(i);
		_CCD.Axis(Y,mid+i) = km_per_pixel()*double(i);
		_CCD.Axis(Y,mid-i) = -km_per_pixel()*double(i);
	}
	for (int z(0);z<_CCD.N(Z);z++) _CCD.Axis(Z,z)=double(colorvector.at(z));
	
	
		
// calculate and output geometry
	_geo = cube(_CCD.N(X), _CCD.N(Y), 8);
	_geo.copyaxis(_CCD);
	
/* _geo is in VIMS backplane format as:

    BAND_BIN_FILTER_NAME = (0"latitude",1"longitude",2"sample_resolution",
      3"line_resolution",4"phase_angle",5"incidence_angle",6"emission_angle",
      7"north_azimuth")  */
	
	
	for (int px(0);px<_geo.N(X);px++) {
		double x(_CCD.Axis(X,px));
		for (int py(0);py<_geo.N(Y);py++) {
			double y(_CCD.Axis(Y,py));
			
			double rho(sqrt(x*x+y*y));
			angle c(angle_asin(rho/2575.));   // hardcoded Titan radius, obv.
			
			angle lat_0(-position().latitude());
			angle lon_0(position().longitude());
			
			angle lat(angle_asin(c.cos()*lat_0.sin() + y/rho*c.sin()*lat_0.cos()));
			if (rho == 0.) lat = angle(angle_asin(c.cos()*lat_0.sin()));   // why was this commented out I wonder
			angle lon(lon_0 + angle_atan2(x*c.sin(),(rho*c.cos()*lat_0.cos()-y*c.sin()*lat_0.sin())));
			
			_geo(px,py,0) =-lat.degrees();
			_geo(px,py,1) = lon.degrees();
			
			geomvector surfacepoint(geomvector::geolatitude(-lat), geomvector::geolongitude(lon), 1.0);
			geomvector Sun(1., 0., 0.);  // shit; hard-coded Sun direction
			
			_geo(px,py,4) = position().unitvector().angulardistancel(Sun).degrees();		// phase
			_geo(px,py,5) = surfacepoint.angulardistancel(Sun).degrees();  	// incidence
			_geo(px,py,6) = surfacepoint.angulardistancel(position().unitvector()).degrees(); // emission
		}
	}

}

colorCCD::colorCCD(geomvector inloc, value a, value b, vector<value> colorvector) : detector(inloc)
{
	if(osuppress<-4) cout <<"\ncreating a colorCCD with ";
	a=a.convert("km");
	b=b.convert("km");
	int pixels(0);
	if (a>b) { //a=, b= square radius of image
		pixels = int(a/b)*2+1;
		_iFOV=angle(double(b)/_position.r(), angle::RAD);
	} else { //a= square radius of image, b=
		pixels = int(b/a)*2+1;
		_iFOV=angle(double(a)/_position.r(), angle::RAD);
	}
	if (osuppress<-4) cout <<pixels<<" pixels\n";
		
	*this = colorCCD(inloc, pixels, _iFOV, colorvector);
}


void colorCCD::write(string filename)
{
	_CCD.write(filename+string(".colorCCD.Jcube"));
	
	geo().write(filename+string("_geo.colorCCD.Jcube"));
}



cube& detector::CCD()
{
	return _CCD;
}

cube& detector::geo()
{
	return _geo;
}

detector& detector::operator/=(double d)
{
	CCD()/=d;
	return *this;
} 

void detector::keyword(string key, string valuestring)
{
	CCD().keyword(key, valuestring);
} 

double colorCCD::km_per_pixel()
{
	return _iFOV.radians()*position().r();
}

double colorCCD::pixelarea_in_square_km()
{
	double side(km_per_pixel());
	return side*side;
}

angle colorCCD::iFOV()
{
	return _iFOV;
}

void colorCCD::detect(photon inphoton)
{	
	geomvector detector_to_photon(inphoton.position()-this->position());
	if (inphoton.amplitude()<0) cout << "\nnegative inphoton!  " << inphoton.amplitude() <<"\n";
	
	// now rotate that d_2_p direction such that the boresight becomes
	// along the x-axis.  Then you can do an easy lookup.

	geomvector rotated;
	if (osuppress<-4) cout << "\nDtoP is initially " << boresight().as(SPHERICALPOLAR) << " or " << boresight() <<"\n";
	rotated = detector_to_photon.euler_ZXZ(
			-(this->boresight().lon()-angle(90.,angle::DEG)),
			-(this->boresight().lat()),
			angle(-90.,angle::DEG));
	
	
	long x(_CCD.Jclosest(X,-rotated.y()));
	long y(_CCD.Jclosest(Y,-rotated.z()));
	long w(_CCD.Jclosest_dumb(Z,inphoton.lambda()));
	
	static bool output(0);
	if (osuppress<-4 || output) cout << "\nDtoP is initially " << boresight().as(SPHERICALPOLAR) << " or " << boresight().as(XYZ) <<"\n";
	if (osuppress<-4 || output) cout << "Rotated now to " << rotated << "\n";
	if (osuppress<-4 || output) cout << "so I'm thinking that the detector-X is " << -rotated.y();
	if (osuppress<-4 || output) cout << " and the detector-y is " << rotated.z() << "\n";
	if (osuppress<-4 || output) cout << "Those correspond to pixels " << x;
	if (osuppress<-4 || output) cout << " and " << y << " at " << w << " (" << inphoton.lambda() << ")\n";
	


	
	float* pixel(&_CCD(x,y,w));
	float newcharge(float(inphoton.amplitude()));
	if (osuppress<-4 || output) cout << "Detector now adding newcharge " << newcharge << "\n";
	
#pragma omp atomic
	*pixel += newcharge;
}

elephant::elephant(colorCCD inCCD, int maxscatterdepth) : colorCCD(inCCD)
{
	_CCD_vv.resize(maxscatterdepth+1);
	for (int i(0);i<=maxscatterdepth;i++) {
		_CCD_vv.at(i) = vector<colorCCD> (pow(2,i),inCCD);
	}
}

void elephant::write(string prefix)
{
	cout << "Writing an elephant\n";
	for (unsigned long i(0);i<_CCD_vv.size();i++) {
		for (unsigned long j(0);j<_CCD_vv.at(i).size();j++) {
			string outfilename(prefix);
			for (unsigned long k(0);k<_CCD_vv.size()-i;k++) outfilename+=string("_");
			bitset<64> thishistory(j);
			outfilename +=thishistory.to_string().substr(64-i,i);
			outfilename += string(".Jcube");
			_CCD_vv.at(i).at(j).write(outfilename);
		}
	}
	
}


void elephant::detect(list<photon> history)
{
	static unsigned long count01(0), count10(0);
	
	if (SRTC::ompdebug) cout << "Thread " << omp_get_thread_num() << " heading into elephant::detect\n"; cout.flush();
	unsigned long numscatters(history.size()-1);
	if (numscatters > _CCD_vv.size()-1) numscatters = _CCD_vv.size()-1;
	
	unsigned long binary_history_code(0);
	long n(0);
	for (list<photon>::iterator i(history.begin());i!=history.end();i++,n++) {
		if (i==history.begin() || n < long(history.size())-(long(_CCD_vv.size())-1)) continue;
		binary_history_code <<= 1;
		binary_history_code += !(i->scattering_zone()->issurface());
	}
	
	if (SRTC::ompdebug) cout << "Thread " << omp_get_thread_num() << " heading into colorCCD::detect\n"; cout.flush();
	_CCD_vv.at(numscatters).at(binary_history_code).detect(history.back());
	if (SRTC::ompdebug) cout << "Thread " << omp_get_thread_num() << " done with colorCCD::detect\n"; cout.flush();
}

elephant& elephant::operator/=(double d)
{
	for (int x(0);x<_CCD_vv.size();x++)
		for (int y(0);y<_CCD_vv.at(x).size();y++)
			_CCD_vv.at(x).at(y) /= d;
	return *this;
}

void elephant::keyword(string key, string valuestring)
{
	for (int x(0);x<_CCD_vv.size();x++)
		for (int y(0);y<_CCD_vv.at(x).size();y++)
			_CCD_vv.at(x).at(y).keyword(key, valuestring);
}

cube& elephant::CCD(unsigned int n_scatters, bitset<64> scatterhistory)
{	
	cout <<"\t\treturning CCD at "<<n_scatters<<" at "<<scatterhistory.to_ulong()<<endl;
	return _CCD_vv.at(n_scatters).at(scatterhistory.to_ulong()).CCD();
}


void elephant::calculate_components()
{	// added smack 2017 DEC 12
	// n = number of scattering events
	// h = history number 
	cube directanswer(_CCD_vv.at(0).at(0).CCD(),0.); 
	cube bluranswer(directanswer), addanswer(directanswer);

	for (int n(0); n<_CCD_vv.size(); n++){
		for (int h(0); h<_CCD_vv.at(n).size(); h++){
			bitset<64> htobit(h);
			if ( htobit[64] == 0 ) {
				directanswer+=_CCD_vv.at(n).at(h).CCD();
			}
			else if (htobit.count()==n)  {
				addanswer+=_CCD_vv.at(n).at(h).CCD();
			}
			else if (htobit[64] == 1) {
				bluranswer+=_CCD_vv.at(n).at(h).CCD();
			}
		}		
	}
	_directCCD= directanswer;
	_blurredCCD= bluranswer;
	_additiveCCD= addanswer;
}

cube& elephant::directCCD()
{	// smack 2017 DEC 11
	return _directCCD;
}

cube& elephant::blurredCCD()
{	// smack 2017 DEC 11
	return _blurredCCD;
}

cube& elephant::additiveCCD()
{	// smack 2017 DEC 11
	return _additiveCCD;
}
