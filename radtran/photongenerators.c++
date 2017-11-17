#include "photongenerators.h"
#include "atmozones.h"

photongenerator::photongenerator(value inlength, long long innumber,
		value indistance, vector<value> inwavelengths) :
		_lengthscale(inlength), _number(innumber), _distance(indistance),
		_wavelengths(inwavelengths)
{
	cube::Jraninit();
}

photongenerator::photongenerator(value inlength, long long innumber,
		value indistance, value inwavelength) :
		_lengthscale(inlength), _number(innumber), _distance(indistance)
{
	cube::Jraninit();
	_wavelengths.resize(0);
	_wavelengths.push_back(inwavelength);
}



vector<value> photongenerator::wavelengths()
{
	return _wavelengths;
}



photongenerator_square::photongenerator_square(value inlength, long long innumber,
		value indistance, vector<value> inwavelengths) : 
		photongenerator(inlength, innumber, indistance, inwavelengths)
{ //all square suns are default set to come in along the x axis in the +x hat direction
	
	
	if (osuppress<0) cout <<"\n----PHOTONGENERATOR_SQUARE CREATED----\n";
	indistance = indistance.convert("km");
	_fieldcenter = geomvector(-double(indistance), 0., 0.);
	
	_makehole=false;
	
	double length_km(inlength.convert("km")); //total length of square radius
	if (osuppress<0) cout <<"length of photon square (km)"<<2*length_km<<endl;
	length_km /= double(innumber); //length in between each photon
	
	_zvector = geomvector(0., 0., length_km);
	_yvector = geomvector(0., length_km, 0.);
	
	_oneside = innumber*2+1;
	
	_photondirection = geomvector(1.,0.,0.);
	if (osuppress<0) cout <<"length between photons (km)\t"<<length_km<<endl;
	if (osuppress<0) cout <<"oneside (# of photon cells along square side)\t"<<_oneside<<endl;
	if (osuppress<0) cout <<"photons per square km\t"<<photons_per_square_km()<<endl;
}

void photongenerator_square::photondirection(geomvector newdirection)
{ // added 2 Sept 2016 by smack in order to be able to move the sun
	_photondirection=newdirection;
}

void photongenerator_square::fieldcenter(geomvector newfieldcenter)
{ // added 2 Sept 2016 by smack in order to be able to move the sun
	_fieldcenter =	newfieldcenter;
}

void photongenerator_square::inc_offaxis(angle inc)
{
	_fieldcenter = geomvector(-double(_distance)*inc.cos(),-double(_distance)*inc.sin(),0.);
	_photondirection = geomvector(1.,1.,0);
	cout <<"\nSun now centered at "<<_fieldcenter.as(XYZ)<<endl;
 
}

photon photongenerator_square::generatephoton(long long id)
{	
	int lambda(id/_oneside/_oneside);
	long long whatsleft(id%(_oneside*_oneside));
	long long y(whatsleft/_oneside-_number);			// row number (truncated)
	long long z(whatsleft%_oneside-_number);        // effectively x, or column number
	
	if (osuppress<0) cout << "Photon #" << id << " computed to be (" << z << ", ";
	if (osuppress<0) cout << y << ", wave=" << lambda << ":" << _wavelengths.at(lambda) << ")\n";
	
	geomvector startinglocation(_fieldcenter+_zvector*double(z)+_yvector*double(y));
	if (osuppress<0) cout << "Fieldcenter = " << _fieldcenter << "\n";
	if (osuppress<0) cout << "zoffset = " << _zvector*double(z) << "\n";
	if (osuppress<0) cout << "yoffset = " << _yvector*double(y) << "\n";
	if (osuppress<0) cout << "starting location = " << startinglocation << "\n";
 	
	double amplitude(1.);
	if (_makehole) if (isinsidehole(startinglocation)) amplitude=0.;
		
	photon answer(startinglocation,_photondirection, amplitude, _wavelengths.at(lambda));
	if (osuppress<0) cout << "starting direction = " << _photondirection << "\n";
	if (osuppress<0) cout << "or " << _photondirection.as(SPHERICAL) << "\n";
	
	answer.set_scattering_zone(&(atmozone::Space()));
	
	return answer;
}

long long photongenerator_square::size()
{
	return _oneside*_oneside*_wavelengths.size();
}

double photongenerator_square::photons_per_square_km()
{
	return 1./(_zvector.z()*_yvector.y());
}


void photongenerator_square::makehole(value length)
{ // set flag to punch a square hole in the photon grid during generatephoton
  // and stores the limits of the hole in a geomvector for ease
 // length= (length in km of z/yside of rectangle hole)/2  
 // added 24 Aug 2016, smack
	
	_makehole=true;
	geomvector answer(_fieldcenter.x(),double(length.convert("km")),double(length.convert("km")));
	_holevector=answer;
}


bool photongenerator_square::isinsidehole(geomvector inloc)
{	// deterines if photon position is inside the photongenerator_square hole
	// added 24 Aug 2016, smack
	bool answer(false);
	
	double zright(_holevector.z()),ytop(_holevector.y());
	
	if (inloc.z()<zright && inloc.z()>-zright){
		if (inloc.y() <ytop && inloc.y()>-ytop){
		//	cout <<"Hole limits are:\n";
		//	cout <<"     \t ----------- \ty ="<<ytop<<endl;
		//	cout <<"z= "<<-zright<<"\t|      0     |\tz="<<zright<<endl;
		//	cout <<"     \t ----------- \ty ="<<-ytop<<endl;
		//	cout <<"Photon position \t"<< inloc <<"\tis within the hole!\n";
			answer=true;
		}
	}
			
	return answer;
}


photongenerator_diffuselamp::photongenerator_diffuselamp(long long innumber, 
		vector<value> inwavelengths, geomvector inposition, geomvector indirection):
		photongenerator(value(0.,"km"), innumber, value(inposition.r(), "km"), inwavelengths),
		_sourcelocation(inposition), _sourcepointing(indirection)
{ 	// for calculating the actual modulation transfer function, this photongenerator
	// emmits photons in a lambertian distribution 
	// smack spring 2017
	
	cout <<"\n---- PHOTONGENERATOR_diffuselamp CREATED ----\n";
	cout <<"centered at "<<inposition.lat()<<" N"<<inposition.lon()<<" E"<<endl;
	
}

photon photongenerator_diffuselamp::generatephoton(long long id)
{	
	// smack 2016 Oct 28
	int lambda=id/_number; 
//	cout << "Photon #" << id << " computed to be";
//	cout <<" wave=" << lambda << ":" << _wavelengths.at(lambda) << ")\n";
	
	 	
	double amplitude(1.);
	
	// passing in both the location and direction in PLANETCENTRIC coordinates
	static phasefunction_lambertian thisphasefunction;
	phasefunction_lambertian* lamb(&thisphasefunction);
	
	photon dummyphoton(_sourcelocation,_sourcepointing,amplitude,_wavelengths.at(lambda));
	// lambertian surface phase function generates photons in lambertian distribution
	// in surface-centered coordinates
	geomvector photondirection_surfacecentric(lamb->randomscatter(dummyphoton));
	
	// translate back to planetocentric coordiates 
	geomvector
			photondirection_planetocentric(photondirection_surfacecentric.transform_from(_sourcelocation));
	photon answer(_sourcelocation,photondirection_planetocentric, amplitude, _wavelengths.at(lambda));
	//cout << "-----------------------\n starting direction = " << answer.direction().as(XYZ) << "\t"<<answer.direction().as(SPHERICAL)<<"\n";
	//cout << "starting position = " << answer.position().as(XYZ) <<"\t"<<answer.position().as(SPHERICAL)<< "\n";
	//cout << "amplitude= "<<answer.amplitude()<<endl<<"-----------------------\n";
	
	//cout <<answer.direction().as(XYZ)<<endl;
	
	//answer.randirection(photondirection_surfacecentric);
	answer.randirection(photondirection_planetocentric);
			
	answer.set_scattering_zone(&(atmozone::Space()));
	
	return answer;
}

long long photongenerator_diffuselamp::size()
{	
	
	return _number*_wavelengths.size();
}


photongenerator_laser::photongenerator_laser(long long innumber, vector<value> inwavelengths, geomvector inposition, 
		geomvector indirection):
		photongenerator(value(0.,"km"), innumber, value(inposition.r(), "km"), inwavelengths),
		_sourcelocation(inposition), _photondirection(indirection)		
{
	cout <<"\n---- PHOTONGENERATOR_LASER CREATED ----\n";
	cout <<"centered at "<<inposition.lat()<<" N"<<inposition.lon()<<" E"<<endl;
	cout <<"with photons traveling in "<<indirection.as(XYZ)<<" direction."<<endl;
}

long long photongenerator_laser::size()
{	
	return _number*_wavelengths.size();
}

photon photongenerator_laser::generatephoton(long long id) 
{
	int lambda=id/_number; 
	double amplitude(1.);
	
	photon answer(_sourcelocation,_photondirection, amplitude, _wavelengths.at(lambda));
	
	answer.set_scattering_zone(&(atmozone::Space()));
	answer.randirection(_photondirection);
	
	return answer;
	
}

photongenerator_multitarget::photongenerator_multitarget(value inside, long long innumber_per_side, 
		value indistance, vector<value> inwavelengths) :
		photongenerator(inside, innumber_per_side, indistance, inwavelengths)
{
}

void photongenerator_multitarget::addfieldcenter(geomvector g)
{
	photongenerator_square newtarget(_lengthscale, _number, _distance, _wavelengths);
	_target_areas.push_back(newtarget);
	_target_areas.back().photondirection(geomvector(-1.,0.,0.));   // hard-coded direction for now
	_target_areas.back().fieldcenter(g);
}

photon photongenerator_multitarget::generatephoton(long long id)
{
	int i(0);
	for (i=0;i<_target_areas.size();i++) {
		if (id >= _target_areas.at(i).size()) id-=_target_areas.at(i).size();
		else break; 
	}
	return _target_areas.at(i).generatephoton(id);
}

long long photongenerator_multitarget::size()
{
	long long answer(0);
	for (int i(0);i<_target_areas.size();i++) answer+=_target_areas.at(i).size();
	return answer;
}

double photongenerator_multitarget::photons_per_square_km()
{
	return _target_areas.front().photons_per_square_km();
}
