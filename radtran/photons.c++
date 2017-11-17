#include "photons.h"
#include "atmolayers.h"

photon::photon(geomvector pos, geomvector direc, double A, value inlambda)
			: _photposition(pos), _photdirection(direc), _amplitude(A)
{ 
	atmoscatterers=0;
	surfscatterers=0;
	
	_lambda_um = double(inlambda.convert("um"));
}	

double photon::generate_random_tau()
{		
	double ran(cube::Jrandom());
	 		
	double thisln = -log(ran);		// "log" is natural log -- "log10" is base 10.
			
	return thisln;
}

geomvector photon::project_km(double d)
{
	double x_projected(position().x() + d*direction().lat().cos()*direction().lon().cos());
	double y_projected(position().y() + d*direction().lat().cos()*direction().lon().sin());
	double z_projected(position().z() + d*direction().lat().sin());
	
	geomvector answer(x_projected, y_projected, z_projected);
	
	return answer;
}

photon photon::scatter(atmozone* to_scattering_zone)
{
	double wave_um(lambda_um());
	double photonlat_deg(position().lat().degrees());
	double photonlon_deg(position().lon().degrees());
	double singlescatteringalbedo(to_scattering_zone->ss_albedo_latdeg_londeg_rkm_waveum(
		photonlat_deg, photonlon_deg, position().r(), wave_um));
	if (osuppress < -4) cout <<to_scattering_zone->layername()<<position().r()<<"\t"<< singlescatteringalbedo<<endl;
	if (osuppress < -4) cout <<"ssa at position: "<< singlescatteringalbedo<<endl;
	
	amplitude(amplitude()*singlescatteringalbedo);  // do the single-scattering albedo FIRST, affecting *this
	
	photon answer(*this);
		
	answer.set_scattering_zone(to_scattering_zone);
	if(osuppress<-3)cout << "in photon::scatter in layer " << to_scattering_zone->toparentlayer()->layername() << "\n";
	if(osuppress<-2)cout << "photon's initial direction before euler = " << direction().as(SPHERICAL) << "\n";//smack 8.apr
	
	phasefunction* ourphasefunc(to_scattering_zone->zonephasefunction());
	
	if (to_scattering_zone->issurface()) {
				 // surface scatter
		// first convert into surface-centered coordinate system with +z as "up", or radially outward from surface
		geomvector direction_in_surfacecoords(direction().transform_to(position()));
		direction(direction_in_surfacecoords);
		
		geomvector scattered_surfacecentric(ourphasefunc->randomscatter(*this));
		
		geomvector scattered_planetocentric(scattered_surfacecentric.transform_from(position()));
		
		answer.direction(scattered_planetocentric);		
	} else {  // atmospheric scatter
	
		// this is the deviation in direction caused by the scatter, in some abstract frame relative to the
		// photon's present direction
		geomvector scatterby(ourphasefunc->randomscatter(*this));
		
		// and then here we rotate that by the direction of the photon to get the final outgoing direction
		answer.direction(geomvector(scatterby.euler_ZXZ(angle(0., angle::DEGREES),
				direction().colatitude(), direction().lon()+angle(90., angle::DEGREES))));
		if(osuppress<-2)cout << "new photon random scatter angle to be eulered is " << scatterby << "\n";
		if(osuppress<-2)cout << "new direction for photon after euler is " << answer.direction().as(XYZ) << "\n"; if(osuppress<-2)cout.flush();
		if(osuppress<-2)cout << "new direction for photon after euler is " << answer.direction().as(SPHERICAL) << "\n"; 
		if(osuppress<-2)cout << "scatter latitude is " << scatterby.latitude() << " or " << scatterby.latitude().as(angle::DEGREES) << "\n";
		if(osuppress<-2)cout.flush();
		
		if (ourphasefunc->keeptrack()) ourphasefunc->add_to_test_distribution(scatterby);

	}
	
	return answer;
}


void photon::ranthetaphi(double theta, double phi)
{
	_ranthetaphi.first=theta;
	_ranthetaphi.second=phi;
}

std::ostream& operator<< (std::ostream& out, photon& p)
{
	out << "Photon &" << &p << ":  amplitude " << p.amplitude() << ". wavel ";
	out << p.lambda() << "\n";
	out << "Position " << p.position() << ", direction " << p.direction();
	
	return out << "\n";
}
