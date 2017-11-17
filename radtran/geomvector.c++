#include "geomvector.h"

geomvector::geomvector(double inx, double iny, double inz) :
	_x(inx), _y(iny), _z(inz), _system(CARTESIAN)
{
}


geomvector::geomvector(geomvector::geolatitude inlat, geomvector::geolongitude inlon, double inr) :
	_lat(inlat), _lon(inlon), _r(inr), _system(SPHERICAL)
{
}

geomvector::geomvector(geomvector::geolongitude inlon, geomvector::geolatitude inlat, double inr) :
	_lat(inlat), _lon(inlon), _r(inr), _system(SPHERICAL)
{
}

geomvector::geomvector(double inr, geomvector::geolatitude inlat, geomvector::geolongitude inlon) :
	_lat(inlat), _lon(inlon), _r(inr), _system(SPHERICAL)
{
}

geomvector::geomvector(double inr, geomvector::geolongitude inlon, geomvector::geolatitude inlat) :
	_lat(inlat), _lon(inlon), _r(inr), _system(SPHERICAL)
{
}

double geomvector::x() const
{
	double answer(0.);
	
	if (_system==CARTESIAN) answer=_x;
	else answer = _r * _lat.cos() * _lon.cos();
	
	return answer;
}

double geomvector::y() const
{
	double answer(0.);
	
	if (_system==CARTESIAN) answer=_y;
	else answer = _r * _lat.cos() * _lon.sin();
	
	return answer;
}

double geomvector::z() const
{
	double answer(0.);
	
	if (_system==CARTESIAN) answer=_z;
	else answer = _r * _lat.sin();
	
	return answer;
}

angle geomvector::lat() const
{
	angle answer(0., angle::DEGREES);
	
	if (_system==SPHERICAL) answer=_lat;
	else answer = angle_asin(_z/r());
	
	return answer; 
}

angle geomvector::latitude() const
{
	return lat();
}

angle geomvector::colatitude() const
{
	angle answer(angle(90., angle::DEGREES) - lat());
	return answer;
}

angle geomvector::lon() const
{
	angle answer(0., angle::DEGREES);
	
	if (_system==SPHERICAL) answer=_lon;
	else answer = angle_atan2(_y,_x);
	
	return answer;
}

angle geomvector::longitude() const
{
	return lon();
}

double geomvector::r() const
{
	double answer(0.);
	
	if (_system==SPHERICAL) answer=_r;
	else {
		answer = sqrt(_x*_x+_y*_y+_z*_z);
	}
	
	return answer;
}

geomvector geomvector::unitvector() const
{
	return geomvector(1., geomvector::geolatitude(lat()), geomvector::geolongitude(lon()));
}

geomvector geomvector::as(coordinatesystem newsys) const
{
	geomvector answer(*this);
			
	return answer.to(newsys);
}

geomvector geomvector::to(coordinatesystem newsys)
{
	if (_system != newsys) {
		if (newsys == SPHERICAL) {
			_r = r();
			_lat = lat();
			_lon = lon();
			_system = SPHERICAL;
		}
		if (newsys == CARTESIAN) {
			_x = x();
			_y = y();
			_z = z();
			_system = CARTESIAN;
		}
	}
	
	return *this;
}

geomvector geomvector::transform_to(const geomvector surfpos) const
{
	return euler_ZXZ(angle(-90., angle::DEGREES)-surfpos.lon(),
			-surfpos.colatitude(), angle(0., angle::RAD));
}

geomvector geomvector::transform_from(const geomvector surfpos) const
{
	return euler_ZXZ(angle(0., angle::RAD),
			surfpos.colatitude(),angle(90., angle::DEGREES)+surfpos.lon());
}

angle geomvector::angulardistance(const geomvector b) const
{
	angle answer(angle_acos(unitvector().dot(b.unitvector())));
	if (!(answer.degrees()<=0.) && !(answer.degrees()>=0.)) 
		answer = angle(0., angle::DEG);
	return answer;
}

angle geomvector::angulardistancel(const geomvector b) const
{
	long double dotproduct(0.);
	dotproduct += (long double)(x()*b.x());  // doing this by hand in long doubles to try to preserve precision
	dotproduct += (long double)(y()*b.y());
	dotproduct += (long double)(z()*b.z());
	angle answer(angle(double(acosl(dotproduct)), angle::RADIANS));
	if (!(answer.degrees()<=0.) && !(answer.degrees()>=0.)) 
		answer = angle(0., angle::DEG);
	return answer;
}

geomvector geomvector::operator-() const
{
	return (*this)*-1.;
}

geomvector geomvector::operator*(double rhs) const
{
	geomvector answer(*this);

	if (_system==SPHERICAL)
		answer = geomvector(r()*rhs, geomvector::geolatitude(lat()), geomvector::geolongitude(lon()));
	else if (_system==CARTESIAN)
		answer = geomvector(x()*rhs, y()*rhs, z()*rhs);

	return answer;
}

geomvector geomvector::operator/(double rhs) const
{
	geomvector answer(*this);
	
	rhs = 1./rhs;

	return answer*rhs;
}

geomvector geomvector::operator+(geomvector rhs) const
{
	geomvector answer(x()+rhs.x(), y()+rhs.y(), z()+rhs.z());

	return answer;
}

geomvector geomvector::operator-(geomvector rhs) const
{
	geomvector answer((*this)+(-rhs));

	return answer;
}

geomvector geomvector::euler_ZXZ(angle z1, angle x2, angle z3) const
{	
	double xx(0.);
	xx = this->x()*(z3.cos()*z1.cos()-x2.cos()*z3.sin()*z1.sin());
	xx-= this->y()*(z3.cos()*z1.sin()+x2.cos()*z1.cos()*z3.sin());
	xx+= this->z()*z3.sin()*x2.sin();
	
	double yy(0.);
	yy = this->x()*(z1.cos()*z3.sin()+z3.cos()*x2.cos()*z1.sin());
	yy+= this->y()*(z3.cos()*x2.cos()*z1.cos()-z3.sin()*z1.sin());
	yy-= this->z()*z3.cos()*x2.sin();
	
	double zz(0.);
	zz = this->x()*x2.sin()*z1.sin();
	zz+= this->y()*z1.cos()*x2.sin();
	zz+= this->z()*x2.cos();
	
	geomvector answer(xx,yy,zz);
	
	return answer;
}

double geomvector::dot(const geomvector rhs) const
{
	double answer(0.);
	answer += x()*rhs.x();
	answer += y()*rhs.y();
	answer += z()*rhs.z();
	return answer;
}

geomvector geomvector::euler(angle z1, angle x2, angle z3) const
{
	return euler_ZXZ(z1,x2,z3);
}

geomvector::geolatitude geomvector::geocolatitude(angle a) 
{
	angle answer(angle(angle::pi/2., angle::RAD)-a);
	return geolatitude(answer);
}

std::ostream& operator<< (std::ostream& out, const geomvector& g)
{
	if (g._system==CARTESIAN)
		out << "(" << g.x() << ", " << g.y() << ", " << g.z() << ")";
	else
		out << "(r=" << g.r() << ", theta=" << g.lat() << ", phi=" << g.lon() << ")"; 
	
	return out;
}
