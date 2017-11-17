#include "Jcrap.h"

#ifndef JPROJECTIONS
#define JPROJECTIONS 1


class projection
{
	public:
		//projection(double slon, double slat) : _subSClon(slon), _subSClat(slat) {} //vestigial
		projection(lonangle slon, latangle slat): _subSClon(slon), _subSClat(slat){}
		projection(latangle slat, lonangle slon): _subSClon(slon),_subSClat(slat){}
		
		virtual pair<double, double> operator()(lonangle, latangle)=0;
		virtual pair<lonangle, latangle> operator()(double, double)=0;
		virtual projtype isprojtype()=0;

		void subSClon(lonangle);
		void subSClat(latangle);
		
		
		
		lonangle subSClon(){return _subSClon;}
		latangle subSClat(){return _subSClat;}
		
		
	protected:
		
		lonangle _subSClon;
		latangle _subSClat;
		~projection(){}
	
		
};


class azimuthalprojection: public projection
{
	public:
		
		azimuthalprojection(cube);
		azimuthalprojection(lonangle, latangle,int, int, int, int, double, angle);
		
		virtual void calcangdist(double)=0;
		virtual projtype isprojtype(){return AZIMUTHAL;}
		virtual double kfactor(lonangle,latangle)=0;
		virtual double kfactor(latangle,lonangle)=0;
		
		pair<double, double> operator()(lonangle, latangle);
		pair<lonangle, latangle> operator()(double, double);
		vector<pair<double, double> > convert_vector_of_angles(vector<pair<lonangle,latangle> >);
		
		void angdist(double);
		void radius(double);
		void xcenter(int);
		void ycenter(int);
		void xsize(int);
		void ysize(int);		
		void Nangle(angle);
		void originatingcyllonconvention(lonconvention);		
		
		double radius(){return _radius;}
		double angdist(){return _angdist;}
		int xcenter(){return _xcenter;}
		int ycenter(){return _ycenter;}
		int xsize(){return _xsize;}
		int ysize(){return _ysize;}
		angle Nangle(){return _Nangle;}
		lonconvention originatingcyllonconvention(){return
				_originatingcyllonconvention;}
		
	protected:
		double _radius;
		double _angdist; //"c" factor
		angle _Nangle;
		double  _xcenter,_ycenter,_xsize,_ysize;
		lonconvention _originatingcyllonconvention;
		~azimuthalprojection(){}
		
	
};


class orthographic: public azimuthalprojection
{
	public: 
		orthographic( lonangle, latangle,int, int, int, int, double, angle);
		orthographic(cube);
		~orthographic(){}
		
		void calcangdist(double);
		projtype isprojtype(){return ORTHO;}	
		double kfactor(lonangle,latangle);
		double kfactor(latangle,lonangle);
				
		using azimuthalprojection::operator();
		
};


class azimuthalstereographic: public azimuthalprojection
{
	public: 
		azimuthalstereographic(lonangle,latangle, int, int, int, int, double, angle);
		azimuthalstereographic(cube);
		void calcangdist(double);
		projtype isprojtype(){return AZSTEREOGR;}	
		double kfactor(lonangle,latangle);
		double kfactor(latangle,lonangle);
				
		using azimuthalprojection::operator();
		
};

class lambertazimuthal: public azimuthalprojection
{
	public: 
		lambertazimuthal(lonangle, latangle, int, int, int, int, double, angle);
		lambertazimuthal(cube);
		void calcangdist(double);
		projtype isprojtype(){return LAMBAZ;}	
		double kfactor(lonangle,latangle);
		double kfactor(latangle,lonangle);
				
		using azimuthalprojection::operator();
		
};

class cylindricalprojection: public projection
{
	public:
		cylindricalprojection(cube);
		pair<double, double> operator()(lonangle,latangle);
		pair<lonangle, latangle> operator()(int, int);
		
		lonconvention cyllonconvention(){return _lonconvention;}
		latconvention cyllatconvention(){return _latconvention;}
		void cyllonconvention(lonconvention);
		void cyllatconvention(latconvention);
		
		projtype isprojtype(){return CYL;}
	protected:
		lonconvention _lonconvention;
		latconvention _latconvention;
};


// Note sure what this is for.... doesn't work with new Angle definition of subSClat and subSClon



// class perspectiveprojection : public projection 
// {
// public:
// 	perspectiveprojection(double p, double slon, double slat) : projection(slon, slat), P(p)  {} // p is distance of the point in sphere radii
// 	
// 
// 	pair<double, double> operator()(float lon, float lat){
// 		pair <double, double> answer; // x, y
// 		
// 		// following http://mathworld.wolfram.com/VerticalPerspectiveProjection.html
// 		double cos_c(sin(subSClat)*sin(lat)+cos(subSClat)*cos(lat)*cos(lon-subSClon));
// 		double k_prime((P-1.)-(P-cos_c));
// 		
// 		answer.first /*X*/ = k_prime * cos(lat)*sin(lon-subSClon);
// 		answer.second/*Y*/ = k_prime *(cos(subSClat)*sin(lat)-sin(subSClat)*cos(lat)*cos(lon-subSClon));
// 		
// 		return answer;
// 	}
// 
// private:
// 	double P;
// };


#endif
