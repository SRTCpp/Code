#include <iostream>
#include <math.h>

#ifndef JANGLE
#define JANGLE 1


enum latconvention {NORTH_POSITIVE,NORTH_NEGATIVE}; //North positive, North negative(jcube default); 2.16 smack
enum lonconvention {CENTER_0,CENTER_NEGATIVE180,CENTER_POSITIVE180}; //center of map at 0 or at -180; 2.16 smack


template <class nTYPE=double>
class Jangle
{
	public:
		enum Jangleunit {RADIANS, DEGREES, RAD=RADIANS, DEG=DEGREES, R=RAD, D=DEG};
		Jangle(nTYPE, Jangleunit);
		Jangle();
	
		nTYPE degrees() const;
		nTYPE radians() const;
		
		Jangle as(Jangleunit) const;
		Jangle& to(Jangleunit);
		
		nTYPE number() const; // returns just the NUMERIC value -- DANGEROUS AVOID IF POSSIBLE
	
		nTYPE sin() const;
		nTYPE cos() const;
		
		Jangle<nTYPE> abs() const;
		
		Jangle  operator+(Jangle);
		Jangle& operator+=(Jangle);
		Jangle  operator-();
		Jangle  operator-(Jangle);
		Jangle& operator-=(Jangle);
		Jangle  operator*(nTYPE);
		Jangle& operator*=(nTYPE);
		
		bool  operator<(Jangle);

		static const nTYPE pi;
			
		template <class vTYPE> friend std::ostream& operator<< (std::ostream&, const Jangle<vTYPE>&);
		template <class vTYPE> friend Jangle angle_atan2(double,double);
		template <class vTYPE> friend Jangle angle_acos(double);
		template <class vTYPE> friend Jangle angle_asin(double);

	private:
		nTYPE _angle;
		Jangleunit _unit;
};

template <class nTYPE>
const nTYPE Jangle<nTYPE>::pi=(3.141592653589793238462643383279502884197);

template <class nTYPE>
Jangle<nTYPE>::Jangle(nTYPE a, Jangle<nTYPE>::Jangleunit u) :
		_angle(a), _unit(u) {}

template <class nTYPE>
Jangle<nTYPE>::Jangle() : _angle(nTYPE(0.)), _unit(Jangle<nTYPE>::RADIANS) {}

template <class nTYPE>
nTYPE Jangle<nTYPE>::degrees() const
{
	nTYPE answer(_angle);
	if (_unit==RADIANS) answer*=180./pi;
	return answer;
}

template <class nTYPE>
nTYPE Jangle<nTYPE>::radians() const
{
	nTYPE answer(_angle);
	if (_unit==DEGREES) answer/=180./pi;
	return answer;
}

template <class vTYPE>
std::ostream& operator<< (std::ostream& out, const Jangle<vTYPE>& a)	
{ 
	out << a._angle;
	if (a._unit==Jangle<vTYPE>::DEGREES) out << char(176);
	return out << "";
}
		
template <class nTYPE>
nTYPE Jangle<nTYPE>::number() const
// returns just the NUMERIC value -- DANGEROUS AVOID IF POSSIBLE
{
	return _angle;
}

template <class nTYPE>
nTYPE Jangle<nTYPE>::sin() const
{
	return ::sin(radians());
}

template <class nTYPE>
nTYPE Jangle<nTYPE>::cos() const
{
	return ::cos(radians());
}

template <class nTYPE>
nTYPE sin(Jangle<nTYPE> a)
{
	return a.sin();
}

template <class nTYPE>
nTYPE cos(Jangle<nTYPE> a)
{
	return a.cos();
}

template <class nTYPE>
Jangle<nTYPE> angle_atan2(nTYPE y,nTYPE x)
{
	Jangle<nTYPE> answer(atan2(y,x), Jangle<nTYPE>::RADIANS);
	return answer;
}

template <class nTYPE>
Jangle<nTYPE> angle_acos(nTYPE a)
{
	Jangle<nTYPE> answer(acos(a), Jangle<nTYPE>::RADIANS);
	return answer;
}

template <class nTYPE>
Jangle<nTYPE> angle_asin(nTYPE a)
{
	Jangle<nTYPE> answer(asin(a), Jangle<nTYPE>::RADIANS);
	return answer;
}

template <class nTYPE>
Jangle<nTYPE> Jangle<nTYPE>::as(Jangle<nTYPE>::Jangleunit u) const
{
	Jangle<nTYPE> answer(*this);
	if (u==Jangle<nTYPE>::DEGREES) answer=Jangle<nTYPE>(degrees(), u);
	else if (u==Jangle<nTYPE>::RADIANS) answer=Jangle<nTYPE>(radians(), u);
	return answer;
}

template <class nTYPE>
Jangle<nTYPE>& Jangle<nTYPE>::to(Jangle<nTYPE>::Jangleunit u)
{
	return (*this)=as(u);
}

template <class nTYPE>
Jangle<nTYPE> Jangle<nTYPE>::operator+(Jangle<nTYPE> rhs)
{
	Jangle<nTYPE> answer(*this);
	return answer+=rhs;
}

template <class nTYPE>
Jangle<nTYPE>& Jangle<nTYPE>::operator+=(Jangle<nTYPE> rhs)
{
	if (rhs._unit == _unit)
		_angle += rhs._angle;
	else 
		(*this)+= rhs.as(_unit);
	return *this;
}

template <class nTYPE>
Jangle<nTYPE> Jangle<nTYPE>::operator-()
{
	_angle *= nTYPE(-1.);
	return *this;
}

template <class nTYPE>
Jangle<nTYPE> Jangle<nTYPE>::operator-(Jangle<nTYPE> rhs)
{
	return (*this)+(-rhs);
}

template <class nTYPE>
Jangle<nTYPE>& Jangle<nTYPE>::operator-=(Jangle<nTYPE> rhs)
{
	return (*this)+=(-rhs);
}

template <class nTYPE>
Jangle<nTYPE> Jangle<nTYPE>::operator*(nTYPE m)
{
	Jangle<nTYPE> answer(*this);
	return answer*=m;
}

template <class nTYPE>
Jangle<nTYPE>& Jangle<nTYPE>::operator*=(nTYPE m)
{
	_angle*=m;
	return *this;
}

template <class nTYPE>
bool Jangle<nTYPE>::operator<(Jangle<nTYPE> r)
{
	return degrees() < r.degrees();
}

template <class nTYPE>
Jangle<nTYPE> Jangle<nTYPE>::abs() const
{
	Jangle<nTYPE> answer(abs(degrees()), Jangle<nTYPE>::DEG);
	return answer;
}

typedef Jangle<double> angle;


class latangle: public angle
{
	public: 
		latangle(angle , latconvention );
		latangle(const latangle& );

		angle north_positive() const;
		angle north_negative() const;
	
		latangle as(latconvention) const;
		latangle& to(latconvention);
		
	private:
		angle _lat;
		latconvention _convention;
};



class lonangle: public angle
{
	public: 
		lonangle(angle , lonconvention);
		lonangle(const lonangle& );
				
		angle center_0() const;
		angle center_negative180() const;
		angle center_positive180() const;
		
		lonangle as(lonconvention) const;
		lonangle& to(lonconvention);
		lonconvention isstoredas() {return _convention;}
		
	private:
		angle _lon;
		lonconvention _convention;
};





#endif
