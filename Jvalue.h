#include <iostream>
#include <math.h>
#include "Junit.h"

string float2str(float, int);

using namespace std;

template <class nTYPE=double> 
class Value {
	public:
		Value(nTYPE, Jcrap::unit/* =" "*/ );
		Value(nTYPE&);
		Value(const Value&);
		Value();
		
		Value operator*(Value);
		Value operator*(nTYPE);
		Value operator/(Value);
		Value operator/(nTYPE);
		Value operator+(Value);
		Value operator-(Value);
		Value operator-(nTYPE);
		Value operator*=(Value);
		Value operator*=(nTYPE);
		Value operator=(nTYPE);
		bool operator>(Value);
		operator nTYPE() const { return number; }
		Jcrap::unit Units() const {return units; }
//		Value pow(Value);
		Value pow(nTYPE);
		
		Value simplify(int);
		Value simplify(Jcrap::unit);
		Value simplify();
		Value convert(Jcrap::unit);
		
		string tostring();

		template <class vTYPE> friend ostream& operator<< (ostream&, const Value<vTYPE>&);
	private:
		nTYPE number;
		Jcrap::unit units;
};


template <class nTYPE>
Value<nTYPE>::Value(nTYPE d, Jcrap::unit u): number(d), units(u) {}
/* 3/7/2001 JB */

template <class nTYPE>
Value<nTYPE>::Value(nTYPE& d): number(d), units("") {}
/* 2017 November 2 JWB */

template <class nTYPE>
Value<nTYPE>::Value(const Value& r)
{
	number = r.number;
//	cout << "Setting number = " << number << "\n";	
	
	units = r.units;
//	cout << "Setting units = \"" << units << "\"\n";
	
//	cout << "new Value = " << *this << "\n";
}

template <class nTYPE>
Value<nTYPE>::Value(): number(0) {units.resize(0);}


template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator*(Value<nTYPE> r)
/* 3/7/2001 JB */
{
	Value<nTYPE> answer;
	answer.number = number*r.number;
	answer.units = units*r.units;
	answer *= nTYPE(answer.units);
	answer.units.removedimensionless();
	return answer;
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator*(nTYPE r)
/* 5/6/2001 JB */
{
	Value<nTYPE> rvalue(r, "");
	return (*this)*rvalue;
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator/(Value<nTYPE> r) 
/* 3/29/2001 JB */
{	
	Value<nTYPE> answer;
	answer.number = number/r.number;
	answer.units = units/r.units;
	answer.number /= nTYPE(answer.units);
	answer.units.removedimensionless();
	return answer;
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator/(nTYPE r) 
/* 3/30/2001 JB */
{
	Value<nTYPE> rvalue(r, "");
	return (*this)/rvalue;
}
		
template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator+(Value<nTYPE> r) 
/* 3/29/2001 JB */
{	
	Value<nTYPE> answer;
	r.simplify(0);
	if (units.diff(r.units).order() != 0) r.convert(units);
	if (units.diff(r.units).order() != 0){		
		std::cout << "ADDITION ERROR in Value operator+.\n";
		cout << "differing unit types in both sides, " << *this << " + "<<r<<"\n";
	}		
//	cout << " unit types in both sides: " << *this << " + "<<r<<"\n";

	answer.number = number+r.number;
	answer.units = units;
	return answer;
}
		
template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator-(Value<nTYPE> r) 
/* 3/30/2001 JB */
{	
	Value<nTYPE> answer;
	r.simplify(0);
	if (units.diff(r.units).order() != 0) r.convert(units);
	if (units.diff(r.units).order() != 0){		
		cout << "ADDITION ERROR in Value operator-.\n";
		cout << "differing unit types in both sides, " << *this << " - "<<r<<"\n";
	}		

	answer.number = number-r.number;
	answer.units = units;
	return answer;
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator-(nTYPE r)
/* 5/6/2001 JB */
{
	return (*this)-Value<nTYPE>(r);
}


template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator=(nTYPE r) { 	
	number=r; units=Jcrap::unit(" ");
	return (*this);
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator*=(Value<nTYPE> r)
/* 3/7/2001 JB 
	Future work:  Make this more efficient.
*/
{
	number *= r.number;
	units  *= r.units;
//	cout << "Value*=Value -- unitmult = \"" << units << "\"\n";
	number *= nTYPE(units);
//	cout << "before dimensionremove: \"" << *this << "\"\n";
	units.removedimensionless();
//	cout << "dimensionlessremoved: \"" << *this << "\"\n";
//	cout << "units.size() = " << units.size() << "\n";
	return (*this);
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::operator*=(nTYPE r)
/* 3/7/2001 JB 
	Future work:  Make this more efficient.
*/
{
	number *= r;
	units.removedimensionless();
	return (*this);
}

template <class nTYPE>
bool Value<nTYPE>::operator> (Value<nTYPE> r)
/* 5/18/2001 JB */
{
	if (units.diff(r.units).order() != 0) r.convert(units);
	if (units.diff(r.units).order() != 0){		
		cout << "COMPARISON ERROR in Value operator>.\n";
		cout << "differing unit types in both sides, " << *this << " - "<<r<<"\n";
		return 0;
	}
	else return ((nTYPE)(*this) > (nTYPE)r);
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::pow(nTYPE e)
/* 3/30/2001 JB */
{
	Value<nTYPE> answer;
	answer.number = ::pow(number, e);
	answer.units = units;
	answer.units.pow(e);
	return answer;
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::simplify(int level)
/* 3/20/2001 JB
	Calls the string version */
{	
	units.simplify(level);
	number*=(nTYPE)units;
	units.removedimensionless();
	return *this;
}
	
template <class nTYPE>
Value<nTYPE> Value<nTYPE>::simplify(Jcrap::unit u)
/* 3/20/2001 JB
	Calls unit.simplify and then takes care of dimensionless multiplicative
	differences */
{
	units.simplify(u);
	number*=(nTYPE)units;
	units.removedimensionless();
	return (*this);
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::simplify()
/* 3/23/2001 JB */
{
	units.simplify();
	number*=(nTYPE)units;
	units.removedimensionless();
	return *this;
}

template <class nTYPE>
Value<nTYPE> Value<nTYPE>::convert(Jcrap::unit aim)
/* 3/23/2001 JB
	Converts to the units in parenthesis as well as it can. */
{
	
	if (units == Jcrap::unit("F")) {
		cout << "The units are just Fahrenheit!  \n";
		number += 459.;
		units = Jcrap::unit("R");
	}	
	if (units == Jcrap::unit("C")) {
		cout << "The units are just Celsius!  \n";
		number += 273.15;
		units = Jcrap::unit("K");
	}
	
	units.convert(aim);
	number*=nTYPE(units);	
	
	if (aim == Jcrap::unit("F")) {
		number -= 459.;
		units = Jcrap::unit("F");
	}	
	if (aim == Jcrap::unit("C")) {
		number -= 273.15;
		units = Jcrap::unit("C");
	}
	
	units.removedimensionless();
	return *this;
}


template <class vTYPE>
ostream& operator<< (ostream& out, const Value<vTYPE>& v)	
{ 
	out << v.number << " ";
	return out << v.Units().tostring();
}


template <class nTYPE>
string Value<nTYPE>::tostring()
{
	string sans(float2str(number, 10) + Units().tostring());
	return sans;
}

typedef Value<double> value;
