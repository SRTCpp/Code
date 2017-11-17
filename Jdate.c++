#include <strstream>
#include <iostream>
#include <iomanip>

#include "Jdate.h"
#include "Jcrap.h"
#include "../../nr/nr.h"


bool operator<(const Jcrap::date& l, const Jcrap::date& r) 
{ return l.JD()<r.JD(); }


date::date(double jd) : _JD_(jd) {}

date::date(string ds) 
{
	cout << "date(string) not yet implimented!\n";
}
 
date::date(int yr, int mo, int da, int hr, int mi, int s)
{
	_JD_ =  (double)julday(yr, mo, da);
	_JD_+= ((double)hr)/24.;
	_JD_+= ((double)mi)/(24.*60.);
	_JD_+= ((double) s)/(24.*3600.);
}

string date::GD() const
{
	string answer;
	double midnightJD=_JD_+0.5;		// nr's JD routines use noon as the start
												// of a JD day, seems like astro uses
												// midnight
	answer.resize(19);
	ostrstream ans(&answer[0], answer.size());

	int yr, mo, da;
	caldat((long)midnightJD, &mo, &da, &yr);	
	ans << setw(4) << yr << "-" << setw(2) << mo << "-" << setw(2) << da;
	int hr, mi;
	double s;
	double residual;
	residual = midnightJD - (double)(int)midnightJD;
	
	hr = (int)(residual*24.);
	ans << "T" << setw(2) << hr << ":";
	residual -= ((double)(hr))/24.;
	
	mi = (int)(residual*24.*60.);
	ans << setw(2) << mi << ":";
	residual -= ((double)(mi))/(24.*60.);
	
	s  = (residual*24.*3600.);
	ans << setw(2) << (int)s;
	
	return answer;
}

string date::MST() const
{
	return ((*this)-7./24.).GD()+string(" MST");
}

string date::CTIO() const
{
	return ((*this)-4./24.).GD()+string(" WST");
}

string date::UTC(double offset) const
{
	return ((*this)+offset/24.).GD()+string(" UTC")+int2str((int)offset);
}

double date::JD() const
{
	return _JD_;
}
