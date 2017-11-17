#include<string>
#include "Jcrap.h"

using namespace std;

class Jcrap::date
{
	public:
			date(double);				// date from JD
			date(string);				// date from a date string
			date(int,int,int,int=0,int=0,int=0);  // yr,mm,dd,hh,mm,ss
			date() {};
			
			friend bool operator<(const date& l, const date& r);
			date operator+=(double d) { _JD_ += d; return *this; }
			date operator+(double d) const { return date(_JD_+d); }
			date operator-(double d) const { return date(_JD_-d); }
			
			string GD() const;
			string MST()const;
			string CTIO()const;
			string UTC(double=0.)const;
			double JD() const;
			
			
	private:
			double _JD_;
};
