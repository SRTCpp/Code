#include "Jcrap.h"
#include <stdlib.h>
#include <fstream>

double str2dbl(string s)
{
	double d; char cs[4000];
	sprintf(cs, "%s", s.c_str());
	istrstream ics(cs);
	ics >> d;
	return d;
}

using namespace Jcrap;

string readuntil(istream& in, string stoppers)
/* created 2/28/2001 JB
	Reads from in until one of the characters in stoppers is reached.  Throws away
	characters from in that match stoppers before the string part. */
{
	string answer; char c;			// result and temporary character
	
	while (in.get(c)  &&  stoppers.find(c)!=string::npos);  // disallow null inputs
	answer+=c;
	while (in.get(c)  &&  stoppers.find(c)==string::npos) answer+=c;
	
	return answer;
}


unit unitelement::unitdatabase;  // initialize static class variable

int unit::defaultsimplifyint(0);


unit::unit(char const* ustr)
/* created 2/24/2001 JB
	Constructor for Junit that takes a string of form "unit^exponent unit^exponent ..."
	
FUTURE WORK:
	allow for W / m^2 um   instead of W m^-2 um^-1*/
		
/* Dammit this unit(char*) unitelement(string) thing is going to kill me yet.
	unit("") is for a null list
	unit(" ") is for a list with only the dimensionless unitelement. */
{
	istrstream instrustr(ustr);	// convert from char* to istrstream
	string indiv;						// where a single unit element string is stored

	if (ustr[0]==' ')
		push_back(unitelement(" "));  // the dimensionless one 
	else while (instrustr >> indiv)
		push_back(unitelement(indiv));
}



unit::unit(double d)
/* created 3/17/2001 JB
	This just allows for a single, unitless with a proper mult */
{	push_back(unitelement(d)); }

unit::unit(unitelement u)
/* created 3/23/2001 JB */
{	push_back(u); }

unitelement::unitelement() { };

unitelement::unitelement(double d)
/* created 3/17/2001 JB
	This creates the dimensionless unitelement with mult=d */
{
	(*this)=unitelement(" ");
	mult=d;
}
	

unitelement::unitelement(string u)
/* created 2/24/2001 JB
	Takes single unit^exponent and creates a Junitelement. 
	Modified to be unitelement creator, 3/7/2001*/
{	
// establish database of units if not already done
	if (unitdatabase.size()==0) readunitfile();
	
// the dimensionless cases
	if (u==" ") { (*this)=(*unitdatabase.begin()); return; }
	if (isdigit(u[0])) { (*this)=unitelement(str2dbl(u)); return; }
	
	
// read in the abbreviation we're looking for and the exponent from the string
	int exploc;	
	float tempexp=1.;		// so the exponent doesn't get stepped on by the new unit
	if ((exploc=u.find("^"))!=u.npos) {
		abbr=u.substr(0,exploc);
		istrstream e(u.substr(exploc+1).c_str());
		e >> tempexp;
	}
	else abbr=u;

// look through the database for this unit, and assign the answer from the database match
	for (list<unitelement>::iterator i=unitdatabase.begin();i!=unitdatabase.end();i++)
		if (abbr==(*i).abbr) { (*this)=*i; break; }
	exp=tempexp;
	subunits.pow(exp);
}

void unitelement::readunitfile()
/* Created 2/24/2001 JB
	Reads in database from "Junits.data".
	put into unitelement class instead of unit 3/7/2001
	
FUTURE WORK:
	Allow for varying locations of Junits.data */
{
	ifstream ufile;// input file stream
	ufile.open("Junits.data");
	if (!ufile.is_open()) {
		ufile.clear();
		ufile.open("/r2d2ZFS/jbarnes/.Junits.data");
		if (!ufile.is_open()) {	
			ufile.clear();
			ufile.open("/r2d2ZFS/jbarnes/.Junits.data");	
			if (!ufile.is_open()){
				cout << "\n\n Can't find Junits data, put file in ~jbarnes/.Junits.data\n";
				exit (1);
			}
		}
	}
	unitelement uel;					// storage for each element before it gets added to database

// Add one hardcoded unit, the null unit or unitless
	unitelement null;
	null.abbr = "";
	null.name = "";
	null.type = "dimensionless";
	null.mult = 1.0;
	null.exp  = 1.0;
	null.specialtype = unit::dimensionless;
	unitdatabase.push_back(null);
			
// until the EOF, read in the units for the database from ufile
	while (ufile >> uel.abbr) {
		uel.exp  = 1.0;
		uel.name = readuntil(ufile, "\t\n");

	// read in the type and set the specialtype enum
		uel.type = readuntil(ufile, "\t\n");
		if (uel.type == "conglomerate") uel.specialtype=unit::conglomerate;
		else if (uel.type == "system") uel.specialtype=unit::system;
		else uel.specialtype=unit::base;
	
	// do the conglomerate recursion if necessary
		if (uel.specialtype==unit::conglomerate    ||    uel.specialtype == unit::system) {
			uel.subunits=unit(readuntil(ufile, "\n").c_str());
			uel.mult=1.0;
		} 
		else ufile >> uel.mult;
		
		
	// add it to the database
		unitdatabase.push_back(uel);
	}
	ufile.close();
	
//	cout << "Database read in: \n\n";
//	for (unit::iterator i(unitdatabase.begin());i!=unitdatabase.end();i++)
//		cout << *i << ", mult = " << i->mult << "\n";
}


unit unit::pow(double exponent)
/* 3/21/2001 JB
	This takes a single unit to a power.*/
{
// Take each unitelement and take it to the exponent power, then recursively do same
// for the subunits.
	for (unit::iterator i=begin();i!=end();i++) {
		i->exp *= exponent;
		(i->subunits).pow(exponent);
	}
	return (*this);
}


unit unit::simplify(void)
/* added 3/27/2001 to use default */
{
	return simplify(defaultsimplifyint);
}
	
	
unit unit::simplify(int level)
/*  All of the simplify(int) stuff moved here from simplify(char *) 3/27/2001 JB */
{
	// level 0 cancellage -- if it matches add the exps and toast the second one
	if (level==0) {
		for (unit::iterator i=begin();i!=end();) {
			unit::iterator iplus;
			i++; iplus=i; i--;
			for (unit::iterator j=iplus;j!=end();){
				if (i->abbr == j->abbr) {
					// add the exponents, or the mults for dimensionless
					if (i->type == "dimensionless") i->mult *= j->mult;
					else i->exp+=j->exp;
//					(*j)=unitelement(" ");
					// delete the redundant Junitelement
					unit::iterator tmp(j);
					tmp++; erase(j); j=tmp;
				}
				else j++;
			}
			if (i->exp == 0) {
				unit::iterator tmp(i);
				tmp++; erase(i); i=tmp;
			} else i++;
		}
	}
	
// level 10 cancel -- match types and leave multiplicative difference
	if (level==10) {
		for (unit::iterator i=begin();i!=end();i++) {
			unit::iterator iplus;
			i++; iplus=i; i--;
			for (unit::iterator j=iplus;j!=end();j++){
				if (i->type == j->type  &&  i->type != "dimensionless") {
					// add the exponents 
					i->exp+=j->exp;
					// figure out the conversion, and put it in as a dimless in place of j
					(*j)=unitelement(::pow(j->mult/i->mult, j->exp));
				}
			}
		}
		simplify(0);  // to put all the mults together
	}	
	
// level 20 cancel -- reduce to fundamentals.  Don't put together if 21.
	if (level==20  ||  level==21) {
		for (unit::iterator i=begin();i!=end();) {
			if (i->type == "conglomerate") {
				//copy the cong's subunits into the main at the end
				splice(end(), i->subunits);
				// delete the cong.
				unit::iterator tmp(i);
				tmp++; erase(i); i=tmp;
			} else i++;
		}
		if (level==20) simplify(10);
	}
	
// level 30 cancel -- build up into larger unit-objects
	if (level==30) {
		unit best(*this);			// will be the best-fit
		double min=expadd();		// the minimum expadd to test for goodness of fit
		unitelement t(" "); t++;// t is the unitelement to be tested
		for (;t.type != "dimensionless";t++) {
			unit cpy(*this);			// the copy of this that can be messed with
			double thistestnum;
			unit thisresult;
			thisresult=cpy.simplify(t);
			thistestnum=thisresult.expadd();
			if (thistestnum < min) { best = thisresult; min = thistestnum; }
		}
		*this = best;
	}
	
	return *this;
}

unit unit::simplify(unit aim)
/*  All of the simplify(unit) stuff moved here from simplify(char *) 3/27/2001 JB */
{		
//	cout << "inside simplify("<<aim<<")\n";
//	cout << "check for simplify(unit):  aim.stype="<<(int)aim.stype()<<" system="<<(int)Jcrap::unit::system<<"\n";
// convert to a unit system
	if (aim.stype() == Jcrap::unit::system) {
		// reduce to fundamentals
		simplify(20);
		// tack the new ones on the front end and convert by doing a simplify(10)
		aim.pow(0.);
		splice(begin(), aim.begin()->subunits);
		simplify(10);
	}
	
// directed cancel, convert to a specific unit
	else {
		// for each unitelement of the aim unit, figure out the optimal exponent
		// to use in *this.  Then divide by that and continue with the next aim unitelement
		unit answer;							// the component of aim in *this
		unit leftover(*this);				// the non-aim *this component
		
		for (unit::iterator i=aim.begin();i!=aim.end();i++){
			double optimalexp=0, min;
			min=order();
			for (double m=-10.;m<=10;m++) {
				double thisdiff;
				thisdiff=(leftover.diff(unit(*i).pow(m))).order();
//				cout << "Checking "<<*i<<" against "<<*this<<": m="<<m<<" thisdiff= " << thisdiff << "\n";
				if (thisdiff < min) { 
					min = thisdiff; 
					optimalexp=m; 
				}
			}
			unit unitpart(*i);
			unitpart.pow(optimalexp);
			answer*=unitpart;
			leftover*=unitpart.pow(-1.).simplify(20);
		}
		(*this) = answer*leftover;
	} 
		
	return *this;
}

/*unit unit::simplify(const char* lvlstr)*/
/* Created 2/28/2001 JB
	Cancels shit out.  If it gets a string, it parses the string as a Junit and tries to
	make the answer look like that (WELL, AT LEAST IT WILL IN THE FUTURE).  If it gets a
	number, that number is assumed to correspond to the simplification level:
		0:   Cancel out any Junitelements that match (except for the exponents)
		10:  Cancel out anything whose type matches, with the resulting conversion difference
			  such that what goes in is equal to what comes out
		20:  Reduce everything to fundamental units (non-conglomerates), with mult difference
		30:  Build up to the largest &/| closest signature conglomerate
		80:  Convert to a specific unit system specified in lvlstr, fundamental
			  to use this feature, call simplify.("mks") or whatever, and as long as that
			  system is in Junits.data it will try to convert to it.
		85:  Convert to lvlstr system, largest
		90:  specially reserved for the string case where it aims for that type */
/*{
// either simplify from simplify(int)
	istrstream lvlstrstr(lvlstr);
	if (isdigit(lvlstr[0])) {
		int level; 
		lvlstrstr >> level;
		simplify(level);
	}
// or simplify(unit)
	else simplify(unit(lvlstr));
	
	return (*this);
}*/
	
		
unit unit::convert(unit aim)
{
	aim.pow(0.);
	splice(begin(), aim);
	simplify(10);
	return *this;
}

unit unit::removedimensionless()
/* 3/20/2001 JB 
	Removes all dimensionless unitelements from the unit */
{
	for (unit::iterator i=begin();i!=end();)
		if (i->type == "dimensionless") {
			unit::iterator nexti;
			i++; nexti=i; i--;	
			erase(i);
			i=nexti;
		} else i++;
	return (*this);
}


unit unit::reorder()
/* 3/21/2001 JB
	Puts the elements of *this into the same order as in the unitdatabase */
{
	sort();				// assumes the < ordering of unitelement::operator< below
	return (*this);
}
	

unit unit::diff(unit r)
/* Created  3/22/2001 JB
	Returns the difference between two units.
	Method:  Multiply  *this * (1/r)*/
{
	return (r.pow(-1.)*(*this));
}
	
double unit::order()
/* Created  3/22/2001 JB
	Adds the absolute value of the fundamental exponents within a unit */
{
	double answer=0;
	unit fundthis(*this);
	fundthis.simplify(20);
	for (unit::iterator i=fundthis.begin();i!=fundthis.end();i++)
		answer += fabs(i->exp);
	return answer;
}

double unit::expadd()
/* 3/23/2001 JB
	Just adds the |exponents|, does not convert to Fund first */
{
	double answer=0;
	for (unit::iterator i=begin();i!=end();i++)
		if (i->stype() != unit::dimensionless  ||  i->mult != 1.) answer += fabs(i->exp);
	return answer;
}
	
void unit::setdefault(int def)
/* 3/27/2001 JB
	sets default simplification level */
{
	defaultsimplifyint=def;
}

unit::builtintype unit::stype()
/* 3/27/2001 JB
	checks all the unitelements and returns a builtintype as appropriate */
{
	if (size()==0) return unit::null;
	if (size()==1) return begin()->stype();
	return unit::conglomerate;
}

bool unit::operator==(unit r)
/* Created 3/22/2001 JB */
{
	unit l(*this);
	l=l.simplify();
	unit r2(r);
	r2=r.simplify();
	l.sort();  r2.sort();
	for (unit::iterator i=l.begin(), j=r2.begin();i!=l.end(), j!=r2.end();i++,j++) 
		if (!(*i == *j)) return 0;
	return 1;
}

unit unit::operator*=(unit r)
/* Created 2/28/2001 JB */
{
//	cout << *this << "  *  " << r << "  =  ";
	splice(end(), r);
//	cout << *this << "\n";
	simplify();
//	cout << "simplified:  \"" << *this << "\"\n";
	return (*this);
}

unit unit::operator*(unit r)
/* Created 2/28/2001 JB */
{
	unit answer(*this);
	answer*=r;
	return answer;	
}

unit unit::operator/(unit r)
/* Created 3/30/2001 JB */
{
	unit answer(*this);
	answer*=r.pow(-1);
	return answer;	
}

unit::operator double()
{
	double mfactor=1.;
	for (unit::iterator i=begin();i!=end();i++)
		if (i->type == "dimensionless") mfactor*=(i->mult);
	
	return mfactor;
}


bool unitelement::operator==(const unitelement& r)
/* Created 2/28/2001 JB */
{
	if (type != "dimensionless") return (abbr==r.abbr);
	else return (mult==r.mult);
}

//unit::~unit() { cout << "Now destructing unit "<<this<<"\n"; }


bool unitelement::operator< (unitelement r)
/* Created 3/21/2001 JB
	returns true if *this is encountered earlier in the unitdatabase than r*/ 
{
	for (list<unitelement>::iterator i=unitdatabase.begin();i!=unitdatabase.end();i++){
		if (abbr   == i->abbr) return 1;
		if (r.abbr == i->abbr) return 0;
	}
	return 0;
}

unit unitelement::operator++(int whatever)
/* 3/23/2001 JB
	Moves to next in unit database.  If end(), wrap around. */
{
	list<unitelement>::iterator i;
	for (i=unitdatabase.begin();i!=unitdatabase.end();i++)
		if (abbr   == i->abbr) break;
	i++;
	if (i==unitdatabase.end()) i = unitdatabase.begin();
	(*this) = (*i);
	return *this;
}
	

ostream& operator<<(ostream& out, const Jcrap::unitelement& j)
/* output override for a Junitelement */
{
	if (j.type=="dimensionless") {
		if (j.mult != 1.0) return out << j.mult;
		else return out << "";
	}
	else if (j.exp==1.) return out << j.abbr;
	return out << j.abbr << "^" << j.exp;
}



ostream& operator<<(ostream& out, const Jcrap::unit j)
/* output override for a Junit 

FUTURE WORK:
	allow for W / m^2 um   instead of W m^-2 um^-1
	possibly alter for formatting in plotutils graph function from Jcube*/
{
	return out << j.tostring();
}


string unit::tostring() const
{
	char answer[10000];
	ostrstream out(answer,10000);
	for (Jcrap::unit::const_iterator i(begin());i!=end();i++) 
		out << *i << " ";
	out << '\0';
	return string(answer);
}
