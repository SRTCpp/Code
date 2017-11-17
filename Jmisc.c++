#include <string>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <iomanip>
#include <iostream>
#include "Jcrap.h"

using namespace std;

extern long idum;

int osuppress=0;

void printpercent(int n, int d)
{
	if (d==0) d=1;
	if (d<100  ||  !(n % (int)(d/100) )){
		int pct;
		if ( (pct=(100*n/d)) >= 10){
			if (pct<100) cout << char(8) << char(8) << char(8) << pct << '%';
			else cout << char(8) << char(8) << char(8) << char(8) << pct << '%';
		}
		else
			cout << char(8) << char(8) << pct << '%';
		cout.flush();
	}
}

string float2str(float f, int n)
{
	string s;
	char cs[20];
	sprintf(cs, "%e", f);
	cs[n+1]= '\0';
	s = cs;
	return s;
}

string double2str(double f, int n)
{
	string s;
	char cs[20];
	sprintf(cs, "%le", f);
	cs[n+1]= '\0';
//	cout << "f=" << f << " translated as " << cs;
	s = cs;
//	cout << ", i.e., " << s << "\n";
	return s;
}


string float2fstr(float f, int n)
{
	string s;
	char cs[20];
	sprintf(cs, "%f", f);
	cs[n+1]= '\0';
	s = cs;
	return s;
}

string double2fstr(double d, int n)
{
	string s;
	char cs[200];
//	sprintf(cs, "%.20lf", f);
	sprintf(cs, "%lf", d);
	cs[n]= '\0';
	s = cs;
	return s;
}
string float2binstr(float f)
{
	string s;
	char cs[20];
	void *dataloc=&f;
	char *cp;
	cp=(char*)dataloc;
	strncpy(cs, cp, 4);
	s = cs;
	return s;
}

float binstr2float(string s)
{
	float f;
	char cs[100];
	strcpy(cs, s.c_str());
	if (cs[0] == '\0') f=0.;
	else {
		void *dataloc=&f;
		char *fp;
		fp=(char*)dataloc;
		strncpy(fp, cs, 4);
	}
	return f;
}

int str2int(string s)
{
	int i;
	char cs[1000];
	strcpy(cs, s.c_str());
	sscanf(cs, "%d", &i);
	return i;
}


unsigned long str2ulong(string s)
{
	unsigned long answer;
	char cs[1000];
	strcpy(cs, s.c_str());
	sscanf(cs, "%lu", &answer);
	return answer;
}

string int2str(int i)
{
	string s;
	char cs[100];
	sprintf(cs, "%d", i);
	s = cs;
	return s;
}

string int2str(int i, int n)
{
	string s;
	char cs[100];
	sprintf(cs, "%d", i);
	cs[n+1]= '\0';
	s = cs;
	return s;
}

string int2str(long long ll)
{
	string s;
	char cs[100];
	sprintf(cs, "%lli", ll);
	s = cs;
	return s;
}

string int2str(long l)
{
	string s;
	char cs[100];
	sprintf(cs, "%ld", l);
	s = cs;
	return s;
}

string lblstr(int lblsize)					//  converts int to length 5 string
{
	string s;
	char cs[200];
	sprintf(cs,"%-5d", lblsize);
	s=cs;
	return s;
}

float str2float(string s)
{
	float f;
	char cs[1000];
	strcpy(cs, s.c_str());
	if (cs[0] == '\0') f=0;
	else sscanf(cs, "%f", &f);
	return f;
}

double str2double(string s)
{
	double d;
	char cs[1000];
	strcpy(cs, s.c_str());
	if (cs[0] == '\0') d=0.;
	else sscanf(cs, "%lf", &d);
	return d;
}



inline long exp(int b, int p)
{
	long ans=1;
	for (int i=0;i<p;i++) ans*=b;
	return ans;
}

inline float exp(float b, int p)
{
	float ans=1;
	for (int i=0;i<p;i++) ans*=b;
	return ans;
}

void raninit()
{
	time_t tp;
	long secs;
	idum = 0-time(&tp);				// idum is a global long
}

string asstring(otype o)
{
	string answer("");
	
	if (o==oX)         answer = "X";
	else if (o==oPS)   answer = "ps";
	else if (o==oGIF)  answer = "gif";
	else if (o==oBW)   answer = "bw";
	else if (o==oCOLOR)answer = "color";
	else if (o==oMOVIE)answer = "mpg";
	else if (o==oCOLORMOVIE)answer = "mpg";
	else if (o==oFIG)  answer = "fig";
	else answer = "none";
	
	return answer;
}
