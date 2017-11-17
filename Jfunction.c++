#include "Jnr.h"

#define MAXIT 60


template <class Datatype>
Datatype JMAX(Datatype x, Datatype y)
{
	if (x > y) return x;
	else return y;
}

template <class Datatype>
Jcrap::function<Datatype>::function() 
{
	jMAX=200; 
	Dimensionality=1;
			
	// these used in bracketminimum() method
	GOLD = 1.618034;
	GLIMIT = 100.;
	TINY = 1.0e-20;
}

#define NRANSI
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

template <class Datatype>
pair<Datatype,Datatype> Jcrap::function<Datatype>::derivative(Datatype x, Datatype h) 
/* After NR dfridr   *** NOT TESTED *** */
{
	pair<Datatype,Datatype> answer;
	Datatype &err(answer.second);
	Datatype &ans(answer.first);
	
	if (h == -1.) h = x/10.;
	if (h ==  0.) h = 1.;
	int i,j;
	Datatype errt,fac,hh;
	
	vector<vector<Datatype> > aa(NTAB+1);
	for (i=0;i<=NTAB;i++) aa.at(i).resize(NTAB+1);
	
	hh=h;
	aa.at(1).at(1)=(the_function(x+hh)-the_function(x-hh))/(2.0*hh);
	err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		aa[1][i]=(the_function(x+hh)-the_function(x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			aa[j][i]=(aa[j-1][i]*fac-aa[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=JMAX(jabs(aa[j][i]-aa[j-1][i]),jabs(aa[j][i]-aa[j-1][i-1]));
			if (errt <= err) {
				err=errt;
				ans=aa[j][i];
			}
		}
		if (jabs(aa[i][i]-aa[i-1][i-1]) >= SAFE*(err)) break;
	}
	return answer;
}

#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE
#undef NRANSI

template <class Datatype>
pair<Datatype,Datatype> Jcrap::function<Datatype>::derivative_dumb(Datatype x, Datatype h) 
/* created 2014 June 10 JWB */
{
	pair<Datatype,Datatype> answer;
	
	if (h == -1.) h = x/10.;
	if (h ==  0.) h = 1.;
	
	Datatype measure1((the_function(x)-the_function(x-h))/h);
	Datatype measure2((the_function(x+h)-the_function(x))/h);
	answer.first = (measure1+measure2)/2.;
	answer.second= abs(measure1-measure2)/2.;
	
	return answer;
}

template <class Datatype>
pair<Datatype,Datatype> Jcrap::function<Datatype>::derivative2nd_dumb(Datatype x, Datatype h) 
/* created 2014 June 10 JWB */
{
	pair<Datatype,Datatype> answer;
	
	if (h == -1.) h = x/10.;
	if (h ==  0.) h = 1.;
	
	Datatype measure1((derivative_dumb(x,h/10.).first-derivative_dumb(x-h,h/10.).first)/h);
	Datatype measure2((derivative_dumb(x+h,h/10.).first-derivative_dumb(x,h/10.).first)/h);
	answer.first = (measure1+measure2)/2.;
	answer.second= abs(measure1-measure2)/2.;
	
	return answer;
}
#define JMAXP (jMAX+1)

template <class Datatype>
Datatype Jcrap::function<Datatype>::integrate_smooth(Datatype a, Datatype b, Datatype EPS,
	int K, string where)
{
	vector<Datatype> s(JMAXP), h(JMAXP+1);
	
	bool debug(0);
	
	h[0]=1.0;
	for (int j=0;j<jMAX;j++) {
		if (debug) cout << j << "th iteration inside Jcrap::function::integrate\n";
		s[j]=trapzd(a,b,j+1);
		if (debug) cout << "\tj=" << j << ", result is " << s[j] << "\n";
		if (j+1 >= K) {
			pair<Datatype, Datatype> answer;
			Datatype &ss(answer.first), &dss(answer.second);
			if (debug) cout << "before Jpolint\n";
//			answer = Jpolint(&h[-1],&s[-1],j+1,0.0);	// use all the points (new way)
//			answer = Jpolint(&h[1],&s[1],j-1,0.0);		// use all but 2 of the points	(new way 5/2003)	
			answer = Jpolint(&h[j-K],&s[j-K],K,0.0);  // its like this in NR.  Why did I change it?
																	// new way is faster, but might use some bum points
			Datatype limit;
			if (ss <= EPS) limit = EPS;    // manhandled to deal with zeroes better
			else if (ss <= sqrt(EPS)) limit = ss*sqrt(EPS);
			else limit = ss*EPS;			
			if (debug) cout << "ss=" << ss << " dss=" << dss << ";  limit = " << limit << "\n";

			if (jabs(dss) <= limit) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	cout << "Too many steps in routine Jcrap::function::integrate_smooth (qromb) from " << where << "\n";
	cout << "Integrating from " << a << " to " << b << " EPS=" << EPS << "\n";
//	cube thisintegration(plot(double(a), double(b), 100));
//	thisintegration.graph();
	return 0.0;
}
#undef JMAXP

#define JMAX 20

template <class Datatype>
Datatype Jcrap::function<Datatype>::integrate_rough(Datatype a, Datatype b, Datatype EPS, int jnk)
// translated to my c++ May 23, 2003 -- originally NR qtrap
// last int is not used, its just there to change between integrate_rough and integrate_smooth freely
{
	int j;
	Datatype s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(a,b,j);
		if (j > 5)
			if (jabs(s-olds) < EPS*jabs(olds) ||
				(s < EPS && olds < EPS)) return s;
		olds=s;
	}
	cout << "Too many steps in routine Jcrap::function::integrate_rough (qtrap)\n";
	cube thisintegration(plot(double(a), double(b), 100));
	thisintegration.graph();
	return s;
}
#undef JMAX



template <class Datatype>
Datatype Jcrap::function<Datatype>::trapzd(Datatype a, Datatype b, int n)
{
	Datatype x,tnm,sum,del;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(the_function(a)+the_function(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
			sum += the_function(x);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}





template <class Datatype>
cube Jcrap::function<Datatype>::plot(Datatype x0, Datatype x1, int points)
{
	cube answer(1,1,points);
	int i=0;
	if (osuppress < 1) cout << "Plotting from scratch --  00%";
	for (Datatype x=x0;i<points;x+=(x1-x0)/(points-1),i++) {
		if (osuppress < 1) printpercent(i, points-1);
		answer(0,0,i) = the_function(x);
		answer.Axis(Z,i) = x;
	}
	if (osuppress < 1) cout << "\n";
	return answer;
}

template <class Datatype>
cube Jcrap::function<Datatype>::plot(const cube& c)
{
	cube answer;
	
	if (Dimensionality <= 1) {
		answer = cube(1,1,c.N(Z));
		cout << "Plotting using cube template --  00%";
		for (int i=0;i<c.N(Z);i++) {
			printpercent(i, c.N(Z)-1);
			answer(0,0,i) = the_function(c.Axis(Z,i));
			answer.Axis(Z,i) = c.Axis(Z,i);
		}
		cout << "\n";
	} else if (Dimensionality <= 3) {
		answer = cube(c.N(X), c.N(Y), c.N(Z));
		cout << "Plotting (2D) --  00%";
		for (int x(0);x<c.N(X);x++)  {
			printpercent(x, c.N(X)-1);
			for (int y(0);y<c.N(Y);y++) {
				for (int z(0);z<c.N(Z);z++) {
					vector<double> point(Dimensionality);
					point[0] = answer.Axis(X,x) = c.Axis(X,x); 
					point[1] = answer.Axis(Y,y) = c.Axis(Y,y);
					if (Dimensionality == 3) point[2] = answer.Axis(Z,z) = c.Axis(Z,z);
					answer(x,y,z) = the_function(point);
				}
			}
		}
	}
	printpercent(99,99);
	cout << "\n";
	return answer;
}


template <class Datatype>
Datatype Jcrap::function<Datatype>::findzero(Datatype x1, Datatype x2, 
		Datatype xacc, string where)
{
	return findzero(x1,x2,xacc,0.,where);
}


template <class Datatype>
Datatype Jcrap::function<Datatype>::findzero(Datatype x1, Datatype x2, 
		Datatype xacc, Datatype offset, string where)
//  after NR zriddr
// added offset 2008 March 18.  Allows to find where the function hits a value "offset"
//     (instead of when it hits zero) without having to create a whole new function that
//		 subtracts it off and then doing findzero on that.
{
	int j;
	Datatype ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

	fl=the_function(x1)-offset;
	fh=the_function(x2)-offset;
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=-1.11e30;		// UNUSED
		for (j=1;j<=MAXIT;j++) {
			xm=0.5*(xl+xh);
			fm=the_function(xm)-offset;
			s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (jabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			fnew=the_function(ans)-offset;
			if (fnew == 0.0) return ans;
			if (Jsign(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (Jsign(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (Jsign(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else cout << "never get here.";
			if (jabs(xh-xl) <= xacc) return ans;
		}
		cout << "Jcrap::function<>::findzero exceed maximum iterations from " <<
				where << "\n";
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		cerr << "root must be bracketed in Jcrap::function<>::findzero from " <<
				where << ".\n";
		return 1.e6;
	}
	return 0.0;
}
#undef MAXIT


#include <math.h>
#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

// Returns an x,y pair corresponding to the minimum value y and where it occurs, x

// ax, bx, cx, increasing abcissas such that answer is between ax and cx, where
// f(bx) < f(ax) and f(cx).  tol is tolerance to which result is desired.

template <class Datatype>
pair<Datatype,Datatype> Jcrap::function<Datatype>::findminimum(Datatype ax, Datatype bx, Datatype cx, Datatype tol)
{
	pair<Datatype,Datatype> answer;
	Datatype f1,f2,x0,x1,x2,x3;

	x0=ax;
	x3=cx;
	if (jabs(cx-bx) > jabs(bx-ax)) {
		x1=bx;
		x2=bx+C*(cx-bx);
	} else {
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=the_function(x1);
	f2=the_function(x2);
	while (jabs(x3-x0) > tol*(jabs(x1)+jabs(x2))) {
		if (f2 < f1) {
			SHFT3(x0,x1,x2,R*x1+C*x3)
			SHFT2(f1,f2,the_function(x2))
		} else {
			SHFT3(x3,x2,x1,R*x2+C*x0)
			SHFT2(f2,f1,the_function(x1))
		}
	}
	if (f1 < f2) {
		answer.first=x1;
		answer.second=f1;
	} else {
		answer.first=x2;
		answer.second=f2;
	}
	return answer;
}

template <class Datatype>
pair<Datatype,Datatype> Jcrap::function<Datatype>::findmaximum(Datatype ax, Datatype bx, Datatype cx, Datatype tol)
{
	pair<Datatype,Datatype> answer;
	Datatype f1,f2,x0,x1,x2,x3;

	x0=ax;
	x3=cx;
	if (jabs(cx-bx) > jabs(bx-ax)) {
		x1=bx;
		x2=bx+C*(cx-bx);
	} else {
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=neg_function(x1);
	f2=neg_function(x2);
	while (jabs(x3-x0) > tol*(jabs(x1)+jabs(x2))) {
		if (f2 < f1) {
			SHFT3(x0,x1,x2,R*x1+C*x3)
			SHFT2(f1,f2,neg_function(x2))
		} else {
			SHFT3(x3,x2,x1,R*x2+C*x0)
			SHFT2(f2,f1,neg_function(x1))
		}
	}	
	if (f1 < f2) {
		answer.first=x1;
		answer.second=-f1;
	} else {
		answer.first=x2;
		answer.second=-f2;
	}
	return answer;
}
#undef C
#undef R
#undef SHFT2
#undef SHFT3


template <class Datatype>
inline void shft2(Datatype &a, Datatype &b, const Datatype c)
{
	a=b;
	b=c;
}
	
template <class Datatype>
inline void shft3(Datatype &a, Datatype &b, Datatype &c, const Datatype d)
{
	a=b;
	b=c;
	c=d;
}

template <class Datatype>
inline void mov3(Datatype &a, Datatype &b, Datatype &c, const Datatype d, const Datatype e,
	const Datatype f)
{
	a=d; b=e; c=f;
}

	
template <class Datatype>
vector<pair<Datatype, Datatype> > Jcrap::function<Datatype>::bracketminimum(Datatype a, Datatype b)
// modified from NR3 mins.h page 491 JWB 2014 June 18
// takes two initial x values for the function, which it then follows downhill
// and returns a size-3 vector containing the x,y values for the points that
// bracket the minimum.
{
//	static Datatype GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;  these now defined in the parent class
	vector<pair<Datatype, Datatype> > answer(3);
	double &ax(answer[0].first);
	double &bx(answer[1].first);
	double &cx(answer[2].first);
	double &fa(answer[0].second);
	double &fb(answer[1].second);
	double &fc(answer[2].second);
	ax=a;  bx=b;
	Datatype fu;
	fa=the_function(answer[0].first);
	fb=the_function(answer[1].first);
	if (fb > fa) {
		swap(ax,bx);
		swap(fb,fa);
	}
	cx=bx+GOLD*(bx-ax);
	fc=the_function(answer[2].first);
	while (fb > fc) {
		Datatype r=(bx-ax)*(fb-fc);
		Datatype q=(bx-cx)*(fb-fa);
		Datatype u=bx-((bx-cx)*q-(bx-ax)*r)/
			(2.0*SIGN(MAX(abs(q-r),TINY),q-r));
		Datatype ulim=bx+GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu=the_function(u);
			if (fu < fc) {
				ax=bx;
				bx=u;
				fa=fb;
				fb=fu;
				return answer;
			} else if (fu > fb) {
				cx=u;
				fc=fu;
				return answer;
			}
			u=cx+GOLD*(cx-bx);
			fu=the_function(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu=the_function(u);
			if (fu < fc) {
				shft3(bx,cx,u,u+GOLD*(u-cx));
				shft3(fb,fc,fu,the_function(u));
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u=ulim;
			fu=the_function(u);
		} else {
			u=cx+GOLD*(cx-bx);
			fu=the_function(u);
		}
		shft3(ax,bx,cx,u);
		shft3(fa,fb,fc,fu);
	}
	return answer;
}
