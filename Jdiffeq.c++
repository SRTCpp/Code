#include <iomanip>
#include <math.h>
#define SIGN(a,b) ((b) >= 0.0 ? abs(a) : -abs(a)) 

template <class Datatype>
Datatype abs(Datatype a)
{
	if (a>=0) return a;
	else return -a;
}

template <class Datatype>
void Jcrap::diffeq::odeint (Datatype x1, Datatype x2)
/*  April 18, 2001 JB
	Call with starting and ending x values to integrate between.
	Modified from nr::odeint */
{
/* added to convert to c++ 4/18/2001 */
	xp.resize(kmax);
	yp.resize(nvar);
	for (int i=0;i<nvar;i++) yp[i].resize(kmax);
	
	cout << "Inside odeint, ystart=" << ystart[0] << "\n";
	
	int nstp,i;
	Datatype xsav,x,hnext,hdid,h;

	x=x1;
	h=SIGN(h1,x2-x1);
	nok = (nbad) = kount = 0;
	for (i=0;i<nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=maxsteps;nstp++) {
		derivs(x,y,dydx);
		for (i=0;i<nvar;i++)
			yscal[i]=abs(y[i])+abs(dydx[i]*h)+tiny;
		if (kmax > 0 && kount < kmax && abs(x-xsav) > abs(dxsav)) {
			xp[kount]=x;
			for (i=0;i<nvar;i++) yp[i][kount]=y[i];		
			kount++;
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		rkqs(y,dydx,x,h,hdid,hnext);
		if (hdid == h) ++nok; else ++nbad;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=0;i<nvar;i++) ystart[i]=y[i];
			if (kmax) {
				xp[kount]=x;
				for (i=0;i<nvar;i++) yp[i][kount]=y[i];
				kount++;
			}
			return;
		}
		if (abs(hnext) <= hmin) {cout <<"Step size too small in odeint"; exit(1); }
		h=hnext;
	}
	cout << "Too many steps in routine odeint"; exit(1);
}

template <class Datatype>
Datatype max(Datatype a, Datatype b)
{
	if (a >= b) return a;
	else return b;
}

template <class Datatype>
Datatype min(Datatype a, Datatype b)
{
	if (a <= b) return a;
	else return b;
}

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

template <class Datatype>
void Jcrap::diffeq::rkqs(vector<Datatype>& y, vector<Datatype>& dydx, Datatype& x, 
		Datatype htry, Datatype& hdid, Datatype& hnext)
{
	int& n(nvar);
	Datatype errmax,h,htemp,xnew;
	vector<Datatype> yerr(nvar),ytemp(nvar);

	h=htry;
	for (;;) {
		this->rkck(x,h,ytemp,yerr);
		errmax=0.0;
		for (int i=0;i<n;i++) errmax=max(errmax,abs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? max(htemp,Datatype(h*0.1)) : min(htemp,Datatype(h*0.1)));
		xnew=x+h;
		if (xnew == x) { cout << "stepsize underflow in rkqs"; exit(1); }
	}
	if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
	else hnext=5.0*h;
	x += (hdid=h);
	for (int i=0;i<n;i++) y[i]=ytemp[i];
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

template <class Datatype>
void Jcrap::diffeq::rkck(Datatype x, Datatype h, vector<Datatype>& yout, vector<Datatype>& yerr)
{
	int i;
	int n=nvar;
	static Datatype a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	Datatype dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	vector<Datatype> ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n);

	for (i=0;i<n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	derivs(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	derivs(x+a3*h,ytemp,ak3);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	derivs(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	derivs(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	derivs(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}



template <class Datatype>
void Jcrap::diffeq<Datatype>::graph(string pustr="")
{
	static int i=0;
	char gname[200];
	sprintf (gname, "graph%i.out", i);
	ofstream out(gname);
	out << setprecision(20) << setw(23);
	for (int i=0;i<kount;i++) {
		out << xp[i]/3.15e7/1e6 << " ";
		for (int j=0;j<1;j++)
			out << yp[j][i] << " ";
		out << "\n";
	}
	out.close();
	
	string command="graph -T X ";
	command += pustr;
	command += " ";
	command += gname;
	command += " &";
	system(command.c_str());
	cout << "command: " << command << "\n";
	i++;
}




template <class Datatype>
Datatype Jcrap::diffeq<Datatype>::qromb(Datatype a, Datatype b)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	pair<Datatype, Datatype> ssP;
	Datatype &ss(ssP.first), &dss(ssP.second);
	vector<Datatype> s(JMAXP),h(JMAXP+1);
	int j;

	h[0]=1.0;
	for (j=0;j<JMAX;j++) {
		s[j]=trapzd(a,b,j);
		if (j >= K-1) {
			ssP = this->polint(&h[j-K],&s[j-K],K,0.0);
			if (abs(dss) <= EPS*abs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	cout << "Too many steps in routine qromb!  j>"<<JMAX<<"\n";
	exit(0);
	return 0.0;
}



#define NRANSI

// returns pair<y, dy>
template <class Datatype>
pair<Datatype, Datatype> Jcrap::diffeq<Datatype>::polint(Datatype *xa, Datatype *ya, int n, Datatype x)
{
	int i,m,ns=0;
	pair<Datatype, Datatype> answer;
	Datatype &y(answer.first), &dy(answer.second);
	Datatype den,dif,dift,ho,hp,w;
	vector<Datatype> c(n), d(n);

	dif=abs(x-xa[0]);
	for (i=0;i<n;i++) {
		if ( (dift=abs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) cout << "Error in routine diffeq::polint\n";
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	return answer;
}
#undef NRANSI

template <class Datatype>
Datatype Jcrap::diffeq<Datatype>::trapzd(Datatype a, Datatype b, int n)
/* Checked out and declared correct JB 5/7/2001 */
{
	Datatype x,tnm,sum,del;
	static Datatype s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(func(a)+func(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += func(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

template <class Datatype>
Datatype Jcrap::diffeq<Datatype>::qtrap(Datatype a, Datatype b)
{
	int j;
	Datatype s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(a,b,j);
		if (j > 5)
			if (abs(s-olds) < EPS*abs(olds) ||
				(s == 0.0 && olds == 0.0)) return s;
		olds=s;
	}
	cout << "Too many steps in routine qtrap, more than "<<JMAX<<"\n";
	return 0.0;
}

#define MAXIT 60

template <class Datatype>
Datatype Jcrap::diffeq<Datatype>::zriddr(Datatype x1, Datatype x2, Datatype xacc)
{
	int j;
	Datatype ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

	fl=auxfunc(x1);
	fh=auxfunc(x2);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=-1.11e30;		// UNUSED
		for (j=1;j<=MAXIT;j++) {
			xm=0.5*(xl+xh);
			fm=auxfunc(xm);
			s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			fnew=auxfunc(ans);
			if (fnew == 0.0) return ans;
			if (SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else cout << "never get here.";
			if (fabs(xh-xl) <= xacc) return ans;
		}
		cout << "Jcrap::diffeq<>::zriddr exceed maximum iterations\n";
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
		cerr << "root must be bracketed in Jcrap::diffeq<>::zriddr.\n";
		return 1e6;
	}
	return 0.0;
}
#undef MAXIT

template <class nTYPE>
unsigned long Jcrap::diffeq<nTYPE>::locate(vector<nTYPE> xx, nTYPE x)
/* Send this a vector with the values, and x, the value you're looking for */
{
	unsigned long n=kount;
	unsigned long ju,jm,jl, j;
	int ascnd;
	
	jl=-1;
	ju=n;
	ascnd=(xx[n-1] >= xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	if (x == xx[0]) j=0;
	else if(x == xx[n-1]) j=n-2;
	else j=jl;
	return j;
}


#define NRANSI

template <class nTYPE>
nTYPE polint(nTYPE xa[], nTYPE ya[], int n, nTYPE x)
{
	int i,m,ns=0;
	nTYPE den,dif,dift,ho,hp,w, y;
	vector<nTYPE> c(n), d(n);

	dif=abs(x-xa[0]);
	for (i=0;i<n;i++) {
		if ( (dift=abs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) cout <<"Error in routine polint\n";
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (2*ns < (n-m) ? c[ns+1] : d[ns--]);
	}
	return y;
}
#undef NRANSI

