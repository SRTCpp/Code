#include<string>
#include<cstdio>
#include<map>
#include<iostream>
#include<math.h>
#include<strstream>
#include<fstream>
#include<utility>
#include<algorithm>
#include<vector>
#include<list>
#include<math.h>
#include"../../nr/nrutilj.h"

#include"Jcrap.h"

void jfourn(vector<float> &data, unsigned long nn[], int ndim, int isign);

long idum=0;

bool is2pow(int t)
{
	int ipow=1;
	if (t<1) return 0;
	while (ipow<=t) {
		if (ipow==t) return 1;
		ipow*=2;
	}
	return 0;
}
	
cube cube::ysigma()
{
	cube mean;
	mean=(*this).plane(Y);
	mean.write("sigxtplane");
	mean/=ny;
	mean.write("sigmamean");
	
	cube sigma(mean.Nx(), 1, mean.Nz(), 0.0);			// create new cube same dims as mean
	for (int x=0;x<nx;x++)
		for (int z=0;z<nz;z++)
			for (int y=0;y<ny;y++){
				float b=(*this)(x,y,z)-mean(x,0,z);
				sigma(x,0,z)+=b*b;
//				cout << b << " ";
			}
	sigma/=ny-1;
	sigma.write("b4sqrt");
	sigma=sigma.squareroot();
	
	return sigma;	
}

cube cube::ymedian(float f)
{
	cube median(nx, 1, nz, 0.0);
	
	for (int x=0;x<nx;x++) median.xaxis[x]=xaxis[x];
	for (int y=0;y<ny;y++) median.yaxis[0]+=yaxis[y];
	median.yaxis[0]/=ny;
	for (int z=0;z<nz;z++) median.zaxis[z]=zaxis[z];
	
	for (int x=0;x<nx;x++)
		for (int z=0;z<nz;z++){
			float yarr[ny+1];					//nr code starts arrays at 1
			for (int y=0;y<ny;y++) yarr[y+1]=(*this)(x,y,z);
			median(x,0,z)=select(int((ny-1)*f+1), ny, yarr);
		}
	
	return median;
}

cube cube::xmedian(float f)
{
	cube median(1, ny, nz, 0.0);
	
	for (int x=0;x<nx;x++) median.xaxis[0]+=xaxis[x];
	median.xaxis[0]/=nx;
	for (int y=0;y<ny;y++) median.yaxis[y]=yaxis[y];
	for (int z=0;z<nz;z++) median.zaxis[z]=zaxis[z];
	
	for (int y=0;y<ny;y++)
		for (int z=0;z<nz;z++){
			float xarr[nx+1];					//nr code starts arrays at 1
			for (int x=0;x<nx;x++) xarr[x+1]=(*this)(x,y,z);
			median(0,y,z)=select(int((nx-1)*f+1), nx, xarr);
		}
	
	return median;
}

cube cube::zmedian(float f)
{
	cube median(nx, ny, 1, 0.0);
	
	for (int z=0;z<nz;z++) median.zaxis[0]+=zaxis[z];
	median.zaxis[0]/=nz;
	for (int y=0;y<ny;y++) median.yaxis[y]=yaxis[y];
	for (int x=0;x<nx;x++) median.xaxis[x]=zaxis[x];
	median.hdr = hdr;
	median.cubetype = cubetype;
	median.byteorder = byteorder;
	
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++){
			float zarr[nz+1];					//nr code starts arrays at 1
			for (int z=0;z<nz;z++) zarr[z+1]=(*this)(x,y,z);
			median(x,y,0)=select(int((nz-1)*f+1), nz, zarr);
		}
	
	return median;
}

cube cube::squareroot()
{
	cube temp(*this);
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				temp(x,y,z)=sqrt((*this)(x,y,z));
	return temp;
}

cube cube::exp()
// output pixel = e^(input pixel)
{
	cube answer(*this);
	for (int x=0;x<answer.N(X);x++)
		for (int y=0;y<answer.N(Y);y++)
			for (int z=0;z<answer.N(Z);z++)
				answer(x,y,z) = expf((*this)(x,y,z));
	return answer;
}

cube cube::cosine()
// output pixel = cos(input pixel)
{
	cube answer(*this);
	for (int x=0;x<answer.N(X);x++)
		for (int y=0;y<answer.N(Y);y++)
			for (int z=0;z<answer.N(Z);z++)
				answer(x,y,z) = ::cos((*this)(x,y,z));
	return answer;
}

cube cube::floatfunc(float f(float))
// output pixel = f(input pixel)
{
	cube answer(*this);
	for (int x=0;x<answer.N(X);x++)
		for (int y=0;y<answer.N(Y);y++)
			for (int z=0;z<answer.N(Z);z++)
				answer(x,y,z) = f((*this)(x,y,z));
	return answer;
}

cube cube::ytrysigma()
{
	cube sigma;
	sigma=(*this).ysigma();
	cube mean;
	mean=((*this).plane(Y,0,ny-1));
	
	cube bad(nx, ny, nz, 0);
	cout << "Trysigma  0%";
	for(int x=0;x<nx;x++){
		printpercent(x,nx-1);
		for(int y=0;y<ny;y++){
			for(int z=0;z<nz;z++){
				cube testbad(1,ny,1,0);
				for(int Y=0;Y<ny;Y++){
					if (y!=Y) testbad(0,Y,0)=(*this)(x,y,z);
					else testbad(0,Y,0)=mean(x,0,z);
					float nusigma=(testbad.ymedian())(0,0,0);
					bad(x,y,z)=sigma(x,0,z)-nusigma;
				}
			}
		}
	}
	return bad;
}

float* cube::polyfit(int degree)
{
/*  polyfit takes as its input the cube of data to be fit.
	It expects this data in the form of a 2xNx2 cube, with each of the N
	points having an x, y (z=0) and dx, dy (z=1).  If the incoming cube
	is 2xNx1, we assume unit standard deviations.  
	oh, by the way, this function can only handle nonzero sigma's in the Y
	coordinate at the moment, don't know when/if x will be implimented.
	This function returns a float array with [0] being the chisq parameter
	and [1]..[degree+1] as the polynomial coefficients in ascending order.  
	*/
	 
	float x[ny+1], y[ny+1], sig[ny+1];
	for (int i=1;i<=ny;i++){
		x[i] = (*this)(0,i-1,0);
		y[i] = (*this)(1,i-1,0);
		if (nz == 2) sig[i] = (*this)(1,i-1,1);
		else sig[i]=1;
		if (!sig[i]) sig[i] = 1;
	}
	
	float* a;
	a = new float[degree+2];
	
	float *u[ny+1], *v[degree+2], w[degree+2];
	float chisq;
	
	for (int i=0;i<ny+1;i++) u[i]=new float[degree+2];
	for (int i=0;i<degree+2;i++) v[i]=new float[degree+2];
	
	jsvdfit (x, y, sig, ny, a, degree+1, u, v, w, &chisq, &fpoly);
	
	a[0]=chisq;
	return a;
}

void Jfpoly(double x, vector<double>& p)
{
	int j;
	int np=p.size()-1;
	
	p[1]=1.0;
	for (j=2;j<=np;j++) p[j]=p[j-1]*x;
}


vector<float> cube::vpolyfit(int degree, float (*f)(float))
/* Takes a cube with one of its dimensions being 1, one being 2(or 1 for stddev=f(data))
	and the other of whatever length and returns a vector of the resulting polynomial fit
	with chisq=[0] and [1]..[degree+1] as the polynomial coefficients in ascending order.
	
	created 11/13/2k JB
	added noise function 10/31/2001
	*/
{
	int n, nsig;
	axis longaxis, sigaxis;
	if (N(X)>N(Y) && N(X)>N(Z)) { 
		longaxis=X; n=N(X);
		if (N(Y)>N(Z)) sigaxis=Y;
		else sigaxis=Z;
	}
	if (N(Y)>N(X) && N(Y)>N(Z)) { 
		longaxis=Y; n=N(Y);
		if (N(X)>N(Z)) sigaxis=X;
		else sigaxis=Z;
	}
	if (N(Z)>N(X) && N(Z)>N(Y)) { 
		longaxis=Z; n=N(Z);
		if (N(X)>N(Y)) sigaxis=X;
		else sigaxis=Y;
	}
	nsig=N(sigaxis);
	if (nsig>2) cout << "ERROR --- vpolyfit not yet ready for more than 1 dimension!!!\n";
			
	vector<double> x(n+1);
	vector<double> y(n+1);
	vector<double> sigma(n+1);
	
	for (int i=1;i<=n;i++){
		x[i]=Axis(longaxis,i-1);
		if (N(sigaxis)==1) {
			if (f==0) sigma[i]=1.;
			else sigma[i] = f((*this)(longaxis, i-1));
		}
		if (longaxis==X && sigaxis==Y) {
			y[i]=(*this)(i-1,0,0);
			if (N(sigaxis)!=1) sigma[i] = (*this)(i-1,1,0);
		}		
		if (longaxis==X && sigaxis==Z) {
			y[i]=(*this)(i-1,0,0);
			if (N(sigaxis)!=1) sigma[i] = (*this)(i-1,0,1);
		}		
		if (longaxis==Y && sigaxis==X) {
			y[i]=(*this)(0,i-1,0);
			if (N(sigaxis)!=1) sigma[i] = (*this)(1,i-1,0);
		}		
		if (longaxis==Y && sigaxis==Z) {
			y[i]=(*this)(0,i-1,0);
			if (N(sigaxis)!=1) sigma[i] = (*this)(0,i-1,1);
		}		
		if (longaxis==Z && sigaxis==X) {
			y[i]=(*this)(0,0,i-1);
			if (N(sigaxis)!=1) sigma[i] = (*this)(1,0,i-1);
		}		
		if (longaxis==Z && sigaxis==Y) {
			y[i]=(*this)(0,0,i-1);
			if (N(sigaxis)!=1) sigma[i] = (*this)(0,1,i-1);
		}
	}
	
	vector<double> a(degree+2);
	
	cout << "Entering Jsvdfit, degree="<<degree<<", a.size()="<<a.size()<<"\n";
	Jsvdfit (x, y, sigma, a, &Jfpoly);
	
	cout << "vpolyfit with a=" << " ";
	for (int i=1;i<a.size();i++) cout <<" a["<<i<<"]="<<a[i] << " ";
	cout << "\n";
	
	vector<float> af(a.size());
	for (int i=0;i<a.size();i++) af[i]=a[i];
	return af;
}

cube cube::cpolyfit(int degree, float (*f)(float))
/* Date unknown.  inspection (and memory leak fix) on Oct 31, 2001 implies
	this cube fits a polynomial to a 1-dimensional cube along any axis (! this is
	back when I was writing axis independant code, its pretty cool!)  It seems
	to then output as a cube the resulting straight line. */
{
	vector<float> a;
	a=(*this).vpolyfit(degree, f);
	return cpoly(a);
}

cube cube::cpoly(vector<float> a)
/* Oct 31, 2001 JB
	Created from the 2nd half of cpolyfit so that vpolyfit can be run outside and
	then the result fed into this.  That way you could, if you wanted, get both
	the coefficients and the graph out.  Or something */		
{	
	int n;
	axis longaxis;
	if (N(X)>N(Y) && N(X)>N(Z)) { longaxis=X; n=N(X); }
	if (N(Y)>N(X) && N(Y)>N(Z)) { longaxis=Y; n=N(Y); }
	if (N(Z)>N(X) && N(Z)>N(Y)) { longaxis=Z; n=N(Z); }
	
	cube result;
	if (longaxis==X) result=cube(n,1,1);
	if (longaxis==Y) result=cube(1,n,1);
	if (longaxis==Z) result=cube(1,1,n);
	
	result.copyaxis(*this,longaxis);
	
// added 10.30.2001 so that lines come out reasonably in graph()
// toasted and replaced with dirinc 2015 February 13 -- CHANGE UNTESTED
/*	vector<double> xvals;
	for (int w=0;w<n;w++) xvals.push_back(result.Axis(longaxis,w));
	xvals.sort();
	
	vector<double>::iterator j(xvals.begin());*/
	result = result.dirinc(longaxis);
	for (int i=0;i<n;i++) {
		if (longaxis==X) result(i,0,0)=polyval(result.Axis(longaxis,i), a);
		if (longaxis==Y) result(0,i,0)=polyval(result.Axis(longaxis,i), a);
		if (longaxis==Z) result(0,0,i)=polyval(result.Axis(longaxis,i), a);
	}
	
	
	return result;
}
		
	



float cube::interp(double x, int d, axis a)
{
/*  interp takes as its input the cube of data to be interpolated.
	It expects this data in the form of a 2xNx1 cube, with each of the N
	points having an x(0,i,0), y(1,i,0).   
	
	Returns the interpolated value
	
	Updated 2012 August 22:  Changed to allow for new double-sized axis values
	
	Updated 2/24/2k:  The most insidious bug surmounted by checking to see if medianfac is
	0, and if it is, trying the average instead of the median.  More than half of my points
	were zero, but it took me 3 hours to figure out that that would make the median zero.
	Go figure.
	
	Updated 2/19/00:  Made more efficient by not creating 2 float[N] arrays for EVERY 
	call by writing Jlocate and then just making the particular part of the array that is
	necessary.  Now it stores the median so that it is calculated only every time the 
	array (or actually the first element of the array) changes.
	
	Updated 2/17/00 to handle wavelength data along any axis, given as 3rd parameter
	of input.  If none given, interp assumes that you mean the to interpolate using
	wavelength data from the longest axis, defaulting to Y if anything funny is going on.
	
	Updated 1/5/00 to take a wavelength calibrated 1xNx1 cube also (either works)
	*/ 
	

	if (a==UNK){
		if (nx > ny && nx > nz) { a=X; }
		if (nz > ny && nz > nx) { a=Z; }
		if (a==UNK) { a=Y; }
	}
	int n(N(a));
	
	vector<float> lf(n);
	for (int i(0);i<n;i++) lf.at(i)=Axis(a,i);
	float *l(&lf.at(0));

	double medianfac=1.;
	if (  keyword("medianverify") != ""  &&
			str2int(keyword("medianverify")) == ((long long)this)  ) {
		medianfac = str2float(keyword("median"));
	} else {
		float xa[n+1], ya[n+1];
	
		for (int i=1;i<=n;i++) xa[i]=l[i-1];	
		if (a==X) for (int i=1;i<=nx;i++) ya[i] = fabs((*this)(i-1,0,0));
		if (nx==2 && a==Y){
			for (int i=1;i<=ny;i++){
				xa[i] = (*this)(0,i-1,0);
				ya[i] = fabs((*this)(1,i-1,0));
			}
		} else if (nx==1 && a==Y) for (int i=1;i<=ny;i++) ya[i] = fabs((*this)(0,i-1,0));		
		if (a==Z) for (int i=1;i<=nz;i++) ya[i] = fabs((*this)(0,0,i-1));
		
		double sum=0;
		for (int i=1;i<=n;i++) sum+=ya[i];
		unsigned long mid=n/2;
		if (!mid) mid=1;
		medianfac = select(mid, n, &ya[0]);
		if (medianfac==0.) {
			medianfac=sum/n;
		}
		keyword("median", float2str(medianfac));
		keyword("medianverify", int2str((long long)this));
		if (medianfac == 0.) { cout << "interp of entirely zero function\n"; return 0.; }
//		cout << "Median calculated as "<<medianfac<<" for "<<(int)this<<"\n";
	} 

	
//	for (int i=1;i<=n;i++) ya[i]/=medianfac;
	
	if (x < (l[0] - 2*(l[1]-l[0])) || x > (l[n-1] + 2*(l[n-1]-l[n-2]))){
		cout << " POSSIBLE ERROR IN INTERP:\n";
		cout << "requested value (" << x << ") is beyond available range (" <<
				l[0] << ".." << l[n-1] << ")\n";
	}
	
	int loc;
	loc=Jlocate(x,a);
	
	if (d%2) loc -= (d+1)/2 - 1;
	else if (d<0) loc += d/2;
	else loc -= d/2 - 1;
	
	if (loc < 0){
		d += loc;
		loc=0;
	}
	
	if (loc+d > n-1) d -= loc+d-n+1;

// 		Creating a simple (and small) array for nr::polint
	float yc[d+2];
//	for (int i=0;i<d+1;i++) xc[i]=l[loc+i];
	if (a==X) for (int i=0;i<d+1;i++) yc[i]=(*this)(loc+i,0,0)/medianfac;
	if (a==Y) for (int i=0;i<d+1;i++) yc[i]=(*this)(0,loc+i,0)/medianfac;
	if (a==Z) for (int i=0;i<d+1;i++) yc[i]=(*this)(0,0,loc+i)/medianfac;
	
	float y, dy;
	polint(&l[loc-1], &yc[-1], d+1, x, &y, &dy);
	y*=medianfac;
	return (y);
}

float cube::interp2d(float c1, float c2, int d) const
/* Created 3/6/2k
	This is a (somewhat?) generalized 2d interpolator.  It determines what dimension to
	interpolate in by which one is only size one, and defaults to XY.  */
{
	// added special bicubic function to hopefully be faster:  2006 January 23
	bool bicubic(0);
	if (d==3)
		bicubic = 1;
	
	
	int dorig=d;
	axis a;
	if (nx==1) a=X;
	if (ny==1) a=Y;
	if (nz==1) a=Z;
	if (a==UNK) a=Z;
	if (N(Uax(a))<d+1) cout << "ERROR -- interp2d finds dimension less than degree+1\n";
	if (N(Dax(a))<d+1) cout << "ERROR -- interp2d finds dimension less than degree+1\n";

	if (c1 < (Axis(Uax(a),0) - 2*(Axis(Uax(a),1)-Axis(Uax(a),0))) || 
		c1 > (Axis(Uax(a),N(Uax(a))-1)) + 2*(Axis(Uax(a),N(Uax(a))-1)-Axis(Uax(a),N(Uax(a))-2))){
		cout << " POSSIBLE ERROR IN INTERP2d:\n";
		cout << "requested value (" << c1 << ") is beyond available range (" <<
				Axis(Uax(a),0) << ".." << Axis(Uax(a),N(Uax(a))-1) << ")\n";
	}
	if (c2 < (Axis(Dax(a),0) - 2*(Axis(Dax(a),1)-Axis(Dax(a),0))) || 
		c2 > (Axis(Dax(a),N(Dax(a))-1)) + 2*(Axis(Dax(a),N(Dax(a))-1)-Axis(Dax(a),N(Dax(a))-2))){
		cout << " POSSIBLE ERROR IN INTERP2d:\n";
		cout << "requested value (" << c1 << ") is beyond available range (" <<
				Axis(Dax(a),0) << ".." << Axis(Dax(a),N(Dax(a))-1) << ")\n";
	}
	// determine location to start the interpsubcube
	int loc;
	loc=Jlocate(c1,Uax(a));
	if (d%2) loc -= (d+1)/2 - 1;
	else loc -= d/2 - 1;
	if (loc < 0){
		d += loc;
		loc=0;
	}	
	if (loc+d > N(Uax(a))-1) d -= loc+d-N(Uax(a))+1;
		
// interpolate values for interpsubcube
	if (!bicubic) {
		cube si(1,1,d+1);  // create new cube for the subcube interpolations
		for (int i=0;i<d+1;i++){
			cube c;
			c=skewer(Dax(a), 0, loc+i);
			if (keyword("median")!="")c.keyword("medianverify",int2str((long long)&c));
			si(0,0,i) = c.interp(c2, dorig);
			si.Zaxis(i) = Axis(Uax(a), loc+i);
			
		}
		return si.interp(c1,d);
	} else {
		
	}
}
	
cube cube::cylindermap(float cx, float cy, float r, float res, methodtype m,
		float e, double rotationangle, double subsclat, double subsclon)
/* Created 3/6/2k
	Maps a disk onto a cylindrical projection.  Right now only works for North at top and
	South at bottom, equator through center L to R.  But this should be enough; to
	generalize just stick a cylinderrotate function that rotates around Euler angles to
	move the result of this function into the result of the funny-angled disk. 
	This function is XYcentric.  Deal.
	e is oblateness see circle
	res is the number of pixels per degree in output image
	*/
/* Added the prerotation thing 2004 November 17 */
/* changed to ortho projection in Curtis Handout 2005 April 21 

   Seems that phi = longitude, theta = latitude*/
{
	cube merc((int)(360*res), (int)(180*res), nz);
	merc.copyaxis(*this, Z);
	
	merc.Xtype(DEGREES);
	for (int x=0;x<merc.N(X);x++) merc.Axis(X, x) = -180.+x/res;
	merc.Ytype(DEGREES);
	for (int y=0;y<merc.N(Y);y++) merc.Axis(Y, y) = -90.+y/res; // SOUTH latitude
	
	osuppress++;
	if (osuppress-1==0) cout << "Performing cylindrical mapping --   00%";
	float phi, theta, PI=3.14159265;
	for (int z=0;z<merc.N(Z);z++) {
//		cube frame;
//		frame=(*this).plane(Z, z);
		for (int x=0;x<merc.N(X);x++) {
			if (osuppress-1==0) printpercent(x+z*merc.N(X),merc.N(Z)*(merc.N(X)-1));
			for (int y=0;y<merc.N(Y);y++) {
				phi=merc.Axis(X, x)*PI/180.;
				theta=merc.Axis(Y, y)*PI/180.;
				if (phi < PI/2 && phi > -PI/2.  &&  theta < PI && theta > -PI) {
					float dx, dy;
					dx=r*cosf(theta)*sinf(phi-subsclon);
					dy=r*(cosf(subsclat)*sinf(theta)-sinf(subsclat)*cosf(theta)*cos(phi-subsclon));
					if (rotationangle) {
						double Dx(dx);
						dx = dx*cos(rotationangle) + dy*sin(rotationangle);
						dy =-Dx*sin(rotationangle) + dy*cos(rotationangle);
					}
					dx /= 1-e;  // oblatness correction
					dx += cx;  // add in centers
					dy += cy; 
					if (dx>0 && dx<N(X)-2  &&  dy>0 && dy<N(Y)-2)
						merc(x,y,z) = (*this)(dx, dy, z, m);
				}
			}
		}
	}
	if (osuppress-1==0) { printpercent(1,1); cout << '\n'; }
	
	merc.Dtype(Dtype());
	osuppress--;
	return merc;
}

float ew0(float inlon)
{
	float answer((360.-inlon));// to E lon
	while (answer > 180.) answer -= 360.;
	while (answer <-180.) answer += 360.;
	return answer;
}

float dontmess(float inlon)
{
	return inlon;
}

cube cube::cylinderrotate(cube rdeg, methodtype m)
/* created 3/9/2k
	Rotates cylindrical map around the z (theta=90deg) axis.
	only works for degrees for the moment
	takes cube 1x1xnz with amounts to rotate.  if 1x1x(n<nz) pads with rdeg(0,0,0)*/
{
	if (rdeg.N(Z)!=nz){
		cube newrdeg(1,1,nz,rdeg(0,0,0));
		for (int z=0;z<rdeg.N(Z);z++) newrdeg(0,0,z)=rdeg(0,0,z);
		rdeg=newrdeg;
	}
	cube r(*this);
	r*=0.;
	
	if (!osuppress) cout << "Cylinderrotating --";
	if (nz==1) osuppress++;
	if (!osuppress) cout << "   00%";
	
	float thisrot;
	LOOPOVERALLPCT { 
		thisrot=rdeg(0,0,z);
		while (Axis(X,x)-thisrot>(2*Axis(X,N(X)-1)-Axis(X,N(X)-2))) thisrot+=360.;
		while (Axis(X,x)-thisrot<Axis(X,0)) thisrot-=360.;
		r(x,y,z)=(*this)(Axis(X,x)-thisrot, y, z, m);
	}
	if (!osuppress) { printpercent(1,1); cout << "\n"; }
	if (nz==1) osuppress--;
	return r;
}

cube cube::cylindereuler(double phi, double theta, double psi, methodtype m, float(*lonconvention)(float))
{
	return this->cylindereuler(*this, phi, theta, psi, m, lonconvention);
}

static double dpr(180./3.1415926535897932384);

cube cube::cylindereuler(const cube& templatecube, double phi, double theta,
									double psi, methodtype m, float (*lonconvention)(float))		
{
	if (lonconvention==0)
		lonconvention = &dontmess;
	
// setting up output cube
	cube eulered(templatecube.N(X), templatecube.N(Y), N(Z), 0.);
	eulered.copyaxis(templatecube, X);
	eulered.copyaxis(templatecube, Y);
	eulered.copyaxis(*this, Z);

// take care of output nuances
	if (!osuppress) cout << "Euler angle rotating --";
	osuppress++;
	if (!(osuppress-1)) cout << "   00%";
	
//	double dpr = 180./3.1415926535897932384;

// for each new euler pixel, calculate back which original pixel corresponds
	for (int x=0;x<eulered.N(X);x++) {
		printpercent(x, eulered.N(X)-1);
		for (int y=0;y<eulered.N(Y);y++) {

			/*
//	first convert from answer's lat/lon to x,y,z with unit radius
			double xE, yE, zE;    // Eulered coords
			double lonE(eulered.Axis(X,x)/dpr), latE(eulered.Axis(Y,y)/dpr);  // Eulered latlon
			xE =-sin(lonE)*cos(latE);
			yE = cos(lonE)*cos(latE);
			zE = sin(latE);
		
// figure out x0, y0, z0 by doing Euler rotation -psi, -theta, -phi (i.e., backwards)
			double a(-psi), b(-theta), c(-phi);
			double x0, y0, z0;
			x0 = xE*cos(c)*cos(a) - xE*cos(b)*sin(a)*sin(c) + yE*cos(c)*sin(a)
					+ yE*cos(b)*cos(a)*sin(c) + zE*sin(c)*sin(b);  // this could be sped-up by precomputing sin & cos, perhaps
			y0 =-xE*sin(c)*cos(a) - xE*cos(b)*sin(a)*cos(c) - yE*sin(c)*sin(a)
					+ yE*cos(b)*cos(a)*cos(c) + zE*cos(c)*sin(b);
			z0 = xE*sin(b)*sin(a) - yE*sin(b)*cos(a) + zE*cos(b);
			
// convert x0, y0, z0 back into lat0, lon0 for interpolation from original cube
			float lat0, lon0;
			lat0 = asin(z0)*dpr;
			lon0 = atan2(-x0,y0)*dpr;*/
			
			pair<double, double> origloc;
			origloc = eulerangles(eulered.Axis(X,x), eulered.Axis(Y,y),
					psi, theta, phi);
			
			double &lat0(origloc.first), &lon0(origloc.second);
			while (lat0 >  90.) lat0 -= 180.;
			while (lat0 < -90.) lat0 += 180.;
//			cout << "Was going to " << lon0 << ", but now going to ";
			lon0 = (*lonconvention)(lon0);
//			cout << lon0 << "\n";
//			while (lon0 > 180.) lon0 -= 360.;
//			while (lon0 <-180.) lon0 += 360.;
			
// loop over Z in the cube and assign all the proper pixels into the new Eulered cube
			for (int z=0;z<N(Z);z++) {
				eulered(x,y,z) = (*this)(lon0, lat0, z, m);
			}
		}
	}
	if (osuppress-1==0) { printpercent(1,1); cout << '\n'; }
	
	osuppress--;
	return eulered;
	
}

vector<double> cube::eulerangles(vector<double> in, vector<angle> e) 
// 2012 November 29
// (1) phi -- clockwise rotation around z-axis, as viewed with z-axis toward you
// (2) theta -- clockwise rotation around x-axis, as viewed with x-axis toward you
// (3) psi -- clockwise rotation around new z-axis
{
//	cout << "in eulerangles(vector)\n";
	double &xin(in.at(0)), &yin(in.at(1)), &zin(in.at(2));    // input coords
	angle &phi(e.at(0)), &theta(e.at(1)), &psi(e.at(2));    // input coords
	
	vector<double> answer(3);
	double &x(answer.at(0)), &y(answer.at(1)), &z(answer.at(2));
	
//	cout << "variables linked\n";
	x = xin*cos(psi)*cos(phi) - xin*cos(theta)*sin(phi)*sin(psi) + yin*cos(psi)*sin(phi)
			+ yin*cos(theta)*cos(phi)*sin(psi) + zin*sin(psi)*sin(theta);  // this could be sped-up by precomputing sin & cos, perhaps
	y =-xin*sin(psi)*cos(phi) - xin*cos(theta)*sin(phi)*cos(psi) - yin*sin(psi)*sin(phi)
			+ yin*cos(theta)*cos(phi)*cos(psi) + zin*cos(psi)*sin(theta);
	z = xin*sin(theta)*sin(phi) - yin*sin(theta)*cos(phi) + zin*cos(theta);
	
//	cout << "returning\n";
	return answer;
}


vector<angle> cube::eulerangles(vector<double> v1, vector<double> v2) 
// 2012 November 30
// this returns the euler angle transformation that would transform v1 into v2.
{
//	if (osuppress < 1) cout << "in vector<angles> eulerangles(vector,vector(" << v2.size() << "))\n";
	if (v2.size()!=3) {
		v2.resize(3);
		v2.at(0)=0.; v2.at(1)=0.; v2.at(2)=0.;
		cout << "resized v2 input vector\n";
	}
	
	double &x1(v1.at(0)), &y1(v1.at(1)), &z1(v1.at(2));    // input coords
	double &x2(v2.at(0)), &y2(v2.at(1)), &z2(v2.at(2));    // input coords
	
	vector<angle> answer;
	answer.push_back(angle(0.,angle::RADIANS));
	answer.push_back(angle(0.,angle::RADIANS));
	answer.push_back(angle(0.,angle::RADIANS));
	angle &phi(answer.at(0)), &theta(answer.at(1)), &psi(answer.at(2));
	
	if (!(x1==0. && y1==0.)) 
		phi = -angle_atan2(x1,y1);
	else phi=angle(0.,angle::RADIANS);
//	cout << "phi=" << phi.as(angle::DEGREES) << ", ";
	
	angle theta1(angle(90.,angle::DEGREES)-angle_atan2(z1,sqrt(x1*x1+y1*y1)));
//	cout << "theta1=" << theta1.as(angle::DEGREES) << ", ";
	angle theta2(angle(90.,angle::DEGREES)-angle_atan2(z2,sqrt(x2*x2+y2*y2)));
//	cout << "theta2=" << theta2.as(angle::DEGREES) << ", ";
	theta = theta2-theta1;
//	cout << "theta=" << theta.as(angle::DEGREES) << ", ";

	if (!(x2==0. && y2==0.)) 
		psi = angle_atan2(x2,y2);
	else psi=angle(0.,angle::RADIANS);
//	cout << "psi=" << psi.as(angle::DEGREES) << ", ";
	
//	cout << "returning\n";
	return answer;
}

vector<double> cube::eulerangles_inverse(vector<double> in, vector<angle> e) 
// 2012 December 3
// this performs the euler angle transformation that transforms v2 back into v1.
{
	vector<angle> e2(e);
	e2.at(0)=-e.at(2);
	e2.at(1)=-e.at(1);
	e2.at(2)=-e.at(0);
	return eulerangles(in, e2);
}
	
	
pair<double, double> cube::eulerangles(double lonE, double latE, double psi,
		double theta, double phi) const
// input and output in degrees.
{
	pair<double, double> answer;
	
//	first convert from answer's lat/lon to x,y,z with unit radius
	vector<double> eulered(3);
	double &xE(eulered.at(0)), &yE(eulered.at(1)), &zE(eulered.at(2));    // Eulered coords
	lonE /= dpr;
	latE /= dpr;  // Eulered latlon
	xE =-sin(lonE)*cos(latE);
	yE = cos(lonE)*cos(latE);
	zE = sin(latE);
		
// figure out x0, y0, z0 by doing Euler rotation -psi, -theta, -phi (i.e., backwards)
	double a(-psi), b(-theta), c(-phi);
	vector<double> original(3);
	double &x0(original.at(0)), &y0(original.at(1)), &z0(original.at(2));
	x0 = xE*cos(c)*cos(a) - xE*cos(b)*sin(a)*sin(c) + yE*cos(c)*sin(a)
			+ yE*cos(b)*cos(a)*sin(c) + zE*sin(c)*sin(b);  // this could be sped-up by precomputing sin & cos, perhaps
	y0 =-xE*sin(c)*cos(a) - xE*cos(b)*sin(a)*cos(c) - yE*sin(c)*sin(a)
			+ yE*cos(b)*cos(a)*cos(c) + zE*cos(c)*sin(b);
	z0 = xE*sin(b)*sin(a) - yE*sin(b)*cos(a) + zE*cos(b);
			
// convert x0, y0, z0 back into lat0, lon0 for interpolation from original cube	
	double& lat0(answer.first), &lon0(answer.second);
	lat0 = asin(z0)*dpr;
	lon0 = atan2(-x0,y0)*dpr;
	
	return answer;
	
}

cube cube::orthoproject(int xsize, int ysize, int centerx, int centery, 
		float radius, methodtype m, float Nangle, float subSClat, float subSClon, projection* Projection) const
/* Created 2005 April 30 for Titan brightspot possible image release.
	needs the size of the output image (xsize, ysize), the center of the
	projected disk (centerx, centery), the radius of the planet in pixels
	(radius), how to interpolate (m), which direction to put as North in the
	resulting image (Nangle), and the sub-spacecraft lat/lon.
*/
{
	Nangle = -Nangle; // rotate clockwise as positive
	subSClat = -subSClat; // incoming latitude is in North Positive!
	int memz(1);
	if (xsize*ysize*N(Z)*4 < maxmemsize) memz=3;
	cube answer(xsize, ysize, N(Z), 0., cubetype, X, Y, Z, memz);
	answer.addheaderfrom(*this);
	answer.addazprojectionhdr("ORTHO",subSClat,subSClon,centerx,centery,radius,Nangle,xsize,ysize);
	
	
	
	for (int i=0;i<N(Z);i++)
		answer.Axis(Z,i) = Axis(Z,i);
	
	float rotx,roty,rho,c,lat,lon;
	int relx,rely;
	if (!osuppress) cout << "Orthoprojecting --  00%";
	for (int x=0;x<xsize;x++) {
		if (!osuppress) printpercent(x, xsize-1);
		for (int y=0;y<ysize;y++) {
			relx=x-centerx; rely=y-centery;
			rho = sqrtf(relx*relx+rely*rely);
			if (Projection==0) {  // orthographic case
				if (rho < radius) {
					rotx = cos(Nangle)*float(relx)+sin(Nangle)*float(rely);
					roty = -sin(Nangle)*float(relx)+cos(Nangle)*float(rely);
					
					// angular distance from subSC point
					c = asinf(rho/radius);
						
					lat = asinf(cos(c)*sin(subSClat)+(roty*sin(c)*cos(subSClat)/rho));
					
					if (subSClat > 3.1415926535897932385/2.)
						lon = subSClon + atan2f(rotx, -roty);
					else if (subSClat <-3.1415897932385/2.)
						lon = subSClon + atan2f(rotx,  roty);
					else
						lon = subSClon + atan2f(rotx*sin(c), 
								(rho*cos(subSClat)*cos(c)-roty*sin(subSClat)*sin(c)));
					lat =  lat * 180. / 3.1415926535897932384; // into cyl format
					lon =  lon * 180. / 3.1415926535897932384;
					while (lon < 180.) lon += 360.;
					while (lon > 180.) lon -= 360.;
							
					// correction for center pixel
					if (x == centerx  &&  y == centery) {
						lat = subSClat* 180. / 3.1415926535;
						lon = subSClon* 180. / 3.1415926535;
					}
					
					if (lon < (*this).Axis(X,0)) lon += 360.;
					if (lon > (*this).Axis(X,N(X)-1)) lon -= 360.;
					
					for (int z=0;z<N(Z);z++) {
						answer(x,y,z) = (*this)(lon, lat, z, m);
					}
					
					
				}
			} else { // non orthographic -- do numerical back-calculation of lat, lon as a function of given x, y
				
			}
		}
		
	}
	if (!osuppress) printpercent(xsize, xsize);
	if (!osuppress) cout << "\n";
	
	return answer;
}


	
cube cube::azimuthalstereographicproject(int xsize, int ysize, int centerx, int centery, float radius, 
	methodtype m, float Nangle, float subSClat, float subSClon) 
/* Created 2016 January (smack). Takes same inputs as orthoproject, but produces 
	an azimuthal stereographic projection, which by nature conserves angles/shapes.
	Name is nuts, I know.
*/
{
	Nangle = -Nangle; // rotate clockwise as positive
	subSClat = -subSClat; // incoming latitude is in North Positive!
	int memz(1);
	if (xsize*ysize*N(Z)*4 < maxmemsize) memz=3;
	
	cube answer(xsize, ysize, N(Z), 0., cubetype, X,Y,Z,memz);
	answer.addheaderfrom(*this);
	answer.addazprojectionhdr("AZSTEREOGR",subSClat,subSClon,centerx,centery,radius,Nangle,xsize,ysize);
		
	for (int i=0;i<N(Z);i++) 
		answer.Axis(Z,i) = Axis(Z,i);
	
	float rotx,roty,rho,c,lat,lon;
	int relx,rely;
			
	for (int x=0; x<xsize;x++){
		for (int y=0; y<ysize;y++){
				
			relx=x-centerx; 
			rely=y-centery;
			
			rho= sqrtf(relx*relx+rely*rely);
			
			if(rho<radius){
					
				c= 2*atanf(rho/(2*radius));
				
				rotx = cos(Nangle)*float(relx)+sin(Nangle)*float(rely);
				roty = -sin(Nangle)*float(relx)+cos(Nangle)*float(rely);
			
				lat= asinf(cos(c)*sin(subSClat)+roty*sin(c)*cos(subSClat)/rho);
				lon= subSClon+ atanf(rotx*sin(c)/(rho*cos(subSClat)*cos(c)-roty*sin(subSClat)*sin(c)));
						
				lat =  lat * 180. / 3.1415926535897932384; // into degrees, for cyl ingestion
				lon =  lon * 180. / 3.1415926535897932384;
				
				if (y<=centery) lon=lon+180;
				
				if (x == centerx  &&  y == centery) {
					lat = subSClat* 180. / 3.1415926535;
					lon = subSClon* 180. / 3.1415926535;
				}						
				
			
				if (lon < (*this).Axis(X,0)) lon += 360.;
				if (lon > (*this).Axis(X,N(X)-1)) lon -= 360.;
				
				for (int z=0;z<N(Z);z++) {
					answer(x,y,z) = (*this)(lon, lat, z, m);
				}	
			}
				
		}	
	}			
			
	return answer;
	
}	


cube cube::lambertazimuthalproject(int xsize, int ysize, int centerx, int centery, float radius, 
	methodtype m, float Nangle, float subSClat, float subSClon)
/* Created 2016 January (smack). Takes same inputs as orthoproject, but produces 
	an lambert azimuthal equal area projection, which, as you probably guessed, 
	conserves area.
*/
{	
	Nangle = -Nangle; // rotate clockwise as positive
	subSClat = -subSClat; // incoming latitude is in North Positive!
	int memz(1);
	if (xsize*ysize*N(Z)*4 < maxmemsize) memz=3;
	
	cube answer(xsize, ysize, N(Z), 0., cubetype, X,Y,Z,memz);
	answer.addheaderfrom(*this);
	answer.addazprojectionhdr("LAMBAZ",subSClat,subSClon,centerx,centery,radius,Nangle,xsize,ysize);
	
	for (int i=0;i<N(Z);i++) 
		answer.Axis(Z,i) = Axis(Z,i);
		
	
	float rotx,roty,rho,c,lat,lon;
	int relx,rely;
			
	for (int x=0; x<xsize;x++){
		for (int y=0; y<ysize;y++){
			relx=x-centerx; 
			rely=y-centery;
			
			rho= sqrtf(relx*relx+rely*rely);
			
			if(rho<radius){
					
				c= 2*asin(rho/(2*radius));
				
				rotx = cos(Nangle)*float(relx)+sin(Nangle)*float(rely);
				roty = -sin(Nangle)*float(relx)+cos(Nangle)*float(rely);
				
				lat= asinf(cos(c)*sin(subSClat)+roty*sin(c)*cos(subSClat)/rho);
				
				if (subSClat > 3.1415926535897932385/2.) //i.e. if center lat == +90 w/ float precision
						lon = subSClon + atan2f(rotx, -roty);
				else if (subSClat <-3.1415897932385/2.) //i.e. if center lat == -90 w/ float precision
						lon = subSClon + atan2f(rotx,  roty);
				else
				 	lon= subSClon+ atanf(rotx*sin(c)/
							(rho*cos(subSClat)*cos(c)-roty*sin(subSClat)*sin(c)));
						
				lat =  lat * 180. / 3.1415926535897932384; // into degrees, for cyl ingestion
				lon =  lon * 180. / 3.1415926535897932384;
				
				
				if (x == centerx  &&  y == centery) {
					lat = subSClat* 180. / 3.1415926535;
					lon = subSClon* 180. / 3.1415926535;
				}						
				
			
				if (lon < (*this).Axis(X,0)) lon += 360.;
				if (lon > (*this).Axis(X,N(X)-1)) lon -= 360.;
				
				
				for (int z=0;z<N(Z);z++) {
					answer(x,y,z) = (*this)(lon, lat, z, m);
				}	
			}
				
		}	
	}		
	
	
	return answer;

}


	    
		 


void fourier(cube hr, cube hi, cube *Hr, cube *Hi, int inv)
{
	cout << "performing fourier transform\n";
//   Creating the array for NR to use
	float H[hr.Nz()+1][hr.Ny()+1][2*hr.Nx()+1];
	for (int z=0;z<hr.Nz();z++){
		for (int y=0;y<hr.Ny();y++){
			for (int x=0;x<hr.Nx();x++){
				H[z+1][y+1][2*x+1] = hr(x,y,z);
				H[z+1][y+1][2*x+2] = hi(x,y,z);
			}
		}
	}
	
	unsigned long nn[4];
	nn[0]=0;
	nn[1]=hr.Nz();
	nn[2]=hr.Ny();
	nn[3]=hr.Nx();
	
	float *pH;
	pH = (float*) H;
	cout << "calling fourn\n";
	fourn(pH, nn, 3, inv);
	
	for (int z=0;z<hr.Nz();z++){
		for (int y=0;y<hr.Ny();y++){
			for (int x=0;x<hr.Nx();x++){
				(*Hr)(x,y,z) = H[z+1][y+1][2*x+1];
				(*Hi)(x,y,z) = H[z+1][y+1][2*x+2];
			}
		}
	}
	return;
}

cube cube::noise(noisetype nt, float sigma) const
{
	cube noisy(*this);
	if (!idum) raninit();
	if (nt==GAUSSIAN)
		for (int x=0;x<nx;x++) {
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++) {
					float thisstddev;
					if (sigma == -1) thisstddev = sqrt((*this)(x,y,z));
					else thisstddev = sigma;
					noisy(x,y,z) = thisstddev*gasdev(&idum);
				}	
			}
	}		
	else if (nt==POISSON){
		for (int x=0;x<nx;x++){
			for (int y=0;y<ny;y++){
				for (int z=0;z<nz;z++){					
					float thisstddev;
					if (sigma == -1) thisstddev = sqrt((*this)(x,y,z));
					else thisstddev = sigma;
					// if this point is small enough use poisson, if its big (and will make
					// poidev crash) use data+gaussian distribution with stddev=sqrt(x), 
					// which poisson(x) approximates for large x
					if (thisstddev < 3e3)
						noisy(x,y,z) = dpoidev((double) thisstddev, &idum);
					else 
						noisy(x,y,z) = (*this)(x,y,z) + thisstddev*gasdev(&idum);
				}
			}
		}
	} else cout << "Noise generation of type " << nt << " has not yet been implimented\n";
	
	return noisy;
}

cube blackbody(cube templatecube, float T, datatype D)
{
	cube bb(templatecube);
	
	if (!osuppress) cout << "Generating T="<<T<<" blackbody curve --    00%";
	
	float c1=3.741774e8;  // in microns^4*kg/s^3
	float c2=14387.69;  // in microns * Kelvin
	for (int x=0;x<bb.Nx();x++)
		for (int y=0;y<bb.Ny();y++)
			for (int z=0;z<bb.Nz();z++){
				if (!osuppress) printpercent(bb.Nx()*bb.Ny()*z, bb.Nx()*bb.Ny()*bb.Nz());
				bb(x,y,z)=c1/pow(bb.Zaxis(z),5)/(expf(c2/(bb.Zaxis(z)*T))-1);
				if (D == WM2) {
					float range;
					if (z == 0)
						range = bb.Axis(Z,1) - bb.Axis(Z,0);
					else if (z == bb.N(Z)-1)
						range = bb.Axis(Z, bb.N(Z)-1) - bb.Axis(Z, bb.N(Z)-2);
					else
						range = (bb.Axis(Z, z+1) - bb.Axis(Z, z-1))/2.;
					bb(x,y,z) *= range;
				}
			}
	bb.Dtype(D);
	if (!osuppress) cout << '\n';
	return bb.stripneg();
}

cube photonblackbody(cube templatecube, float T, datatype D)
{
	cube bb(blackbody(templatecube, T, D));
	
	for (int z=0;z<bb.N(Z);z++) {
		double divfactor=1.98644521e-25/bb.Axis(Z,z);
		for (int y=0;y<bb.N(Y);y++) {
			for (int x=0;x<bb.N(X);x++) {
				bb(x,y,z) /= divfactor;
			}
		}
	}
	
	return bb;
}
		
float cube::min() const
/* finds minimum value in cube
	11/14/2k JB */
// altered to call mincoord() 2014 June 3
{
/*	float m;
	m=(*this)(0,0,0);
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				if ((*this)(x,y,z)<m) m=(*this)(x,y,z);
	return m;*/
	return mincoord().first;
}

float cube::max() const
{
	return maxcoord().first;
}

pair<float, vector<int> > cube::maxcoord() const
{
	float m=(*this)(0,0,0);
	int maxx(0),maxy(0),maxz(0);
	for (int x=0;x<N(X);x++)
		for (int y=0;y<N(Y);y++)
			for (int z=0;z<N(Z);z++) 
				if ((*this)(x,y,z) > m) {
					m=(*this)(x,y,z);
					maxx=x;
					maxy=y;
					maxz=z;
				}
	
	pair<float, vector<int> > answer;
	
	answer.first = m;
	answer.second.resize(3);
	answer.second[0] = maxx;
	answer.second[1] = maxy;
	answer.second[2] = maxz;
	
	return answer;
}

pair<float, vector<int> > cube::mincoord() const
// created 2014 June 3 JWB based on previously existing maxcoord()
{
	float m=(*this)(0,0,0);
	int minx(0),miny(0),minz(0);
	for (int x=0;x<N(X);x++)
		for (int y=0;y<N(Y);y++)
			for (int z=0;z<N(Z);z++) 
				if ((*this)(x,y,z) < m) {
					m=(*this)(x,y,z);
					minx=x;
					miny=y;
					minz=z;
				}
	
	pair<float, vector<int> > answer;
	
	answer.first = m;
	answer.second.resize(3);
	answer.second[0] = minx;
	answer.second[1] = miny;
	answer.second[2] = minz;
	
	return answer;
}

cube cube::max(axis a)
/* created 3/10/2k but already due for obsolescence in Jcube2.0 so just basic */
{
	cube m(1,1,nz);
	for (int z=0;z<nz;z++) { m(0,0,z)=(*this)(0,0,z); m.Axis(Z,z)=0.; }
	loopoverall if ((*this)(x,y,z) > m(0,0,z)) {
		m(0,0,z)=(*this)(x,y,z);
		if (a==X) m.Axis(Z,z)=x;
		if (a==Y) m.Axis(Z,z)=y;
	}
	return m;
}

double cube::maxaxis(axis A) const
// added 2010 March 1 for jroberts by jbarnes
{
	double answer(Axis(A,0));
	
	for (int i(0);i<N(A);i++)
		if (Axis(A,i) > answer) answer = Axis(A,i);
	
	return answer;
}

double cube::minaxis(axis A) const
// added 2010 March 1 for jroberts by jbarnes
{
	double answer(Axis(A,0));
	
	for (int i(0);i<N(A);i++)
		if (Axis(A,i) < answer) answer = Axis(A,i);
	
	return answer;
}

double cube::meanaxis(axis A) const
// added 2010 March 1 by jbarnes
{
	double answer(Axis(A,0));
	
	for (int i(1);i<N(A);i++)
		answer += Axis(A,i);

	return answer/double(N(A));
}


cube cube::centroid() const
// Written 2010 January 22 JWB -- Jcube is 10 years old!  That shit is scary.
{
	cube answer(1,1,1,0.);
	
	double xx(0.), yy(0.), zz(0.), totalsum(sum());
	
	for (int x(0);x<N(X);x++)
		for (int y(0);y<N(Y);y++)
			for (int z(0);z<N(Z);z++) {
				xx += Axis(X,x)*(*this)(x,y,z);
				yy += Axis(Y,y)*(*this)(x,y,z);
				zz += Axis(Z,z)*(*this)(x,y,z);
			}
	
	answer.Axis(X,0) = xx / totalsum;
	answer.Axis(Y,0) = yy / totalsum;
	answer.Axis(Z,0) = zz / totalsum;

	return answer;
}


cube cube::abscontinuum(int interval)
/*   Overly simple program to generate the continuum of an absorption spectrum by taking
the maximum of every (interval) pixels and interpolating between them.  Assumes wavelength
data is in Z dimension for now
	should barf big-time if fed heavily negative data
	Created 2/22/2k
*/
{
	if (Ztype()==ARRAYVAL) { cout << "Z axis not wavelength in abscontinuum\n"; }	
	
// take the maximum over each interval range and put that point in the continuum bin
	cube maxes(Nx(), Ny(), Nz()/interval);
	maxes.hdr=hdr;
	maxes.keyword("Nz", int2str(maxes.Nz()) );
	for (int cz=0;cz<maxes.Nz();cz++) {
		printpercent(cz,maxes.Nz()-1);
		for (int z=cz*interval; z < FMIN((cz+1)*interval,Nz());z++) {
			if ( fabs((*this)(0,0,z)) > maxes(0,0,cz) ) {
				maxes(0,0,cz)   = (*this)(0,0,z);
				maxes.Zaxis(cz) = Zaxis(z);
			}
		}
	}
	
	return maxes;
}
			
cube cube::dopplercorrelate(cube& g, double vmax, double res)
/* Created 2/24/2k
	Determines correlation as a function of doppler velocities between +vmax and -vmax
	at resolution res (all in km/s). Returns cube of correlation as a function of v.
	*/
{
	cube& f=*this;			// alias *this so I don't have to use it all the time.
	
	osuppress++;
	if (Ztype() > 4  || Ztype()==0  || g.Ztype() > 4 || g.Ztype() == 0)
		cout << "ERROR -- Dopplercorrelate called with incorrect Z axis type -- Must be wavelength";
	if (f.Nx()!=g.Nx() || f.Ny()!=g.Ny())
		cout << "ERROR -- dopplercorrelate called with incompatible cubes (nx,ny different)";
	f=f.dirinc();
	g=g.dirinc();
	
	if (!(osuppress-1)) cout << "Computing dopplercorrelation.   ";

	
// determine overlap region to use when integrating up.  Region = min..max
	cube temp;
	double min=f.Zaxis(0), max=f.Zaxis(f.Nz()-1);
	temp=g.dopplershift(vmax);
	if (temp.Zaxis(0) > min) min = temp.Zaxis(0);
	if (temp.Zaxis(temp.Nz()-1) < max) max = temp.Zaxis(temp.Nz()-1);
	temp=g.dopplershift(-vmax);
	if (temp.Zaxis(0) > min) min = temp.Zaxis(0);
	if (temp.Zaxis(temp.Nz()-1) < max) max = temp.Zaxis(temp.Nz()-1);
	
// create the cube that will store the results of the multiplication at each v
// use resolution of f and range of min..max
	cube mult;
	
// the big loop over radial velocities
	cube corr(nx, ny, (int)(2*vmax/res));
	corr.Ztype(KMPERS);
	cube fsm, gsm;
	fsm=f(-1,-1,-1,-1,min,max);
	for (int i=0;i<corr.Nz();i++) {
		if (!(osuppress-1)) printpercent(i,corr.Nz()-1);
		float v = -vmax + i*res;
		corr.Zaxis(i) = v;
		gsm=g.dopplershift(v);
		gsm=gsm.resample(fsm);		
		mult = fsm * gsm;
		mult = mult.plane(Z);
		corr.insert(mult, 0, 0, i);
	}
	
	osuppress--;
	return corr;
}
		
cube cube::circle(float cx, float cy, float r, float e)
/* created 3/6/2k
	Just makes a filled in circle centered at cx, cy of radius r in each XY plane of the
	template cube 
	e is the planetary oblateness, (Req-Rpol)/Req, r is the polar radius.*/
{
	cube c;
	c=*this;
	c*=0;
	
	e=(1.-e);
	
	if (!osuppress) cout << "Creating circle ";
	if (nx*ny*nz<1e6) osuppress++;
	if (!osuppress) cout << "--  00%";
	loopoverall {
		if (!osuppress) printpercent(x*ny*nz+y*nz+z, nx*ny*nz);
		if (sqrt((x-cx)*(x-cx)*e*e+(y-cy)*(y-cy)) < r) c(x,y,z)=1.;
	}
	char b=8;
	if (!osuppress) cout << b<<b<<b<<b<<"100%\n";
	if (nx*ny*nz<1e6) osuppress--;
	return c;
}

cube cube::realFFT(int isign)
/* pads with zeroes and then runs FFT 
created 11/5/2k*/
{
	int testx;
	testx=134;
	cube blanx((*this)*0);
	cube padded;
	padded=(*this).blocksz(blanx);
	cube result;
	cout << "Performing real fourier transform:\n";
	result=padded.FFT(isign);
	return result;
}

cube cube::FFT(int isign)
/* A pretty slow and dirty fourier transform routine.  Takes in a 3-d cube with
	dims all multiples of 2 (If they aren't, the last value is padded w/0).  
	The 1st half (in z) is suppsed to be real,and the second imaginary.
	isign=1  foward transform
	isign=-1 inverse transform
written 11/5/2k  

	output cube(N(X),N(Y), 2N(Z)); first z plane= reals, second z plane= imaginaries		
		
future work:  Change this to a Jcube native routine without the conversion to NR
arrays overhead
*/
{
	int bx, by, bz, dimN, ictr;
	dimN=3; ictr=1;
	long unsigned int nn[4];
	if ((nx%2) && nx!=1) cout << "FFT called with odd array size, crash.\n";
	if ((ny%2) && ny!=1) cout << "FFT called with odd array size, crash.\n";
	if ((nz%4) && nz!=2) cout << "Total errorage, FFT called with (nz%2)!=0.  Need imaginary components!\n";
	
	bx=nx;
	by=ny; 
	bz=nz/2; 
	
	if (nx==1) dimN--; else { 
		nn[ictr]=nx; 
		if (!is2pow(nx)) cout << "X dimension not power of 2!!!\n"; 
		ictr++; }
	if (ny==1) dimN--; else { 
		nn[ictr]=ny; 
		if (!is2pow(ny)) cout << "Y dimension not power of 2!!!\n"; 
		ictr++; }
	if (nz==2) dimN--; else {
		nn[ictr]=nz/2; 
		if (!is2pow(nz)) cout << "Z dimension not power of 2!!!\n"; 
		ictr++; }
	
	vector<float> nrdata((int)(nx*ny*nz+1));
	vector<float>::iterator nrptr;
	nrptr=nrdata.begin();
	nrptr++;

	for (int z=0;z<nz/2;z++) {
		for (int y=0;y<ny;y++)  {
			for (int x=0;x<nx;x++) {
				nrdata[2*(nx*ny*z+nx*y+x)+2]=(*this)(x,y,z);
				nrptr++;
				nrdata[2*(nx*ny*z+nx*y+x)+1+2]=(*this)(x,y,z+nz/2);
				nrptr++;
			}
		}
	}
	
	cout << "Performing " << dimN << " dimensional FFT on data dims("<<nx<<","<<ny<<","<<nz<<")\n";
	jfourn(nrdata, nn, dimN, isign);
	
	nrptr=nrdata.begin();
	nrptr++;
	cube result(nx, ny, nz);
	for (int z=0;z<nz/2;z++) {
		for (int y=0;y<ny;y++)  {
			for (int x=0;x<nx;x++) {
				result(x,y,z)=nrdata[2*(nx*ny*z+nx*y+x)+2];
				nrptr++;
				result(x,y,z+nz/2)=nrdata[2*(nx*ny*z+nx*y+x)+1+2];
				nrptr++;
				if (z<=nz/4) {
					result.Axis(Z,z) = z;
					result.Axis(Z,z+nz/2) = z+1;
				} else {
					result.Axis(Z,z) = -(nz/2-z);
					result.Axis(Z,z+nz/2) = -(nz/2-z)+1;
					
				}
			}
		}
	}
	
	if (isign==1) result/=nx*ny*nz/2;
//	cout << "dividing by " << nx*ny*nz/2 << "\n";
//	result.copyaxis(*this, X);
//	result.copyaxis(*this, Y);
//	result.copyaxis(*this, Z);
	return result;
}

cube cube::wavelet(int direction, int ncoeff)
// wavelet driver, based on wt1 from NR
// ncoeff should be 4, 12, or 20.
// only works in 1 dimension for now -- if the other N(A)'s are not 1, it will NOT work!
{
	cube answer(*this);

	if (N(Z) < 4) return answer;
	
	pwtset(ncoeff);
	
	if (direction>=0)
		for (unsigned long nn(N());nn>=4;nn>>=1) pwt(&(*(answer.Databegin())),nn,direction);
	else
		for (unsigned long nn(4);nn<=N();nn<<=1) pwt(&(*(answer.Databegin())),nn,direction);
	
	return answer;
}


cube crap;

void funcs(float xval, float afunc[], int ma){ afunc[1]=crap(0,0,(int)xval); }

cube cube::templatefit(cube templatecube, bool sqrtnoise)
/* Fits *this with the template as a least-squares thing. 
	Created 11/14/2k JB
*/
{
	float k;
	vector<float> x(nx*ny*nz+1), y(nx*ny*nz+1), s(nx*ny*nz+1), a(2);
	vector<int> ia(2);
	float **covar, chisq;
	
	covar=matrix(1,nx*ny*nz,1, nx*ny*nz);
	
	int i=1;
	for (int xi=0;xi<nx;xi++) {
		for (int yi=0;yi<ny;yi++) {
			for (int zi=0;zi<nz;zi++) {
				x[i]=i;
				y[i]=(*this)(xi,yi,zi);
				if (sqrtnoise) s[i]=sqrtf((*this)(xi,yi,zi));
				else s[i]=1.;
			}
		}
	}
	ia[1]=1; ia[0]=0;
	a[1]=1;

	crap=templatecube;
	lfit(&x[0], &y[0], &s[0], nx*ny*nz, &a[0], &ia[0], 1, covar, &chisq, funcs);
	
	cout << "fit to a=" << a[1] << "\n";
	cube ans;
	ans=templatecube*a[1];
	return ans;
}
	
cube cube::HPfilter(int minfrq)
/* Fourier transforms cube, sets the minfrq lowest frequency components to 0, transforms
	back.
	9/14/2k JB 
	Kludgy and hurried, only does Z
	*/
{
	cube FTed;
	FTed=(*this)-(*this).cpolyfit(1);
	FTed=FTed.realFFT(1);
	for (int i=0;i<minfrq;i++) FTed(0,0,i)=0;
	FTed=FTed.FFT(-1);
	FTed=FTed(0,0,0,0,0,FTed.N(Z)/2-1);
	FTed.copyaxis(*this, Z);
	
	return FTed;
}


cube cube::abs()
/* Jason BArnes 10/26/2001
	returns absolute value of cube */
{
	cube answer(*this);
	for (cube::iterator i(answer.begin());i!=answer.end();i++)
		if ((*i) < 0.) (*i) *= -1.;
	return answer;
}


float cube::sum() const
/*  JB 10/26/2001
	returns sum of all values in cube */
{	
	float answer(0.);
//	for (cube::const_iterator i(begin());i!=end();i++)
//		answer += *i;
	for (int x=0;x<N(X);x++)
		for (int y=0;y<N(Y);y++)
			for (int z=0;z<N(Z);z++)
				answer += (*this)(x,y,z);
	return answer;
}

float cube::median(float f) 
{
	return median(f, UNK, UNK, UNK)(0,0,0);
}

cube cube::median(axis A, axis B, axis C) {return median(0.5, A, B, C);}

cube cube::median(float f, axis A, axis B, axis C)
/* JB 10/26/2001 */
// Severely modified to handle arbitrary dimension medianing, 2010 March 2
{
	if (A==UNK && B==UNK && C==UNK) {
		A=X;
		B=Y;
		C=Z;
	}
	
	cube answer;
	if (memdims==3) {
	// set up the answer cube
		long ax(N(X)), ay(N(Y)), az(N(Z));
		if (A==X || B==X || C==X) ax=1;
		if (A==Y || B==Y || C==Y) ay=1;
		if (A==Z || B==Z || C==Z) az=1;
		answer = cube(ax, ay, az);
		
	// finesse the axes for the answer cube
		if (A!=UNK) answer.Axis(A,0) = meanaxis(A);
		if (B!=UNK) answer.Axis(B,0) = meanaxis(B);
		if (C!=UNK) answer.Axis(C,0) = meanaxis(C);
		if (answer.N(X)>1) answer.copyaxis(*this, X);
		if (answer.N(Y)>1) answer.copyaxis(*this, Y);
		if (answer.N(Z)>1) answer.copyaxis(*this, Z);
		
	// figure out how many floats will have to be in each median
		long totalN;
		if (C != UNK) totalN = N();
		else if (B != UNK) totalN = N(A)*N(B);
		else totalN = N(A);
		
	// we'll each set of them all in this new vector, v.
		vector<float> v(totalN);
		
	// assign a median for each pixel in the answer cube.
		for (long x(0);x<answer.N(X);x++)
			for (long y(0);y<answer.N(Y);y++)
				for (long z(0);z<answer.N(Z);z++) {
				
				// assign the values in v
					long n(0);
					if (C!=UNK) // 3-d median 
						for (long ox(0);ox<N(X);ox++)
							for (long oy(0);oy<N(Y);oy++)
								for (long oz(0);oz<N(Z);oz++,n++)
									v.at(n) = (*this)(ox,oy,oz);
					else if (B!=UNK) // 2-d median
						for (long a(0);a<N(A);a++)
							for (long b(0);b<N(B);b++,n++) {
								axis A3(otheraxis(A,B));
								long a3;	
								if (A3==X) a3=x;
								if (A3==Y) a3=y;
								if (A3==Z) a3=z;
								v.at(n) = (*this)(A, a, B, b, A3, a3);
							}
					else if (A!=UNK) // it better not . . .  1-d median
						for (long a(0);a<N(A);a++,n++) {
							axis A1=Uax(A);
							axis A2=Dax(A);
							long a1, a2;
							if (A1==X) a1=x;
							if (A1==Y) a1=y;
							if (A1==Z) a1=z;
							if (A2==X) a2=x;
							if (A2==Y) a2=y;
							if (A2==Z) a2=z;
							v.at(n) = (*this)(A, a, A1, a1, A2, a2);
						}
								
			
			
				// sort and take the median
					sort(v.begin(), v.end());
					answer(x,y,z) = v.at(long((v.size()-1)*f));
				}
	} else {
		cout << "This cube is just too big to median in memory.\n";
	}
	return answer;
}

cube cube::sigma(axis A, axis B, axis C)
/* JB 10/26/2001 */
// Severely modified to handle arbitrary dimension medianing, 2010 March 2
{
	if (A==UNK && B==UNK && C==UNK) {
		A=X;
		B=Y;
		C=Z;
	}
	
	cube answer;
	if (memdims==3) {
	// set up the answer cube
		int ax(N(X)), ay(N(Y)), az(N(Z));
		if (A==X || B==X || C==X) ax=1;
		if (A==Y || B==Y || C==Y) ay=1;
		if (A==Z || B==Z || C==Z) az=1;
		answer = cube(ax, ay, az);
		
	// finesse the axes for the answer cube
		if (A!=UNK) answer.Axis(A,0) = meanaxis(A);
		if (B!=UNK) answer.Axis(B,0) = meanaxis(B);
		if (C!=UNK) answer.Axis(C,0) = meanaxis(C);
		if (answer.N(X)>1) answer.copyaxis(*this, X);
		if (answer.N(Y)>1) answer.copyaxis(*this, Y);
		if (answer.N(Z)>1) answer.copyaxis(*this, Z);
		
	// figure out how many floats will have to be in each median
		int totalN;
		if (C != UNK) totalN = N();
		else if (B != UNK) totalN = N(A)*N(B);
		else totalN = N(A);
		
	// we'll each set of them all in this new vector, v.
//		vector<float> v(totalN);
		
		if (osuppress<1) cout << "Calculating Standard Deviation:   00%";
	// assign a sigma for each pixel in the answer cube.
		for (int x(0);x<answer.N(X);x++)
			for (int y(0);y<answer.N(Y);y++)
				for (int z(0);z<answer.N(Z);z++) {
					if (osuppress<1) printpercent(z+N(Z)*y+N(Z)*N(Y)*x, answer.N());
				
				// assign the values in v
					cube subdata;
					if (C!=UNK) // 3-d median 
						subdata = *this;
					else if (B!=UNK) {// 2-d median
						axis A3(otheraxis(A,B));
						int a3;	
						if (A3==X) a3=x;
						if (A3==Y) a3=y;
						if (A3==Z) a3=z;
						osuppress++;
						subdata = plane(A3,a3);
						osuppress--;
					} else if (A!=UNK) {
						axis A1=Uax(A);
						axis A2=Dax(A);
						int a1, a2;
						if (A1==X) a1=x;
						if (A1==Y) a1=y;
						if (A1==Z) a1=z;
						if (A2==X) a2=x;
						if (A2==Y) a2=y;
						if (A2==Z) a2=z;
						subdata = skewer(A, a1, a2);
					}
								
				// calculate std dev
					float thismean(subdata.mean());
					subdata -= thismean;
					subdata*=subdata;
					float thissigma(subdata.sum());
					thissigma/=float(subdata.N());
					thissigma = sqrt(thissigma);
					answer(x,y,z) = thissigma;
					
				}
		if (osuppress<1) printpercent(answer.N(), answer.N());
		if (osuppress<1) cout << "\n";
	} else {
		cout << "This cube is just too big to median in memory.\n";
	}
	return answer;
}


float cube::mean() const
/* JB 10/26/2001 */
{
	float answer;
	answer=sum()/N();
	return answer;
}

/*
float cube::mode()
/* JB 10/26/2001 */
/*
{
	hash_map<float, unsigned long> occurancetable;
	for (cube::iterator i(begin()) ; i!=end() ; i++) {
		(occurancetable[*i])++;
	}
	
	unsigned long max=0;
	float maxval;
//	for (hash_map<float, unsigned long>::iterator i(occurancetable.begin());
//		  i!=occurancetable.end();i++) {
//		if (i->second > max) { max=i->second; maxval=i->first; }
//	}
	return maxval;
}
*/


string cube::info() const
/* JB 10/26/2001 */
{
	char answer[3000];
	ostrstream sout(&answer[0],3000);
	
	sout << (string)"Cube info for cube title " << hdr["title"] ;
	sout << (string)"type " << cubetype << '\n';
	sout << (string)"dimensions:  ("<<N(X)<<","<<N(Y)<<","<<N(Z)<<")\n";
	cout << "about to min/max\n"; cout.flush();
	sout << (string)"    min=" << min() << (string)"      max=" << max() << '\n';
	cout << "about to sum/mean\n"; cout.flush();
	sout << (string)"    sum=" << sum() << (string)"  mean=" << mean();
//	cout << "about to median\n"; cout.flush();
//	sout << (string)"    median="<<median() << (string)"  mode=" << /*mode() << */'\n';
	sout << '\0';
	
	cout << "infodone"; cout.flush();
	return (string)answer;
}


cube cube::generateaxis(axis A, float cmin, float cmax)
// Nov 2, 2001 JB
// returns the given cube with (axis) given values between min and max
{
	cube answer(*this);
	for (int a=0;a<N(A);a++) {
		answer.Axis(A,a) = cmin+(double)a*(cmax-cmin)/(double)N(A);
	}
	
	return answer;
}

cube cube::histogram(int nbins,float hmin, float hmax, cube weightcube)
/* 11/2/2001 JB   -- for unweighted
	Returns histogram with nbins 
	modified for weighting 2008 Jan 6 JWB.
	*/
{	
	cube histogram(1,1,nbins,0.0);
	float cmin((*this).min()), cmax((*this).max());
	if (hmin==-1.) hmin=cmin+((cmax-cmin)/(nbins-1))/2.;
	if (hmax==-1.) hmax=cmax-((cmax-cmin)/(nbins-1))/2.;
	float binsize=fabs(hmax-hmin)/(nbins-1);
	
	histogram=histogram.generateaxis(Z,hmin,hmax);
	
	cube::iterator w(weightcube.begin());
	for (cube::iterator i(begin());i!=end();i++) {
		double weight(1.);
		if (weightcube.N()==N()) weight = *w;
		for (int j=0;j<histogram.N(Z);j++) {
			if (*i  > histogram.Axis(Z,j)-binsize/2.   &&
				 *i  < histogram.Axis(Z,j)+binsize/2.) histogram(0,0,j)+=1.0*weight;
			if (*i == histogram.Axis(Z,j)-binsize/2.   ||
				 *i == histogram.Axis(Z,j)+binsize/2.) histogram(0,0,j)+=0.5*weight;
		}
		if (weightcube.N()==N()) w++;
	} 
	
	return histogram;
}

cube cube::cumulativehistogram(int nbins, float hmin, float hmax)
/* JB 11/2/2001
	Same as prev but cumulative */	
{
	cube answer;
	answer = histogram(nbins,hmin,hmax);
	
	for (int i=1;i<answer.N(Z);i++) answer(0,0,i) = answer(0,0,i)+answer(0,0,i-1);
	
	return answer;
	
}

cube cube::simplegradient()
/*smack Sept 2016 
	takes simple gradient in x,y
	calculates composite and composite direction */
{	
	cube answer(N(X),N(Y),1);
	
	
	// take gradient	
	
	double mean(0.);
	int count(0);
	cout <<"means:\n";
	for (int x(1); x<N(X)-1; x++){
		
		for (int y(1); y<N(Y)-1; y++){
			if ( (*this)(x-1,y,0)==-1e5 || (*this)(x,y-1,0) == -1e5 ||
					(*this)(x+1,y,0)==-1e5 || (*this)(x,y+1,0) ==-1e5 ){		
				answer(x,y,0)=0;
			}
			else{
				double delx(0.),dely(0.),theta(0.),magnitude(0.);	
				delx=((*this)(x+1,y,0)-(*this)(x-1,y,0))/2;
				dely=((*this)(x,y+1,0)-(*this)(x,y-1,0))/2;
				theta=atan(delx/dely);
				magnitude=sqrt(delx*delx+dely*dely);
				if (magnitude == magnitude){
					answer(x,y,0)=magnitude;
					mean+=magnitude;
					count++;
				}
			}
		}		
	}
	mean/=count;
	cout <<"\n"<<mean<<endl;
	double sigma(0.);
	
	for (int x(0); x<answer.N(X);x++){
		for (int y(0); y<answer.N(Y);y++){
			double thisvalue(answer(x,y,0)-mean);
			sigma+=thisvalue*thisvalue;
		}
	}		
	
	sigma/=count;
	sigma=sqrt(sigma);
	cout <<"\n\n--> FOUND SIGMA! "<<sigma<<endl;
	cout <<"--> mean: "<<mean<<endl<<endl;
	
	
	cube sigmacube(answer);
	sigmacube=answer.sigma(X,Y);
	cout <<"\n\n--> cube sigma: "<< sigmacube(0,0,0)<<endl<<endl;
	
	for (int x(0); x<answer.N(X);x++){
		for (int y(0); y<answer.N(Y);y++){
			
			if (answer(x,y,0) < mean+3*sigma )
				answer(x,y,0)=0.; 
			
		}				
	}
	
	
	// check all slopes, throw out values outside of 2sigma, make second plane the error estimate?
	
 	
	
	return answer;
	
}
