#include <string>
#include <map>
#include <stdlib.h>
#include "../../nr/nrutilj.h"

#include "Jcube.h"


cube& cube::operator=(const cube& r)
{
	if (debug) cout << "copying cube (operator=()) ID(" << &r << ") to ID(" << this << ")!\n";
	
	deallocate();
	
	if (this == &r) return *this;
	
	copymetadatafrom(r);       
	nx  = r.nx;
	ny  = r.ny;
	nz  = r.nz;

	memdims=r.memdims;	
	
	if (psidecube   != 0) {
		if (debug)cout << "deleting sidecube\n"; 
		delete psidecube; 
		psidecube=0;
	}
//	if (r.psidecube   != 0) psidecube   = new cube(r.sidecube());
	if (pbottomcube != 0) {
		if (debug)cout << "deleting bottomcube\n"; 
		delete pbottomcube; 
		pbottomcube=0;
	}
//	if (r.pbottomcube != 0) pbottomcube = new cube(r.bottomcube());
	if (pbackcube   != 0) {
		if (debug)cout << "deleting backcube\n"; 
		delete pbackcube; 
		pbackcube=0;}
//	if (r.pbackcube   != 0) pbackcube   = new cube(r.backcube());
	
	xaxis = new double[nx];
	yaxis = new double[ny];
	zaxis = new double[nz];
	for (int x=0;x<nx;x++) xaxis[x] = r.xaxis[x];
	for (int y=0;y<ny;y++) yaxis[y] = r.yaxis[y];
	for (int z=0;z<nz;z++) zaxis[z] = r.zaxis[z];
	
	// recursive.  Woah!
	if (r.Side()!=0) psidecube = new cube(r.sidecube());
	if (r.Bottom()!=0) pbottomcube = new cube (r.bottomcube());
	if (r.Back()!=0) pbackcube = new cube(r.backcube());
			
	createarrays(30.0);
			
	if (memdims==3)
		Data = r.Data;
	else {
		if (memdims<=2)
			if (r.current[storageSLOW]>=0 && r.current[storageSLOW]<N(storageSLOW)) r.pageout();
		else if (memdims<=1)
			if (r.current[storageMEDIUM]>=0 && r.current[storageMEDIUM]<N(storageMEDIUM)) r.pageout();
		else if (memdims==0)
			if (r.current[storageFAST]>=0 && r.current[storageFAST]<N(storageFAST)) r.pageout();
		
		
		cubefilename=string("Jcubetmp.")
				+ int2str((long long)(this))
				+ string(".Jcube");
		
		// I would rather do a r.write(cubefilename); here, but the compiler
		// bitches about const correctness, so I'm using this workaround
		string command("cp ");
		command += r.cubefilename;
		command += " ";
		command += cubefilename;
		cout << command << "\n";
		system(command.c_str());
		
		read(cubefilename,memdims);
		hdr = r.hdr;
		for (int x=0;x<nx;x++) Axis(X,x) = r.Axis(X,x);
		for (int y=0;y<ny;y++) Axis(Y,y) = r.Axis(Y,y);
		for (int z=0;z<nz;z++) Axis(Z,z) = r.Axis(Z,z);
		
		current[X]=-1;
		current[Y]=-1;
		current[Z]=-1;
	}
	
	return *this;
}

bool cube::operator==(cube r)
{
	if (nx != r.nx) return 0;
	if (ny != r.ny) return 0;
	if (nz != r.nz) return 0;
	for (int z=0;z<nz;z++){
		for (int y=0;y<ny;y++){
			for (int x=0;x<nx;x++){
				if ((*this)(x,y,z) != r(x,y,z)) return 0;
			}
		}
	}
	return 1;
}

bool cube::operator==(float r)
{
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				if ((*this)(x,y,z) != r)
					return 0;
	return 1;
}

float& cube::operator()(long x, int y, int z) {return operator()(x,long(y),long(z));}
float& cube::operator()(int x, long y, int z) {return operator()(long(x),y,long(z));}
float& cube::operator()(int x, int y, long z) {return operator()(long(x),long(y),z);}
float& cube::operator()(long x, long y, int z) {return operator()(x,y,long(z));}
float& cube::operator()(long x, int y, long z) {return operator()(x,long(y),z);}
float& cube::operator()(int x, long y, long z) {return operator()(long(x),y,z);}
float& cube::operator()(int x, int y, int z) {return operator()(long(x),long(y),long(z));}

float& cube::operator()(long x, int y, int z) const {return operator()(x,long(y),long(z));}
float& cube::operator()(int x, long y, int z) const {return operator()(long(x),y,long(z));}
float& cube::operator()(int x, int y, long z) const {return operator()(long(x),long(y),z);}
float& cube::operator()(long x, long y, int z) const {return operator()(x,y,long(z));}
float& cube::operator()(long x, int y, long z) const {return operator()(x,long(y),z);}
float& cube::operator()(int x, long y, long z) const {return operator()(long(x),y,z);}
float& cube::operator()(int x, int y, int z) const {return operator()(long(x),long(y),long(z));}

float& cube::operator()(long x, long y, long z)
{
	int s,m,f;	
	if      (storageSLOW==X) s=x;
	else if (storageSLOW==Y) s=y;
	else if (storageSLOW==Z) s=z;
	
	if      (storageMEDIUM==X) m=x;
	else if (storageMEDIUM==Y) m=y;
	else if (storageMEDIUM==Z) m=z;
	
	if      (storageFAST==X) f=x;
	else if (storageFAST==Y) f=y;
	else if (storageFAST==Z) f=z;

// have this version call the const version to return the relavant float
	return storageorderaccess(s,m,f);
}

float& cube::storageorderaccess(long s, long m, long f) 
{
	if (memdims != 3) {
			if (memdims==0 && ((current[storageFAST]>=0 && current[storageFAST]!=f) ||
									(current[storageMEDIUM]>=0 && current[storageMEDIUM]!=m) ||
									(current[storageSLOW]>=0 && current[storageSLOW]!=s))) 
				pageout();
			if (memdims==1 && ((current[storageMEDIUM]>=0 && current[storageMEDIUM]!=m) ||
									(current[storageSLOW]>=0 && current[storageSLOW]!=s))) 
				pageout();
			if (memdims==2 && current[storageSLOW]>=0 && current[storageSLOW]!=s) 
				pageout();
		}
	// have this version call the const version to return the relevant float
	return const_cast<float&>(static_cast<const cube&>(*this).storageorderaccess(s,m,f));
}

float& cube::operator()(long x, long y, long z) const
{
	if (x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz) {
		int s,m,f;  // slow, medium, and fast
		
		if      (storageSLOW==X) s=x;
		else if (storageSLOW==Y) s=y;
		else if (storageSLOW==Z) s=z;
		else cout << "No storageSLOW axis!  Internal error.\n";
		
		if      (storageMEDIUM==X) m=x;
		else if (storageMEDIUM==Y) m=y;
		else if (storageMEDIUM==Z) m=z;
		else cout << "No storageMEDIUM axis!  Internal error.\n";
		
		if      (storageFAST==X) f=x;
		else if (storageFAST==Y) f=y;
		else if (storageFAST==Z) f=z;
		else cout << "No storageFAST axis!  Internal error.\n";	
		
		return storageorderaccess(s,m,f);
	}
	else{
		cout << "OUT OF RANGE ERROR:  ("<<x<<", "<<y<<", "<<z<<") requested";
		cout << "FROM CUBE DIM ("<<nx<<", "<<ny<<", "<<nz<<")\n";
	}		
		
}

float& cube::storageorderaccess(long s, long m, long f) const
// split off from operator() for faster runtime operation, 2006/04/05 JB
{		
	if (Boundschecking) {
		if (!(s>=0 && m>=0 && f>=0  &&  s<N(storageSLOW) && m<N(storageMEDIUM) && f<N(storageFAST))) {
			cout << "OUT OF RANGE ERROR:  ("<<s<<storageSLOW<<", "<<m<<storageMEDIUM<<", "<<f<<storageFAST<<") requested";
			cout << "FROM CUBE DIM ("<<N(storageSLOW)<<", "<<N(storageMEDIUM)<<", "<<N(storageFAST)<<")\n";
			s = m = f = 0;
		}	
	}
	if (memdims==3) {
		try {
			size_t location(size_t(s)*size_t(N(storageFAST))*size_t(N(storageMEDIUM)));
			location += size_t(m)*size_t(N(storageFAST)) + size_t(f);
/*			cout << "s="<<s<<";m="<<m<<";f="<<f<<".  ";
			cout << "New location=" << location;
			cout << " -- old location would have been ";
			cout << s*N(storageFAST)*N(storageMEDIUM)+m*N(storageFAST)+f << "\n";*/
			return Data.at(location);
		}
		catch (const std::exception& error)
		{
			cerr << "Exception!!  ERROR accessing " << storageFAST << "=" << f << ", ";
			cerr << storageMEDIUM << "=" << m << ", " << storageSLOW << "=" << s;
			cerr << " from cube with dimensions " << N(storageFAST) << ", ";
			cerr << N(storageMEDIUM) << ", " << N(storageSLOW) << " with memdims=";
			cerr << memdims << "\n";
		}
	} else if (memdims==2) {
		if (s != current[storageSLOW]) {
			// page in new currentx
			long long pagesize((long long)(N(storageMEDIUM)) * (long long)(N(storageFAST)));
			current[storageSLOW] = s;
//			cout << "Paging in " << storageSLOW << "=" << current[storageSLOW] << "\n"; cout.flush();
			ifstream fi;
			fi.open(cubefilename.c_str());
			fi.seekg(begindata+(long long)(s)*pagesize*sizeof(float));
			for (long long i=0;i<pagesize;i++)
				Jcrap::getfloat(&fi, (void*)(&(Data[i])));
			fi.close();
		}
		return Data[m*N(storageFAST)+f];
	} else if (memdims==1) {
		if (s!=current[storageSLOW] || m!=current[storageMEDIUM]) {
			// page in new currentx
			long long pagesize((long long)N(storageFAST));
			current[storageSLOW] = s;
			current[storageMEDIUM] = m;
//			cout << "Paging in " << storageSLOW << "=" << current[storageSLOW] << ", ";
//			cout << storageMEDIUM << "=" << current[storageMEDIUM] << "\n"; cout.flush();
			ifstream fi;
			fi.open(cubefilename.c_str());
			fi.seekg(begindata+((long long)(s)*N(storageMEDIUM)+(long long)(m))
					*pagesize*sizeof(float));
			for (long long i=0;i<pagesize;i++)
				Jcrap::getfloat(&fi, (void*)(&(Data[i])));
			fi.close();
		}
		return Data[f];
	} else if (memdims==0) {
		if (s!=current[storageSLOW] || m!=current[storageMEDIUM] || f!=current[storageFAST]) {
			current[storageSLOW] = s;
			current[storageMEDIUM] = m;
			current[storageFAST] = f;
							
			ifstream fi;
			fi.open(cubefilename.c_str());
			fi.seekg(begindata+((long long)(s)*N(storageMEDIUM)*N(storageFAST)+
					              (long long)(m)*N(storageFAST)+(long long)(f))*sizeof(float));
			Jcrap::getfloat(&fi, (void*)(&(Data[0])));
			fi.close();
		}
		return Data[0];
	}
}	
		
void cube::pageout() const
// created 2006 Feb 14, pages changes in memory to disk for a bigcube.
// this used to be embedded within operator() above, but I found that I
// also needed to use it in operator=.
{
	long long pagesize, offset;
	if (memdims==2){
		pagesize=((long long)(N(storageMEDIUM)) * (long long)(N(storageFAST)));
		offset=current[storageSLOW]*pagesize;
	}
	if (memdims==1){
		pagesize=(N(storageFAST));
		offset=(current[storageSLOW]*N(storageMEDIUM)+current[storageMEDIUM])*pagesize;		
//		cout << "Paging out " << storageSLOW << "=" << current[storageSLOW] << ", ";
//		cout << storageMEDIUM << "=" << current[storageMEDIUM] << "\n"; cout.flush();
	}
	if (memdims==0){
		pagesize=1;	
//		cout << "Paging out " << storageSLOW << "=" << current[storageSLOW] << ", ";
//		cout << storageMEDIUM << "=" << current[storageMEDIUM] << ", ";
//		cout << storageFAST << "=" << current[storageFAST] << "\n"; cout.flush();
		offset=current[storageSLOW]*N(storageMEDIUM)*N(storageFAST)+
				 current[storageMEDIUM]*N(storageFAST) + current[storageFAST];
	}
	fstream fo;
	fo.open(cubefilename.c_str(), ios::in|ios::out|ios::binary|ios::ate);
	fo.seekp(begindata+offset*sizeof(float), ios::beg);
	
	for (long long i=0;i<pagesize;i++) {
//		cout << "Writing to disk " << storageSLOW << "=" << current[storageSLOW] << ", ";
//		cout << storageMEDIUM << "=" << current[storageMEDIUM] << ", ";
//		cout << storageFAST << "=" << i << " = ";
//		cout << Data.at(i) << " or " << storageorderaccess(current[storageSLOW], current[storageMEDIUM], i) << "\n";
		
		putfloat(&fo,Data.at(i));
	}
	fo.close();
}


float& cube::operator()(axis a1, long n1, axis a2, long n2, axis a3, long n3)
/* JB 10.31.2001
	for axis independance 
	This looks pretty slow as implimented now, work on someday maybe if I'm really bored*/
{
	long x(0), y(0), z(0);
	
	if (a1 == X) x=n1;
	if (a1 == Y) y=n1;
	if (a1 == Z) z=n1;
	
	if (a2 == X) x=n2;
	if (a2 == Y) y=n2;
	if (a2 == Z) z=n2;
	
	if (a3 == X) x=n3;
	if (a3 == Y) y=n3;
	if (a3 == Z) z=n3;
	
	return (*this)(x,y,z);
}

float& cube::operator()(axis a1, int n1, axis a2, int n2, axis a3, int n3)
{
	return operator()(a1, long(n1), a2, long(n2), a3, long(n3));
}

float& cube::operator()(axis a1, long n1, axis a2, long n2, axis a3, long n3) const
// added const 2006 March 26
{
	long x(0), y(0), z(0);
	
	if (a1 == X) x=n1;
	if (a1 == Y) y=n1;
	if (a1 == Z) z=n1;
	
	if (a2 == X) x=n2;
	if (a2 == Y) y=n2;
	if (a2 == Z) z=n2;
	
	if (a3 == X) x=n3;
	if (a3 == Y) y=n3;
	if (a3 == Z) z=n3;
	
	return (*this)(x,y,z);
}

float& cube::operator()(axis a1, int n1, axis a2, int n2, axis a3, int n3) const
{
	return operator()(a1, long(n1), a2, long(n2), a3, long(n3));
}

float P(float x) 
{
	float answer(x);
	if (x<=0) answer = 0.;
	return answer;
}

float R(float x)
// as defined at http://astronomy.swin.edu.au/~pbourke/colour/bicubic/
{
	float answer(P(x+2)*P(x+2)*P(x+2) -
			4.*P(x+1)*P(x+1)*P(x+1) +
			6.*P(x)*P(x)*P(x) -
			4.*P(x-1)*P(x-1)*P(x-1));
	return answer/6.;
}


float& cube::operator()(double x, double y, int z)
/* Created 3/6/2k
	Tried to allow this to be more easily generalized by introduction of interp2d. */
// changed this to return references if nearest neighbor 2016 August 18 JWB
{
	static float null(0.);
	int ix, iy;
	ix=Jlocate(x, X);
	if (ix!=nx-1  &&  x-Xaxis(ix) > (Xaxis(ix+1)-Xaxis(ix))/2) ix=ix+1;
	iy=Jlocate(y, Y);
	if (iy!=ny-1  &&  y-Axis(Y,iy) > (Axis(Y,iy+1)-Axis(Y,iy))/2) iy=iy+1;
	if (ix>=nx-1) {
		if (ix <= nx+1) ix = nx-1;
		else {
			cout << "ERROR -- operator()(double, double, int, mFAST) edge effects not yet coded\n";
			cout << "cube(" << x << ", " << y << ", " << z << ") from (";
			cout << N(X) << ", " << N(Y) << ", " << N(Z) << ")  ix=" << ix << "\n";
			return null=0.;
		}
	}
		if (iy>=ny-1){
		if (iy <= ny+1) iy = ny-1;
		else {
			cout << "ERROR -- operator()(double, double, int, mFAST) edge effects not yet coded\n";
			cout << "cube(" << x << ", " << y << ", " << z << ") from (";
			cout << N(X) << ", " << N(Y) << ", " << N(Z) << ")  ix=" << ix << ")\n";
			return null=0.;
		}
	}		
	return (*this)(ix,iy,z);
}

float cube::operator()(double x, double y, int z, methodtype m) const
/* Added a const version of this. */
{
	if (m==mFAST) {
		static float null(0.);
		int ix, iy;
		ix=Jlocate(x, X);
		if (ix!=nx-1  &&  x-Axis(X,ix) > (Axis(X,ix+1)-Axis(X,ix))/2) ix=ix+1;
		iy=Jlocate(y, Y);
		if (iy!=ny-1  &&  y-Axis(Y,iy) > (Axis(Y,iy+1)-Axis(Y,iy))/2) iy=iy+1;
		if (ix>=nx-1) {
			if (ix <= nx+1) ix = nx-1;
			else {
				cout << "ERROR -- operator()(double, double, int, mFAST) edge effects not yet coded\n";
				cout << "cube(" << x << ", " << y << ", " << z << ") from (";
				cout << N(X) << ", " << N(Y) << ", " << N(Z) << ")  ix=" << ix << "\n";
				return null=0.;
			}
		}
		if (iy>=ny-1){
			if (iy <= ny+1) iy = ny-1;
			else {
				cout << "ERROR -- operator()(double, double, int, mFAST) edge effects not yet coded\n";
				cout << "cube(" << x << ", " << y << ", " << z << ") from (";
				cout << N(X) << ", " << N(Y) << ", " << N(Z) << ")  ix=" << ix << ")\n";
				return null=0.;
			}
		}		
		return (*this)(ix,iy,z);
	}
	if (m==mSLOW || m==mBICUBIC) {
		cube c;
		const cube* u;
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(double, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		if (z==0 && nz==1) u=this;
		else { c=plane(Z, z); u=&c; }
		return (*u).interp2d(x, y, 2);
	}
	if (m==mLINEAR) {  // there's a memory leak somewhere in this section.
		float ans(0.);
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(double, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		int lx(Jlocate(x, X)+1), ly(Jlocate(y, Y)+1);  // low or left-end coord
		if (lx >= N(X) || lx < 1) return ans;
		if (ly >= N(Y) || ly < 1) return ans;
		float xfrac=(x-Axis(X,lx-1)) / (Axis(X,lx)-Axis(X,lx-1));
		float yfrac=(y-Axis(Y,ly-1)) / (Axis(Y,ly)-Axis(Y,ly-1));
		float xtop((1.-xfrac)*(*this)(lx-1,ly-1,z) + xfrac*(*this)(lx,ly-1,z));
		float xbot((1.-xfrac)*(*this)(lx-1,ly  ,z) + xfrac*(*this)(lx,ly  ,z));
		ans= float((1.-yfrac) * xtop               + yfrac * xbot);
		if (!(ans > 0  ||  ans <= 0)) {
			cout << "operator(double, double, int) returns non-number\n";
			cout << "x = " << x << ", y = " << y << ", z = " << z << "\n";
			cout << "lx = " << lx << ", ly = " << ly << ", or " << Axis(X,lx) << ", " << Axis(Y,ly) << "\n";
			cout << "Xvalues:  " << (*this)(lx-1,ly-1,z) << ", ";
			cout << xtop;
			cout << ", " << (*this)(lx,ly-1,z) << "\n";
			cout << "          " << (*this)(lx-1,ly,z) << ", ";
			cout << xbot;
			cout << ", " << (*this)(lx,ly-1,z) << "\n";
			cout << "xfrac = " << xfrac << ", yfrac = " << yfrac << "\n";
			cout << "ans = " << ans << ", setting equal to zero\n";
			ans = 0.;
		}
		
		return ans;
	}
	
	// remnants of non-const version inserted here for posterity 2016 August 18 JWB
	if (m==mSLOW || m==mBICUBIC) {  // done right 2006 May 7
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(double, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		
		float answer(0.);
		int lx(Jlocate(x, X)-1), ly(Jlocate(y, Y)-1);  // low or left-end coord
		if (lx >= N(X)-1 || lx < 0) return answer;
		if (ly >= N(Y)-1 || ly < 0) return answer;
//		lx--; ly--;
//		if (lx < 0) lx = 0;
//		if (lx > N(X)-4) lx = N(X)-4;
//		if (ly < 0) ly = 0;
//		if (ly > N(X)-4) ly = N(X)-4;
		float dx(x-float(lx)), dy(y-float(ly));
		
//		cout << "in bicubic -- requested = " << x << "," << y << "\n";
//		cout << "\tlx=" << lx << ",ly=" << ly << "\n";
//		cout << "\tdx=" << dx << ",dy=" << dy << "\n";
		
		for (int m(-1);m<3;m++) {
			for (int n(-1);n<3;n++) {
				int usex(lx+m), usey(ly+n);
				if (usex<0) usex=0;
				if (usex>N(X)-1) usex=N(X)-1;
				if (usey<0) usey=0;
				if (usey>N(Y)-1) usey=N(Y)-1;
				answer += (*this)(usex,usey,z)*R(float(m)-dx)*R(dy-float(n));
		// R defined just before this function.
			}
		}
		
		return answer;
	}
	if (m==mLANCZOS) {
		static const float rparam(1.);
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(double, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		float answer(0.);
		int lx(Jlocate(x, X)-1), ly(Jlocate(y, Y)-1);  // low or left-end coord
		if (lx >= N(X)-1 || lx < 0) return answer;
		if (ly >= N(Y)-1 || ly < 0) return answer;
		
		int xmin(lx-int(rparam)), xmax(lx+1+int(rparam));
		int ymin(ly-int(rparam)), ymax(ly+1+int(rparam));
		if (xmin<0) xmin=0;
		if (xmax>N(X)-1) xmax=N(X)-1;
		if (ymin<0) ymin=0;
		if (ymin>N(Y)-1) ymax=N(Y)-1;

		float normalizer(0.);		
		for (int nx(xmin);nx<=xmax;nx++) {
			for (int ny(ymin);ny<=ymax;ny++) {
				static const float pi(3.1415926535897932384);
				static const float pi3(pi/3.);
				float dist(sqrtf((float(nx)-x)*(float(nx)-x) +
						(float(ny)-y)*(float(ny)-y)));
				if (dist > rparam) continue;
				
				float coeff(sin(pi*dist/rparam)/(pi*dist/rparam) *
						sin(pi3*dist/rparam)/(pi3*dist/rparam));
				normalizer+=coeff;
				answer += coeff * (*this)(nx,ny,z);
			
			}
		}
		answer /= normalizer;
		return answer;
	}
	
	if (m==mLINEAR) {  // updated, works well 2006 May 7
		float answer(0.);
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(double, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		int lx(Jlocate(x, X)+1), ly(Jlocate(y, Y)+1);  // low or left-end coord
		//cout<<lx<<"\t"<<ly<<"\n";
		if (lx >= N(X) || lx < 1) return answer;
		if (ly >= N(Y) || ly < 1) return answer;
		float xfrac=(x-Axis(X,lx-1)) / (Axis(X,lx)-Axis(X,lx-1));
		float yfrac=(y-Axis(Y,ly-1)) / (Axis(Y,ly)-Axis(Y,ly-1));
		float xtop((1.-xfrac)*(*this)(lx-1,ly-1,z) + xfrac*(*this)(lx,ly-1,z));
		float xbot((1.-xfrac)*(*this)(lx-1,ly  ,z) + xfrac*(*this)(lx,ly  ,z));
		answer=float((1.-yfrac) * xtop               + yfrac * xbot);
		if (!(answer > 0  ||  answer <= 0)) {
			cout << "operator(double, double, int) returns non-number\n";
			cout << "x = " << x << ", y = " << y << ", z = " << z << "\n";
			cout << "lx = " << lx << ", ly = " << ly << ", or " << Axis(X,lx) << ", " << Axis(Y,ly) << "\n";
			cout << "Xvalues:  " << (*this)(lx-1,ly-1,z) << ", ";
			cout << xtop;
			cout << ", " << (*this)(lx,ly-1,z) << "\n";
			cout << "          " << (*this)(lx-1,ly,z) << ", ";
			cout << xbot;
			cout << ", " << (*this)(lx,ly-1,z) << "\n";
			cout << "xfrac = " << xfrac << ", yfrac = " << yfrac << "\n";
			cout << "*ans = " << answer << ", setting equal to zero\n";
			answer = 0.;
		}
		
		return answer;
	}
}

float& cube::operator()(double x, int y, int z)
// split from double, int, int, methodtype 2016 August 18 JWB
{
	int ix;
	ix=Jlocate(x, X);
	if (ix<nx-1  &&  x-Xaxis(ix) > (Xaxis(ix+1)-Xaxis(ix))/2) ix++;
	if (ix>nx-1) 
		cout << "ERROR -- operator()(double, double, int, mFAST) edge effects not yet coded\n";
	return (*this)(ix, y, z);
}

float cube::operator()(double x, int y, int z, methodtype m) const
/* Created 3/9/2k*/
// altered for const and parallelization 2016 August 18 JWB
{
	float answer(0.);
	if (m==mFAST) {
		int ix;
		ix=Jlocate(x, X);
		if (ix<nx-1  &&  x-Axis(X,ix) > (Axis(X,ix+1)-Axis(X,ix))/2) ix++;
		if (ix>nx-1) 
			cout << "ERROR -- operator()(double, double, int, mFAST) edge effects not yet coded\n";
		answer = (*this)(ix, y, z);
	}
	if (m==mSLOW) {
		cube c;
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(double, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		c=skewer(X, y, z);
		answer = c.interp(x);
	}
	return answer;
}

float& cube::operator()(int x, double y, int z)
// split from int, double, int, methodtype 2016 August 18 JWB
{		
	int iy;
	iy=Jlocate(y, Y);
	if (iy<ny-1  &&  y-Axis(Y,iy) > (Axis(Y,iy+1)-Axis(Y,iy))/2) iy++;
	if (iy>ny-1) 
		cout << "ERROR -- operator()(int, double, int, mFAST) edge effects not yet coded\n";
	return (*this)(x, iy, z);
}

float cube::operator()(int x, double y, int z, methodtype m) const
/* Created 3/9/2k*/
// rejiggered and split off the reference version (above) 2016 August 18 JWB
{
	float answer(0.);
	if (m==mFAST) {
		int iy;
		iy=Jlocate(y, Y);
		if (iy<ny-1  &&  y-Axis(Y,iy) > (Axis(Y,iy+1)-Axis(Y,iy))/2) iy++;
		if (iy>ny-1) 
			cout << "ERROR -- operator()(int, double, int, mFAST) edge effects not yet coded\n";
		answer = (*this)(x, iy, z);
	}
	if (m==mSLOW) {
		cube c;
		if (z > nz-1  ||  z<0  ||  x>nx-1 || x<0){ 
			cout << "ERROR -- operator(int, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		c=skewer(Y, z, x);
		float *ans;
		ans = new float;
		*ans =c.interp(y);
		answer = *ans;
	}
	if (m==mLINEAR) {  // there's a memory leak somewhere in this section.
		static float ans(0.);
		if (z > nz-1  ||  z<0){ 
			cout << "ERROR -- operator(int, double, int) z("<<z;
			cout << " out of range 0.."<<nz-1<<"\n";
		}
		if (x > nx-1  ||  x<0){ 
			cout << "ERROR -- operator(int, double, int) x("<<x;
			cout << " out of range 0.."<<nx-1<<"\n";
		}
		int ly(Jlocate(y, Y));  // low or left-end coord
		if (ly == N(Y)-1)  {ly = N(Y)-2;}
		if (ly > N(Y) || ly < 0) return ans;
		float yfrac=(y-Axis(Y,ly)) / (Axis(Y,ly+1)-Axis(Y,ly));
	
		ans= float((1.-yfrac) * (*this)(x,ly,z)     + yfrac * (*this)(x,ly+1,z));
		
		if (!(ans > 0  ||  ans <= 0)) {
			cout << "operator(int, double, int) returns non-number\n";
			cout << "x = " << x << ", y = " << y << ", z = " << z << "\n";
			cout << ", ly = " << ly << ", or " << ", " << Axis(Y,ly) << "\n";
			cout << "Xvalues:  " << (*this)(x,ly,z) << ", ";
			cout << (*this)(x,ly,z);
			cout << ", " << (*this)(x,ly,z) << "\n";
			cout << "          " << (*this)(x,ly+1,z) << ", ";
			cout << (*this)(x,ly+1,z);
			cout << ", " << (*this)(x,ly,z) << "\n";
			cout << ", yfrac = " << yfrac << "\n";
			cout << "ans = " << ans << ", setting equal to zero\n";
			ans = 0.;
		}
		
		answer = ans;
	}
	return answer;
}

float& cube::operator()(int x, int y, double z)
// broken off of (int, int, double, methodtype) for parallelization to 
// lose the static, 2016 August 18 JWB
{		
	int iz;
	iz=Jlocate(z, Z);
	if (iz<nz-1  &&  z-Axis(Z,iz) > (Axis(Z,iz+1)-Axis(Z,iz))/2) iz++;
	if (iz>nz-1) 
		cout << "ERROR -- & operator()(int, int, double) edge effects not yet coded\n";
	return (*this)(x, y, iz);
}

float cube::operator()(int x, int y, double z, methodtype m) const
/* Created 11/3/2k -- modified 2016 August 18 JWB*/
{
	if (m==mFAST) {
		int iz;
		iz=Jlocate(z, Z);
		if (iz<nz-1  &&  z-Axis(Z,iz) > (Axis(Z,iz+1)-Axis(Z,iz))/2) iz++;
		if (iz>nz-1) 
			cout << "ERROR -- operator()(int, int, double, mFAST) edge effects not yet coded\n";
		return (*this)(x, y, iz);
	}
	if (m==mSLOW || m==mLINEAR) {
		float answer((*this)(x,y,0));
		if (N(Z)>1) {
			int iz;
			iz=Jlocate(z, Z);
			if (z > Axis(Z,N(Z)-1)) answer = (*this)(x,y,N(Z)-1);
			if (z < Axis(Z,0))      answer = (*this)(x,y,0);
			float diff = (*this)(x, y, iz+1) - (*this)(x, y, iz);
			float dist = Axis(Z,iz+1) - Axis(Z,iz);
			
			answer = ((z-Axis(Z,  iz))/dist) * (*this)(x, y, iz+1)  +
					   ((Axis(Z,iz+1)-z)/dist) * (*this)(x, y, iz);
		}
		return answer;
		
	}
}

cube cube::operator()(double x1, double x2, double y1, double y2, int z1, int z2) const
{
	int xmin, xmax, ymin, ymax;
	xmin=Jlocate(x1, X);
	xmax=Jlocate(x2, X)+1;
	if (xmax>=nx) xmax=nx-1;
	ymin=Jlocate(y1, Y);
	ymax=Jlocate(y2, Y)+1;
	if (ymax>=ny) ymax=ny-1;
	
	cube sc;
	sc = (*this)(xmin, xmax, ymin, ymax, z1, z2);
	return sc;
}

cube cube::operator()(double x1, double x2, int y1, int y2, int z1, int z2) const
{
	int xmin, xmax;
	xmin=Jlocate(x1, X);
	xmax=Jlocate(x2, X)+1;
	if (xmax>=nx) xmax=nx-1;
	
	cube sc;
	sc = (*this)(xmin, xmax, y1, y2, z1, z2);
	return sc;
}

cube cube::operator()(int x1, int x2, double y1, double y2, int z1, int z2) const
// Created 2013 September 30 JWB
{
	int ymin, ymax;
	ymin=Jlocate(y1, Y);
	ymax=Jlocate(y2, Y)+1;
	if (ymax>=N(Y)) ymax=N(Y)-1;
	
	return (*this)(x1, x2, ymin, ymax, z1, z2);
}

cube cube::operator()(int x1, int x2, int y1, int y2, double z1, double z2) const
{
	int zmin=0, zmax=nz-1;
	if (Axis(Z,0) < Axis(Z,1)) {
		while(Axis(Z,zmin) < z1) zmin++;
		while(Axis(Z,zmax) > z2) zmax--;
	}
	else {
		cout << "Axis(Z,zmin) = " << Axis(Z,zmin) << " and z2 is " << z2 << "\n";
		while(Axis(Z,zmin) > z2) zmin++;
		while(Axis(Z,zmax) < z1) zmax--;
	}
	if (z1==z2) zmin=zmax=Jlocate(z1, Z);
	if (osuppress<=0) cout << "zrange "<<z1<<" to "<<z2<<" interpreted as "<<zmin<<" to "<<zmax<<"\n";
	
	return (*this)(x1, x2, y1, y2, zmin, zmax);
}
	

cube cube::operator()(int x1, int x2, int y1, int y2, int z1, int z2) const
{
	if (x1==-1) x1=0;
	if (x2==-1) x2=nx-1;
	if (y1==-1) y1=0;
	if (y2==-1) y2=ny-1;
	if (z1==-1) z1=0;
	if (z2==-1) z2=nz-1;
	 
	
	if (x1<0 || x2>=nx || y1<0 || y2>=ny || z1<0 || z2>=nz || x2<0 || y2<0 || z2<0){
		cout << "Error in subcube() operator: subcube larger than source\n";
		cout << "Attempt to obtain (" << x1 << ".." << x2 << "," << y1 << "..";
		cout << y2 << "," << z1 << ".." << z2 << ") from total size (";
		cout << nx << "," << ny << "," << nz << ").\n";
	}
	
	cube subcube(x2-x1+1, y2-y1+1, z2-z1+1, 0.0, cubetype);
	subcube.copymetadatafrom(*this);
	
//	cout << "Subcube created ("<<x2-x1+1<<","<<y2-y1+1<<","<<z2-z1+1<<"\n";
	
	for (int x=x1;x<=x2;x++) subcube.Axis(X,x-x1)=Axis(X,x);
	for (int y=y1;y<=y2;y++) subcube.Axis(Y,y-y1)=Axis(Y,y);
	for (int z=z1;z<=z2;z++) subcube.Axis(Z,z-z1)=Axis(Z,z);
	
	subcube.hdr["Nx"]=int2str(subcube.Nx());
	subcube.hdr["Ny"]=int2str(subcube.Ny());
	subcube.hdr["Nz"]=int2str(subcube.Nz());
	
	
		
//	cout << "axes copied\n";
	
	for (int x=x1;x<=x2;x++)
		for (int y=y1;y<=y2;y++)
			for (int z=z1;z<=z2;z++)
				subcube(x-x1,y-y1,z-z1)=(*this)(x,y,z);

	return subcube;
}

cube cube::operator()(axis a1, int a1start, int a1end,
							 axis a2, int a2start, int a2end,
							 axis a3, int a3start, int a3end) const
// JB 2005 Jaunary 6 -- generalized subcube operator()
{
// figure out non-subcubed axes
	if (a3 == UNK  &&  a2 != UNK)
		if (a2 == Uax(a1)) a3 = Dax(a1);
		else a3 = Uax(a1);
	if (a3 == UNK  &&  a2 == UNK) {
		a2 = Uax(a1);
		a3 = Dax(a1);
	}
	
// convert default values to usable form
	if (a1start == -1) a1start = 0;
	if (a1end   == -1) a1end   = N(a1)-1;
	if (a2start == -1) a2start = 0;
	if (a2end   == -1) a2end   = N(a2)-1;
	if (a3start == -1) a3start = 0;
	if (a3end   == -1) a3end   = N(a3)-1;
	 
	
// check for pathological parameters
	if (a1start<0 || a1end>=N(a1)    || 
		 a2start<0 || a2end>=N(a2)    || 
		 a3start<0 || a3end>=N(a3)    || 
		 a1end<0 || a2end<0 || a3end<0){
		cout << "Error in subcube() operator: subcube larger than source\n";
		cout << "Attempt to obtain (" << a1 << "=" << a1start << ".." << a1end;
		cout << "," << a2 << "=" << a2start << ".." << a2end;
		cout << "," << a3 << "=" << a3start << ".." << a3end << ") from total size (";
		cout << nx << "," << ny << "," << nz << ").\n";
	}
	
// create the output cube
	cube subcube(a1, a1end-a1start+1, a2, a2end-a2start+1, a3, a3end-a3start+1,
			0.0, cubetype);
	
// copy the axis values across	
	for (long x=a1start;x<=a1end;x++) subcube.Axis(a1,x-a1start)=Axis(a1,x);
	for (long y=a2start;y<=a2end;y++) subcube.Axis(a2,y-a2start)=Axis(a2,y);
	for (long z=a3start;z<=a3end;z++) subcube.Axis(a3,z-a3start)=Axis(a3,z);
	
// copy over other information.
	subcube.hdr=hdr;
	subcube.byteorder=byteorder;
	subcube.hdr["Nx"]=int2str(subcube.N(X));
	subcube.hdr["Ny"]=int2str(subcube.N(Y));
	subcube.hdr["Nz"]=int2str(subcube.N(Z));
	
// now, finally copy over all the data
	for (long x=a1start;x<=a1end;x++)
		for (long y=a2start;y<=a2end;y++)
			for (long z=a3start;z<=a3end;z++)
				subcube(a1,x-a1start, a2,y-a2start, a3,z-a3start) =
						(*this)(a1,x, a2,y, a3,z);

	return subcube;
}

cube cube::operator+(cube r)
/*  Fucked with, hopefully improved, 2/21/2k.  Just redid things so that this,
	one of the first operators created, was more like the newer ones.  Also added
	interpolative addition section for Z 
*/
{
	if (nx == r.nx  &&  ny == r.ny  &&  nz == r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		cube sum(*this);
		if (!osuppress) cout <<"Summing . . .     00%";
		for (int z=0;z<nz;z++){
			if (!osuppress) printpercent(z, nz-1);
			for (int y=0;y<ny;y++)
				for (int x=0;x<nx;x++)
					sum(x,y,z) += r(x,y,z);
		}
		if (!osuppress) printpercent(nz, nz);
		if (!osuppress) cout << '\n';
		return sum;
	}		
	else if (nx == r.nx  &&  ny == r.ny  &&  nz == r.nz){
		cube sum(*this);
		if (!osuppress) cout <<"Summing . . .     00%";
		for (int z=0;z<nz;z++){
			if (!osuppress) printpercent(z, nz-1);
			for (int y=0;y<ny;y++)
				for (int x=0;x<nx;x++)
					sum(x,y,z) += r(x,y,z);
		}
		if (!osuppress) printpercent(nz, nz);
		if (!osuppress) cout << '\n';
		return sum;
	}	
	else if (r.nx == 1  &&  ny == r.ny  &&  nz == r.nz){
		cube sum(*this);
		if (!osuppress) cout <<"Summing . . .     00%";
		for (int z=0;z<nz;z++){
			if (!osuppress) printpercent(z, nz-1);
			for (int y=0;y<ny;y++)
				for (int x=0;x<nx;x++)
					sum(x,y,z) += r(0,y,z);
		}
		if (!osuppress) printpercent(nz, nz);
		if (!osuppress) cout << '\n';
		return sum;
	}	
	else if (Ztype()==r.Ztype()  &&  nx==r.Nx()&&ny==r.Ny() && Ztype()>0 && Ztype()<=4) {
		if (!osuppress) cout << "Performing interpolative addition.\n";
		cube c1, c2;
		interpolativeZ(r, c1, c2);
		
		cube sum;
		sum=c1+c2;
		return sum;
	}
	else {
		cerr << "INCOMPATIBLE CUBES TO BE SUMMED";
		cerr << nx << " " << ny << " " << nz << " + ";
                cerr << r.nx << " " << r.ny << " " << r.nz << " + ";
		exit(1);
	}
}

cube cube::operator-(cube r) const
{
	cube difference(*this);
	
	if (nx==r.nx && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					difference(x,y,z)=(*this)(x,y,z)-r(x,0,z);
	}
	else if (nx==r.nx && ny==r.ny && nz==r.nz/* && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)*/){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					difference(x,y,z)=(*this)(x,y,z)-r(x,y,z);
	}	
	else if (r.nx==1 && r.ny==1 && nz==r.nz){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					difference(x,y,z)=(*this)(x,y,z)-r(0,0,z);
	}	
	else if (nx==r.nx && ny==r.ny && r.nz==1){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					difference(x,y,z)=(*this)(x,y,z)-r(x,y,0);
	}
	else if (r.nx==1 && ny==r.ny && nz==r.nz){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					difference(x,y,z)=(*this)(x,y,z)-r(0,y,z);
	}
	else if (Ztype()==r.Ztype()  &&  nx==r.Nx()&&ny==r.Ny() && Ztype()>0 && Ztype()<=4) {
		cout << "Performing interpolative subtraction -- BROKEN BY CONST 2006/June30.\n";
		cout << "It was ancient and shitty, anyway.\n";
		cube c1, c2;
//		interpolativeZ(r, c1, c2);
		
		cube diff;
		diff=c1-c2;
		return diff;
	}
	else{
		cout << "PART OF SUBTRACTION FUNCTION REQUESTED IS NOT YET IMPLIMENTED:";
		cout << " ("<<nx<<","<<ny<<","<<nz<<")-("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
	}

	return difference;
}

cube cube::operator*(cube r) const
/* Created 2/22/2k -- Unbelivably this didn't exist before */
{
	cube product(*this);
	
// (x,#,z)*(x,1,z)
	if (nx==r.nx && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					product(x,y,z)=(*this)(x,y,z)*r(x,0,z);
	}
// (x,y,z)*(x,y,z)
	else if (nx==r.nx && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					product(x,y,z)=(*this)(x,y,z)*r(x,y,z);
	}	
// (#,y,z)*(1,y,z)
	else if (ny==r.ny && nz==r.nz && r.nx==1 && Axis(Y,0)==r.Axis(Y,0) && Axis(Y,ny-1)==r.Axis(Y,ny-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					product(x,y,z)=(*this)(x,y,z)*r(0,y,z);
	}
// (x,y,#)*(x,y,1)
	else if (ny==r.ny && nx==r.nx && r.nz==1){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					product(x,y,z)=(*this)(x,y,z)*r(x,y,0);
	}
// (#,#,z)*(1,1,z)
	else if (r.nx==1 && r.ny==1 && nz==r.nz){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					product(x,y,z)=(*this)(x,y,z)*r(0,0,z);
	}
	else if (Ztype()==r.Ztype()  &&  nx==r.Nx()&&ny==r.Ny() && Ztype()>0 &&	Ztype()<=7) {
		cout << "Performing interpolative multiplication -- BROKEN BY CONST 2006 June.\n";
		cout << "Boring conversation anyway.\n";
		cube c1, c2;
//		interpolativeZ(r, c1, c2);
		
		cube product;
		product=c1*c2;
		return product;
	}		
	else{
		cout << "PART OF * FUNCTION REQUESTED IS NOT YET IMPLIMENTED:";
		cout << " ("<<nx<<","<<ny<<","<<nz<<")*("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
	}

	return product;
}

cube cube::operator/(cube r)
{
	cube quotient(*this);
		
	if (r.nx==1 && r.ny==1 && r.nz==nz) {
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					quotient(x,y,z)=(*this)(x,y,z)/r(0,0,z); }
	else if (nx==r.nx && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					quotient(x,y,z)=(*this)(x,y,z)/r(x,0,z);
	} else if (nx==r.nx && ny==r.ny && r.nz==1){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					quotient(x,y,z)=(*this)(x,y,z)/r(x,y,0);
	} else if (nx==r.nx && ny==r.ny && nz==r.nz) {
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					quotient(x,y,z)=(*this)(x,y,z)/r(x,y,z);
	} else if (r.nx==1 && ny==r.ny && nz==r.nz) {
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					quotient(x,y,z)=(*this)(x,y,z)/r(0,y,z);
	} else if (nx==r.nx && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int y=0;y<ny;y++)
			for (int x=0;x<nx;x++)
				for (int z=0;z<nz;z++)
					quotient(x,y,z)=(*this)(x,y,z)/r(x,y,z);
	} else if (Ztype()==r.Ztype()  &&  nx==r.Nx()&&ny==r.Ny() && Ztype()>0 && Ztype()<=4) {
		cout << "Performing interpolative division.\n";
		cube c1, c2;
		interpolativeZ(r, c1, c2);
		
		cube q;
		q=c1/c2;
		return q;
	} else{
		cout << "PART OF DIVISION FUNCTION REQUESTED IS NOT YET IMPLIMENTED:";
		cout << " ("<<nx<<","<<ny<<","<<nz<<")/("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
		if (Ztype()==r.Ztype()){
			cout << "(" << Axis(Z,0) << ".." << Axis(Z,nz-1) << ")/(";
			cout << r.Axis(Z,0) << ".." << (r.Axis(Z,r.Nz())-1) << ")\n";
		}
	}

	return quotient;
}

cube cube::operator+=(float r)
{
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				(*this)(x,y,z)+=r;
	return *this;
}

cube cube::operator-=(float r)
{
	for (cube::iterator i(begin());i!=end();i++) *i -= r;
	return *this;
}		

cube cube::operator*=(float r)
{
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				(*this)(x,y,z)*=r;
	return *this;
}		

cube cube::operator/=(float d)
{
	if (d==0.0) cerr << "DIVISION BY ZERO!\n";
	if (!osuppress) cout << "Dividing . . .  0%";
	cout.flush();
	for (int x=0;x<nx;x++){
		if (!osuppress) printpercent(x,(nx-1));
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				(*this)(x,y,z)/=d;
	}
	return *this;
}

cube cube::operator*(float r) const
{
	cube p(*this);
	for (cube::iterator i=p.begin();i!=p.end();i++)
		(*i)*=r;
	return p;
}

cube cube::operator/(float r)
{
	cube q(*this);
	for (cube::iterator i=q.begin();i!=q.end();i++)
		(*i)/=r;
	return q;
}

cube cube::operator+(float r)
{
	cube s(*this);
	for (cube::iterator i=s.begin();i!=s.end();i++)
		(*i)+=r;
	return s;
}

cube cube::operator-(float r)
{
	cube d(*this);
	for (cube::iterator i=d.begin();i!=d.end();i++)
		(*i)-=r;
	return d;
}

cube cube::operator+=(cube r)
{
	if (r.N(X)==N(X) && r.N(Y)==N(Y) && r.N(Z)==N(Z)) {
		for (cube::iterator a(begin()), b(r.begin());a!=end();a++,b++)
			*a += *b;
	}
	else if (r.nx==1 && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)+=r(0,y,z);
	} 
	else if (r.nx==1 && r.nz==1 && ny==r.ny && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)+=r(0,y,0);
	} 
	else if (r.nx==1 && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)+=r(0,0,z);
	} 
	else if (r.nx==nx && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)+=r(x,0,z);
	} 
	else{
		cout << "PART OF += NOT YET IMPLIMENTED--("<<nx<<","<<ny<<","<<nz<<")";
		cout << " += ("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
	}
	return *this;
}

cube cube::operator-=(cube& r)
{
	if (nx==r.nx && ny==r.ny && nz==r.nz) {
		for (cube::iterator i=begin(),j=r.begin();i!=end();i++,j++)
			*i -= *j;
	}
	else if (r.nx==1 && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)-=r(0,y,z);
	} 
	else if (r.nx==1 && r.nz==1 && ny==r.ny && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)-=r(0,y,0);
	} 
	else if (r.nx==1 && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)-=r(0,0,z);
	} 
	else if (nx==r.nx && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)-=r(x,0,z);
	}
	else{
		cout << "PART OF -= NOT YET IMPLIMENTED--("<<nx<<","<<ny<<","<<nz<<")";
		cout << " -= ("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
	}
	return *this;
}

cube cube::operator*=(cube &r)
{
	if (r.nx==1 && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)*=r(0,y,z);
	} else if (r.nx==1 && ny==r.ny && r.nz==1 && Axis(Y,0)==r.Axis(Y,0) && Axis(Y,ny-1)==r.Axis(Y,ny-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)*=r(0,y,0);
	} else if (nx==r.nx && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)*=r(x,0,z);
	} else if (nx==r.nx && ny==r.ny && r.nz==1 && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)*=r(x,y,0);
	} else if (nx==r.nx && r.ny==1 && r.nz==1 && Xaxis(0)==r.Xaxis(0) && Xaxis(nx-1)==r.Xaxis(nx-1)){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)*=r(x,0,0);
	} else if (r.nx==1 && r.ny==1 && nz==r.nz){
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)*=r(0,0,z);
	} else
		cout << "PART OF *= NOT YET IMPLIMENTED("	
			  <<nx<<","<<ny<<","<<nz<<")*=("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
	return *this;
}

cube cube::operator/=(cube &r)
{
	if (r.nx==1 && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		cout << "Dividing cube by ytplane\n";		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)/=r(0,y,z);
	} 	
	else if (r.nz==1 && ny==r.ny && nx==r.nx){
		cout << "Dividing cube by xyplane\n";		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)/=r(x,y,0);
	} 
	else if (r.ny==1 && nx==r.nx && nz==r.nz){
		cout << "Dividing cube by xzplane\n";		
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)/=r(x,0,z);
	} 
	else if (r.nx==1 && ny==r.ny && r.nz==1  && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1)){
		cout << "Dividing cube by Y spatium\n";
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)/=r(0,y,0);
	} 
	else if (nx==r.nx && ny==r.ny && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1) ){
		cout << "Dividing cube by cube\n";
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)/=r(x,y,z);
	} 	
	else if (r.nx==1 && r.ny==1 && nz==r.nz && Axis(Z,0)==r.Axis(Z,0) && Axis(Z,nz-1)==r.Axis(Z,nz-1) ){
		cout << "Dividing cube by spectrum\n";
		for (int x=0;x<nx;x++)
			for (int y=0;y<ny;y++)
				for (int z=0;z<nz;z++)
					(*this)(x,y,z)/=r(0,0,z);
	}		
	else
		cout << "PART OF /= NOT YET IMPLIMENTED("
				<<nx<<","<<ny<<","<<nz<<")/=("<<r.nx<<","<<r.ny<<","<<r.nz<<")\n";
	return *this;
}

cube cube::pow(double exp)
/* Created 1/25/2002 while at the 61" by JB */
{
	cube answer(*this);
	for (cube::iterator i(answer.begin());i!=answer.end();i++) 
		*i = ::pow(*i, exp);
	return answer;
}

cube cube::log(double base)
/* Created 2008 Feb 11 by JB */
{
	cube answer(*this);
	for (cube::iterator i(answer.begin());i!=answer.end();i++) 
		*i = ::log(*i)/::log(base);
	return answer;
}
void cube::interpolativeZ(cube& r, cube &c1, cube &c2)
/* Created 2/21/2k

	This routine changes c1 and c2 to be compatible versions of *this and r 
	respectively.
	The length of c1 and c2 is equal to that of the overlap region of *this and r, 
	and the resolution of c1 and c2 is equal to the greater of the two input cubes.
	
	Frigging vaporized the first version, I hope this one works as well. */
{
	if (Axis(Z,0) > Axis(Z,nz-1)) *this=dirinc();
	if (r.Axis(Z,0) > r.Axis(Z,r.Nz()-1)) r=r.dirinc();
	
// Determine overlap region
	double lo, hi;
	lo = FMAX(Axis(Z,0), r.Axis(Z,0));
	hi = FMIN(Axis(Z,nz-1), r.Axis(Z,r.Nz()-1));
	cout << "overlap = " << lo << " to " << hi << " \n";
	
// Chop edges (slower than it used to be in prevaporized version)
	c1 = (*this)(0,0,0,0,lo,hi);
	c2 = r(0,0,0,0,lo,hi);
	
// resample
	if (c2.Nz() > c1.Nz()) c1=c1.resample(c2);
	else c2=c2.resample(c1);
}

void cube::interpolativeX(cube& r, cube &c1, cube &c2)
/* Created 1/24/2k+1

	This routine changes c1 and c2 to be compatible versions of *this and r 
	respectively.
	The length of c1 and c2 is equal to that of the overlap region of *this and r, 
	and the resolution of c1 and c2 is equal to the greater of the two input cubes.
	
	Frigging vaporized the first version, I hope this one works as well. */
{
	if (Xaxis(0) > Xaxis(nx-1)) *this=dirinc();
	if (r.Xaxis(0) > r.Xaxis(r.Nx()-1)) r=r.dirinc();
	
// Determine overlap region
	double lo, hi;
	lo = FMAX(Xaxis(0), r.Xaxis(0));
	hi = FMIN(Xaxis(nx-1), r.Xaxis(r.Nx()-1));
	cout << "overlap = " << lo << " to " << hi << " \n";
	
// Chop edges (slower than it used to be in prevaporized version)
	c1 = (*this)(lo,hi,-1,-1,-1,-1);
	c2 = r(lo,hi,-1,-1,-1,-1);
	
// resample
	if (c2.Nx() > c1.Nx()) c1=c1.resample(c2);
	else c2=c2.resample(c1);
}
