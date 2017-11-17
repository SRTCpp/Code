#include <algorithm>
#include "Jcrap.h"

cube cube::plane(axis a, int v1, int v2) const
/* Created 3/6/2k */
{
	if (a==UNK) cout << "ERROR -- plane called with axis UNK\n";
	if (v1==-1){
		v1=0;
		v2=N(a)-1;
	}
	if (v2==-1) v2=v1;
	if(v1<0 || v1 > v2 || v2 < v1 || v2 > N(a)-1){ 
		cout << "ERROR -- plane("<<a<<") called with v1"<<v1<<" and/or v2("<<v2;
		cout << ") invalid, nz=" << nz << "\n";
	}
	
// establish size of new cube
	int px, py, pz, xit=1, yit=1, zit=1;
	if (a==X) { px=1; py=ny; pz=nz; xit=0; }
	if (a==Y) { px=nx; py=1; pz=nz; yit=0; }
	if (a==Z) { px=nx; py=ny; pz=1; zit=0; }
	cube p(px, py, pz, 0.0, cubetype, a, Uax(a), Dax(a));
	
	
	
// set up new axes
	p.Xtype(Xtype());
	p.Ytype(Ytype());
	p.Ztype(Ztype());
	for (int x=0;x<p.Nx();x++) p.xaxis[x]=xaxis[x];
	for (int y=0;y<p.Ny();y++) p.yaxis[y]=yaxis[y];
	for (int z=0;z<p.Nz();z++) p.zaxis[z]=zaxis[z];
	p.Axis(a,0) = 0.;
	for (int i=v1;i<=v2;i++) p.Axis(a, 0)+=Axis(a,i);
	p.Axis(a,0) /= v2-v1+1;
	
// determine coordinates to loop over and loop
	int xl=0, xh=nx-1, yl=0, yh=ny-1, zl=0, zh=nz-1;
	if (a==X) { xl=v1; xh=v2; }
	if (a==Y) { yl=v1; yh=v2; }
	if (a==Z) { zl=v1; zh=v2; }
	for (int x=xl;x<=xh;x++)
		for (int y=yl;y<=yh;y++)
			for (int z=zl;z<=zh;z++)
				p(xit*x, yit*y, zit*z)+=(*this)(x,y,z);
	if (!osuppress) cout << "plane complete -- " << v1 << " to " << v2 << " averaged in "<<a<<".\n";
	return p;
}
			
			
cube cube::plane(axis a1, axis a2, int v1, int v2)
/* Created 3/6/2k as a generalization of xyplane, ytplanes.  Useful?  you decide. */
{
	axis a;
	if (a1==X && a2==Y  ||  a2==X && a1==Y) a=Z;
	if (a1==X && a2==Z  ||  a2==X && a1==Z) a=Y;
	if (a1==Y && a2==Z  ||  a2==Y && a1==Z) a=X;
	cube p;
	p=plane(a, v1, v2);
	return p;
}

cube cube::plane(axis a, double f1, double f2) const
{
	if (f2==-1) f2=f1;
	return plane(a, int(Jlocate(f1,a)), int(Jlocate(f2,a)));
}

cube cube::xt2xy()
{
	cube xy(nx, nz, ny);
	
	for (int x=0;x<nx;x++) xy.xaxis[x]=xaxis[x];
	for (int y=0;y<nz;y++) xy.yaxis[y]=zaxis[y];
	for (int z=0;z<ny;z++) xy.zaxis[z]=yaxis[z];
	
	xy.hdr=hdr;
	xy.Ytype(Ztype());
	xy.Ztype(Ytype());
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				xy(x,z,y)=(*this)(x,y,z);
	return xy;
}

cube cube::yt2xy()
{
	cube xy(nz, ny, nx);
	xy.hdr=hdr;
	xy.cubetype=cubetype;
	
	for (int x=0;x<nz;x++) xy.xaxis[x]=zaxis[x];
	for (int y=0;y<ny;y++) xy.yaxis[y]=yaxis[y];
	for (int z=0;z<nx;z++) xy.zaxis[z]=xaxis[z];
	xy.Xtype(Ztype());
	xy.Ztype(Xtype());	
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				xy(z,y,x)=(*this)(x,y,z);
	return xy;
}

cube cube::yt2xt()
{
	cube xt(ny, nx, nz);
	
	for (int x=0;x<ny;x++) xt.xaxis[x]=yaxis[x];
	for (int y=0;y<nx;y++) xt.yaxis[y]=xaxis[y];
	for (int z=0;z<nz;z++) xt.zaxis[z]=zaxis[z];
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				xt(y,x,z)=(*this)(x,y,z);
	return xt;
}		

cube cube::xt2yt()
{
	cube xt(ny, nx, nz);
	
	for (int x=0;x<ny;x++) xt.xaxis[x]=yaxis[x];
	for (int y=0;y<nx;y++) xt.yaxis[y]=xaxis[y];
	for (int z=0;z<nz;z++) xt.zaxis[z]=zaxis[z];
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				xt(y,x,z)=(*this)(x,y,z);
	return xt;
}

cube cube::skewer(axis a, long c1, long c2) const
/* Created 3/6/2k
	Axis generalized skewer.  c1, c2 are the Uax and Dax coordinate respectively. */
{
	long sx=1, sy=1, sz=1;
	if (a==UNK) a=Z;
	if (a==X) sx=nx;
	if (a==Y) sy=ny;
	if (a==Z) sz=nz;

	cube s(sx, sy, sz);
	s.copymetadatafrom(*this);
	
	s.copyaxis(*this, a);
	s.Axis(Uax(a), 0) = Axis(Uax(a), c1);
	s.Axis(Dax(a), 0) = Axis(Dax(a), c2);
	
	for (long i(0);i<N(a);i++)
		s(a,i,Uax(a),long(0),Dax(a),long(0)) = (*this)(a,i,Uax(a),c1,Dax(a),c2);
	
	return s;
}
	
cube cube::chunky(int y1,int y2)
{
	if(y1==-1) y1=0;
	if(y2==-1) y2=ny-1;
	
	cout << "Chunking\n";
	
	cube chunked(nx, y2-y1+1, nz, 0.0);
	
	for (int x=0;x<nx;x++) chunked.xaxis[x]=xaxis[x];
	for (int y=y1;y<=y2;y++) chunked.yaxis[y-y1]=yaxis[y];
	for (int z=0;z<nz;z++) chunked.zaxis[z]=zaxis[z];
	
	for(int x=0;x<nx;x++)
		for(int z=0;z<nz;z++)
			for(int y=y1;y<=y2;y++)
				chunked(x,y-y1,z)=(*this)(x,y,z);
	return chunked;
}

cube cube::chunkz(int t1, int t2)
{
	if(t1==-1) t1 = 0;
	if(t2==-1) t2=nz-1;
	
	if (osuppress<1) cout << "Chunking in t\n";
	
	cube chunked(nx, ny, t2-t1+1, 0.0);
	
	for (int x=0;x<nx;x++) chunked.xaxis[x]=xaxis[x];
	for (int y=0;y<ny;y++) chunked.yaxis[y]=yaxis[y];
	for (int z=t1;z<=t2;z++) chunked.zaxis[z-t1]=zaxis[z];
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=t1;z<=t2;z++)
				chunked(x,y,z-t1)=(*this)(x,y,z);
	return chunked;
}

cube cube::chunkx(int x1, int x2)
{
	if(x1==-1) x1 = 0;
	if(x2==-1) x2=nx-1;
	
	if (osuppress < 1)cout << "Chunking in x from " << x1 << " to " << x2 << "\n";
	
	cube chunked(x2-x1+1, ny, nz, 0.0);
	
	for (int x=x1;x<=x2;x++) chunked.xaxis[x-x1]=xaxis[x];
	for (int y=0;y<ny;y++) chunked.yaxis[y]=yaxis[y];
	for (int z=0;z<nz;z++) chunked.zaxis[z]=zaxis[z];
	
	for(int z=0;z<nz;z++)
		for(int y=0;y<ny;y++)
			for(int x=x1;x<=x2;x++)
				chunked(x-x1,y,z)=(*this)(x,y,z);
	return chunked;
}

cube cube::chunk(axis A, long a0, long a1)
{
	if(a0==-1) a0=0;
	if(a1==-1) a1=N(A)-1;
	
	if (osuppress < 1) 
		cout << "Chunking in " << A << " from " << a0 << " to " << a1 << " --  00%";
	
	cube chunked(A,a1-a0+1, Uax(A), N(Uax(A)), Dax(A), N(Dax(A)), 0.);
	
	for (long i=a0;i<=a1;i++) chunked.Axis(A,i-a0)=Axis(A,i);
	for (long i=0;i<N(Uax(A));i++) chunked.Axis(Uax(A),i)=Axis(Uax(A),i);
	for (long i=0;i<N(Dax(A));i++) chunked.Axis(Dax(A),i)=Axis(Dax(A),i);
	
	chunked.hdr=hdr;
	
	for(long a=0;a<chunked.N(A);a++) {
		if (!osuppress) printpercent(a,chunked.N(A)-1);
		for(long c=0;c<chunked.N(Dax(A));c++)
			for(long b=0;b<chunked.N(Uax(A));b++)
				chunked(A,a,Uax(A),b,Dax(A),c)=(*this)(A,a+a0,Uax(A),b,Dax(A),c);
	}
	if (!osuppress) printpercent(12,12);
	if (!osuppress) cout << "\n";
	return chunked;
}

cube cube::chunk(axis A, double f0, double f1) // removed const due to locate dependencies
{
	long a0, a1;
	a0 = Jlocate(f0,A);
	a1 = Jlocate(f1,A)-1;
	if (osuppress < 1) cout << "Chunking in " << A << " between " << f0 << "(" << a0 << ") and ";
	if (osuppress < 1) cout << f1 << "(" << a1 << ")\n";
	
	return chunk(A, a0, a1);
}

cube cube::blocksz(cube b)
{
	return blocks(Z,b);
}
			
cube cube::blocksx(cube b)
/* Created by modifying blocksz 9/18/2k JB */
{
	if (nz != b.N(Z)  ||  ny != b.N(Y)) { 
		cout << "ERROR -- blocksx cubes of different z,y dimension\n";
		cout << "   tried to block z=" << b.N(Z) << ",y=" << b.N(Y) << " to ";
		cout << " z=" << nz << ",y=" << ny << "\n";
	}
	
	cube wall(nx+b.Nx(),ny,nz);
	wall.hdr=hdr;
	
	for (int x=0;x<nx;x++) wall.Xaxis(x) = Xaxis(x);
	for (int y=0;y<ny;y++) wall.Yaxis(y) = Yaxis(y);	
	for (int z=0;z<nz;z++) wall.Zaxis(z) = Zaxis(z);
	for (int x=nx;x<wall.Nx();x++) wall.Xaxis(x) = b.Xaxis(x-nx);
	for (int z=0;z<nz;z++){
		for (int y=0;y<ny;y++){
			for (int x=0;x<nx;x++) wall(x,y,z) = (*this)(x,y,z);
			for (int x=nx;x<wall.Nx();x++) wall(x,y,z) = b(x-nx,y,z);
		}
	}
	return wall;
}

cube cube::blocks(axis A, const cube &b)
// Generalized 2005 October 12 JB
{
	if (N(Uax(A)) != b.N(Uax(A))  ||  N(Dax(A)) != b.N(Dax(A))) {
		cout << "ERROR -- blocks("<<A<<") cubes of different Uax,Dax dimensions\n";
	}
	
	cube wall(A, N(A)+b.N(A), Uax(A), N(Uax(A)), Dax(A), N(Dax(A)));
	wall.hdr = hdr;
	wall.copyaxis(b, Uax(A));
	wall.copyaxis(b, Dax(A));
	
	for (long x=0;x<nx;x++) wall.Xaxis(x) = Xaxis(x);
	for (long y=0;y<ny;y++) wall.Yaxis(y) = Yaxis(y);	
	for (long z=0;z<nz;z++) wall.Zaxis(z) = Zaxis(z);
	for (long i=N(A);i<wall.N(A);i++) wall.Axis(A, i) = b.Axis(A, i-N(A));
	for (long i=0;i<N(Uax(A));i++) {
		for (long j=0;j<N(Dax(A));j++) {
			for (long k=0;k<N(A);k++) wall(Uax(A), i, Dax(A), j, A, k) = 
				(*this)(Uax(A), i, Dax(A), j, A, k);
			for (long k=N(A);k<wall.N(A);k++) wall(Uax(A), i, Dax(A), j, A, k) = 
				b(Uax(A), i, Dax(A), j, A, k-N(A));
		}
	}
	return wall;
}

cube cube::insert(cube& s, int x1, int y1, int z1)
{
	if (!osuppress) {
		cout << "inserting cube into " << x1 << " " << y1 << " " << z1 << "\n";
		cout << "(losing wavelength and axis scaling information)\n";
	}
	for (int x=0;x<s.Nx();x++)
		for (int y=0;y<s.Ny();y++)
			for (int z=0;z<s.Nz();z++)
				(*this)(x+x1,y+y1,z+z1) = s(x,y,z);
	return *this;
}

cube cube::resample(cube& templatecube, methodtype m)
{
	if (!osuppress) cout << "custom resampling in 3D:    00%";
	if (N(Z) != templatecube.N(Z)) cout << "RESAMPLE IN THREE AXES DOES NOT WORK FOR Z RIGHT NOW!!!!\n";
	cube answer(templatecube);
	for (int x(0);x<answer.N(X);x++) {
		if (!osuppress) printpercent(x,answer.N(X)-1);
		for (int y(0);y<answer.N(Y);y++) {
			for (int z(0);z<answer.N(Z);z++) { 
				answer(x,y,z) = (*this)(double(templatecube.Axis(X,x)),
						                  double(templatecube.Axis(Y,y)),
						                  z, m);
//				if (answer(x,y,z) != 0.) {
//					cout << "interped " << double(templatecube.Axis(X,x)) << "W ";
//					cout << double(templatecube.Axis(Y,y)) << "S as " << answer(x,y,z) << "\n";
//				}
				
			}
		}
	}
	if (!osuppress) printpercent(1,1);
	if (!osuppress) cout << "\n";
	return answer;
}

cube cube::resample(axis A, cube templatecube, interpmethod imeth)
/*		Created 2/19/2k
	Resamples data according to a predetermined new wavelength scale, set in the
	templatecube.  Only works in Z for now. */
/*		Added FluxConserve 1/25/01 */
/*    FluxConserve does not treat the edges of each integration region properly.
	Fix this in future */
/* INTEGRATE provision added 11/05/2003 for noisier data that gets upset with
   INTERPOLATE */
/* added the Axis A part so as to allow for 3d interpolation, 2006 May 22 */
{
	if (!osuppress) cout << "custom resampling in " << A << ":    00%";
	if (A!=Z) cout << "RESAMPLE IN ONE AXIS ONLY WORKS FOR Z RIGHT NOW!!!!\n";
	
	cube re(A, templatecube.N(A), Uax(A), N(Uax(A)), Dax(A), N(Dax(A)));
	re.hdr=hdr;
	re.copyaxis(templatecube,A);
	
	cube step2(*this);
	

	if (imeth==FLUXCONSERVE || imeth==INTEGRATE) {
		for (int x=0;x<re.nx;x++){
			for (int y=0;y<re.ny;y++){
				for (int z=0;z<re.nz;z++) {
					float lo, hi;
					if (!osuppress) printpercent(re.nx*re.ny*z,re.nx*re.ny*re.nz);
					if (z==0) {
						lo=templatecube.Axis(A,0)-
								(templatecube.Axis(A,1)-templatecube.Axis(A,0))/2.;
						hi=(templatecube.Axis(A,1)+templatecube.Axis(A,0))/2.;
					}
					else if (z==re.nz-1) {
						lo=(templatecube.Zaxis(z)+templatecube.Zaxis(z-1))/2.;
						hi=templatecube.Zaxis(z)+
								(templatecube.Zaxis(z)-templatecube.Zaxis(z-1))/2.;
					}
					else {
						lo=(templatecube.Zaxis(z)+templatecube.Zaxis(z-1))/2.;
						hi=(templatecube.Zaxis(z)+templatecube.Zaxis(z+1))/2.;
					}
					
					int loi, hii;
					loi=Jlocate(lo,Z);
					hii=Jlocate(hi,Z);
					int n(hii-loi);
					for (int zz=loi;zz<hii;zz++) re(x,y,z)+= (*this)(x,y,zz);	
					if (imeth==INTEGRATE) {
						if (n) re(x,y,z) /= (float)n;
						// If there's nothing to integrate, i.e. a gap in the data, interpolate
						else re(x,y,z)=step2.interp(re.Zaxis(z), 1, Z);
					}
				}
			}
		}
	}
	
	if (imeth==INTERPOLATE){
		for (int x=0;x<re.nx;x++){
			//cout << "N x:\t"<<re.nx<<endl;
			for (int y=0;y<re.ny;y++){
				//cout << "N y:\t"<<re.ny<<endl;
				for (int z=0;z<re.nz;z++) {
					//cout << "N z:\t"<<re.nz<<endl;
					//cout << "re.Zaxis(z):\t"<<re.Zaxis(z)<<endl;
					if (!osuppress) printpercent(re.nx*re.ny*z,re.nx*re.ny*re.nz);
					re(x,y,z)=step2.interp(re.Zaxis(z), 3, Z);
				}
			}
		}
	}
		
	if (!osuppress) cout << char(8)<<char(8)<<char(8)<<"100%\n";
	
	re.hdr=templatecube.hdr;
	re.hdr["Nx"]=int2str(re.Nx());
	re.hdr["Ny"]=int2str(re.Ny());
	re.hdr["Nz"]=int2str(re.Nz());
	return re;
}


int imin(int x1, int x2) 
{
	if (x1<x2) return x1;
	else return x2;
}

cube cube::resample(int xres, int yres, int zres)
/*  Created 2/23/2k.
	This version of resample coadds every res pixels in each dimension into a
	smaller cube for positive Nres, and replicates pixels to make images
	larger for negarive Nres. */
{
	float fxres=1.,fyres=1.,fzres=1.;
	if (xres > 0) fxres=xres; else fxres=-1./(float)xres;
	if (yres > 0) fyres=yres; else fyres=-1./(float)yres;
	if (zres > 0) fzres=zres; else fzres=-1./(float)zres;
	
	int snx=(int)((float)nx/fxres);
	int sny=(int)((float)ny/fyres);
	int snz=(int)((float)nz/fzres);
	if (nx % xres   &&   xres < 0) snx++;
	if (ny % yres   &&   yres < 0) sny++;
	if (nz % zres   &&   zres < 0) snz++;
	cube coadded(snx, sny, snz, 0.0);
	coadded.hdr=hdr;
	coadded.cubetype=cubetype;
	coadded.byteorder=byteorder;
	
	if (osuppress < 0) cout << "coadding/replicating ("<<nx<<","<<ny<<","<<nz<<") by "<<xres<<"x, "<<yres<<"y, ";
	if (osuppress < 0) cout << zres<<"z  to ("<<snx<<","<<sny<<","<<snz<<")\n";
	
// sx is x in subcube, x is in the normal cube
	
// Here let the wavelength values for the new cube be the average of the old ones.
	for (int sx=0;sx<snx;sx++){
		if (xres>0) {  // coadding case
			coadded.Xaxis(sx) = 0.;
			for ( int x=sx*xres; x<imin((sx+1)*xres, nx);x++) coadded.Xaxis(sx) += Xaxis(x);
			if (sx==snx-1 && nx%xres) coadded.Xaxis(sx) /= nx%xres;
			else coadded.Xaxis(sx) /= xres;
		} else { // replicating case -- need to improve someday
			coadded.Xaxis(sx) = Xaxis((int)((float)sx*fxres));
		}
	}
	for (int sy=0;sy<sny;sy++){
		if (yres>0) {
			coadded.Yaxis(sy) = 0.;
			for ( int y=sy*yres; y<imin((sy+1)*yres, ny);y++) coadded.Yaxis(sy) += Yaxis(y);
			if (sy==sny-1 && ny%yres) coadded.Yaxis(sy) /= ny%yres;
			else coadded.Yaxis(sy) /= yres;
		} else {
			coadded.Yaxis(sy) = Yaxis((int)((float)sy*fyres));
		}	
	}	
	for (int sz=0;sz<snz;sz++){
		if (zres>0) {
			coadded.Zaxis(sz) = 0.;
			for ( int z=sz*zres; z<imin((sz+1)*zres, nz);z++) coadded.Zaxis(sz) += Zaxis(z);
			if (sz==snz-1 && nz%zres) coadded.Zaxis(sz) /= nz%zres;
			else coadded.Zaxis(sz) /= zres;
		} else {
			coadded.Zaxis(sz) = Zaxis((int)((float)sz*fzres));
		}
	}	
	
// Perform coadding/replication into new cube
	int qcfx(-(xres<-4));  // quantized correction factor. 
	int qcfy(-(yres<-4));  // b/c of finite differences in floats?  who knows!
	int qcfz(-(zres<-4));  // empirically determined.
	
	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				for (int sx=(int)(x/fxres) ; sx<((float)x/fxres)+1./fxres+qcfx ; sx++)
					for (int sy=(int)(y/fyres) ; sy<((float)y/fyres)+1./fyres+qcfy ; sy++)
						for (int sz=(int)(z/fzres) ; sz<((float)z/fzres)+1./fzres+qcfz ; sz++)
							coadded(sx,sy,sz) += (*this)(x,y,z);
	
// correct for leftover bins
	if (nx%xres) {
		float correctionfactor;
		correctionfactor=(float)xres/float(nx%xres);
		for (int z=0;z<snz;z++)
			for (int y=0;y<sny;y++)
				coadded(snx-1,y,z) *= correctionfactor;
	}
	if (ny%yres) {
		float correctionfactor;
		correctionfactor=float(yres)/float(ny%yres);
		for (int x=0;x<snx;x++)
			for (int z=0;z<snz;z++)
				coadded(x,sny-1,z) *= correctionfactor;
	}
	if (nz%zres) {
		float correctionfactor;
		correctionfactor=float(zres)/float(nz%zres);
		for (int y=0;y<sny;y++)
			for (int x=0;x<snx;x++)
				coadded(x,y,snz-1) *= correctionfactor;
	}
			
	
//	coadded=coadded(0,snx-1, 0,sny-1, 0, snz-2);
	return coadded;
}

cube cube::resample(double xlo, double xhi, double xres, double ylo,
						  double yhi, double yres, double zlo,
						  double zhi, double zres, bool fluxcnsrv)
{
	if (osuppress < 0) cout << "resampling by interpolation\n";
	int xdir=1,ydir=1,zdir=1;				// dir of 1 means increasing, -1 means decreasing
	if (nx>1  &&  xaxis[0] > xaxis[1]) xdir=-1;
	if (ny>1  &&  yaxis[0] > yaxis[1]) ydir=-1;
	if (nz>1  &&  zaxis[0] > zaxis[1]) zdir=-1;
//  setting default values
	if (xlo==-1. && xdir>0) xlo = xaxis[0];
	if (xlo==-1. && xdir<0) xlo = xaxis[nx-1];
	if (xhi==-1. && xdir>0) xhi = xaxis[nx-1];
	if (xhi==-1. && xdir<0) xhi = xaxis[0];
	if (ylo==-1. && ydir>0) ylo = yaxis[0];
	if (ylo==-1. && ydir<0) ylo = yaxis[ny-1];
	if (yhi==-1. && ydir>0) yhi = yaxis[ny-1];
	if (yhi==-1. && ydir<0) yhi = yaxis[0];
	if (zlo==-1. && zdir>0) zlo = zaxis[0];
	if (zlo==-1. && zdir<0) zlo = zaxis[nz-1];
	if (zhi==-1. && zdir>0) zhi = zaxis[nz-1];
	if (zhi==-1. && zdir<0) zhi = zaxis[0];
	if (xres==-1.)xres= 0;
	if (yres==-1.)yres= 0;
	if (zres==-1.)zres= 0;
	
	
	int rx, ry, rz;
	bool xsub, ysub, zsub;
	if (xres) rx = (int) (((xhi-xlo)+xres)/xres);
	else { rx=nx; xsub=0; }
	if (yres) ry = (int) (((yhi-ylo)+yres)/yres);
	else { ry=ny; ysub=0; }
	if (zres) rz = (int) (((zhi-zlo)+zres)/zres);
	else { rz=nz; zsub=0; }
	
	if (osuppress < 0) cout << "Old cube size          " << nx << "," << ny << "," << nz << "\n";
	if (osuppress < 0) cout << "Creating new cube size " << rx << "," << ry << "," << rz << "\n";
	if (osuppress < 0) cout << "New wavelength range ("<<xlo<<"-"<<xhi<<","<<ylo<<"-"<<yhi<<",";
	if (osuppress < 0) cout <<								 zlo<<"-"<<zhi<<")\n";
	cube re(rx,ry,rz);
	for (int x=0;x<rx;x++){
		if (xres) re.Axis(X,x) = xlo + x*xres;
		else re.Axis(X,x) = xaxis[x];
	}
	for (int y=0;y<ry;y++){
		if (yres) re.Axis(Y,y) = ylo + y*yres;
		else re.Axis(Y,y) = yaxis[y];
	}
	for (int z=0;z<rz;z++){
		if (zres) re.Axis(Z,z) = zlo + z*zres;
		else re.Axis(Z,z) = zaxis[z];
	}
	
	if (!fluxcnsrv){
		if (osuppress < 0) cout << "Proceeding with step 1\n";
		cube step1(rx,ny,nz);
		step1.copyaxis(*this);
		int order;
		order = (nx <= 4) ? nx-1 : 3;
		if (osuppress < 0) cout << "entering loop\n";
		
		if (xres){
			for (int z=0;z<nz;z++){
				printpercent(z,nz);
				for (int y=0;y<ny;y++){
					cube thisx;	
					if (osuppress < 0) cout << "subsampling (0.." << nx-1 << "," << y << "," << z << ")";
					if (osuppress < 0) cout << (*this)(0,0,z) << "\n";
					thisx=(*this)(0,nx-1,y,y,z,z);
					if (osuppress < 0) cout << "subsampled\n";
					if (osuppress < 0) cout << "axis copied\n";
					thisx=thisx.xt2yt();
					for (int x=0;x<rx;x++){
						step1(x,y,z) = thisx.interp(re.Xaxis(x), order);
					}
				}
			}
		} else
			step1=*this;
			
//		cout << "step2\n";
		cube step2(rx,ry,nz);
		step2.copyaxis(re);
		
		if (yres){
			for (int x=0;x<rx;x++){
				for (int y=0;y<ry;y++){
					for (int z=0;z<nz;z++){
						cube thisy;
						thisy=step1(x,x,0,ny-1,z,z);
						step2(x,y,z) = thisy.interp(re.Yaxis(y));
					}
				}
			}
		} else
			step2=step1;
						
		if (osuppress < 0) cout << "last step\n  00%";
		for (int x=0;x<rx;x++){
			for (int y=0;y<ry;y++){
				cube thisz;
				thisz=step2(x,x,y,y,0,-1);
//				thisz=thisz.xt2xy();
				for (int z=0;z<rz;z++){
					printpercent(rx*ry*z,rx*ry*rz);					
					if (osuppress < 0) cout << re.Zaxis(z) << ", "; cout.flush();
//					re(x,y,z) = thisz.interp(re.Zaxis(z));
					re(x,y,z) = thisz(0,0,re.Zaxis(z), mLINEAR);
 					if (osuppress < 0) cout << re(x,y,z) << "\n";
				}
			}
		}
	}
	
	if (fluxcnsrv){
		cout << "resampling by integrating --  00%";
		
		int flx;
		if (Dtype()==FLUX || Dtype()==FLUXCGS || Dtype()==FLUXCGSWAVENUM) flx=1;
		
		int lh,ll;
		float ff, ffnet=0;									// flux factor (to conserve flux)
		float begw, endw;										// beginning and end of window
		cube step2(*this);									// to be coadded together
		
		for (int x=0;x<rx;x++){
			for (int y=0;y<ry;y++){
				cube thisz;
				thisz=step2(x,x,y,y,0,-1);
				for (int z=0;z<rz;z++){
					printpercent(rx*ry*z,rx*ry*rz);
					
					if (z == 0) begw=re.Zaxis(0); /* - (re.Zaxis(1)-re.Zaxis(0))/2.;*/
					else begw=(re.Zaxis(z) + re.Zaxis(z-1)) / 2.;
					if (z == rz-1) endw=re.Zaxis(rz-1);    /*+(re.Zaxis(nz-1)+re.Zaxis(nz-2))/2.;*/
					else endw=(re.Zaxis(z) + re.Zaxis(z+1)) / 2.;
					
					ll = thisz.Jlocate(begw, Z);
					lh = thisz.Jlocate(endw, Z);
					
//					cout << "Integrating from "<<begw<<"("<<ll<<") to "<<endw<<"("<<lh<<")\n";
					// left end
					re(x,y,z)=(thisz.interp(begw,1)+thisz(0,0,ll+1))/2.;
//					cout << "starting re() = " << re(x,y,z) << "\n";
//					cout << "("<<thisz.interp(begw,1)<<"+"<<thisz(0,0,ll+1)<<")/2.\n";
					ff = (thisz.Zaxis(ll+1)-begw)/(thisz.Zaxis(ll+1)-thisz.Zaxis(ll));
					re(x,y,z)*=ff;
					ffnet=ff;
//					cout << "after left bit, value="<<re(x,y,z)<<" ff="<<ff<<"\n";
					
					// middle parts
					for (int l=ll+1;l<lh;l++){
						re(x,y,z)+=(thisz(0,0,l)+thisz(0,0,l+1))/2.;
						ffnet++;
//						cout << "after a middle, value="<<re(x,y,z)<<" (ff=1)	ffnet="<<ffnet<<"\n";
					}
					
					// right end
					ff = -(endw-thisz.Zaxis(lh)) / (thisz.Zaxis(lh+1)-thisz.Zaxis(lh));
					re(x,y,z)+=(thisz.interp(endw,1)+thisz(0,0,lh))/2. * ff;
					ffnet+=ff;
//					cout << "after right bit, value="<<re(x,y,z)<<" ff="<<ff<<" ffnet="<<ffnet<<"\n";
					
					if (flx) re(x,y,z)/=ffnet;
				}
			}
		}
	}
		
		
		
		
	re.hdr=hdr;
	re.keyword("Nx", int2str(re.Nx()));
	re.keyword("Ny", int2str(re.Ny()));
	re.keyword("Nz", int2str(re.Nz()));
	cout << "re.nx = " << re.Nx() << " re.ny = " << re.Ny() << " re.nz = " <<re.Nz()<<"\n";
	
	return re;
}
	
cube cube::resample(axis a, float binsize)
/* JB 2005 January 6
	Resamples by BINNING, taking the average value in each bin of size binsize*/		
{
	if (!osuppress) cout << "Resampling by binning\n";
	
	axis a1(Uax(a)); 		// for a=Z, this becomes X
	axis a2(Dax(a));		// for a=Z, this is Y
	
	cube original(this->dirinc(a)); // make it easy to work with
	cout << "orignal info():  " << original.info();
	
// create the target cube, a little big:  faster to just subcube it later than to
// grow it as we go
	cube resampled(a1, N(a1), a2, N(a2), a, 
			int((original.Axis(a,N(a)-1)-original.Axis(a,0))/binsize + 2),
			0., cubetype);
	for (long m=0;m<N(a1);m++) resampled.Axis(a1,m) = Axis(a1,m);
	for (long m=0;m<N(a2);m++) resampled.Axis(a2,m) = Axis(a2,m);
	
// start before the beginning, count in bins until past the end.
	float start(float(long(original.Axis(a,0)/binsize - 3.)) * binsize);
	float end(float(long(original.Axis(a,N(a)-1)/binsize + 2.)) * binsize);
	long i=0;    // keeps track of where we are in the orignal cube.
	long j=0;		// keeps track of where we are in the resampled cube
	for (float x0=start;x0<=end, i<N(a);x0+=binsize) { // x0 is where this bin starts
		long n=0; // keeps track of how many pixels we're coadding this bin.
		while (i<N(a)  &&  original.Axis(a,i) < x0+binsize) {
			n++; 
			double sum=0;
			for (long z1=0;z1<N(a1);z1++)
				for (long z2=0;z2<N(a2);z2++) {
					resampled(a1,z1, a2,z2, a,j) += original(a1,z1, a2,z2, a,i);
				}
			i++;
		}
		if (n) {  // only send a value if there are points in this bin
			for (long z1=0;z1<N(a1);z1++)
				for (long z2=0;z2<N(a2);z2++)
					resampled(a1,z1, a2,z2, a,j) /= n;
			resampled.Axis(a,j) = x0 + binsize/2.;
			j++;
		}
		
	}
	cout << "before subcube, " << resampled.info();
	resampled = resampled(a, 0, j-1, Uax(a), -1, -1, Dax(a), -1, -1);
	cout << "after subcube, " << resampled.info();
	
	return resampled;			
}


cube cube::convert(axis a, axistype t, double param, double param2)
// Altered for new doubles in axis values, 2012 August 22
{	
	axistype f;
	int n;
	double *A;
	if (a==X){
		A=&Axis(X,0);
		n=nx;
		f=Xtype();
	}
	if (a==Y){
		A=&Axis(Y,0);
		n=ny;
		f=Ytype();
	}
	if (a==Z){
		A=&Axis(Z,0);
		n=nz;
		f=Ztype();
	}
	
	cout << "Converting axis "<<a<<" from "<<f<<" to "<<t<<"\n";
	
//  STEP 1:  convert whatever it is to meters.

	double factor=1;
	if (f==MICRONS) factor=1e-6;
	if (f==CM) factor=.01; 
	if (f==ARRAYVAL) 
			cout << "Assuming that axis with ARRAYVAL is in meters (bad assumption probably)\n";
	for (int m=0;m<n;m++) {
		if (f!=WAVENUM) A[m]*=factor;
		else A[m]=.01/A[m];
	}
	
//  STEP 2:  convert to target
	
	factor=1;
	if (t==MICRONS) factor=1e6;
	if (t==CM) factor=100.;
	if (t==ARRAYVAL)
		cout << "Assuming that axis with ARRAYVAL is in meters (bad assumption probably)\n";
	for (int m=0;m<n;m++) {
		if (t!=WAVENUM) A[m]*=factor;
		else A[m]=.01/A[m];
	}
	if (a==X) Xtype(t);
	if (a==Y) Ytype(t);
	if (a==Z) Ztype(t);
	return *this;
}
	
cube cube::convert(datatype t, double param, double param2)    
												//  How param is used is different based
												//	 on the conversion that is being performed.
												//   --> PHOTONS        param = effective telescope
												//										  area
												//								param2 = integration time
{
	double area, itime, factor, lambda;
	cout << "Converting data from type " << Dtype() << " to type " << t << "\n";
//	if (Dtype() == FLUX  &&  t == PHOTONS){
	if (t == PHOTONS) {
		area=param;
		if (area==-1) area=1;
		itime=param2;
		if (itime==-1) itime=1;
	}
		
	axis a=UNK;

	if (Xtype() == MICRONS || Xtype() == WAVENUM) a=X;
	if (Ytype() == MICRONS || Ytype() == WAVENUM) a=Y;
	if (Ztype() == MICRONS || Ztype() == WAVENUM) a=Z;
	if (a == UNK) {
		cout << "No suitable wavelength axis, convert failure.\n";
		return *this;
	}
		
	cout << "Using axis " << a << " for wavelength data\n";
	double nrg, bandpass;
	double b1, b2;

	if ((Dtype() == FLUXCGS || Dtype() == FLUXCGSWAVENUM || Dtype() == FLUXCGSHZ) && t == FLUX) 
		b1 = 1e-3;     // 1e-7 W / (erg/s)   10^4 cm^2 / m^2
	
//   Ok so this isn't elegant but it worx
	if (a==X){		
		for (int x=0;x<nx;x++) {
			lambda = Xaxis(x);
			if (Dtype() == FLUXCGS && t == FLUX) b2 = 1e-4;
			if (Dtype() == FLUXCGSWAVENUM && (t == FLUX || t==FLUXCGS)){
				double wn0, wn1, um0, um1, dum;	// convert wavelength range
				wn0=lambda-.5;							// from per wavenum to per um
				wn1=lambda+.5;
				um0=(1./wn0)*1e4;
				um1=(1./wn1)*1e4;
				dum=fabs(um0-um1);
				b2=1./dum;
				if (t==FLUXCGS) b2*=1e4;
			}
			if (Dtype() == FLUXCGSHZ && (t == FLUX || t == FLUXCGS)){
				double fq0, fq1, um0, um1, dum;	// convert wavelength range
				fq0=3e14/(lambda)-.5;							// from per hz to per um
				fq1=3e14/(lambda)+.5;
				um0=(3e14/fq0);
				um1=(3e14/fq1);
				dum=fabs(um0-um1);
				b2=1./dum;
				if (t==FLUXCGS) b2*=1e4;
			}	
				
			if (Dtype() == FLUX) nrg = 6.626e-34 * 3e8 / (1e-6 * lambda);   
																//  E = h*nu = h*c/lambda
			if (Dtype() == FLUXCGSWAVENUM) nrg = 6.626e-27 * 3e10 * lambda;
			if (Dtype() == FLUXCGSHZ) nrg = 6.626e-34 * 3e8 / (1e-6 * lambda);   
			if (t == PHOTONS){
				if (x==0) bandpass = fabs(Xaxis(1) - Xaxis(0));
				else if (x==nx-1) bandpass == fabs(Xaxis(nx-1) - Xaxis(nx-2));
				else bandpass = fabs(Xaxis(x+1) - Xaxis(x-1)) / 2;
			}
			if (Dtype() == FLUX || Dtype() == FLUXCGSHZ) bandpass *= 1e-6;
			if (t == PHOTONS) factor = area*bandpass*itime/nrg;
			if (t == FLUX) factor=b1*b2;
			for (int y=0;y<ny;y++) {
				for (int z=0;z<nz;z++){
					(*this)(x,y,z) *= factor;
				}
			}
		}
	}
	if (a==Y){		
		for (int y=0;y<ny;y++) {
			lambda = Yaxis(y);
			if (Dtype() == FLUXCGS && t == FLUX) b2 = 1e-4;
			if (Dtype() == FLUXCGSWAVENUM && (t == FLUX || t==FLUXCGS)){
				double wn0, wn1, um0, um1, dum;	// convert wavelength range
				wn0=lambda-.5;							// from per wavenum to per um
				wn1=lambda+.5;
				um0=(1./wn0)*1e4;
				um1=(1./wn1)*1e4;
				dum=fabs(um0-um1);
				b2=1./dum;
				if (t==FLUXCGS) b2*=1e4;
			}			
			if (Dtype() == FLUXCGSHZ && (t == FLUX || t == FLUXCGS)){
				double fq0, fq1, um0, um1, dum;	// convert wavelength range
				fq0=3e14/(lambda)-.5;							// from per hz to per um
				fq1=3e14/(lambda)+.5;
				um0=(3e14/fq0);
				um1=(3e14/fq1);
				dum=fabs(um0-um1);
				b2=1./dum;
				if (t==FLUXCGS) b2*=1e4;
			}

			if (Dtype() == FLUX) nrg = 6.626e-34 * 3e8 / (1e-6 * lambda);   
																	//  E = h*nu = h*c/lambda
			if (Dtype() == FLUXCGSWAVENUM) nrg = 6.626e-27 * 3e10 * lambda;
			if (Dtype() == FLUXCGSHZ) nrg = 6.626e-34 * 3e8 / (1e-6 * lambda); 
			if (t == PHOTONS){
				if (y==0) bandpass = fabs(Yaxis(1) - Yaxis(0));
				else if (y==ny-1) bandpass == fabs(Yaxis(ny-1) - Yaxis(ny-2));
				else bandpass = fabs(Yaxis(y+1) + Yaxis(y-1)) / 2.;
			}
			if (Dtype() == FLUX || Dtype() == FLUXCGSHZ) bandpass *= 1e-6;
			if (t == PHOTONS) factor = area*bandpass*itime/nrg;	
			if (t == FLUX) factor=b1*b2;
			for (int x=0;x<nx;x++){
				for (int z=0;z<nz;z++){
					(*this)(x,y,z) *= factor;
				}
			}
		}
	}
	if (a==Z){		
		for (int z=0;z<nz;z++) {
			lambda = Zaxis(z);
			if (Dtype() == FLUXCGS && t == FLUX) b2 = 1e-4;
			if (Dtype() == FLUXCGSWAVENUM && (t == FLUX || t==FLUXCGS)){
				double wn0, wn1, um0, um1, dum;	// convert wavelength range
				wn0=lambda-.5;							// from per wavenum to per um
				wn1=lambda+.5;
				um0=(1./wn0)*1e4;
				um1=(1./wn1)*1e4;
				dum=fabs(um0-um1);
				b2=1./dum;
				if (t==FLUXCGS) b2*=1e4;
			}			
			if (Dtype() == FLUXCGSHZ && (t == FLUX || t == FLUXCGS)){
				double fq0, fq1, um0, um1, dum;	// convert wavelength range
				fq0=3e14/(lambda)-.5;							// from per hz to per um
				fq1=3e14/(lambda)+.5;
				um0=(3e14/fq0);
				um1=(3e14/fq1);
				dum=fabs(um0-um1);
				b2=1./dum;
				if (t==FLUXCGS) b2*=1e4;
			}
			if (lambda==0) cout << "WAVELENGTH = 0 in convert!\n";
			if (Dtype() == FLUX) nrg = 6.626e-34 * 3e8 / (1e-6 * lambda);   
																//  E = h*nu = h*c/lambda
			if (Dtype() == FLUXCGSWAVENUM) nrg = 6.626e-27 * 3e10 * lambda;
			if (Dtype() == FLUXCGSHZ) nrg = 6.626e-34 * 3e8 / (1e-6 * lambda); 
			if (t == PHOTONS){
				if (z==0) bandpass = fabs(Zaxis(1) - Zaxis(0));
				else if (z==nz-1) bandpass = fabs(Zaxis(nz-1) - Zaxis(nz-2));
				else bandpass = fabs(Zaxis(z+1) - Zaxis(z-1)) / 2.;
			}				
//			if (Dtype() == FLUX || Dtype() == FLUXCGSHZ) bandpass *= 1e-6;
			
			if (t == PHOTONS) {
				factor = area*bandpass*itime/nrg;
//				cout << area << "(area)*" << bandpass<<"(bandpass)*"<<itime<<"(inttime)/";
//				cout << nrg << "(energy)="<<factor<<"(factor) for l="<<lambda<<"\n";
			}
			if (t == FLUX) factor=b1*b2;	
			for (int y=0;y<ny;y++) {
				for (int x=0;x<nx;x++){
					(*this)(x,y,z) *= factor;
				}
			}
		}
	}
	if (Dtype() == FLUXCGS || Dtype() == FLUXCGSWAVENUM && t == FLUX) convert(a,MICRONS);
	Dtype(t);
	cout << "Conversion successful\n";
	return *this;
}

cube cube::tovicar()
{
	if(cubetype==vicar)return *this;
	
	cube v(nx, ny, nz, 0.0, vicar);
	
	for (int x=0;x<nx;x++) v.xaxis[x]=xaxis[x];
	for (int y=0;y<ny;y++) v.yaxis[y]=yaxis[y];
	for (int z=0;z<nz;z++) v.zaxis[z]=zaxis[z];
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				v(x,y,z)=(*this)(x,y,z);
	return v;
}

cube cube::tofits()
{
	if(cubetype==fits)return *this;
	
	cube v(nx, ny, nz, 0.0, fits);
	
	for (int x=0;x<nx;x++) v.xaxis[x]=xaxis[x];
	for (int y=0;y<ny;y++) v.yaxis[y]=yaxis[y];
	for (int z=0;z<nz;z++) v.zaxis[z]=zaxis[z];
	v.hdr=hdr;
	v.byteorder = byteorder;
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				v(x,y,z)=(*this)(x,y,z);
	return v;
}

cube cube::tobsqraw()
{
	if(cubetype==bsqraw) return *this;
	
	cube r(nx,ny,nz,0.0,bsqraw);
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				r(x,y,z)=(*this)(x,y,z);
	return r;
}

cube cube::totext()
{
	if(cubetype==text)return *this;
	
	cube t(nx, ny, nz, 0.0, text);
	
	for (int x=0;x<nx;x++) t.xaxis[x]=xaxis[x];
	for (int y=0;y<ny;y++) t.yaxis[y]=yaxis[y];
	for (int z=0;z<nz;z++) t.zaxis[z]=zaxis[z];
	
	for(int x=0;x<nx;x++)
		for(int y=0;y<ny;y++)
			for(int z=0;z<nz;z++)
				t(x,y,z)=(*this)(x,y,z);
	return t;
}

cube cube::toJcube(int tolambda)			// tolambda's 3 least significant bits
													// signify whether to do this in the 3 axes
													// 2^0 = x   2^1 = y   2^2 = z
{
	if (cubetype==Jcube1) return *this;
	
	int ournx=nx, ourny=ny, ournz=nz;
	int doz=0, doy=0, dox=0;
	if (tolambda/4.>=1) {
		doz=1;
		ournx-=1;
		tolambda-=4;
	}
	if (tolambda/2.>=1) {
		doy=1;
		if (!doz)ournx-=1;
		tolambda-=2;
	}
	if (tolambda==1){
		dox=1;
		ourny-=1;
		tolambda=0;
	}
	cube j(ournx, ourny, ournz);
	j.hdr = hdr;
	if (tolambda==-1){
	
		for (int x=0;x<nx;x++) j.xaxis[x]=xaxis[x];
		for (int y=0;y<ny;y++) j.yaxis[y]=yaxis[y];
		for (int z=0;z<nz;z++) j.zaxis[z]=zaxis[z];
	
		for(int x=0;x<nx;x++)
			for(int y=0;y<ny;y++)
				for(int z=0;z<nz;z++)
					j(x,y,z)=(*this)(x,y,z);
	} else {
		
		if (dox) for (int x=0;x<nx;x++) j.xaxis[x]=(*this)(x,0,0);
		else for (int x=0;x<nx;x++) j.xaxis[x]=xaxis[x];
		if (doy) for (int y=0;y<ny;y++) j.yaxis[y]=(*this)(0,y,0);
		else for (int y=0;y<ny;y++) j.yaxis[y]=yaxis[y];
		if (doz) for (int z=0;z<nz;z++) j.zaxis[z]=(*this)(0,0,z);
		else for (int z=0;z<nz;z++) j.zaxis[z]=zaxis[z];
		
		cout << "dox,doy,doz "<<dox<<','<<doy<<','<<doz<<'\n';		
		cout << "nx ,ny ,nz  "<<nx <<','<<ny <<','<<nz <<'\n';
		for(int x=(doz==1 || doy==1);x<nx;x++)
			for(int y=dox;y<ny;y++)
				for(int z=0;z<nz;z++)
					j(x-(doz==1 || doy==1),y-dox,z)=(*this)(x,y,z);
	}
	
	return j;
}

cube cube::totif()
{
	cube t(*this);
	
	t.cubetype=tif;
	return t;
}

cube cube::astype(filetype f) const
{
	cube answer(*this);
	answer.cubetype=f;
	return answer;
}

cube& cube::totype(filetype f)
{
	cubetype = f;
	return *this;
}

cube cube::asstorageorder(axis S, axis M, axis F) const
{
	cube answer(nx,ny,nz,0.,cubetype,S,M,F,0);
	answer.copymetadatafrom(*this);
	answer.storageSLOW=S;
	answer.storageMEDIUM=M;
	answer.storageFAST=F;
	answer.copyaxis(*this,X);
	answer.copyaxis(*this,Y);
	answer.copyaxis(*this,Z);
	
	cout << "Converting storage order to " << S << M << F << ":   00%";
	for (long s=0;s<N(storageSLOW);s++) {
		printpercent(s,N(storageSLOW)-1);
		for (long m=0;m<N(storageMEDIUM);m++) {
			for (long f=0;f<N(storageFAST);f++) {
				answer(storageSLOW,s,storageMEDIUM,m,storageFAST,f) = 
						(*this).storageorderaccess(s,m,f);
			}
		}
	}		
	printpercent(1,1);
	cout << "\n";
	
	return answer;
}

cube& cube::tostorageorder(axis S, axis M, axis F)
{
	(*this) = asstorageorder(S,M,F);   // could do this more elegantly in future
	return *this;
}

cube cube::txt2spec()
/* JB 10/31/2001
	converts annoying txt files 2,3x1xn into 1x1xn with proper axis values.  Works 
	for anything with a 2(3, 1, and n, not just x,y,z (i.e. 1x2,3xn works too)
			(axis independant)
	*/
{
	axis longaxis(Z), valuedaxis(X);  // defaults so that it doesn't totally bomb
	int ln;
	if (N(X) == 2  ||  N(X) == 3) valuedaxis=X;
	if (N(Y) == 2  ||  N(Y) == 3) valuedaxis=Y;
	if (N(Z) == 2  ||  N(Z) == 3) valuedaxis=Z;
	if (N(X) > 3) { longaxis=X; ln=N(X); }
	if (N(Y) > 3) { longaxis=Y; ln=N(Y); }
	if (N(Z) > 3) { longaxis=Z; ln=N(Z); }
	
	cout << "creating result cube -- longaxis=" << longaxis << ", ln=" << ln << "\n";
	cout << "                      valuedaxis=" << valuedaxis << ", ";
	cout << " n = " << N(valuedaxis)-1 << "\n";
	
	cube result(longaxis, ln, valuedaxis, N(valuedaxis)-1);
	cout << "Result cube made\n";
	for (int i=0;i<ln;i++) {
		result.Axis(longaxis,i) = (*this)(longaxis, i);
		result(longaxis, i) = (*this)(longaxis, i, valuedaxis, 1);
		if (N(valuedaxis) == 3)
			result(longaxis, i, valuedaxis, 1) = (*this)(longaxis, i, valuedaxis, 2);
			
	}
	
	return result;	
}


cube cube::byteorderinvert()
{
	cout << "Inverting byte order --   00%";
	char tmp, *pbyte;
	cube answer(*this);
	
	for (cube::iterator i=begin();i!=end();i++) {
		pbyte=(char*)(void*)&(*i);
		
		tmp = pbyte[0];
		pbyte[0] = pbyte[3];
		pbyte[3] = tmp;
		
		tmp = pbyte[1];
		pbyte[1] = pbyte[2];
		pbyte[2] = tmp;
	}
	cout << "\n";
	
	if (answer.byteorder==bigendian) answer.byteorder=littleendian;
	else answer.byteorder=bigendian;
	return answer;
	
}

void cube::setendian(endian e)
{
	byteorder=e;
}

cube cube::tolittleendian()
{
	cube l(*this);
	if (byteorder==littleendian) return *this;
	else {
		l=byteorderinvert();
		return l;
	}
}


cube cube::dirinc()
{
	int xdir=1,ydir=1,zdir=1;				// dir of 1 means increasing, -1 means decreasing
	for (int x=0;x<N(X)-1;x++) {
		if (xdir==1  &&  xaxis[x] > xaxis[x+1]) xdir=-1;
		if (xdir==-1 &&  xaxis[x] < xaxis[x+1]) { xdir=0; break; }
	}
	for (int y=0;y<N(Y)-1;y++) {
		if (ydir==1  &&  yaxis[y] > yaxis[y+1]) ydir=-1;
		if (ydir==-1 &&  yaxis[y] < yaxis[y+1]) { ydir=0; break; }
	}	
	for (int z=0;z<N(Z)-1;z++) {
		if (zdir==1  &&  zaxis[z] > zaxis[z+1]) zdir=-1;
		if (zdir==-1 &&  zaxis[z] < zaxis[z+1]) { zdir=0; break; }
	}
	if (xdir==1 && ydir==1 && zdir==1) return (*this);

	
	cube inc(*this);
	int ox,oy,oz;
	cout << "Changing direction to increasing:   00%";
	if (xdir*ydir*zdir) {  // if the whole thing seems decreasing, just reverse
		for (int x=0;x<nx;x++){
			if (xdir>=0) ox=x; else ox=nx-1-x;
			inc.Xaxis(x)=Xaxis(ox);
			
			for (int y=0;y<ny;y++){
				if (ydir>=0) oy=y; else oy=ny-1-y;
				inc.Yaxis(y)=Yaxis(oy);
				
				for (int z=0;z<nz;z++){
					printpercent(nx*ny*z,nx*ny*nz);
					if (zdir>=0) oz=z; else oz=nz-1-z;
					inc.Zaxis(z)=Zaxis(oz);
					
					inc(x,y,z)=(*this)(ox,oy,oz);
				}
			}
		}
	}
	
	for (int n=0;n<3;n++) {
		axis a; int dir;
		if (n==0) { a=X; dir=xdir; }
		if (n==1) { a=Y; dir=ydir; }
		if (n==2) { a=Z; dir=zdir; }
		axis a1=Uax(a);
		axis a2=Dax(a);
		if (!dir) {
//			osuppress++;
			inc = inc.dirinc(a);
//			cout << "Bubblesorting axis " << a << "\n";
		// bubblesort!
/*			for (int i=1,imax=1;i<N(a);i++) {
				if (inc.Axis(a,i-1) > inc.Axis(a,i)) {
					float tmp(inc.Axis(a,i));
					inc.Axis(a,i)=inc.Axis(a,i-1);
					inc.Axis(a,i-1)=tmp;
					for (int x1=0;x1<N(a1);x1++) {
						for (int x2=0;x2<N(a2);x2++) {
							tmp = inc(a1, x1, a2, x2, a, i);
							inc(a1, x1, a2, x2, a, i)  = inc(a1, x1, a2, x2, a, i-1);
							inc(a1, x1, a2, x2, a, i-1)= tmp;
						}
					}
					if (i>imax) { imax=i; printpercent(imax,N(a)-1); }
					i=i-2;
					if (i<0) i=0;
				}
				
			}
			printpercent(1,1);*/
		}
		 
	}
	
	cout << '\n';
	return inc;
}

cube cube::dirinc(axis a)
// JB 2005 January 6 -- dirinc just one axis.
// code pulled out of dirinc(), and that now calls this.		
{
	cube inc(*this);
	
	if (!osuppress)
		cout << "Changing axis " << a << " to increasing.  dirinc("<<a<<") --  00%";
	axis a1=Uax(a);
	axis a2=Dax(a);
	long imax(0);  // keeps track of percent done
	for (long i=1,imax=1;i<N(a);i++) {
		if (inc.Axis(a,i-1) > inc.Axis(a,i)) {
			float tmp(inc.Axis(a,i));
			inc.Axis(a,i)=inc.Axis(a,i-1);
			inc.Axis(a,i-1)=tmp;
			for (long x1=0;x1<N(a1);x1++) {
				for (long x2=0;x2<N(a2);x2++) {
					tmp = inc(a1, x1, a2, x2, a, i);
					inc(a1, x1, a2, x2, a, i)  = inc(a1, x1, a2, x2, a, i-1);
					inc(a1, x1, a2, x2, a, i-1)= tmp;
				}
			}
			if (i>imax) { imax=i; if (!osuppress) printpercent(imax,N(a)-1); }
			i=i-2;
			if (i<0) i=0;
		}
		
	}
	if (!osuppress) printpercent(1,1);
	if (!osuppress) cout << "\n";

	return inc;
}
	
	
// removed in favor of Jlocate, 2012 August 22
/*
int cube::locate(double v, axis A)						
		// returns the pixel number whose wavelength is closest to v 
		// but still below v on axis A
{
	float *arrstart;
	unsigned long n, j;
	if (A==X) { n=nx; arrstart=&(Xaxisarray())[-1]; }			// -1 because NR uses
	if (A==Y) { n=ny; arrstart=&(Yaxisarray())[-1]; }			// arrays that start at 1
	if (A==Z) { n=nz; arrstart=&(Zaxisarray())[-1]; }
	
	::locate(arrstart, n, (float)v, &j);							// Numerical recipes locate
	
	return (int)(j);
}*/

cube cube::dopplershift(float velocity, axis a)
/*  Created 2/23/2k, only works for Z at this point 
	Velocity assumed to be in km/s.
*/
{
	cube shifted(*this);
	
	if (a!=Z) cout << "ERROR -- Dopplershift only works in Z at this point \n";
	
	float c=3e5;
	float dlambda;
	dlambda=sqrtf( (1+velocity/c) / (1-velocity/c) );
	for (int i=0;i<nz;i++) shifted.Zaxis(i) = shifted.Zaxis(i) * dlambda;
	if (Dtype() == 2  ||  (Dtype()>=4 && Dtype()<=6) )
		loopoverall shifted(x,y,z) /= dlambda;
	
	return shifted;
} 

/*cube cube::copyaxis(float xarr[], float yarr[], float zarr[])
{
	for (int x=0;x<nx;x++) xaxis[x]=xarr[x];
	for (int y=0;y<ny;y++) yaxis[y]=yarr[y];
	for (int z=0;z<nz;z++) zaxis[z]=zarr[z];
	return *this;
}*/  // phased out 2012 August 22

void cube::copyaxis(const cube& src, axis A)
/* created 2/23/2k
	copies axis a from cube src to cube *this 
	fixed 3/6/2k
	Gutted and redone for const correctness 2006/03/27
	*/
{	
	for (int i=0;i<N(A);i++) Axis(A,i) = src.Axis(A,i);
	
	return;
}
	
	
void cube::copyaxis(const cube& src)
/* created 2012 August 22	*/
{	
	copyaxis(src,X);
	copyaxis(src,Y);
	copyaxis(src,Z);
	
	return;
}

cube cube::stripneg(float defaultval)
/* created 2/23/2k   Struggling with whether or not to have this return a cube
or just operate on the cube in place, decided this way . . . 
	Changed my mind to return a cube, 2004 December 3.  Also changed this to 
		striplt.
*/
{
	return striplt(0., defaultval);
}

cube cube::striplt(float limit, float defaultval)
// Created 2004 December 3
{
	cube answer(*this);
	
	loopoverall
		if ((*this)(x,y,z) < limit) answer(x,y,z) = defaultval;

	return answer;	
}		

cube cube::striple(float limit, float defaultval)
// Created 2004 December 3
{
	cube answer(*this);
	
	loopoverall
		if ((*this)(x,y,z) <= limit) answer(x,y,z) = defaultval;

	return answer;	
}	

cube cube::stripgt(float limit, float defaultval)
/* Created 2/23/2k 
	Modified 2004 December 3*/
{
	cube answer(*this);
	
	loopoverall
		if ((*this)(x,y,z) > limit) answer(x,y,z) = defaultval;

	return answer;	
}	



cube cube::grow(cube templatecube, double padno)
/* Pads a cube to dimensions of templatecube.  Fills in padded spaces with padno
	Created 1/26/01 JB */
{
	cube grown;
	grown=grow(templatecube.nx, templatecube.ny, templatecube.nz, padno);
	return grown;
}

cube cube::grow(int newx, int newy, int newz, double padno)
/* Pads a cube to the given dimensions.  Fills in padded spaces with padno
	Created 1/26/01 JB */
{
	if (newx==-1) newx==nx;
	if (newy==-1) newy==ny;
	if (newz==-1) newz==nz;
	cube grown(newx, newy, newz, padno);
	grown.hdr = hdr;
	grown.cubetype=cubetype;
	grown.byteorder=byteorder;

	for (int x=0;x<nx;x++)
		for (int y=0;y<ny;y++)
			for (int z=0;z<nz;z++)
				grown(x,y,z)=(*this)(x,y,z);
			
	return grown;
}


cube cube::shift(int dx, int dy, int dz, float repl)
/* JB 10/26/2001
	moves the cube around inside by dx, dy, dz.  The elements that go off the edge are
	lost forever, those that are new are zeroes. */	
{
	cube answer(N(X),N(Y),N(Z),repl);
	
	for (int x=0;x<answer.N(X);x++) {
		for (int y=0;y<answer.N(Y);y++) {
			for (int z=0;z<answer.N(Z);z++) {
				int oldx(x-dx),oldy(y-dy),oldz(z-dz);
				if (oldx>=0 && oldx<N(X)  &&
					 oldy>=0 && oldy<N(Y)  &&
					 oldz>=0 && oldz<N(Z))
					answer(x,y,z) = (*this)(oldx,oldy,oldz);
			}
		}
	}
	return answer;
}

cube cube::spliceout(axis A, long i0, long i1)
//  effectively removes one or a set of planes from a cube
// created 2005 October 12
{
	if (i1==-1) i1 = i0;  // -1 is the default, coding for removing just one plane
	cube answer;
	
	if (debug>0) cout << "Splicing " << i0 << " to " << i1 << " out of cube with " << A;
	if (debug>0) cout << "=" << N(A) << "\n";  cout.flush();
	
	if (i0 == 0)
		answer = (*this)(A, i1+1, N(A)-1, Uax(A), -1, -1, Dax(A), -1, -1);
	else if (i1 == N(A)-1)
		answer = (*this)(A, 0, i0-1, Uax(A), -1, -1, Dax(A), -1, -1);
	else {
		answer = (*this)(A, 0, i0-1, Uax(A), -1, -1, Dax(A), -1, -1);
		answer = answer.blocks(A,(*this)(A, i1+1, N(A)-1, Uax(A), -1, -1,
				Dax(A), -1, -1));
	}
	if (debug>0) cout << "  done splicing.  N("<<A<<") = " << N(A) << "\n";
	
	return answer;
}

cube cube::spliceout(axis A, int i0, int i1)
{
	return spliceout(A, long(i0), long(i1));
}
		
cube cube::spliceout(axis A, double d0, double d1)
// created 2016 March 3 JWB
{	
	if (d1 < d0) {
		double temp(d0);
		d0=d1;
		d1=temp;
	}
	
	long i0(Jlocate(d0, A));
	long i1(Jlocate(d1, A));
	
	return spliceout(A, i0, i1);
}


cube cube::smartselectivespliceout(bool (*evaluator)(cube&, axis, int))
//  Creates a result cube that leaves out those channels where evaluator returns 0
// created 2011 May 4
{
	vector<vector<bool> > keeps(3);
	
	keeps.at(X).resize(N(X));
	for (int x(0);x<keeps.at(X).size();x++)
		keeps.at(X).at(x) = (*evaluator)(*this, X, x);
	
	keeps.at(Y).resize(N(Y));
	for (int y(0);y<keeps.at(Y).size();y++)
		keeps.at(Y).at(y) = (*evaluator)(*this, Y, y);
	
	keeps.at(Z).resize(N(Z));
	if (osuppress < 1) cout << "Smartselectivespliceout in Z:   00%";
	for (int z(0);z<keeps.at(Z).size();z++) {
		if (osuppress < 1) printpercent(z, keeps.at(Z).size()-1);
		keeps.at(Z).at(z) = (*evaluator)(*this, Z, z);
	}
	if (osuppress < 1) printpercent(keeps.at(Z).size()-1, keeps.at(Z).size()-1);
	if (osuppress < 1) cout << "\n";
		
		
	cube answer(N(X)-int(count(keeps.at(X).begin(), keeps.at(X).end(), 0)), 
					N(Y)-int(count(keeps.at(Y).begin(), keeps.at(Y).end(), 0)),
					N(Z)-int(count(keeps.at(Z).begin(), keeps.at(Z).end(), 0)));
	answer.copymetadatafrom(*this);
	
	if (osuppress < 1) cout << "count in X:" << answer.N(X) << "  Y:" << answer.N(Y) << "  Z:" << answer.N(Z) << "\n";
	
	for (axis A(X); ; A=Uax(A)) {
		for (int oldi(0), newi(0);oldi<N(A);oldi++) {
//			cout << "A=" << A << ", oldi=" << oldi << ", newi=" << newi << ", keeps=";
//			cout << keeps.at(A).at(oldi) << "\n";
			if (keeps.at(A).at(oldi) == true) {
				answer.Axis(A, newi) = Axis(A, oldi);
				newi++;
			} 
		}
		if (A==Z) break;
	}
	
	for (int oldx(0), newx(0);oldx<N(X);oldx++) {
		if (keeps.at(X).at(oldx)) { 
			for (int oldy(0), newy(0);oldy<N(Y);oldy++) {
				if (keeps.at(Y).at(oldy)) {
					for (int oldz(0), newz(0);oldz<N(Z);oldz++) {
						if (keeps.at(Z).at(oldz)) {
							answer(newx,newy,newz) = (*this)(oldx,oldy,oldz);
							newz++;
						}
					}
					newy++;
				}
			}
			newx++;
		}
	}
	
	return answer;
}

void cube::headerwipe() { hdr.clear(); }

vector<double> Jcube2nr(cube incube)
{
	int xx(1), yy(1), zz(1);
	
	while (xx < incube.N(X)) xx *= 2;
	while (yy < incube.N(Y)) yy *= 2;  // pad all of the dimensions to be multiples of 2 in size
	while (zz < incube.N(Z)) zz *= 2;
	
	vector<double> answer(xx*yy*zz);
	
	int xoffset((xx-incube.N(X))/2);
	int yoffset((yy-incube.N(Y))/2);
	int zoffset((zz-incube.N(Z))/2);
	
	for (int x(0);x<incube.N(X);x++)
		for (int y(0);y<incube.N(Y);y++)
			for (int z(0);z<incube.N(Z);z++)
				if (x>=xoffset && x<xoffset+incube.N(X) && y>=yoffset && y<yoffset+incube.N(Y) && z>=zoffset && z<zoffset+incube.N(Z))
					answer.at(x*yy*zz+y*zz+z) = incube(x-xoffset,y-yoffset,z-zoffset);
				else
					answer.at(x*yy*zz+y*zz+z) = 0.;
	
	return answer;
}


cube cube::convolved_with(cube r)
{
	cube answer(*this);
		
	if (1. || (N()*r.N() < 1e8)) {  // explicit convolution -- slow for large cubes
	
		int xradius((r.N(X)-1)/2);
		int yradius((r.N(Y)-1)/2);
		int zradius((r.N(Z)-1)/2);
		answer *= 0.;
	

		if (osuppress<1) cout << "Convolving ("<<N(X)<<","<<N(Y)<<","<<N(Z)<<")**("<<r.N(X)<<","<<r.N(Y)<<","<<N(Z)<<")";
		if (osuppress<1) cout << " explicitly:   00%";
		for (int x(0);x<N(X);x++)  {
			if (osuppress<1) printpercent(x,N(X)-1);
			for (int y(0);y<N(Y);y++) 
				for (int z(0);z<N(Z);z++) 
					for (int cx(-xradius);cx<=xradius;cx++) 
						for (int cy(-yradius);cy<=yradius;cy++) 
							for (int cz(-zradius);cz<=zradius;cz++) 
								if (x+cx>=0 && x+cx<N(X) && y+cy>=0 && y+cy<N(Y) && z+cz>=0 && z+cz<N(Z))
									if (cx+xradius>=0 && cy+yradius>=0 && cz+zradius>=0 && cx+xradius<r.N(X) &&
											cy+yradius<r.N(Y) && cz+zradius<r.N(Z))  // shouldn't need this, but segfaulting
										answer(x+cx,y+cy,z+cz) += (*this)(x,y,z)*r(cx+xradius,cy+yradius,cz+zradius);
		}
		if (osuppress<1) printpercent(100, 100);
		if (osuppress<1) cout << "\n";	
	} else {  // convolution in fourier space
		// move both arrays into appropriately-sized and -ordered vector<double>s for NR routines
		vector<double> left(Jcube2nr(*this));
		vector<double> right(Jcube2nr(r));   // left unfinished.
	}
				
	
	cout << "returning from convolution\n";
	return answer;
}




cube cube::FFTconvolve(cube response, int deconvolve)
{
	/* 2016.11 smack
		runs, but isn't giving sensical answer. It seems like the padding funcitonality of realFFT()
		doesn't work as I'd expect it to.
	*/
	cube answer(*this);	
	
	// need to pad the incoming cube first
	
	
	
	// first 3 steps of nr's convlv p.645 (put response in array, pad with zeroes, and take FFT) are taken care of by cube::realFFT(int)
	cube s(realFFT(1)),r(response.realFFT(1));
	cube ans(s);
	// s and r are cubes of N(X), N(Y), 2N(Z) their input cubes.
	
	//nr recipes ans= signal
	//nr recipes temp=r
	
	//"no2 = n>>1"?! 
	int no2=N(X)*N(Y)/2; 
		
	if (deconvolve != 1){
		//multiply FFTs to convolve
		for (int z(2); z<ans.N(Z); z+=2){
			for (int x(1); x<ans.N(X);x+=2){
				for (int y(1); y<ans.N(Y);y+=2){
					double thisans(s(x,y,z));
					ans(x,y,z)= (s(x,y,z)*r(x,y,z)-s(x,y,z+1)*r(x,y,z+1))/no2;
					ans(x,y,z+1)=(s(x,y,z+1)*r(x,y,z)+thisans*r(x,y,z+1))/no2;
				}	
			}	
		} 
		
		ans(0,0,0)=s(0,0,0)*r(0,0,0)/no2;
		ans(0,0,1)=s(0,0,1)*r(0,0,1)/no2;		
	
	} else {
		// divide FFTs to deconvolve
		for (int z(2); z<ans.N(Z); z+=2){
			for (int x(1); x<ans.N(X);x+=2){
				for (int y(1); y<ans.N(Y);y+=2){
					double thisans(s(x,y,z));
					double mag(sqrt(r(x,y,z))+sqrt(r(x,y,z+1)));
					 
					ans(x,y,z)= (s(x,y,z)*r(x,y,z)+s(x,y,z+1)*r(x,y,z+1))/mag/no2;
					ans(x,y,z+1)=(s(x,y,z+1)*r(x,y,z)-thisans*r(x,y,z+1))/mag/no2;
				}	
			}	
		}
		
		ans(0,0,0)=s(0,0,0)/r(0,0,0)/no2;
		ans(0,0,1)=s(0,0,1)/r(0,0,1)/no2;		
	 
	} 
	
	answer=ans.realFFT(-1); //inverse transform back to original cube domain
	return answer;	
}

cube cube::rotated(double d, axis A, double lc, double hc) const
// created 2012 November 28 JWB
// rotates the cube by an angle d in degrees clockwise
// around axis A (Z by default)
// around a center at (lc,hc).  This center bit is in x,y if A==Z, yz if A==X,
// and xz if A==Y.
{
	double rad(d/180.*Jcrap::pi);  // clockwise rotation in radians
	
	axis B(UNK),C(UNK);
	if (A==Z) {B=X; C=Y;}
	if (A==X) {B=Y; C=Z;}
	if (A==Y) {B=X; C=Z;}
	
	if (lc==-1) lc=double((N(B)-1.)/2.);
	if (hc==-1) hc=double((N(C)-1.)/2.);
	
// If the B&C axis values are interesting, the rotation is meaningless.
// so make 'original' with just normal pixel numbers for Axis values
	cube original(*this);
	for (long b(0);b<N(B);b++) original.Axis(B,b)=double(b);
	for (long c(0);c<N(C);c++) original.Axis(C,c)=double(c);
	
	cube answer(original*0.);
	if (osuppress<=0) cout << "Rotating by " << d << " --  00%";
	for (long a(0);a<N(A);a++) {
		for (long b(0);b<N(B);b++) {
			if (osuppress<=0) printpercent(N(B)*N(C)*a+b*N(C), N());
			for (long c(0);c<N(C);c++) {
				double bwrtc(double(b)-lc);
				double cwrtc(double(c)-hc);
				double r(sqrt(bwrtc*bwrtc+cwrtc*cwrtc));
				double theta(atan2(-cwrtc/r,bwrtc/r));
	
				// here's the rotation
				theta += rad;
				
				double newb(r*cos(theta)+lc);
				double newc(r*sin(theta)+hc);
				
				if (newb < -0.5) continue;
				if (newb > double(N(B))-0.5) continue;
				if (newc < -0.5) continue;
				if (newc > double(N(C))-0.5) continue;
				
				// not pretty -- doesn't interpolate or anything yet
				answer(A,a,B,b,C,c) = original(A,a,B,long(newb+0.5),C,long(newc+0.5));
			}
		}
	}
	if (osuppress<=0) {printpercent(N(), N());  cout << "\n";}
	
	return answer;
}
