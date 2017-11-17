#include "Jcrap.h"

float mad(cube& reference, cube& b, int dx, int dy, int maxx, int maxy)
{
// compute difference of reference and this xyplane for this offset
	cube diff;
	diff = reference(maxx,b.N(X)-maxx-1,maxy,b.N(Y)-maxy-1,0,0);
	{	cube c(b(maxx+dx,b.N(X)-maxx-1+dx,maxy+dy,b.N(Y)-maxy-1+dy,0,0));
		diff -= c;
	}
	diff=diff.abs();
	
	float ad(diff.sum());	
	cout << dx << ","<<dy<<","<<0<<">  " << ad << "\n";
	return ad;
}

//cube cube::madroughgradient(int zref=0,int rought

cube cube::madroughdumb(int zref,int roughness,int maxx,int maxy)
// 3/4/2002 JB -- does quick and dirty rough MAD by coadding roughness x roughness pixels.
{
	if (maxx == -1) maxx = 10;
	if (maxy == -1) maxy = 10;
	
// resample and MAD the small cube in the dumb way.
	cube roughed(resample(roughness,roughness,1));
	vector<pair<int,int> > xyoffsets;
	xyoffsets = roughed.madraster(zref,maxx,maxy);
	
// fix the offsets, shift, and output
	cube answer(N(X),N(Y),0);  answer.hdr=hdr;  answer.cubetype=cubetype;  answer.byteorder=byteorder;
	for (vector<pair<int,int> >::iterator i=xyoffsets.begin();i!=xyoffsets.end();i++) {
		i->first *= roughness;
		i->second*= roughness;
	}
	for (int z=0;z<N(Z);z++)
		answer=answer.blocksz(this->plane(Z,z).shift(-xyoffsets[z].first,-xyoffsets[z].second));
	cout << "Exiting madroughdumb\n";
	
	return answer;
}

vector<pair<int,int> > cube::madraster(int zref,int maxx,int maxy)
// imported from maddumb 3/4/2002
/* Created 10/26/2001 JB 
	performs minimum absolute difference on the cube (assumes xy is image plane)
	zref = the plane that all are to be shifted to
	maxx = max x shift (default = 3% nx)
	maxy = max y shift (default = 3% ny) 
	
	This version returns a vector of x,y offsets between the reference
	and each individual plane(Z).
	*/
{
	if (maxx == -1) maxx = (int)(.03*N(X));
	if (maxy == -1) maxy = (int)(.03*N(Y));
	
	vector<pair<int,int> > xyoffsets(N(Z));
	cube reference=(*this).plane(Z,zref,zref);
	cout << reference.info();
	cout << plane(Z,1).info();
	for (int z=0;z<N(Z);z++){
		int deltax(0), deltay(0);
		double minad(1.e99);		// hopefully impossibly high values
		if (z!=zref) {
			cube roamer(plane(Z,z));
			for (int dx=-maxx;dx<=maxx;dx++) {
				for (int dy=-maxy;dy<=maxy;dy++) {
					double ad;
					ad = mad(reference, roamer, dx, dy, maxx, maxy);
					if (ad < minad) { minad=ad; deltax=dx; deltay=dy; }
				}
			}
		}
		
	// put shifts into answer cube
	xyoffsets[z].first = deltax;
	xyoffsets[z].second= deltay;
	cout << "For z="<<z<<", xoff="<<deltax<<", yoff="<<deltay<<"\n";
	}
	
	cout << "MADraster complete!\n";
	return xyoffsets;
}


cube cube::maddumb(int zref,int maxx,int maxy)
// code pushed out into madraster 3/4/2002
{
	vector<pair<int,int> > xyoffsets;
	xyoffsets = madraster(zref,maxx,maxy);
	
// perform the shifts and put them into the answer cube
	cube answer(N(X),N(Y),0);
	answer.hdr = hdr;
	answer.cubetype = cubetype;
	answer.byteorder= byteorder;
	
	for (int z=0;z<N(Z);z++)
		answer.blocksz(plane(Z,z).shift(xyoffsets[z].first,xyoffsets[z].second));
	return answer;
}

cube cube::madgradient(cube searchspace, int zref,int maxx,int maxy){
	if (maxx == -1) maxx = (int)(.03*N(X));
	if (maxy == -1) maxy = (int)(.03*N(Y));
	
	cube answer(N(X),N(Y),0);
	answer.hdr = hdr;
	answer.cubetype = cubetype;
	answer.byteorder = byteorder;
	

	
	cube reference=(*this).plane(Z,zref,zref);	
	
	for (int z=0;z<N(Z);z++){
		int deltax(0), deltay(0);
		double minad(1.e99);		// hopefully an impossibly high value
		cube roamer(plane(Z,z));
		int dx=0, dy=0;
		if (z!=zref) {
			int xarea=(searchspace.N(X)-1)/2, yarea=(searchspace.N(Y)-1)/2;
			cube matrix(searchspace);
			cout << "starting matrix = " << matrix;
			while (dx>-maxx && dx < maxx  &&  dy>-maxy && dy<maxy) {
				// set t
				for (cube::iterator i=searchspace.begin(),j=matrix.begin();i!=searchspace.end();i++,j++)
					if (*i == -2) *j = 1.e99;
				
				// for the values around the current dx,dy, get the MADs
				for (int tryx=-xarea;tryx<=xarea;tryx++) {
					for (int tryy=-yarea;tryy<=yarea;tryy++) {
						if (matrix(tryx+xarea,tryy+yarea,0) == -1.   &&
							 dx+tryx <= maxx   &&   dx+tryx >= -maxx &&
							 dy+tryy <= maxy   &&   dy+tryy >= -maxy) 
							matrix(tryx+xarea,tryy+yarea,0) = mad(reference,roamer,dx+tryx,dy+tryy,maxx,maxy);
					}
				}
				cout << "MAD matrix: \n" << matrix;
				
				// find out which around here is the lowest.
				int xdir=0,ydir=0;
				float minsofar(1e99);//matrix(xarea+1,yarea+1,0));
				for (int x=-xarea;x<=xarea;x++) {
					for (int y=-yarea;y<=yarea;y++) {
						if (matrix(x+xarea,y+yarea,0) < minsofar) {
							minsofar = matrix(x+xarea,y+yarea,0);
							xdir=x; ydir=y;
						}
					}
				}
				
				// go ahead and move by xdir and ydir (unless we've found the local min)
				if (xdir==0 && ydir==0) break; // (else)
				dx+=xdir; dy+=ydir;
				matrix=matrix.shift(-xdir,-ydir,0,-1.);
				cout << "moving by " << xdir << "," << ydir << "\n";
			}
		}

	// put shifted xyplane into answer cube
		answer = answer.blocksz(roamer.shift(-dx,-dy));
	}
	
	return answer;
}

cube cube::madgradient(int zref,int maxx,int maxy, int radius)
/* Created 10/26/2001 JB 
	performs minimum absolute difference on the cube (assumes xy is image plane)
	zref = the plane that all are to be shifted to
	maxx = max x shift (default = 3% nx)
	maxy = max y shift (default = 3% ny)
	This version uses the gradient method.  Will we run into local
	minimum/maximum problems?  Who knows!
	 */
// this method became a special case of madgradient(cube searchspace) 3/4/2002
{
	cube searchspace(radius*2+1,radius*2+1,1,-1.);
	return madgradient(searchspace,zref,maxx,maxy);
}


cube cube::clean(cube badmap)
/* 2/25/2002 JB
	This method returns a cube cleaned based on the cube badmap.
	Each pixel of the badmap should be '0' if the pixel is good,
	and, if converted to an integer byte, the '1' set if the bad
	pixel should be interpolated in X, the '2' for Y, and the '4'
	for Z.*/
{
	cube cleaned(*this);
	
	for (int x=0;x<N(X);x++) {  
		for (int y=0;y<N(Y);y++) {
			for (int z=0;z<N(Z);z++) {
				int badx, bady, badz;
				if (badmap.N(X) == 1) badx = 0; else badx = x;
				if (badmap.N(Y) == 1) bady = 0; else bady = y;
				if (badmap.N(Z) == 1) badz = 0; else badz = z;
				
				unsigned char diminterp;
				diminterp = (unsigned char) badmap(badx,bady,badz);
				if (diminterp & 0x01)
					cleaned(x,y,z) = (cleaned(x-1,y,z)+cleaned(x+1,y,z))/2.;
				if (diminterp & 0x02)
					cleaned(x,y,z) = (cleaned(x,y-1,z)+cleaned(x,y+1,z))/2.;
				if (diminterp & 0x04)
					cleaned(x,y,z) = (cleaned(x,y,z-1)+cleaned(x,y,z+1))/2.;
			}
		}
	}
	
	return cleaned;
}

cube cube::boxcarfilter(int xadd, int yadd, int zadd)
/*  Boxcar filters such that if xadd = 1 then only one pixel is added, if 2 then
    3 pixels, 3 then 5 pixels, etc.
	 
	 Added 2004-12-01 JB */
{
	cube answer((*this)*0.);
	if (!osuppress) cout << "Boxcar filtering " << xadd << "x" << yadd << "x" << zadd;
	if (!osuppress) cout << " --  00%";
	for (int x=0;x<N(X);x++) {
		if (!osuppress) printpercent(x, N(X)-1);
		for (int y=0;y<N(Y);y++) {
			for (int z=0;z<N(Z);z++) {
				int count=0;
				for (int bx=-(xadd-1);bx<=xadd-1;bx++) {
					for (int by=-(yadd-1);by<=yadd-1;by++) {
						for (int bz=-(zadd-1);bz<=zadd-1;bz++) {
							if (x+bx >= 0  &&  x+bx < N(X)   &&
								 y+by >= 0  &&  y+by < N(Y)   &&
								 z+bz >= 0  &&  z+bz < N(Z)) {
								count++;
								answer(x, y, z) += (*this)(x+bx, y+by, z+bz);
							}
						}
					}
				}
				answer(x, y, z) /= (double)count;
			}
		}
	}
	if (!osuppress) printpercent(nz, nz);
	if (!osuppress) cout << "\n";
	return answer;
}


cube cube::stripNaN(float newval)
{
	float v;
	cube answer(*this);
	loopoverall {
		v = (*this)(x,y,z);
		if (!(v<0.) && !(v==0.) && !(v>0.))
			answer(x,y,z) = newval;
	}
	return answer;
		
}

cube cube::stripinf(float newval)
{
	float v;
	cube answer(*this);
	loopoverall {
		v = (*this)(x,y,z);
		if ((v<-1.e30) || v>1.e30)
			answer(x,y,z) = newval;
	}
	return answer;
		
}
