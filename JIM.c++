#include "Jcube.h"
#include "string.h"
#include "Magick++.h"
#include <vector>

using namespace Magick;


void cube::IMread(const char *fname)
{
	cout << "IMread\n";
	Image in;
	in.read(fname);
	
	nz=0;
	
	*this = cube(in.baseColumns(), in.baseRows(), 3);
	
	Color c;
	float *pfl;
	pfl = &Data[0];
	for (int x=0;x<nx;x++){
		printpercent(x,nx);
		for (int y=0;y<ny;y++){
			c = in.pixelColor(x,y);
			ColorRGB *p;
			p=(ColorRGB*)&c;
			
			*pfl=(float)p->red();
			pfl++;
			*pfl=(float)p->green();
			pfl++;
			*pfl=(float)p->blue();
			pfl++;
		}
	}
	printpercent (nx,nx);
	cout << "\n";
}

void cube::IMwrite(const char *infname, float minval, float maxval, float (*s)(float))
{
	cout << "In IMwrite\n";
	char fname[1000];
	strcpy (fname, infname);
	if (cubetype==jpg && strcmp(&fname[strlen(fname)-3], "jpg")) strcat(fname, ".jpg");
	if (cubetype==gif && strcmp(&fname[strlen(fname)-3], "gif")) strcat(fname, ".gif");
	if (cubetype==tif && strcmp(&fname[strlen(fname)-3], "tif")) strcat(fname, ".tif");

	
	
	int rz,gz,bz;
	if (N(Z) == 3) { rz=0; gz=1; bz=2; }
	else {rz=gz=bz=0;}
	
	Geometry dims(nx,ny);
	ColorRGB white(1.,1.,1.);
	
	vector<Magick::Image> frames;
	int numberofframes(1);
	if (cubetype==gif && N(Z)>3) numberofframes=N(Z);
	else numberofframes=1;
	
	for (int i(0);i<numberofframes;i++) {
		cout << "creating Image\n";
		Image out(dims, white);
		if (cubetype==tif32) out.depth(32);
		else out.depth(8);
	
		if (minval==-1) minval=Jselect((long unsigned int)(nx*ny*nz*.1));
		if (maxval==-1) maxval=Jselect((long unsigned int)(nx*ny*nz*.9));
		cout << "Saving scaled between " << minval << " and " << maxval << "\n";
		cout << "Converting for output --  00%";	
		for (int x=0;x<nx;x++){
			printpercent (x, nx-1);
			for (int y=0;y<ny;y++){
				float r, g, b;
				
				if (cubetype == tif32) {
					r = (*this)(x,y,rz+i);
					g = (*this)(x,y,gz+i);
					b = (*this)(x,y,bz+i);
				} else {
					if (s) r=(s((*this)(x,y,rz+i))-s(minval))/(s(maxval)-s(minval));
					else r=(((*this)(x,y,rz+i))-(minval))/((maxval)-(minval));
					if (r<0.) r=0.;
					if (r>1.) r=1.;
					
					if (s) g=(s((*this)(x,y,gz+i))-s(minval))/(s(maxval)-s(minval));
					else g=(((*this)(x,y,gz+i))-(minval))/((maxval)-(minval));
					if (g<0.) g=0.;
					if (g>1.) g=1.;
					
					if (s) b=(s((*this)(x,y,bz+i))-s(minval))/(s(maxval)-s(minval));
					else b=(((*this)(x,y,bz+i))-(minval))/((maxval)-(minval));
					if (b<0.) b=0.;
					if (b>1.) b=1.;
				}
				
				ColorRGB p (r,g,b);
				out.pixelColor(x,y, p);
			}
		}
		
		if (this->keyword("blacktransparent").size())
			out.transparent(ColorRGB(0.,0.,0.));
		
		if (cubetype!=gif && N(Z)<=3) {
			if (cubetype==jpg) out.quality(100);
			cout << "about to write " << fname << "\n";
			out.write(string(fname));
		} else {
			out.animationDelay(30.);
			frames.push_back(out);
		}
	}
	
	if (frames.size() && cubetype==gif && N(Z)>3)
		Magick::writeImages(frames.begin(), frames.end(), fname);
}

float nullfunc(float a) {return a;}

void cube::display(otype dest, float minval, float maxval, 
		float (*s)(float), float delay)
/* Created 3/3/2k.  Converts to imagemagick format to display on the screen, scaling
between minval and maxval, or if default, the 10% and 90% levels in the image */
{
	if (minval==-1) minval=Jselect((long unsigned int)(nx*ny*nz*.1));
	if (maxval==-1) maxval=Jselect((long unsigned int)(nx*ny*nz*.9));
	cout << "Displaying scaled between " << minval << " and " << maxval << "\n";
	
	Geometry dims(nx,ny);
	ColorRGB white(1.,1.,1.);
	Color clear(1., 1., 1., 1.);
	
	Image out(dims, clear);
	
	if (dest==oCOLOR){
		ColorRGB p;
		float i;
		cout << "creating color IMage (" << nx << ", " << ny << ")\n";
		for (int x=0;x<nx;x++){
			for (int y=0;y<ny;y++){
				if ((*this)(x,y,0) != -1.) {  // if it is -1, leave transparent
					if (s) i=(s((*this)(x,y,0))-s(minval))/(s(maxval)-s(minval));
					else i=(((*this)(x,y,0))-(minval))/((maxval)-(minval));
					if (i<=0.0001) i=0.0001;
					if (i>1.) i=1.;
					p.red(i);
				}
				
				if ((*this)(x,y,1) != -1.) {  // if it is -1, leave transparent
					if (s) i=(s((*this)(x,y,1))-s(minval))/(s(maxval)-s(minval));
					else i=(((*this)(x,y,1))-(minval))/((maxval)-(minval));
					if (i<=0.0001) i=0.0001;
					if (i>1.) i=1;
					p.green(i);	
				}
						
				if ((*this)(x,y,2) != -1.) {  // if it is -1, leave transparent	
					if (s) i=(s((*this)(x,y,2))-s(minval))/(s(maxval)-s(minval));
					else i=(((*this)(x,y,2))-(minval))/((maxval)-(minval));
					if (i<=0.0001) i=0.0001;
					if (i>1.) i=1.;
					p.blue(i);
				}
				
				Color c(p);
				out.pixelColor(x,y,c);
			}
		}
		cout << "displaying\n";
		out.display();
	}
	if (dest==oBW){
		ColorRGB p;
		float i;
		cout << "creating IMage\n";
		for (int x=0;x<nx;x++){
			for (int y=0;y<ny;y++){
				if (s) i = (s((*this)(x,y,0))-s(minval))/(s(maxval)-s(minval));
				else i=((*this)(x,y,0)-minval)/(maxval-minval);
				if (i<=0.0001) i=0.0001;
				if (i>1.) i=1.;
				
				p.red(i);
				p.green(i);
				p.blue(i);
				Color c(p);
				out.pixelColor(x,y,p);
			}
		}
		cout << "displaying\n";
		out.display();
	}
	
	if (dest==oMOVIE || dest==oCOLORMOVIE){
		vector<Image> movie;
		out.animationDelay((unsigned int)(100*delay));
		out.animationIterations(0);
		ColorRGB p;
		float i;
		
		int addz = 1;
		if (dest == oCOLORMOVIE) addz=3;
		
		cout << "preparing movie --   00%";
		for (int z=0;z<nz;z+=addz){
			movie.push_back(out);
			for (int x=0;x<nx;x++){
				printpercent(z*nx+x, nx*nz);
				for (int y=0;y<ny;y++){
					if (s) i=(s((*this)(x,y,z))-s(minval))/(s(maxval)-s(minval));
					else i=(((*this)(x,y,z))-(minval))/((maxval)-(minval));
					if (i<=0.0001) i=0.0001;
					if (i>=0.9999) i=0.9999;
					if (dest == oMOVIE) {
						p.red(i);
						p.green(i);
						p.blue(i);
					}
					if (dest == oCOLORMOVIE) {
						p.red(i);					
						if (s) i=(s((*this)(x,y,z+1))-s(minval))/(s(maxval)-s(minval));
						else i=(((*this)(x,y,z+1))-(minval))/((maxval)-(minval));	
						if (i<=0.0001) i=0.0001;
						if (i>=0.9999) i=0.9999;
						p.green(i);					
						if (s) i=(s((*this)(x,y,z+2))-s(minval))/(s(maxval)-s(minval));
						else i=(((*this)(x,y,z+2))-(minval))/((maxval)-(minval));			
						if (i<=0.0001) i=0.0001;
						if (i>=0.9999) i=0.9999;
						p.blue(i);
					}
					Color c(p);
					movie[z/addz].pixelColor(x,y,c);
				}
			}
		}
		cout << "\n";
		animateImages(movie.begin(), movie.end());
		writeImages(movie.begin(), movie.end(), "movie.mpg");
	}
}
