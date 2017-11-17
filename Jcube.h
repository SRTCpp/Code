#include<string>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<map>
#include<vector>
#include<math.h>
//#include <isis3/Cube.h>
#include "../../nr/nr.h"

#include "Jangle.h"

#define EGCS 1

#define loopoverall for(int x=0;x<nx;x++) for(int y=0;y<ny;y++) for(int z=0;z<nz;z++)
#define LOOPOVERALLPCT for(int z=0;z<nz;z++,printpercent(z,nz-1)) for(int y=0;y<ny;y++) for(int x=0;x<nx;x++) 

using namespace std;

enum filetype {vicar,isis,isis3,text,Jcube1,jpg,gif,tif,fits,Jcube1a,
		Jcube1b,bsqraw,graphout,tif32};
ostream& operator<<(ostream& out, ::filetype& f);

enum datatype {NONE,DN,FLUX,PHOTONS,FLUXCGS,FLUXCGSWAVENUM,FLUXCGSHZ,ARB,KM,
               PHOTONFLUX,WM2};  
				// FLUX assumed to be Watts per m^2 per micron
				//	FLUXCGS assumed to be erg/cm^2/s/cm
				// FLUXCGSWAVENUM erg/cm^2/s/cm^-1
				// FLUXCGSHZ erg/cm^2/s/hz
				// PHOTONFLUX is photons/sec/micron
				// WM2 is watts/m^2 straight up
enum axistype {ARRAYVAL,MICRONS,WAVENUM,CM,M,KMPERS,DEGREES,TIME,HZ};
enum noisetype{UNIFORM, GAUSSIAN, POISSON};
enum methodtype{mFAST, mSLOW, mLINEAR, mQUADRATIC, mBICUBIC, mLANCZOS, 
	mNN=mFAST};
enum interpmethod{INTERPOLATE, FLUXCONSERVE, INTEGRATE};
enum endian{bigendian,littleendian};

enum axis {X=0, Y=1, Z=2, UNK=3, T=2, XY=2, XZ=1, YZ=0, XT=1, YT=0};


axis Uax(axis a);			//  (Jcube) returns next hightest axis
axis Dax(axis a);  		//  (Jcube) returns next lowest axis
axis otheraxis(axis a, axis b); // (Jcube) returns the one axis that is not sent.

enum otype {oX, oPS, oGIF, oBW, oCOLOR, oMOVIE, oCOLORMOVIE, oFIG};
enum projtype {CYL,ORTHO,AZSTEREOGR,LAMBAZ,AZIMUTHAL}; //projection types

class projection; //2.16 smack
class azimuthalprojection; //2.16 smack


extern int osuppress;
extern long idum;

float nullfunc(float a);
inline long exp(int, int);
inline float exp(float, int);
inline float polyval(float x, vector<float> a)
/* 11/13/2k JB  What's with the other ver of this, I don't think it should work*/
{
	int degree=a.size()-2;
	float ans=a[1];
	for (int j=degree+1;j>1;j--) ans+=pow(x,j-1)*a[j];
	return ans;
}
void printpercent(int n, int d);				// Jmisc

int str2int(string s);
unsigned long str2ulong(string s);
float str2float(string s);
double str2double(string s);
string float2str(float f, int n=10);						// Jmisc
string float2fstr(float f, int n=6);						// Jmisc
string double2str(double, int n=19);						// Jmisc
string double2fstr(double, int n=15);						// Jmisc
string int2str(int i);
string int2str(int i, int n);
string int2str(long long);						// Jmisc
string int2str(long);							// Jmisc
string lblstr(int lblsize);					// Jmisc
float binstr2float(string s);					// Jmisc
string float2binstr(float f);					// Jmisc
string asstring(otype);							// Jmisc
string projtype2str(projtype);				// Jprojections

namespace Jcrap {
	class pixelnode;
	void getfloat(ifstream*, void *);  // Jio
	inline void getbigendint(ifstream* &, void *); // Jio
}

class cube {
	public:
		/* creators and annihilators*/
		cube(string, int=-1);												// creates from filename
		cube(const cube&, int=-1);												// creates from another cube
		cube(long x=1, long y=1, long z=1, float initval=0.0, filetype type=Jcube1b, 
				axis=X, axis=Y, axis=Z, int=-1);
		cube(axis, int, axis=UNK, int=1, axis=UNK, int=1, float=0., filetype=Jcube1b,
				axis=X, axis=Y, axis=Z, int=-1);
		explicit cube(float);
		~cube();
		void deallocate();
		
		/* cube generators */
		static cube binary2cube(const char [], int, int, int, axis=X, axis=Y, axis=Z); // (Jio) turns a binary file into a cube.
		cube cpoly(vector<float> a);					// (Jfunc) generates a polynomial based on *this and a[]
		cube generateaxis(axis, float, float);		// (Jfunc) returns the given cube with (axis) given values between min and max
		cube histogram(int=10,float=-1.,float=-1.,cube=cube(1,1,1));// (Jfunc) generates histogram
		cube cumulativehistogram(int=10,float=-1.,float=-1.);// (Jfunc)
		
		/* cube information functions */
		long Nx() const;
		long Ny() const;
		long Nz() const;
		long N(axis a) const;
		long long N() const;										// (Jcube) total n = nx*ny*nz
		long Side() const;
		long Bottom() const;
		long Back() const;
		datatype Dtype() const;
		axistype Xtype() const;
		axistype Ytype() const;
		axistype Ztype() const;
		string& keyword(string);
		string keyword(string) const;
//		int locate(double, axis);  // removed in favor of Jlocate, 2012 August 22
		double res(axis a);								// (Jcube)  returns median resolution in axis a
		float min() const;										// (Jfunc)  returns minimum value in cube
		float max() const;										// (Jfunc)  returns max value in cube
		pair<float, vector<int> > maxcoord() const;     // (Jfunc) returns max and the x,y,z coords where that max is
		pair<float, vector<int> > mincoord() const;     // (Jfunc) returns min and the x,y,z coords where that max is
		double maxaxis(axis) const;								// (Jfunc) returns max of a given axis value
		double minaxis(axis) const;								// (Jfunc) min of an axis 
		double meanaxis(axis) const;							// (Jfunc) avg axis value
		float sum() const;								// (Jfunc)  sum of all values in cube
		float median(float=0.5);									
		cube median(float, axis, axis=UNK, axis=UNK);// (Jfunc)  median
		cube median(axis, axis=UNK, axis=UNK); // (Jfunc) for convenience; assumes 0.5 median.
		cube xmedian(float=0.5);
		cube ymedian(float=0.5);							//returns xt of median
		cube zmedian(float=0.5);
		cube ysigma();											//takes std dev in y, returns	an xtplane with the sigmas
		cube ytrysigma();
		cube sigma(axis=UNK, axis=UNK, axis=UNK);		// (Jfunc) returns standard deviation in axis; default is XYZ
		float mean() const;										// (Jfunc)  avg
		float mode();										// (Jfunc)  mode
		cube centroid() const;							// (Jfunc) returns 1x1x1 cube with Axis values equal to centroid
		string info() const;										// (Jfunc)  returns some interesting info
		pair<double,double> Javevar();		// (Jnr)    Jason's version of NR variance and avg generator
		int dimensionsstored() { return memdims; }		
		
		// access
		vector<float>::iterator Databegin();

//		typedef vector<float>::iterator iterator;
		class iterator : public std::iterator_traits<vector<float>::iterator >
// ho boy.  This is the hard part.  2006 Feb
		{
			public:
				iterator();
				iterator(vector<float>::iterator);
				iterator(long long, cube*);
				friend bool operator==(const cube::iterator&, const cube::iterator&);
				friend bool operator!=(const cube::iterator&, const cube::iterator&);
				
				float* operator=(iterator);
				float* operator=(void*);
				float* operator=(vector<float>::iterator);	
				float& operator*();
				float* operator->();
				cube::iterator operator++();
				cube::iterator operator++(int);
				
			private:
				vector<float>::iterator curr;
				long long i;
				cube *thiscube;
		};
		
		iterator begin();
		iterator end();
		

		map<string,string,less<string> >& header() const {return hdr;}
				
		/* operations that alter the cube given */
		cube& operator=(const cube&);		// Joperators.c++
		bool operator==(cube);
		bool operator==(float);
		cube operator+=(cube);
		cube operator-=(cube&);
		cube operator*=(cube&);
		cube operator/=(cube&);
		cube operator+=(float);
		cube operator-=(float);
		cube operator*=(float);
		cube operator/=(float);
		float& operator()(long, int, int);  // added for long transition, 2013 Apr
		float& operator()(int, long, int);
		float& operator()(int, int, long);
		float& operator()(long, long, int);
		float& operator()(long, int, long);
		float& operator()(int, long, long);
		float& operator()(int, int, int);
		float& operator()(long, int, int) const;
		float& operator()(int, long, int) const;
		float& operator()(int, int, long) const;
		float& operator()(long, long, int) const;
		float& operator()(long, int, long) const;
		float& operator()(int, long, long) const;
		float& operator()(int, int, int) const;		
		float& operator()(long x, long y, long z);
		float& storageorderaccess(long, long, long);
		float& operator()(long x, long y, long z) const;
		float& storageorderaccess(long, long, long) const;
		float& operator()(axis, long, axis=UNK, long=0, axis=UNK, long=0); // (Joperators)
		float& operator()(axis, long, axis=UNK, long=0, axis=UNK, long=0) const; // (Joperators)
		float& operator()(axis, int, axis=UNK, int=0, axis=UNK, int=0); // (Joperators)
		float& operator()(axis, int, axis=UNK, int=0, axis=UNK, int=0) const; // (Joperators)
		
		/* NR methods */
		long Jlocate(double, axis a=UNK) const;				//Jnr
		long Jclosest(axis a, double) const;				//Jnr (despite no NR)
		long Jclosest_dumb(axis a, double) const;				//Jnr (despite no NR)
		float Jselect(unsigned long k);				// (Jnr)
		void Jsvdfit(vector<double>&, vector<double>&, vector<double>&, vector<double>&, void (*funcs)(double,vector<double>&));
		void Jsvdcmp(vector<vector<double> >&, vector<double>&, vector<vector<double> >&);
		vector<double> Jsvbksb(vector<vector<double> >&, vector<double>&, vector<vector<double> >&, vector<double>&);
		pair<cube, vector<double> > Jjacobi();  // (Jnr) Find eigenvalues and eigenvectors
		pair<cube, vector<int> > Jludcmp();
		vector<float> Jlubksb(pair<cube,vector<int> >, vector<float>);
		cube matrixinvert(); // in Jnr
		cube Jgaussj(cube);
		
		/* functions that alter the cube given */
		double& Xaxis(int x);
		double& Yaxis(int y);
		double& Zaxis(int z);
		double& Axis(axis a, int v);
		double Axis(axis a, int v) const;
		cube insert(cube&, int x1=0, int y1=0, int z1=0);
		void copyaxis(const cube& src, axis a);              			// (manip) modern version of the classic, copies axis a from src to *this.
		void copyaxis(const cube& src);              // (manip) copies ALL axes
		void Dtype(datatype);
		void Xtype(axistype);
		void Ytype(axistype);
		void Ztype(axistype);
		void keyword(string, string);	
		void addazprojectionhdr(string, double, double,int , int, int, float,int,int); //Jcube

		void headerwipe();
		void history_add();
		string& history();
		cube convert(axis, axistype, double param=-1, double param2=-1);	//  (manip)
		cube convert(datatype, double param=-1, double param2=-1);   		//  (manip)
		cube& sidecube(); 
		cube& backcube(); 
		cube& bottomcube();
		cube& sidecube() const; 
		cube& backcube() const; 
		cube& bottomcube() const;
		void vimswaveshift(); //Jcube.c++
					
		/* operators that do not alter the cubes they use */
		cube operator()(int, int, int, int, int, int) const;
		cube operator()(axis, int, int, axis=UNK, int=-1, int=-1, 
				axis=UNK, int=-1, int=-1) const;  // (operators) generalized subcube
		cube operator()(double, double, int, int, int, int) const;
		cube operator()(int, int, double, double, int, int) const;
		cube operator()(int, int, int, int, double, double) const;
		cube operator()(double, double, double, double, int, int) const;
		float& operator()(double, double, int);
		float operator()(double, double, int, methodtype=mFAST) const;
		float& operator()(double, int, int);
		float operator()(double, int, int, methodtype) const;
		float& operator()(int, double, int);
		float operator()(int, double, int, methodtype) const;
		float& operator()(int, int, double);
		float operator()(int, int, double, methodtype) const;
		cube operator+(cube);
		cube operator-(cube) const;
		cube operator*(cube) const;
		cube operator/(cube);
		cube operator+(float);
		cube operator-(float);
		cube operator*(float) const;
		cube operator/(float);
		cube pow(double);															// Joperators
		cube log(double = ::exp(1.));                                 
		void interpolativeZ(cube& r, cube& c1,cube &c2);
		void interpolativeX(cube& r, cube& c1,cube &c2);
		void copymetadatafrom(const cube& r);			// (Jcube)
		void copyheaderfrom(const cube& r);			// (Jcube)
		void addheaderfrom(const cube& r);			// (Jcube)
		cube convolved_with(cube);   // (manip)
		cube FFTconvolve(cube, int); //(manip) 2016.11.21 smack:: NOT YET WORKING
		vector<pair<double,double> > vimswaveshift_cub(); //Jcube.c++
		double finddate();
		void setwavelengthaxisfromhdr();
		double findexptime(); //only grabs IR exposure time!! Jcube.c++

		
		/* functions that don't alter the cube given */
		void write(const char outfile[1000],float=-1,float=-1,float (*s)(float)=0);										// Jio.c++
		void write(string,float=-1,float=-1,float (*s)(float)=0);										// Jio.c++
		void graph(const char* arg="",otype dest=oX);							// Jio
		void display(otype o=oBW, float = -1, float = -1, 
				float(*s)(float)=0, float=1.);														// JIM	
		cube tovicar();
		cube totext();
		cube tobsqraw();
		cube toJcube(int tolambda=-1);
		cube totif();
		cube tofits();
		cube toisis();
		cube astype(filetype) const;					// (Jmanip) makes new cube from old, with new type
		cube& totype(filetype);						// (Jmanip) alters *this cube
		filetype type() const;						
		cube asstorageorder(axis,axis,axis) const;					// (Jmanip) -- make this const someday
		cube& tostorageorder(axis,axis,axis);					// (Jmanip) 
		cube byteorderinvert();								//  (Jmanip)
		cube tolittleendian();									//  (Jmanip) converts to littleendian
		void setendian(endian);
		cube chunkx(int x1=-1, int x2=-1);	//Do not use these.
		cube chunky(int y1=-1, int y2=-1);  //Do not use these.
		cube chunkz(int t1=-1, int t2=-1); //Do not use these.
		cube chunk(axis, long a0, long a1=-1);
		cube chunk(axis a, int a0, int a1=-1) {return chunk(a, long(a0), long(a1));}
		cube chunk(axis, double, double) ;
		cube blocksx(cube b);
		cube blocksz(cube b);
		cube blocks(axis, const cube&);
		cube grow(int=-1, int=-1, int=-1, double=0.0);							// Jmanip
		cube grow(cube, double=0.0);													//Jmanip
		cube spliceout(axis, int, int i1=-1);				// (Jmanip)
		cube spliceout(axis, long, long i1=-1);				// (Jmanip)
		cube spliceout(axis, double, double);				// (Jmanip)
		cube smartselectivespliceout(bool (*evaluator)(cube&, axis, int));				// (Jmanip) -- evaluator returns 1 to keep, 0 to splice out
		cube skewer(axis a, long c1, long c2) const;			//  (Jmanip)
		cube plane(axis, axis, int=-1, int=-1);		//  (Jmanip)  Generalized plane to take xyplane, ytplane, etc's place
		cube plane(axis, int=-1, int=-1) const;				//  (Jmanip)  XYplane = Zplane
		cube plane(axis, double, double=-1) const;
		cube xt2xy();											//rotates the cube front face	to bottom face
		cube yt2xy();
		cube yt2xt();
		cube xt2yt();
		cube rotated(double d, axis A=Z, double lc=-1, double hc=-1) const;	// (Jmanip) rotates the cube around axis by d degrees clockwise around center 'x' 'y' (lower,higher)
		cube squareroot();
		cube cosine();    // (Jfunc.c++) cos of every input pixel is the output
		cube exp();    // (Jfunc) exponential, i.e. output pixel = e^input pixel
		cube floatfunc(float f(float));
		cube abs();												// (Jfunc) absolute value 
		cube resample(axis, cube, interpmethod=INTERPOLATE);
		cube resample(cube&, methodtype=mNN);
		cube resample(int xres, int yres=1, int zres=1); // (manip) coadd version of resample
		cube resample(double xlo, double xhi, double xres, double ylo=-1,
						  double yhi=-1, double yres=0, double zlo=-1,
						  double zhi=-1, double zres=0, bool fluxcnsrv=0); // (manip) interp version
		cube resample(axis, float); // (manip) resample by bins
		cube txt2spec();										// (Jmanip) converts annoying txt files 2/3x1xn into 1/2x1xn with proper axis values.  Works for anything with a 2, 1, and n, not just x,y,z (i.e. 1x2xn works too)
		cube noise(noisetype=GAUSSIAN,float=-1) const;		// (func) Generates noise with <n>=data
		cube dirinc();											// (manip) returns cube with monotonically increasing wavelength values
		cube dirinc(axis);									// (manip) just dirinc() axis a.
		cube abscontinuum(int interval=10);				// (func)  Generates continuum based on maximum values within a range (interval) and interpolates between
		cube dopplershift(float velocity, axis a=Z); // (manip) shifts by velocity
		cube dopplercorrelate(cube& g, double vmax, double res);   
																	// (Jfunc)  Computes the correlation of two cubes as a function of dopplershift to +/-vmax (in km) and resolution res (km) 
		cube shift(int dx, int dy=0, int dz=0,float repl=0.);		// (Jmanip) shift the image by dx,dy,dz pixels, the ones that go off the edge are lost and those that are new become zeroes.
		cube FFT(int=1);										// (Jfunc)  Takes fourier transform of input, assumed to be 0..nx/2-1 real, nx/2..nx-1 imaginary.  Returns in same format.
		cube realFFT(int=1);									// (Jfunc)  Takes fourier transform of *this, returns cube (2Nx,2Ny, 2Nz), with real components first, (0..nx-1), imaginary next (nx,2nx-1)
		cube wavelet(int=1, int=20);									// (Jfunc)  Wavelet transform, per NR
		pair<cube,pair<int,double> > Jperiod(float=4., float=1.);
		pair<cube,pair<int,double> > Jperiod_parallel(float=4., float=1.);
		cube max(axis);
		cube cylindermap(float, float, float, float =1., methodtype=mFAST, float=0., double =0., double=0., double=0.);  //  (Jfunc)
		cube cylinderrotate(cube, methodtype=mFAST); // (Jfunc)
		pair<double, double> eulerangles(double, double, double, double, double) const; // (Jfunc) -- the guts of the euler angle calculation
		static vector<double> eulerangles(vector<double>, vector<angle>); // (Jfunc) -- the guts of the euler angle calculation
		static vector<double> eulerangles_inverse(vector<double>, vector<angle>); // (Jfunc) -- backwards euler rotation
		static vector<angle> eulerangles(vector<double>, vector<double> = vector<double>());  // (Jfunc) determine the euler angle rotation to get from one vector to another
		cube cylindereuler(double, double, double, methodtype=mFAST,float (*f)(float)=0); // (Jfunc) Rotates cyl with Euler angles
		cube cylindereuler(const cube&, double, double, double, methodtype=mFAST,float(*f)(float)=0); // (Jfunc) Rotates cyl with Euler angles
		cube simplegradient(); // (Jfunc)
		cube azproject(azimuthalprojection*, methodtype=mFAST); //Jprojections.c++; 2.16 smack

		
		cube orthoproject(int, int, int, int, float, methodtype=mFAST, float=0.,
				float=0., float=0., projection* =0) const; // (Jfunc)
		cube azimuthalstereographicproject(int, int, int, int, float, methodtype=mFAST, float=0., float=0., float=0.) ;// (Jfunc), azimuthal stereographic conformal projection (smack 1-2016)
		cube lambertazimuthalproject(int, int, int, int, float, methodtype=mFAST, float=0., float=0., float=0.) ;// (Jfunc), lambert azimuthal equal area projection (smack 1-2016)
		
		
		cube circle(float, float, float, float=0.);	// (Jfunc)
		vector<float> vpolyfit(int,float (*f)(float)=0);// (Jfunc) returns vector with polynomial fit elements.
		cube cpolyfit(int, float (*f)(float)=0);		// (Jfunc) returns cube of polynomial fit curve
		cube templatefit(cube, bool=0);					// (Jfunc) returns least-squares fit of k*template to *this
		cube HPfilter(int minfrq=5);						// (Jfunc) returns high pass filtered data
		bool static boundschecking() {return Boundschecking;}
		void static boundschecking(bool b) {Boundschecking=b;}
		
		
	/* Some processing functions (not necessarily exhaustive) */
		// MAD
		cube maddumb(int zref=0,int maxx=-1,int maxy=-1);	 // (Jproc) minimum absolute difference
		cube madroughdumb(int=0,int=10,int=-1,int=-1);// (Jproc)
		vector<pair<int,int> > madraster(int=0,int=-1,int=-1);			// (Jfunc) maddumb that returns offsets
		cube madgradient(int zref=0,int maxx=-1,int maxy=-1,int=1);  // (Jproc) minimum absolute difference using gradient
		cube madgradient(cube,int zref=0,int maxx=-1,int maxy=-1);  // (Jproc) minimum absolute difference using gradient

		pair<cube, vector<double> > findprincipalcomponents() const;  // (Jalgorithms)
		cube principalcomponentprojection(cube) const;  // (Jalgorithms)  takes a cube containing PC vectors and reprojects *this into that space.
		cube principalcomponentdeprojection(cube) const;
		
		Jcrap::pixelnode classify(methodtype=mFAST);
		pair<cube,cube>& pixelclassify(Jcrap::pixelnode, cube,unsigned int);
		pair<cube,cube>& pixelclassify(Jcrap::pixelnode, cube&);
		
		// some cleanup functions
		cube clean(cube);							// (Jproc) cleans telescope data
		cube stripNaN(float=0.);			// (Jproc) get rit of NaNs, set to new double value.
		cube stripinf(float=0.);
		cube stripneg(float =0.0);												//  (manip)
		cube striplt(float limit, float defaultval=0.);					//  (manip)
		cube striple(float limit, float defaultval=0.);					//  (manip)
		cube stripgt(float limit, float defaultval=0.);					//  (manip)
		
		cube boxcarfilter(int xadd, int yadd=1, int zadd=1);  //Jproc
		
		
		
		
		/* some basic drawing capabilities */
		cube draw_circle(int,int,float,float=1., float=1., float=1.);
		cube draw_ellipse(int,int,float,float,float,float=1., float=1., float=1., bool=false);
		cube draw_line(int,int, int,int, float=1., float=1., float=1.); //(Jdraw)
		cube draw_filled_polygon(vector<pair<lonangle, latangle> >, float=1., float=1., float=1.); //(Jdraw); smack 3.16
		cube draw_filled_polygon(vector<pair<double, double> >, float=1., float=1., float=1.); //(Jdraw)
		void draw_filled_polygon_on(vector<pair<double, double> >, float=1., float=1., float=1.); //(Jdraw)
		static bool is_inside_polygon(vector<pair<double, double> >, pair<double, double>); // (Jdraw)
		cube autocrop(int=0, int=0, int=0);
		
		/* multiple cube functions */
		float* polyfit(int degree);
		float interp(double, int d=3, axis a=UNK);
		float interp2d(float, float, int d=3) const;
		
		/* subroutines associated with I/O */
		void createarrays(float=0.0);									// (Jcube)
		void IMread(const char []);							//JIM.c++
		void IMwrite(const char [],float=-1,float=-1,float (*s)(float)=0);							//Jim
		void FITSwrite(const char []) const;
		void FITShdrread(const char []);
		void FITSreaddata(const char []);
		void determinetype(ifstream *, const char []);	//Jio.c++
		int getlabelsize(ifstream *);
		void read(string, int=-1);
		void readheader(ifstream *, const char []);
		void readdata(ifstream *, const char []);
		float getfloat(ifstream *);
		float getbigendfloat(ifstream *);
		short getshort(ifstream *);
		long get32int(ifstream *);
		int pad(string&, int);
		void putheader(ofstream *);
		void writedata(ofstream *);
		void putfloat(ofstream*, const float&) const;
		void putfloat(fstream*, const float&) const;
		inline void putshort(ofstream* &, short&);
		static void Jraninit(unsigned long long=0);
		static double Jrandom();
		static unsigned long long Jrandom_ull();
		static unsigned long long Jrandom_ull(unsigned long long&);
		
	private:   // some functions internal to Jcube
		void pageout() const;							// Joperators -- copies changes in memory to disk for bigcubes.

//		cube copyaxis(float [], float [], float []);  // (manip) -- legacy.  Phase out when you get time.
		
		int choosememdims() const;
		
	private:
		mutable map<string,string,less<string> > hdr;
		mutable vector<float> Data;   // totally lame, I know, should split "Data" into Data and mutable DataPage for bigcubes
		double *xaxis, *yaxis, *zaxis;
		long nx, ny, nz;
		axis storageSLOW, storageMEDIUM, storageFAST;   // added 2006/03/22 for storage order flexibility
		cube *pbackcube, *pbottomcube, *psidecube;
		filetype cubetype;
		endian byteorder;
		
//		Isis::Cube isis3cube;
		
		static bool Boundschecking;  // for use in operator()
		static vector<unsigned long long> _random_v;
				
		// for big cubes & swapping
//		int bigcube;
		int memdims;
		mutable vector<long long> current;
		long long begindata;
		string cubefilename;
	public:
		static long long maxmemsize;
		static bool debug;
};

//  Functions for interaction with the STL
ostream& operator<<(ostream& out, const cube& c);			//Jio
ostream& operator<<(ostream& out, axis a);			//Jio
ostream& operator<<(ostream& out, datatype d);		//Jio
ostream& operator<<(ostream& out, axistype t);		//Jio

//  Supplementary functions
void raninit();												// (Jmisc) initializes random number generator ran1

//  Functions that don't really use a single cube
void fourier(cube hr, cube hi, cube *Hr, cube *Hi, int inv=1);
cube blackbody(cube,float,datatype=FLUX);  // Jfunc
cube photonblackbody(cube,float,datatype=FLUX);  // Jfunc




template<class T> struct printpairendl : public unary_function<T, void>
{
	printpairendl(ofstream *foutin) : fout(foutin) {}
	void operator() (T x) {
		(*fout) << x.first << ' ' << x.second << "\n";
	}
	ofstream *fout;
};

template<class T> struct printpair : public unary_function<T, void>
{
	printpair(ofstream *foutin) : fout(foutin) {}
	void operator() (T x) {
		(*fout) << x.first << ' ' << x.second << "\n";
	}
	ofstream *fout;
};

