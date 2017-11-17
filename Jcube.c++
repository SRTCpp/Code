#include<string>
#include<cstring>
#include<cstdio>
#include<map>
#include<iostream>
#include<math.h>
#include<strstream>
#include<fstream>
#include<utility> 
#include<algorithm>
#include<time.h>

#include<gdal_priv.h>


#include"Jcube.h"
 

long long cube::maxmemsize(1024LL*1024LL*1024LL*8LL);  // about 8 GB

bool cube::debug(0);
bool cube::Boundschecking(false);
int graphno;

vector<unsigned long long> cube::_random_v(0);

float byteorderswap (float a)
{
   union u {float vi; unsigned char c[sizeof(float)];}; 
   union v {float ni; unsigned char d[sizeof(float)];};
   union u un; 
   union v vn; 
   un.vi = a; 
   vn.d[0]=un.c[3]; 
   vn.d[1]=un.c[2]; 
   vn.d[2]=un.c[1]; 
   vn.d[3]=un.c[0]; 
   return (vn.ni); 
}

cube::cube(string fname, int stordims) :
	psidecube(0), pbottomcube(0), pbackcube(0)
{
	Boundschecking=false;
	cout << "creating cube from file:  " << fname << " -- ";
	cubefilename=fname;
	read(fname, stordims);

}

void cube::read(string fname, int stordims) 
// split from cube::cube(string) 2006 February
{
	memdims = stordims;
	
	ifstream fin(fname.c_str());
	if (!fin.is_open()) cout << "NO SUCH FILE " << fname << " OR OTHER FILE ERROR\n";;
	determinetype(&fin, fname.c_str());
	if (!osuppress) cout << "CUBETYPE is " << cubetype << "\n";
	if(cubetype==gif || cubetype==jpg || cubetype==tif){
		fin.close();
		IMread(fname.c_str());
	}  else {

		if (osuppress < -2) cout << "Going to readheader\n";
		readheader(&fin, fname.c_str());
		if (osuppress < -2) cout << "header read in successfully\n";
		byteorder=littleendian; 				// unless otherwise indicated . . .
		if(cubetype==vicar){
			nx=str2ulong(hdr["NS"]);
			ny=str2ulong(hdr["NL"]);
			nz=str2ulong(hdr["NB"]);
			if (hdr["INTFMT"]=="'HIGH'") byteorder=bigendian;
		}
		else if(cubetype==isis){
			if (hdr["AXIS_NAME"]=="(SAMPLE,BAND,LINE)" || hdr["AXES_NAME"]=="(SAMPLE,BAND,LINE)")
				sscanf( (hdr["CORE_ITEMS"]).c_str(), "(%ld,%ld,%ld)", &nx, &nz, &ny);
			if (hdr["AXIS_NAME"]=="(SAMPLE, BAND, LINE)" || hdr["AXES_NAME"]=="(SAMPLE,BAND,LINE)")
				sscanf( (hdr["CORE_ITEMS"]).c_str(), "(%ld,%ld,%ld)", &nx, &nz, &ny);
			if (hdr["AXIS_NAME"]=="(SAMPLE,LINE,BAND)" || hdr["AXES_NAME"]=="(SAMPLE,LINE,BAND)")
				sscanf( (hdr["CORE_ITEMS"]).c_str(), "(%ld,%ld,%ld)", &nx, &ny, &nz);
			long sideplanes, backplanes, bottomplanes;
			sscanf( (hdr["SUFFIX_ITEMS"]).c_str(), "(%ld,%ld,%ld", &sideplanes,
					&backplanes, &bottomplanes);
			if (osuppress < -2) cout << "Nx=" << nx << ", Ny=" << ny << ", Nz=" << nz << "\n";
/*			nx+=sideplanes;
			ny+=bottomplanes;
			nz+=backplanes;*/
			if (psidecube != 0) {delete psidecube; psidecube=0;}
			if (sideplanes) psidecube = new cube(sideplanes, ny, nz);
			if (pbottomcube != 0) {delete pbottomcube; pbottomcube=0;}
			if (bottomplanes) pbottomcube = new cube(nx, bottomplanes, nz);
			if (pbackcube != 0) {delete pbackcube; pbackcube=0;}
			if (backplanes) pbackcube = new cube(nx, ny, backplanes);
			if (hdr["CORE_ITEM_TYPE"].find("SUN_REAL")!=string::npos) byteorder=bigendian;
			
					
		} else if (cubetype==isis3) {  // GDAL: added 2016 August 9 JWB
			GDALDataset *inDataset;
			GDALAllRegister();
			inDataset = (GDALDataset *) GDALOpen(fname.c_str(), GA_ReadOnly);
			if (inDataset == NULL)
				cout << "Attempted GDAL readin of cube fails to open file\n";
			else {
				nz = inDataset->GetRasterCount();
				cout << "read in data with " << inDataset->GetRasterCount() << " bands\n";
				GDALRasterBand *xyTestBand;
				xyTestBand = inDataset->GetRasterBand(1);
				nx = xyTestBand->GetXSize();
				ny = xyTestBand->GetYSize();
				
				cout << "x is apparently " << nx << "; y is " << ny << "\n";
			}
			GDALClose(inDataset);
		} else if(cubetype==Jcube1 || cubetype==Jcube1a || cubetype==Jcube1b){
			nx=str2ulong(hdr["Nx"]);
			ny=str2ulong(hdr["Ny"]);
			nz=str2ulong(hdr["Nz"]);
		}	
		else if(cubetype==fits);
		
		if (memdims < 0  ||  memdims > 3) memdims = choosememdims();
		
		current.resize(3);
		current[X] = -1;
		current[Y] = -1;
		current[Z] = -1;
		
		readdata(&fin, fname.c_str());
		if (cubetype==isis) setwavelengthaxisfromhdr();
		
		if (!osuppress) cout << "read in " << nx << ", " << ny << ", "<< nz << ", memdims=" << memdims << "\n";
		
		if (byteorder==bigendian) {
			cout << "Swapping byte order.\n";
			for (cube::iterator i=begin();i!=end();i++)
				*i = byteorderswap(*i);
		}

	}
}

int cube::choosememdims() const
{
//	cout << "in choosememdims, storageSLOW = " << storageSLOW;
	int answer(3);
	if (N(storageFAST) > maxmemsize)
		answer = 0;
	else if (N(storageFAST) * N(storageSLOW) > maxmemsize)
		answer = 1;
	else if (N() > maxmemsize) 
		answer = 2;
	else 
		answer = 3;
	
//	cout << " memdims chosen to be " << answer << "\n";
	
	return answer;
}

cube::cube(long x, long y, long z, float initval, filetype type, axis a100, axis a10, axis a1, int dimnum) :
	memdims(dimnum),
	psidecube(0), pbottomcube(0), pbackcube(0)
{
	if(debug) cout << "creating cube ID(" << this << ") dim " << x << ", " << y << ", ";
	if(debug) cout << z << " from scratch\n";
	
	nx=x;
	ny=y;
	nz=z;
	cubetype=type;
	byteorder=littleendian;
	Boundschecking=false;

	storageSLOW   = a100;
	storageMEDIUM =  a10;
	storageFAST   =   a1;
	
	xaxis = new double[nx];
	yaxis = new double[ny];
	zaxis = new double[nz];
	for(int x=0;x<nx;x++) xaxis[x]=double(x);		//setting up the axis calibration
	for(int y=0;y<ny;y++) yaxis[y]=double(y);		//labels.
	for(int z=0;z<nz;z++) zaxis[z]=double(z);
	
	if(type==vicar){
		hdr["N1"]=lblstr(nx);					// L8 + 2
		hdr["N2"]=lblstr(ny);					// L8 + 2
		hdr["N3"]=lblstr(nz);					// L8 + 2
		hdr["NL"]=hdr["N2"];						//same as above for these 3
		hdr["NS"]=hdr["N1"];
		hdr["NB"]=hdr["N3"];
		hdr["DIM"]="3";							// L5 + 2
		hdr["FORMAT"]="'REAL'";					// L13 + 2
		hdr["TYPE"]="'IMAGE'";					// L12 + 2
		hdr["EOL"]="0";							// L5 + 2
		hdr["NBB"]="0";							// L5 + 2
		hdr["NLB"]="0";							// L5 + 2
		hdr["ORG"]="'BSQ'";						// L9 + 2
		hdr["TASK"]="'BarnesProgram'";		// L20 + 2
		hdr["RECSIZE"]=lblstr(4*nx);			// L13 + 2
		hdr["BUFSIZE"]=hdr["RECSIZE"];		// L13 + 2
	}	
	if (type==Jcube1){
		time_t tp;
		hdr["Nx"]=int2str(nx);
		hdr["Ny"]=int2str(ny);
		hdr["Nz"]=int2str(nz);
		hdr["Jcubever"]="1.0";
		hdr["Xaxistype"]="0";
		hdr["Yaxistype"]="0";
		hdr["Zaxistype"]="0";
		hdr["datatype"]="0";
		time(&tp);
		char timeasc[80];
		strcpy(timeasc, ctime(&tp));
		for (int i=0;i<strlen(timeasc);i++) if (timeasc[i]=='\n') timeasc[i]=' ';		
		hdr["Created"]=timeasc;
	}	
	if (type==Jcube1a){
		time_t tp;
		hdr["Nx"]=int2str(nx);
		hdr["Ny"]=int2str(ny);
		hdr["Nz"]=int2str(nz);
		hdr["Jcubever"]="1.1";
		hdr["Xaxistype"]="0";
		hdr["Yaxistype"]="0";
		hdr["Zaxistype"]="0";
		hdr["datatype"]="0";
		time(&tp);
		char timeasc[80];
		strcpy(timeasc, ctime(&tp));
		for (int i=0;i<strlen(timeasc);i++) if (timeasc[i]=='\n') timeasc[i]=' ';		
		hdr["Created"]=timeasc;
	}	
	if (type==Jcube1b){
		time_t tp;
		hdr["Nx"]=int2str(nx);
		hdr["Ny"]=int2str(ny);
		hdr["Nz"]=int2str(nz);
		hdr["Jcubever"]="1.2";
		hdr["Xaxistype"]="0";
		hdr["Yaxistype"]="0";
		hdr["Zaxistype"]="0";
		hdr["datatype"]="0";
		time(&tp);
		char timeasc[80];
		strcpy(timeasc, ctime(&tp));
		for (int i=0;i<strlen(timeasc);i++) if (timeasc[i]=='\n') timeasc[i]=' ';		
		hdr["Created"]=timeasc;
	}
	
	if (memdims == -1) memdims = choosememdims();
	
	if ((cubetype==Jcube1b || cubetype==Jcube1a || cubetype==isis || cubetype==Jcube1
			) && memdims!=3) {
		current.resize(3);
		current[X]=-1LL;
		current[Y]=-1LL;
		current[Z]=-1LL;
		
		if (!osuppress) cout << "this is a big cube ("<<nx<<"x"<<ny<<"x"<<nz<<"), memdims=" << memdims << " ";
		cubefilename=string("Jcubetmp.") + 
				int2str((long long)(this)) +
				string(".Jcube");
		
		ofstream fo;
		fo.open(cubefilename.c_str(), ios::trunc);
		
		putheader(&fo);
			
		for (int x=0;x<nx;x++) putfloat(&fo, xaxis[x]);
		for (int y=0;y<ny;y++) putfloat(&fo, yaxis[y]);
		for (int z=0;z<nz;z++) putfloat(&fo, zaxis[z]);
		begindata = fo.tellp();
		if (debug) cout << "begindata = " << begindata << "\n";
		for (long long i=0;i<(long long)(N(X))*(long long)(N(Y))*(long long)(N(Z));i++)
			putfloat(&fo, initval);
		if (debug) cout << "disking as " << cubefilename << "\n"; cout.flush();
		
		fo.close();
	} else
		memdims = 3;
	
	if (debug) cout << "about to initialize arrays with memdims=" << memdims << "\n";
	createarrays(initval);
	if (debug) cout << "arrays initialized successfully\n";
}

void cube::createarrays(float initval)
// Altered 2012 July 22 to try to address allocation issue for large cubes
// added 'newsize' in size_t space instead of int.
{	
	size_t newsize;
	
	if (memdims==3)
		newsize = size_t(nx)*size_t(ny)*size_t(nz);
	else if (memdims==2)
		newsize = size_t(N(storageMEDIUM))*size_t(N(storageFAST));
	else if (memdims==1)
		newsize = size_t(N(storageFAST));
	else if (memdims==0)
		newsize = 1;
	
	Data.resize(newsize, initval);
}
	
double cube::findexptime()
{
	// Only finds IR exp time (that's how keyword works?)
	double answer(0.);

	string a("EXPOSURE_DURATION");
	string irtimestring(keyword(a));
	
	sscanf(irtimestring.c_str(), "%lf",&answer);

	return answer;
}

double findshift(double date){
	
	double answer(0.),shift(0.);
	bool found(false);
	string inshiftcal("/vims/vimswavecal_simple.txt");
	ifstream inshift(inshiftcal.c_str());
	pair< double, double> thisshift;
	vector< pair<double, double> > shiftvector;
		
	getline(inshift,inshiftcal);
	while(!inshift.eof()){
		inshift >> thisshift.first >> thisshift.second;
		//cout << thisshift.first <<"\t"<<thisshift.second<<endl;	
		shiftvector.push_back(thisshift);
		getline(inshift,inshiftcal);
	}
	
	for (int l(0); l<shiftvector.size();l++){
		if (found == false){
			if (date == shiftvector.at(l).first) {
				shift = shiftvector.at(l).second;
				found == true;
			}
		}
	}
	
	answer=shift;
	return answer;
	
}

double cube::finddate(){
		double answer(0.), year(0.), doy(0.), date(0.), shift(0.); 
	
		string a("START_TIME"),datey,dated;
	
		string datestr (keyword(a));	
		//cout <<"Date string from header\t"<<datestr<<endl;
		datestr = datestr.substr(0,datestr.find("T"));
	//	cout <<"Date string cropped \t"<<datestr<<endl;
			
		if(type() !=isis){
			datey=datestr.substr(3,4);
		//	cout <<"Date year\t"<<datey<<endl;
			dated=datestr.substr(datestr.find("-")+1,datestr.length());
		//	cout <<"Date day\t"<<dated<<endl;
				
		}	
		else if (type() == isis){
			
			datey=datestr.substr(1,4);
		//	cout <<"Date year\t"<<datey<<endl;
			dated=datestr.substr(6,datestr.find("T"));
		//	cout <<"Date day\t"<<dated<<endl;
			
		}
		
		sscanf(datey.c_str(), "%lf",&year);
		sscanf(dated.c_str(),"%lf",&doy);
		//cout <<"Date:\t"<<date<<endl;
		date = year+(1./365.25)*doy;	
		cout <<"Date #:\t"<<date<<endl;
		answer = date;
		return answer;
}

void cube::setwavelengthaxisfromhdr()
{ // created: 16 Aug 2016, smack 
 //  set Z axis to wavelength values stored in the header file
 // replaces the need to import wavelength.cal.txt
// AND with newcalpipeline cubs, no need to apply the old wavelength shift proxy
 
	cout <<"\nGoing to wavelength axis from header\n";	
	// add in functionality to check if band bin center exists
	string wavelengthstring(keyword("BAND_BIN_CENTER"));
	wavelengthstring.erase(0,wavelengthstring.find("(")+1);
	int pos(0);
	vector<double> wavelengthvector;
			
		
	while ((pos =wavelengthstring.find(",")) != string::npos){
		string substring;
				
		substring=wavelengthstring.substr(0,pos);
				
		wavelengthvector.push_back(str2double(substring));
		wavelengthstring.erase(0,pos+1);
	}
	wavelengthvector.push_back(str2double(wavelengthstring.substr(0,wavelengthstring.find(")"))));
	
	for (int z=0; z<N(Z); z++) Axis(Z,z)=wavelengthvector.at(z);
			
	
	cout <<"\nSet wavelength axis from header\n";
}

void cube::vimswaveshift(){
	
	if (type()==isis) "Use cub-specific waveshift function instead.";
	else{
		double date(finddate()); 
		
		
		double longdate = round(trunc(date*100)/10);
		//cout<<"Longdate:\t"<<longdate<<endl;
		date = (longdate-fmod(longdate,5))/10;
		//cout <<"Date round to multpile:\t"<<date<<endl;
	
	
	
		double shift=findshift(date);
		cout <<"SHIFT::\t"<<shift<<endl;
 		
		for (int z(0);z<N(Z);z++){
			//cout<<"Preshift:\t"<<Axis(Z,z)<<"\t";
			Axis(Z,z)=Axis(Z,z) +shift;
			//cout<<"Postshift:\t"<<Axis(Z,z)<<endl;
		}
	}

}


vector<pair<double,double> > cube::vimswaveshift_cub(){
	/*output is vector of plane#, shiftedwavelength,
	NOTE: This doens't depend on I/F at all-- the cube fed in is 
	only used to calculate the date. */
	vector<pair<double,double> > answer;
	
	if (type()!=isis) "Use cyl-specific waveshift function instead.";
	
	else{
		
		if (osuppress >=2) cout<<"In vimswaveshift_cub method\n";
		
		double date(finddate()); 
				
		double longdate = round(trunc(date*100)/10);
		//cout<<"Longdate:\t"<<longdate<<endl;
		date = (longdate-fmod(longdate,5))/10;
		//cout <<"Date round to multpile:\t"<<date<<endl;
	
		double shift(findshift(date));
		if (osuppress >=2) cout <<"SHIFT::\t"<<shift<<endl;
 		
		string inwavesfile;
		inwavesfile="/vims/wavelength.cal.txt";
	
		ifstream inwaves(inwavesfile.c_str());
		pair<double, double> thisplanewave;
		vector< pair<double, double> > planewavevector;
	
		getline(inwaves,inwavesfile);
		int plane(0);
		while(!inwaves.eof()){
			double thiswave(0.);
			thisplanewave.first=plane;							
			inwaves >> thiswave;
			thisplanewave.second = thiswave+shift;
			//cout << thisplanewave.first <<"\t"<<thisplanewave.second<<endl;	
			planewavevector.push_back(thisplanewave);
			getline(inwaves,inwavesfile);
			plane=plane+1;
		}
		
		answer=planewavevector;
	}
	
	return answer;
}


void cube::copymetadatafrom(const cube& r)
{
	hdr = r.hdr;
	cubetype = r.cubetype;
	byteorder = r.byteorder;
	
	storageSLOW   = r.storageSLOW;
	storageMEDIUM = r.storageMEDIUM;
	storageFAST   = r.storageFAST;
}

void cube::copyheaderfrom(const cube& r)
{
	hdr = r.hdr;
}

void cube::addheaderfrom(const cube& r)
{
	for (map<string, string,less<string> >::iterator i(r.hdr.begin());i!=r.hdr.end();i++)
		hdr[i->first]=i->second;
}

cube::cube(axis a1, int n1, axis a2, int n2, axis a3, int n3, float initval,
		filetype ftype, axis a100, axis a010, axis a001, int stormems)
/* JB 10/31/2001  how scary I'm messing with cube constructors . . .
	Allows for easier axis-independant cube generation */		
{
	int x(1), y(1), z(1);
	
	if (a1 == X) x=n1;
	if (a1 == Y) y=n1;
	if (a1 == Z) z=n1;
	
	if (a2 == X) x=n2;
	if (a2 == Y) y=n2;
	if (a2 == Z) z=n2;
	
	if (a3 == X) x=n3;
	if (a3 == Y) y=n3;
	if (a3 == Z) z=n3;
	
	nx=ny=nz=0;	
	psidecube=pbackcube=pbottomcube=0;
	if (debug) cout << "calling fixed-axis cube() constructor\n";
	*this=cube(x,y,z,initval,ftype,a100,a010,a001,stormems);	 
}


cube::cube(const cube& in, int inmem) :
	psidecube(0), pbottomcube(0), pbackcube(0)
{
	if(debug) cout << "In cube(cube) constructor\n";
	nx=ny=nz=0;
	memdims=inmem;
	cubefilename="";
	*this=in;
}

cube::cube(float v) :
	psidecube(0), pbottomcube(0), pbackcube(0)
{
	cube* pa;
	pa=new cube(1,1,1,v);
	*this=*pa;
}

cube::~cube()
{
	if(debug) cout << "destructing cube of ID(" << this << ")\n";
	
	deallocate();
}

void cube::deallocate()
// split off from ~cube 2006 Feb 14 to consolidate this function from here and op=
{
	if (debug) cout << "Deallocating " << this << ":  axes . . ."; if(debug)cout.flush();
	if (nz) delete []xaxis;
	if (nz) delete []yaxis;
	if (nz) delete []zaxis;
	
	
	if (debug) cout <<"  sidecubes . . ."; if(debug)cout.flush();
	if (psidecube != 0) {delete psidecube; psidecube=0;}
	if (pbottomcube != 0) {delete pbottomcube; pbottomcube=0;}
	if (pbackcube != 0) {delete pbackcube; pbackcube=0;}
	
	if (debug) cout << "  tempfiles . . ."; if(debug)cout.flush();
	if (memdims!=3) {
		string command("rm ");
		command += cubefilename;
		if (cubefilename.find("Jcubetmp") == 0)
			system(command.c_str());
	}		
	if (debug) cout << "  done.\n";
}

long cube::Nx() const {return (*this).nx;}	
long cube::Ny() const {return (*this).ny;}
long cube::Nz() const {return (*this).nz;}
long long cube::N() const {return (long long)(nx)*(long long)(ny)*(long long)(nz);}
long cube::N(axis a) const /*created 3/6/2k for axis generalization */ {
	if (a==X) return nx;
	if (a==Y) return ny;
	if (a==Z) return nz;
	if (a==UNK) cout << "N(axis) error, given UNK";
	return 1;
}
long cube::Side()  const {if (psidecube!=0) return (*this).sidecube().N(X); else return 0;}
long cube::Bottom()const {if (pbottomcube!=0)return(*this).bottomcube().N(Y); else return 0;}
long cube::Back()  const {if (pbackcube!=0) return (*this).backcube().N(Z); else return 0;}
double& cube::Xaxis(int x) { 
	if (x>=0 && x<nx) return xaxis[x]; 
	else cout << "Error in Xaxis(), tried to get "<<x<<"from cube xdimensionsize "<<nx<<"\n";
}
double& cube::Yaxis(int y) { 
	if (y>=0 && y<ny) return yaxis[y]; 
	else cout << "Error in Yaxis(), tried to get "<<y<<"from cube ydimensionsize "<<ny<<"\n";
}
double& cube::Zaxis(int z) {
	if (z>=0 && z<nz) return zaxis[z];
	else cout << "Error in Zaxis(), tried to get "<<z<<"from cube zdimensionsize "<<nz<<"\n";
}
double& cube::Axis(axis a, int v) /*created 3/6/2k*/ {
	if (a==UNK){ cout << "ERROR -- in axis, given UNK\n"; a=X; }
	if (a==X) return Xaxis(v);
	if (a==Y) return Yaxis(v);
	if (a==Z) return Zaxis(v);
}
double cube::Axis(axis a, int v) const /*created 2006 Feb 14*/ {
	if (a==UNK){ cout << "ERROR -- in axis, given UNK\n"; a=X; }
	if (a==X) return xaxis[v];
	if (a==Y) return yaxis[v];
	if (a==Z) return zaxis[v];
}
//double* cube::Xaxisarray() {return xaxis;}  // removed these 2012 August 22
//double* cube::Yaxisarray() {return yaxis;}  // they're ugly and unnecessary
//double* cube::Zaxisarray() {return zaxis;}
datatype cube::Dtype() const {return (datatype)str2int(hdr["datatype"]);}
axistype cube::Xtype() const {return (axistype)str2int(hdr["Xaxistype"]);}
axistype cube::Ytype() const {return (axistype)str2int(hdr["Yaxistype"]);}
axistype cube::Ztype() const {return (axistype)str2int(hdr["Zaxistype"]);}
void cube::Dtype(datatype t) {hdr["datatype"]=int2str((int)t);}
void cube::Xtype(axistype t) {hdr["Xaxistype"]=int2str((int)t);}
void cube::Ytype(axistype t) {hdr["Yaxistype"]=int2str((int)t);}
void cube::Ztype(axistype t) {hdr["Zaxistype"]=int2str((int)t);}
void cube::keyword(string a, string b) {hdr[a]=b;}
string& cube::keyword(string a) {return hdr[a];}
string cube::keyword(string a) const {return hdr[a];}


void cube::addazprojectionhdr(string thistype, double subSClat_deg, double subSClon_deg,
		int centerx_pix, int centery_pix, int radius_pix, float Nangle, int xsize, int ysize)
{
	hdr["cylindrical_map"]= "no";
	hdr["azimuthal_map"]=thistype;
	hdr["subSClat_Nneg"]=double2fstr(subSClat_deg);
	hdr["subSClon"]=double2fstr(subSClon_deg);
	hdr["centerx"]=int2str(centerx_pix);
	hdr["centery"]=int2str(centery_pix);
	hdr["radius_in_pixels"]=int2str(radius_pix);
	hdr["Nangle_Clockwisepos"]=float2str(Nangle);
		
}


double cube::res(axis a)
{
	(*this)=dirinc();
	long n;
	double *arr;
	if (a==UNK) a=Z;
	if (a==X) { n=nx; arr=&Axis(X,0); }
	if (a==Y) { n=ny; arr=&Axis(Y,0); }
	if (a==Z) { n=nz; arr=&Axis(Z,0); }

	cube diffs(n-1,1,1);
	for (long i=0;i<n-1;i++)
		diffs(i,0,0) = arr[i+1]-arr[i];
	cube m;
	m=diffs.xmedian();
	return m(0,0,0);
}

axis Uax(axis a) { // returns next highest axis
	axis b;
	if (a==X) b=Y;
	if (a==Y) b=Z;
	if (a==Z) b=X;
	if (a==UNK) b=X;
	return b;
}

axis Dax(axis a) { // returns next lowest axis
	axis b;
	if (a==X) b=Z;
	if (a==Y) b=X;
	if (a==Z) b=Y;
	if (a==UNK) b=Y;
	return b;
}

axis otheraxis(axis a, axis b) 
{
	axis answer(UNK);
	
	if (a==UNK || b==UNK) answer = UNK;
	else answer = axis(3-int(a)-int(b));  // cute, huh?
	
	return answer;
}

// added 12/12/2k --
vector<float>::iterator cube::Databegin() { return Data.begin(); }
		
cube::iterator::iterator() : curr(0), i(0), thiscube(0) { }

cube::iterator::iterator(vector<float>::iterator rit) :
	curr(rit), i(0), thiscube(0) {}

cube::iterator::iterator(long long n, cube* inc) :
	curr(0), i(n), thiscube(inc) {}

bool operator==(const cube::iterator& l, const cube::iterator& r)
{
	if (l.thiscube && r.thiscube)
		return (l.i == r.i);
	else if (l.thiscube == 0  &&  r.thiscube == 0)
		return l.curr == r.curr;
	return 0;
}	
	
bool operator!=(const cube::iterator& l, const cube::iterator& r)
{
	if (l.thiscube && r.thiscube)
		return l.i != r.i;
	else if (l.thiscube == 0  &&  r.thiscube == 0)
		return l.curr != r.curr;
	return 0;
}

float* cube::iterator::operator=(iterator r)
{
	curr=r.curr;
	thiscube = r.thiscube;
	i = r.i;
}

float& cube::iterator::operator*() 
{
	if (thiscube==0)
		return *curr;	 
	if (thiscube) {
		int s,m,f;
		s = i/(thiscube->N(thiscube->storageMEDIUM)*thiscube->N(thiscube->storageFAST));
		m = i%(thiscube->N(thiscube->storageMEDIUM)*thiscube->N(thiscube->storageFAST))/thiscube->N(thiscube->storageFAST);
		f = i%(thiscube->N(thiscube->storageMEDIUM)*thiscube->N(thiscube->storageFAST))%thiscube->N(thiscube->storageFAST);
		return thiscube->storageorderaccess(s,m,f);
	}
}

float* cube::iterator::operator->() 
{ 
	return &(**this);
}

cube::iterator cube::iterator::operator++()
{
	if (thiscube==0)
		++curr;
	if (thiscube)  // disk-stored cube
		++i;
	return *this;
}

cube::iterator cube::iterator::operator++(int)
{
	cube::iterator tmp(*this);
	++(*this);
	return tmp;
}

cube::iterator cube::begin() 
{
	if (memdims==3)
		return cube::iterator(Data.begin());
	else
		return cube::iterator(0, this);
}

cube::iterator cube::end()   
{
	if (memdims==3)
		return cube::iterator(Data.end());
	else
		return cube::iterator(N(), this);
}

string & cube::history()
{
  return keyword("history_record");
}

cube& cube::sidecube()
{
	if (psidecube != 0) return *psidecube;
	else {
		psidecube = new cube(0, ny, nz);
		return *psidecube;
	}
}

cube& cube::backcube()
{	
	if (pbackcube != 0) return *pbackcube;
	else {
		pbackcube = new cube(nx, ny, 0);
		return *pbackcube;
	}
}

cube& cube::bottomcube()
{	
	if (pbottomcube != 0) return *pbottomcube;
	else {
		pbottomcube = new cube(nx, 0, nz);
		return *pbottomcube;
	}
}

cube& cube::sidecube() const
{
	if (psidecube != 0) return *psidecube;
}

cube& cube::backcube() const
{	
	if (pbackcube != 0) return *pbackcube;
}

cube& cube::bottomcube() const
{	
	if (pbottomcube != 0) return *pbottomcube;
}

filetype cube::type() const
{
	return cubetype;
}


