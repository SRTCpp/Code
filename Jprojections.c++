//#include "Jcrap.h"
#include "Jprojections.h"

string projtype2str(projtype p) 
{
	string s;
	if (p == CYL) s="CYL";
	else if (p == ORTHO) s="ORTHO";
	else if (p == AZSTEREOGR) s="AZSTEREOGR";
	else if (p == LAMBAZ) s="LAMBAZ";
	
	return s;
}



azimuthalprojection::azimuthalprojection(cube c)
			:projection(lonangle(angle(0.,angle::DEG),CENTER_NEGATIVE180),latangle(angle(0.,angle::DEG),NORTH_POSITIVE))
{
	if (c.keyword("cylindrical_map") == " no" | c.keyword("cylindrical_map") == "no"){
		
		angle readlat(str2double(c.keyword("subSClat_Nneg")), angle::DEG);
		angle readlon(str2double(c.keyword("subSClon")), angle::DEG);
		angle readNangle(str2double(c.keyword("Nangle_Clockwisepos")), angle::DEG);
		
		Nangle(readNangle);
		
		if (c.keyword("lonconvention")== " CENTER_NEGATIVE180" | c.keyword("lonconvention")== "CENTER_NEGATIVE180") {
			subSClon(lonangle(readlon,CENTER_NEGATIVE180));
			originatingcyllonconvention(CENTER_NEGATIVE180);
		} else if (c.keyword("lonconvention")== " CENTER_0" | c.keyword("lonconvention")== "CENTER_0") {
			subSClon(lonangle(readlon,CENTER_0));
			originatingcyllonconvention(CENTER_0);
		} else if (c.keyword("lonconvention")== " CENTER_POSITIVE180" | c.keyword("lonconvention")== "CENTER_POSITIVE180") {
			subSClon(lonangle(readlon,CENTER_POSITIVE180));
		} 
		else if (c.Axis(X,0)<-180.) originatingcyllonconvention(CENTER_NEGATIVE180);
		else if (c.Axis(X,c.N(X)-1)>0) originatingcyllonconvention(CENTER_0);
			
				
		subSClat(latangle(readlat,NORTH_NEGATIVE));
		
		radius(str2double(c.keyword("radius_in_pixels")));
		xcenter(str2int(c.keyword("centerx")));
		ycenter(str2int(c.keyword("centery")));
		xsize(str2int(c.keyword("Nx")));
		ysize(str2int(c.keyword("Ny")));
		
		
	}
	else if (c.keyword("cylindrical_map")==" yes") cout <<"------ Nonazimuthal cube detected. Cannot instantiate azimuthalproject object\n";
	else cout <<"Initializing empty azimuthal object.\n";
}	




 azimuthalprojection::azimuthalprojection(lonangle inSClon, latangle inSClat, int sizex, int sizey,
		 int centerx, int centery, double inrad, angle inNangle)
		:projection(inSClon,inSClat)
{
	xsize(sizex);
	ysize(sizey);
	xcenter(centerx);
	ycenter(centery);
	radius(inrad);
	
} 


void azimuthalprojection::angdist(double c)
{
	_angdist=c;	
}


void projection::subSClon(lonangle a)
{
	_subSClon=a;
}
void projection::subSClat(latangle a)
{
	_subSClat=a;
}
void azimuthalprojection::Nangle(angle a)
{
	_Nangle=a;
}
void azimuthalprojection::radius(double d)
{
	_radius=d;
}
void azimuthalprojection::xcenter(int i)
{
	_xcenter=i;
}
void azimuthalprojection::ycenter(int i)
{
	_ycenter=i;
}
void azimuthalprojection::ysize(int i)
{
	_ysize=i;
}
void azimuthalprojection::xsize(int i)
{
	_xsize=i;
}

void azimuthalprojection::originatingcyllonconvention(lonconvention l)
{
	_originatingcyllonconvention=l;	
}


pair<double, double> azimuthalprojection::operator()(lonangle lon,latangle lat)
{
	pair<double, double> answer;
	answer.first=0.;
	answer.second=0.;
	double x(0.),y(0.),rotx(0.),roty(0.),rely(0.),relx(0.);

	
	
	x = kfactor(lon,lat)*cos(lat.north_positive())*sin(lon.center_negative180()-subSClon().center_negative180());
	y = kfactor(lon,lat)*(cos(subSClat().north_positive())*sin(lat.north_positive())-sin(subSClat().north_positive())*cos(lat.north_positive())*cos(lon.center_negative180()-subSClon().center_negative180()));
		 
	rotx = x;
	roty = -y;
		
	answer.first=rotx+xcenter();
	answer.second=roty+ycenter();

	
	
	return answer;	
	
}


pair<lonangle, latangle> azimuthalprojection::operator()(double inx, double iny)
{
	
	//cout <<"\t"<<subSClon().center_0().degrees()<<" "<<subSClat().north_positive().degrees()<<endl;
	latangle thislat(angle(90.,angle::DEG),NORTH_POSITIVE);
	lonangle thislon(angle(0.,angle::DEG),CENTER_0);
	//subSClon(thislon);
	//subSClat(thislat);
	
	pair<lonangle, latangle> lonlatanswer(thislon,thislat);
	
	double relx,rely;
	double rotx,roty, rho;
	relx=inx-double(xcenter());
	rely=iny-double(ycenter());
	

	rho = sqrtf(relx*relx+rely*rely);
	if (rho < radius()) {
			
		angle lat,lon;
		angle rotangle(-90.,angle::DEG);
		
		//rotx = cos(Nangle())*relx+sin(Nangle())*rely;
		//roty = -sin(Nangle())*relx+cos(Nangle())*rely;
	
		rotx=relx;
		roty=-rely; 
		
				
		// angular distance from subSC point
		calcangdist(rho);
						
		lat = angle_asin(cos(angdist())*sin(subSClat().north_positive())+(roty*sin(angdist())*cos(subSClat().north_positive())/rho));
				
		if (subSClat().north_positive().radians() > 3.1415926535897932385/2.){
			lon = subSClon().center_0() + angle(atan2(rotx, -roty), angle::RAD);
			//cout <<"LON calc 1:\t"<<lon.degrees()<<endl;
		}
		else if (subSClat().north_positive().radians() <-3.1415897932385/2.){
			lon = subSClon().center_0() + angle(atan2(rotx,  roty), angle::RAD);
			//cout <<"LON calc 2:\t"<<lon.degrees()<<endl;
		
		}
		else {
			lon = subSClon().center_0() + angle(atan2(rotx*sin(angdist()), 
					(rho*cos(subSClat().north_positive())*cos(angdist())-roty*sin(subSClat().north_positive())*sin(angdist()))), angle::RAD);
			//cout <<"LON calc 3:\t"<<lon.degrees()<<endl;
		
		}
		
		//while (lon.degrees() < 180.) lon += angle(360.,angle::DEG);
		//while (lon.degrees() > 180.) lon -= angle(360.,angle::DEG);
		
		//cout <<lon.degrees()<<"\n";
				
		// correction for center pixel
		if (inx == xcenter()  &&  iny == ycenter()) {
			lonlatanswer.second = subSClat();
			lonlatanswer.first= subSClon();
		}
		else{
			lonlatanswer.second=latangle(lat,NORTH_POSITIVE);
			lonlatanswer.first=lonangle(lon,CENTER_0);
		}	
	
	} 
	
	else {
		lonlatanswer.second=latangle(angle(1000.,angle::DEG),NORTH_POSITIVE);
		lonlatanswer.first=lonangle(angle(1000.,angle::DEG),CENTER_0);
	}			
	
	return lonlatanswer; 
					
}


vector<pair<double, double> > azimuthalprojection::convert_vector_of_angles(vector<pair<lonangle,latangle> > anglepoints)
{
	
	vector<pair<double, double> > doublepoints;
	
	for (int i(1);i<anglepoints.size();i++) {
		pair<double, double> thisxypair;
		thisxypair=(*this)(anglepoints.at(i).first,anglepoints.at(i).second);
		doublepoints.push_back(thisxypair);					
	}
	
	return doublepoints;
	
}


orthographic::orthographic(lonangle inlon, latangle inlat, int sizex, int sizey, int centerx, int centery, double inrad, angle inNangle)
		: azimuthalprojection(inlon, inlat,sizex,sizey, centerx, centery, inrad, inNangle)
{
}

orthographic::orthographic(cube c)
			:azimuthalprojection(c)
{}

double orthographic::kfactor(lonangle, latangle)
{	

	return (1.*radius());		
}
double orthographic::kfactor(latangle, lonangle)
{

	return (1.*radius());	
	
}

void orthographic::calcangdist(double rho)
{
	angdist(asinf(rho/radius()));

}


azimuthalstereographic::azimuthalstereographic(lonangle inlon, latangle inlat, int sizex, int sizey, int centerx, int centery, double inrad, angle inNangle)
		: azimuthalprojection(inlon, inlat,sizex,sizey, centerx, centery, inrad, inNangle)
{}

azimuthalstereographic::azimuthalstereographic(cube c)
			:azimuthalprojection(c)
{}


void azimuthalstereographic::calcangdist(double rho)
{
	angdist(2*atanf(rho/(2*radius())));

}

double azimuthalstereographic::kfactor(lonangle lon, latangle lat)
{	
	double answer;

	answer = 2*radius() / (1+sin(subSClat().north_positive())*sin(lat.north_positive())+cos(subSClat().north_positive())*cos(lat.north_positive())*cos(lon.center_0()-subSClon().center_0()));
	
	return (answer);		
}

double azimuthalstereographic::kfactor(latangle lat, lonangle lon)
{
	return(kfactor(lon,lat));
	
}


lambertazimuthal::lambertazimuthal(lonangle inlon, latangle inlat,  int sizex, int sizey, int centerx, int centery, double inrad, angle inNangle)
		: azimuthalprojection(inlon, inlat,sizex,sizey, centerx, centery, inrad, inNangle)
{}

lambertazimuthal::lambertazimuthal(cube c)
			:azimuthalprojection(c)
{ cout <<"Creating lambertazimuthal from cube!"<<endl;}

void lambertazimuthal::calcangdist(double rho)
{
	angdist(2*asin(rho/(2*radius())));
	
}

double lambertazimuthal::kfactor(lonangle lon, latangle lat)
{	
	double answer;

	answer = 2/(1+sin(subSClat().north_positive())*sin(lat.north_positive())+cos(subSClat().north_positive())*cos(lat.north_positive())*cos(lon.center_0()-subSClon().center_0()));
	
	return (answer);		
}

double lambertazimuthal::kfactor(latangle lat, lonangle lon)
{
	return(kfactor(lon,lat));
	
}


cylindricalprojection::cylindricalprojection(cube c)
		:projection(lonangle(angle(0.,angle::DEG),CENTER_NEGATIVE180),latangle(angle(0.,angle::DEG),NORTH_POSITIVE))
{
	if (c.keyword("cylindrical_map")=="yes"){
	
		//set lonconvention,latconvention	
		
	}
	else cout <<"-----Non-cylindrical map detected. Cannot create cylindrical projection object.\n";
}

void cylindricalprojection::cyllonconvention(lonconvention l)
{
	_lonconvention=l;	
}

void cylindricalprojection::cyllatconvention(latconvention l)
{
	_latconvention=l;
}


pair<double, double> cylindricalprojection::operator()(lonangle ,latangle )
{
	pair<double, double> answer;
	answer.first=0.;
	answer.second=0.;
	
	// NOT YET WRITTEN
	
	
	return answer;	
	
}


pair<lonangle, latangle> cylindricalprojection::operator()(int , int )
{
	latangle thislat(angle(0,angle::DEG),NORTH_POSITIVE);
	lonangle thislon(angle(0,angle::DEG),CENTER_POSITIVE180);
	pair<lonangle, latangle> lonlatanswer(thislon,thislat);
	
	// NOT YET WRITTEN

	return lonlatanswer; 
					
}



cube cube::azproject(azimuthalprojection* p, methodtype m)
		//azimuthal projection, added 2-2016 by smack
{ 
	cout<<flush;
	cout <<"In the azproject function"<<endl;
	
	int memz(1);
	if (p->xsize()*p->ysize()*N(Z)*4 < maxmemsize) memz=3;
	cube answer(p->xsize(), p->ysize(), N(Z), 0., cubetype, X, Y, Z, memz);
	
	for (int i=0;i<N(Z);i++)
		answer.Axis(Z,i) = Axis(Z,i);
	
	string thistype;
	if (p->isprojtype() == CYL) thistype="CYL";
	else if (p->isprojtype() == ORTHO) thistype="ORTHO";
	else if (p->isprojtype() == AZSTEREOGR) thistype= "AZSTEREOGR";
	else if (p->isprojtype() == LAMBAZ) thistype="LAMBAZ"; 
					
	answer.copyheaderfrom(*this);
	double hdrsubSClon(0.); 
	
	if (keyword("lonconvention")== " CENTER_NEGATIVE180") hdrsubSClon=p->subSClon().center_negative180().degrees();
	else if (keyword("lonconvention")== " CENTER_POSITIVE180") hdrsubSClon=p->subSClon().center_positive180().degrees();
	else if (keyword("lonconvention")== " CENTER_0") hdrsubSClon=p->subSClon().center_0().degrees();
	
	answer.addazprojectionhdr(thistype,
			p->subSClat().north_negative().degrees(),hdrsubSClon,			
			p->xcenter(),p->ycenter(),p->radius(),
			-1*p->Nangle().degrees(),
			p->xsize(),p->ysize());
	
	
	
	if (!osuppress) cout << "Azimuthal projecting --  00%";
	for (int x=0;x<p->xsize();x++) {
		if (!osuppress) printpercent(x, p->xsize()-1);
		for (int y=0;y<p->ysize();y++) {
			
			latangle thislat(angle(0,angle::DEG),NORTH_POSITIVE);
			lonangle thislon(angle(0,angle::DEG),CENTER_NEGATIVE180);
	
			pair<lonangle, latangle> thislonlat(thislon,thislat);
	
			thislonlat=(*p)(x,y);
			if (thislonlat.first.number() != fabs(1000.) & thislonlat.second.number() != fabs(1000.)){
				for (int z=0;z<N(Z);z++) {
					if ((*this).Axis(X,0)<-180){ 
						answer(x,y,z) = (*this)(thislonlat.first.center_negative180().degrees(), thislonlat.second.north_negative().degrees(), z, m);
					}
					else if ((*this).Axis(X,N(X)-1)>0) {
						answer(x,y,z) = (*this)(thislonlat.first.center_0().degrees(), thislonlat.second.north_negative().degrees(), z, m);
					}			
				}
			}
				
		}
		
	}
	
	
	return answer;
}


