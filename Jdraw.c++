#include <math.h>

#include "Jcrap.h"
#include "Jprojections.h"

cube cube::draw_circle(int x,int y,float radius, float r, float g, float b)
{
	cube answer(*this);
	for (int cx=(int)(x-radius); cx<=(int)(x+radius) && cx<N(X); cx++) {
		if (cx<0) cx=0;
		double cy2;
		int dcy;
		cy2= radius*radius-(double)((x-cx)*(x-cx));
		if (cy2>0) dcy = (int)sqrt(cy2);
		else dcy=0;

		if (y+dcy <= N(Y)-1.  &&  y+dcy >= 0.) {
			answer(cx,y+dcy,0) = r;
			if (answer.N(Z) > 1) {
				answer(cx,y+dcy,1) = g;
				answer(cx,y+dcy,2) = b;
			}
		}
		if (y-dcy <= N(Y)-1.  &&  y-dcy >= 0.) {
			answer(cx,y-dcy,0) = r;
			if (answer.N(Z) > 1) {
				answer(cx,y-dcy,1) = g;
				answer(cx,y-dcy,2) = b;
			}
		}
	}	
	for (int cy=(int)(y-radius); cy<=(int)(y+radius) && cy<N(Y); cy++) {
		if (cy < 0) cy=0;
		double cx2;
		int dcx;
		cx2= radius*radius-(double)((y-cy)*(y-cy));
		if (cx2>0) dcx = (int)sqrt(cx2);
		else dcx=0;
		if (x+dcx <= N(X)-1.  &&  x+dcx >= 0.) {
			answer(x+dcx,cy,0) = r;
			if (answer.N(Z) > 1) {
				answer(x+dcx,cy,1) = g;
				answer(x+dcx,cy,2) = b;
			}
		}
		if (x-dcx <= N(X)-1.  &&  x-dcx >= 0.) {
			answer(x-dcx,cy,0) = r;
			if (answer.N(Z) > 1) {
				answer(x-dcx,cy,1) = g;
				answer(x-dcx,cy,2) = b;
			}
		}
	}
	return answer;
}

cube cube::draw_line(int x1,int y1, int x2,int y2, float r, float g, float b)
{
	cube answer(*this);
	
	if (x1>x2) {
		int tmp(x1);
		x1=x2;
		x2=tmp;
		
		tmp=y1;
		y1=y2;
		y2=tmp;
	}
	double slope(double(y2-y1)/double(x2-x1));
//	cout << "xslope=" << slope << "\n";
	for (int x=x1;x!=x2;x++) {
		int y = int(y1+float(x-x1)*slope);
		if (x>=0 && x<N(X) && y>=0 && y<N(Y)) {
			answer(x,y,0) = r;
			if (answer.N(Z) > 1) {
				answer(x,y,1) = g;
				answer(x,y,2) = b;
			}
		}
	}	
	
	if (y1>y2) {
		int tmp(x1);
		x1=x2;
		x2=tmp;
		
		tmp=y1;
		y1=y2;
		y2=tmp;
	}	
//	slope = double((x2-x1)/(y2-y1));
	slope = 1./slope;
//	cout << "yslope=" << slope << "\n";
	for (int y=y1;y!=y2;y++) {
		int x = int(x1+float(y-y1)*slope);
		if (x>=0 && x<N(X) && y>=0 && y<N(Y)) {
			answer(x,y,0) = r;
			if (answer.N(Z) > 1) {
				answer(x,y,1) = g;
				answer(x,y,2) = b;
			}
		}
	}
	
	return answer;
}

cube cube::draw_ellipse(int x,int y,float a /*semimajor*/, float b/*semiminor*/, 
		float orientation, float red, float grn, float blu, bool filled)
{
	cube answer(*this);
	float &t(orientation);
	
	for (int cy=(int)(-a); cy<=(int)(+a); cy++) {
		if (y+cy < 0) cy=-y;
		if (y+cy > N(Y)-1) continue;

		double cxp, cxm, det;
		det = a*a*sin(t)*sin(t)+b*b*cos(t)*cos(t)-cy*cy;
		cxp=( -cy*sin(t)*cos(t)*(a*a-b*b) + a*b*sqrt(det) )
				                   /
				(a*a*sin(t)*sin(t)+b*b*cos(t)*cos(t));
		cxm=( -cy*sin(t)*cos(t)*(a*a-b*b) - a*b*sqrt(det) )
				                   /
				(a*a*sin(t)*sin(t)+b*b*cos(t)*cos(t));
//		cout << "y=" << y+cy << ", cxp=" << cxp << ", cxm=" << cxm;
//		cout << ", det = " << det << "\n";
		
		if (filled) for (int cx(cxm); cx<=cxp; cx++) {
			if (x+cx<0) cx=-x;
			if (x+cx>N(X)-1) continue;
			
			answer(x+cx,y+cy,0) = red;
			if (answer.N(Z)>=2) answer(x+cx,y+cy,1) = grn;
			if (answer.N(Z)>=3) answer(x+cx,y+cy,2) = blu;
		} else {
			if (x+cxp>=0 && x+cxp<N(X)) answer(x+cxp,y+cy,0) = red;
			if (x+cxm>=0 && x+cxp<N(X)) answer(x+cxm,y+cy,0) = red;
		}
	}
	return answer;
	
}

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define INSIDE 1
#define OUTSIDE 0

bool cube::is_inside_polygon(vector<pair<double, double> > polygon, pair<double, double> p)
{
  int counter = 0;
  int i, N(polygon.size());
  double xinters;
  pair<double, double> p1,p2;

  p1 = polygon[0];
  for (i=1;i<=N;i++) {
    p2 = polygon[i % N];
    if (p.second > MIN(p1.second,p2.second)) {
      if (p.second <= MAX(p1.second,p2.second)) {
        if (p.first <= MAX(p1.first,p2.first)) {
          if (p1.second != p2.second) {
            xinters = (p.second-p1.second)*(p2.first-p1.first)/(p2.second-p1.second)+p1.first;
            if (p1.first == p2.first || p.first <= xinters)
              counter++;
          }
        }
      }
    }
    p1 = p2;
  }

  if (counter % 2 == 0)
    return(OUTSIDE);
  else
    return(INSIDE);
}

cube cube::draw_filled_polygon(vector<pair<lonangle,latangle> > points, float r, float g, float b)
{
	cube answer(*this);
	vector<pair<double, double> > points_dbl;
	if (keyword("cylindrical_map")==" yes"){
		
		for (int i(0); i<points.size(); i++){			
			pair<double, double> thispair;
			if (keyword("lonconvention") == " CENTER_NEGATIVE180") 
				thispair.first=points.at(i).first.center_negative180().degrees();
			else if (keyword("lonconvention") == " CENTER_0")
				thispair.first=points.at(i).first.center_0().degrees();		
			thispair.second=points.at(i).second.north_negative().degrees();
			points_dbl.push_back(thispair);
		}
	}
	
	else if (keyword("cylindrical_map")!=" yes"){
		
		if (keyword("azimuthal_map")==" ORTHO"){
			orthographic thisortho(*this);
			points_dbl=thisortho.convert_vector_of_angles(points);
			
		}
		if (keyword("azimuthal_map")==" AZSTEREOGR"){
			azimuthalstereographic thisazstereogr(*this);
			points_dbl=thisazstereogr.convert_vector_of_angles(points);
		}
		if (keyword("azimuthal_map")==" LAMBAZ"){
			lambertazimuthal thislambaz(*this);
			points_dbl=thislambaz.convert_vector_of_angles(points);
		
		}
		
	
	}
	
	answer.draw_filled_polygon_on(points_dbl,r,g,b);
	
}

cube cube::draw_filled_polygon(vector<pair<double, double> > points, float r, float g, float b)
{
	cube answer(*this);
	
	answer.draw_filled_polygon_on(points,r,g,b);
	
	return answer;
}

void cube::draw_filled_polygon_on(vector<pair<double, double> > points, float r, float g, float b)
{	
	double maxx(points.at(0).first), minx(points.at(0).first);
	double maxy(points.at(0).second),miny(points.at(0).second);
	for (int i(1);i<points.size();i++) {
		if (points.at(i).first > maxx) maxx = points.at(i).first;
		if (points.at(i).first < minx) minx = points.at(i).first;
		if (points.at(i).second> maxy) maxy = points.at(i).second;
		if (points.at(i).second< miny) miny = points.at(i).second;
	}
	cout << "Minx=" << minx << ", Maxx=" << maxx << ", Miny=" << miny << ", Maxy=" << maxy << "\n";
	
	for (int x(0);x<N(X);x++) {
		if (Axis(X,x) >= minx  &&  Axis(X,x) <= maxx) {
			for (int y(0);y<N(Y);y++) {
				if (Axis(Y,y) >= miny  &&  Axis(Y,y) <= maxy) {
					pair<double, double> testpoint;
					testpoint.first=Axis(X,x);
					testpoint.second=Axis(Y,y);
					if (is_inside_polygon(points,testpoint)) {
						(*this)(x,y,0) = r;
						if (N(Z) >= 2) (*this)(x,y,1) = g;
						if (N(Z) >= 3) (*this)(x,y,2) = b;
					}
				} 
			}
		}
	}	
}

cube cube::autocrop(int xmargin, int ymargin, int zmargin)
{
	osuppress++;
	int xmin, xmax, ymin, ymax, zmin, zmax;
	float slicesum(0.);
	for (xmin=0;slicesum==0.;xmin++) { 
		if (xmin == N(X)) {
			cout << "Autocropped empty cube to zero\n"; cout.flush();
			return cube(1,1,1);
		}
		slicesum = (*this).plane(X,xmin).sum();
	}
	slicesum=0.;
	for (xmax=N(X)-1;slicesum == 0.;xmax--)
		slicesum = (*this).plane(X,xmax).sum();
	slicesum=0.;
	for (ymin=0;slicesum == 0.;ymin++)
		slicesum = (*this).plane(Y,ymin).sum();
	slicesum=0.;
	for (ymax=N(Y)-1;slicesum == 0.;ymax--)
		slicesum = (*this).plane(Y,ymax).sum();
	slicesum=0.;
	for (zmin=0;slicesum == 0.;zmin++)
		slicesum = (*this).plane(Z,zmin).sum();
	slicesum=0.;
	for (zmax=(N(Z)-1);slicesum == 0.;zmax--)
		slicesum = (*this).plane(Z,zmax).sum();
	
	xmin--; xmax++; ymin--; ymax++; zmin--; zmax++;
	
	xmin -= xmargin;
	if (xmin < 0) xmin = 0;
	xmax += xmargin;
	if (xmax > N(X)-1) xmax = N(X)-1;
	ymin -= ymargin;
	if (ymin < 0) ymin = 0;
	ymax += ymargin;
	if (ymax > N(Y)-1) ymax = N(Y)-1;
	zmin -= zmargin;
	if (zmin < 0) zmin = 0;
	zmax += zmargin;
	if (zmax > N(Z)-1) zmax = N(Z)-1;
	
	if (osuppress<1) {
		cout << "Cropping " << xmin << " " << xmax << "   " << ymin << " ";
		cout << ymax << "   " << zmin << " " << zmax << "\n";
	}
	
	osuppress--;
	
	return (*this)(xmin,xmax,ymin,ymax,zmin,zmax);
}

double Jcrap::polygonarea(vector<pair<double, double> > polygon)
{
	double answer(0.);
	int j(polygon.size()-1);
	
	for (int i(0);i<polygon.size();i++) {
	// from http://www.mathopenref.com/coordpolygonarea2.html
		answer += (polygon.at(j).first+polygon.at(i).first) * 
				(polygon.at(j).second-polygon.at(i).second);
		j = i; // j is now previous vertex
	}
	
	return fabs(answer/2.);
}
