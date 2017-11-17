#include "photontraverse.h"

segment::segment(double instart, atmozone* inzone, photontraverse* intraverse) :
	_start(instart), _tozone(inzone), _totraverse(intraverse)
{
}

double segment::Start()
{
	return _start;
}

atmozone* segment::Zone()
{
	return _tozone;
}

photontraverse* segment::Traverse()
{
	return _totraverse;
}

double segment::the_function(double d)
{
	double answer(0.);
	
	double height_km(Traverse()->Photon().project_km(d).r());
//	cout << "In segment::the_function for distance " << d << " km, now";
//	cout << " looking up the dTau/dz at height " << height_km << " km.  ";
	
	answer = Zone()->extinction_haze_km_um(height_km, Traverse()->Photon().lambda_um());
//	cout << "got " << answer << "\n"; 
	
	return answer;
}


photontraverse::photontraverse(photon& inphoton, atmosphere& inatmosphere) :
		_thisphoton(inphoton)
{
	
// Step 4a:  calculate distance to shell-crossing points
	// Step 4a(i):  calculate impact parameter, i.e. minimum radius r_min
					//	 from d_min, the distance to the minimum radius from the
					//  starting point
	
	// Equation 7 from Young & Trouille unpublished SRTC paper
	double d_min(inphoton.position().x()*inphoton.direction().lat().cos()*inphoton.direction().lon().cos() +
					 inphoton.position().y()*inphoton.direction().lat().cos()*inphoton.direction().lon().sin() +
					 inphoton.position().z()*inphoton.direction().lat().sin());
	d_min *= -1.; // not sure what's wrong with Eliot's equation here . . .
	
	// Equation 6 from Y&T to get x_min, y_min, z_min
	geomvector minpoint(inphoton.project_km(d_min));

	// Equation 8 from Y&T to get r_min
	double r_min(minpoint.r());		

	// Step 4a(ii):  determine crossing points, assign segments
	
	// misses planet entirely
	if (osuppress<-1)cout << "Atmosphere outer edge is " << inatmosphere.outeredge_km() << "\n";
	if (r_min > inatmosphere.outeredge_km()) {
		_segments.push_back(segment(0., &atmozone::Space(), this));  
		if (osuppress<-1)cout << "Nice: this photon totally misses the atmosphere.  Moving on!\n";
	} else {
		
	// inbound segments
		for (atmosphere::reverse_iterator ri(inatmosphere.rbegin());ri!=inatmosphere.rend();ri++) {
			double r_outeredge(ri->top_km());
			if (r_outeredge < r_min) {  // chord doesn't drop down this far
				//double-up on the bottom segment then quit and move on to egress
				_segments.push_back(segment(d_min, (--ri)->front(), this));  // starting with one zone/layer
				break;
			}
			double d_outeredge(d_min - sqrt(r_outeredge*r_outeredge-r_min*r_min));
			segment thisseg(d_outeredge, (ri->front()), this);
			_segments.push_back(thisseg);
		}
		
	// outbound segments
		if (!_segments.back().Zone()->issurface()  ||  d_min<0.) {  // if last segment ends in atmosphere or surface is behind us
			for (atmosphere::iterator i(++inatmosphere.begin());i!=inatmosphere.end();i++) {
				double r_inneredge(i->bottom_km());
				if (r_min > r_inneredge) continue;  // doesn't get down this far, or this is the re-layer.

				double d_inneredge(d_min + sqrt(r_inneredge*r_inneredge-r_min*r_min));
				_segments.push_back(segment(d_inneredge, (i->front()), this));
			}
			double r_outeredge(inatmosphere.outeredge_km());
			double d_outeredge(d_min + sqrt(r_outeredge*r_outeredge-r_min*r_min));
			
			segment thisseg(d_outeredge, &atmozone::Space(), this);
			_segments.push_back(thisseg);
			
		}
	}
}


double photontraverse::the_function(double d)
{
	double answer(0.);
	
	cube traverseplot(1,1,1);
	
	list<segment>::iterator i(Segments().begin());
	if (d<0.) { 
		if (osuppress<-4) cout << "OMG why do you hate me by giving me a negative d value (" << d << ") in photontraverse::the_function!?\n";
		if (osuppress<-4) cout << "This segment is in a layer named " << i->Zone()->layername() << "\n";
	} else {
		for (list<segment>::iterator i=Segments().begin();i!=Segments().end();i++) {
			list<segment>::iterator next(i);
			next++;
			if (next!=Segments().end() && d>i->Start()) {
				if (next->Start() > 0.) {  // if next->start is less than zero, just ignore this segment entirely					
					// setting d0
					double d0(0.);  // if this segment's start is LESS than 0, then start at 0.
					if (i->Start() > 0.) d0 = i->Start();  // otherwise start at Start()
					
					// setting d1
					double d1(d);
					if (d > next->Start()) d1 = next->Start();
					
					answer += (*i).integrate(d0, d1);
				}
			} 
		}
	}
	
	return answer;
}

photon photontraverse::Photon()
{
	return _thisphoton;
}

photontraverse_extinction::photontraverse_extinction(photon& _p, atmosphere& _a) : photontraverse(_p, _a)
{
}

double photontraverse_extinction::the_function(double d)
{
	double answer(-1.);
	
	for (list<segment>::iterator i=Segments().begin();i!=Segments().end();i++) {
		list<segment>::iterator next(i);
		next++;
		if (i->Start() <= d  &&  next->Start() > d) answer = i->the_function(d);
	}
	
	return answer;
}
