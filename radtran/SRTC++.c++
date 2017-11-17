#include <omp.h>
#include <stdlib.h>
#include "SRTC++.h"

int SRTC::debug(0);
int SRTC::ompdebug(0);

SRTC::SRTC(atmosphere inatmo, photongenerator* inphotons, vector<detector*> indetectors) :
	_atmosphere(inatmo), 
	_tophotongenerator(inphotons),
	_detectorlist(indetectors)
{
	_iterationlimit = 10000;
	numericalprecision = 1.e-6;
	
	for (atmosphere::iterator i(_atmosphere.begin());i!=_atmosphere.end();i++) {
		for (atmolayer::iterator j(i->begin());j!=i->end();j++) {
			(**j).zonephasefunction()->check_that_phase_function_integrates(_tophotongenerator->wavelengths());
			for (int w(0);w<_tophotongenerator->wavelengths().size();w++) {
				(**j).zonephasefunction()->precalculate_CDF(double(_tophotongenerator->wavelengths().at(w).convert("um")));
			}
		}
	}
}

void SRTC::run(int numthreads)
{ 
	if (numthreads > 0) omp_set_num_threads(numthreads);  // otherwise just default to number of cores on the machine	

	long long done(0);
//	osuppress = -2;
	
#pragma omp parallel for schedule(dynamic)
	for (long long n=0;n<tophotongenerator()->size();n++) {	
		if (osuppress<1 && n==0) cout << "Planning to run " << tophotongenerator()->size() << " photons on " << omp_get_num_threads() << " processors\n";
		if (osuppress<1 && n==0) cout << "SRTC Progress:  00%";
		if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting new photon\n"; cout.flush();
		if (osuppress<0) cout << "thread " << omp_get_thread_num() << "/" << omp_get_num_threads();
		if (osuppress<0) cout << " -- Starting photon number " << n << " ...\n";
		photon thisphoton(tophotongenerator()->generatephoton(n));
		
		if (osuppress<-1) cout << "thisphoton.position() = " << thisphoton.position() << "\n";
		if (osuppress<-1) cout << "thisphoton.direction() = " << thisphoton.direction() << "\n";
		
		
		list<photon> photon_history;
		
		
		// iterate through steps 3-6 in Eliot's SRTC paper for each photon
		while (thisphoton.amplitude()) {
			photon_history.push_back(thisphoton);
			
			if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting Step 3\n"; cout.flush();
		// Step 3:  Generate a random Tau excursion for each remaining photon
			double go_to_tau(photon::generate_random_tau());
							
			if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting Step 4\n"; cout.flush();
		// Step 4:  Translate that Tau excursion into meatspace						
			photontraverse trek(thisphoton, Atmosphere());
			if (osuppress<-2) cout << "Testing.  " << trek.Segments().size() << " segments produced.\n";
			
			unsigned int n(0);
		
			
			if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting Step 4 case 0\n"; cout.flush();
			// Case 0:  no intersection
			atmozone *to_scattering_zone(0);
			photon photon_at_scatter(thisphoton);  // we will need these on the next step
			photon scatteredphoton(thisphoton);
			
			if (trek.Segments().front().Zone() == &atmozone::Space()) {
				photon_at_scatter.amplitude(0.);
				scatteredphoton.amplitude(0.);
			} else {
				double totaltau(trek.the_function(2.*_atmosphere.outeredge_km())); //*2 to make sure we catch all the atmosphere 
				// 3 cases:
				// case 1:  sails off into space never to be heard from again
				if (go_to_tau > totaltau  &&
						trek.Segments().back().Zone()==&(atmozone::Space())) {
					photon_at_scatter.amplitude(0.);
					scatteredphoton.amplitude(0.);
					if (osuppress<-1)cout << "This photon sails off into space because gototau of ";
					if (osuppress<-1)cout << go_to_tau << " is greater than totaltau of " << totaltau;
					if (osuppress<-1)cout << "\n";
				} else {						
					if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting Step 4 interesting case.\n"; cout.flush();
					// calculating the location of the scatter, needed for either scatter type
					double scatter_d_km(0.);
					atmolayer *to_scatterlayer(0);
					
					if (go_to_tau > totaltau) {
				// case 2:  surface scatter
						
						if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " step 4 surface scatter\n"; cout.flush();		
						if (osuppress<-1)cout << "This photon crashes into the surface because gototau of ";
						if (osuppress<-1)cout << go_to_tau << " is greater than totaltau of " << totaltau;
						if (osuppress<-1)cout << "\n";
						
						scatter_d_km = trek.Segments().back().Start();
						photon_at_scatter.position(thisphoton.project_km(scatter_d_km));
						
						to_scatterlayer = &(Atmosphere().front());
					} else {
				// case 3:  atmospheric scatter
						if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " step 4 atmosphere scatter\n"; cout.flush();
						if (osuppress<-1)cout << "This photon scatters in the atmosphere because gototau of ";
						if (osuppress<-1)cout << go_to_tau << " is less than totaltau of " << totaltau;
						if (osuppress<-1)cout << "\n";
						
						scatter_d_km = trek.findzero(0., 2.*_atmosphere.outeredge_km(), numericalprecision*1.e3, go_to_tau, 
								string("atmospheric scatter findzero"));
						if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " step 4 atmosphere scatter after findzero\n"; cout.flush();
						
						
						photon_at_scatter.position(thisphoton.project_km(scatter_d_km));
						if (osuppress<-3) cout << "photon history size is " << photon_history.size() << ".  x is " << photon_at_scatter.position().x() << "\n";
 						
						for (atmosphere::iterator i(Atmosphere().begin()++);
								i!=Atmosphere().end();i++) {
							if (photon_at_scatter.position().r()>=i->bottom_km()  &&
								 photon_at_scatter.position().r()<=i->top_km()) {
								to_scatterlayer = &(*i);
								break;
							}
						}
						
						if (to_scatterlayer==0) {
							if (osuppress<-2) {
								cout << "Uhoh -- can't find an appropriate layer in which to scatter!\n";					
								cube traverseplot(trek.plot(0., 4000., 300));
								traverseplot.keyword("title", string("Traverse the_function"));
								traverseplot.graph();
								cout << "Looking for Tau = " << go_to_tau << " with total tau of " << totaltau << "\n";
							}							
							if (photon_at_scatter.position().r()<=Atmosphere().front().top_km())
								to_scatterlayer = &(Atmosphere().front());
							else
								to_scatterlayer = &(Atmosphere().back());
						}
					}
					if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " end of Step 4 interesting case.\n"; cout.flush();
					
					double photonlat_deg(photon_at_scatter.position().lat().degrees());
					double photonlon_deg(photon_at_scatter.position().lon().degrees());
						
					for (atmolayer::reverse_iterator i(to_scatterlayer->rbegin());
								i!=to_scatterlayer->rend(); i++) {
						if (photonlon_deg>=(**i).left_lon().degrees()  &&
							 photonlon_deg<=(**i).right_lon().degrees() &&
							 photonlat_deg>=(**i).top_lat().degrees()   &&
							 photonlat_deg<=(**i).bottom_lat().degrees()) {
							to_scattering_zone = *i;
							break;
						}
					}
				}
				
				if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting Step 5\n"; cout.flush();
			// Step 5:  Execute the scattering
				if (photon_at_scatter.amplitude() && to_scattering_zone) {
					if (osuppress<-1)cout << "Now executing a scattering for this photon.\n";
					if (osuppress<-1)cout << " to_scattering_zone = " << to_scattering_zone << "\n";
					if (osuppress<-1) if (to_scattering_zone)
						              cout << "  Scattering zone is in layer " << to_scattering_zone->layername() << "\n";
					if (osuppress<-1)cout << "Photon scattering at location " << photon_at_scatter.position() << "\n";	

					
			
					
					
					scatteredphoton = photon_at_scatter; // important -- photon::scatter() modifies the calling photon's direction *and* amplitude
					scatteredphoton = scatteredphoton.scatter(to_scattering_zone);   // this scatter updates the amplitude based on the single-scattering albedo of the haze in this zone
					photon_at_scatter.amplitude(scatteredphoton.amplitude());  // hence why here we update rather than multiply by the SSA twice
					photon_at_scatter.set_scattering_zone(scatteredphoton.scattering_zone());					
					
					if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " starting Step 5a\n"; cout.flush();
				// Step 5a:  Update the virtual CCD based on the scatter
					if (photon_at_scatter.amplitude()) {
						for (vector<detector*>::iterator i(Detectors().begin());i<Detectors().end();i++) {
							if (osuppress < -3) cout << "\n\nDetector " << (*i)->name() << "\n";
							angle phase_to_detector, incidence, emission;
							double phase_function_to_detector;
							
							geomvector from_photon_to_detector((*i)->position()-photon_at_scatter.position());
							phase_to_detector = angle_acos( photon_at_scatter.direction().unitvector()
									.dot(from_photon_to_detector.unitvector()) );
							
							if (osuppress < -3) cout << "position vector: " << photon_at_scatter.position().as(POLAR) << "\n";
							if (osuppress < -3) cout << "direction vector: " <<photon_at_scatter.direction().as(POLAR) << "\n";
							if (osuppress < -3) cout << "direction vector: " <<photon_at_scatter.direction().as(CARTESIAN) << "\n";
							if (osuppress < -3) cout << "to_detector vector: " << from_photon_to_detector.as(POLAR) << "\n";
							if (osuppress < -3) cout << "to_detector vector: " << from_photon_to_detector.as(CARTESIAN) << "\n";
							
							if (to_scattering_zone->issurface()) {
								
								
								phase_to_detector = angle(180., angle::DEG)-phase_to_detector;  // surface phase means something different from atmospheric phase
																	// for surface, 0 phase is backscatter.
																	// for atmosphere, 0 phase is forwardscatter.
								if (osuppress < -3) cout << "Calculated surface detector scattering phase to be " << phase_to_detector.as(angle::DEG) << "\n";
								// surface scatter to detector
								geomvector incident(photon_at_scatter.direction().transform_to(photon_at_scatter.position()));
								incidence = angle(90.,angle::DEG)-angle_asin(-(incident.unitvector().z()));
								if (osuppress < -3) cout << "Calculated surface detector scattering incidence to be " << incidence.as(angle::DEG) << "\n";
								
								geomvector emitted(from_photon_to_detector.unitvector().transform_to(photon_at_scatter.position()));
								emission = angle(90.,angle::DEG)-angle_asin(emitted.z());
								if (osuppress < -3) cout << "Calculated surface detector scattering emission to be " << emission.as(angle::DEG) << "\n";
							} else { 
								// atmospheric scatter to detector
								incidence = emission = angle(0.,angle::RAD);
							}
							
							if (emission.as(angle::DEG).number() > 90.)
								phase_function_to_detector = 0.;
							else{
								phase_function_to_detector = to_scattering_zone->zonephasefunction()->normalized_phasefunction(
										photon_at_scatter.lambda_um(), phase_to_detector, incidence, emission);
							}
							
							photon scattered_to_detector(photon_at_scatter);
							scattered_to_detector.direction(from_photon_to_detector);
							scattered_to_detector.amplitude(photon_at_scatter.amplitude()*phase_function_to_detector);

							photontraverse outbound_trek(scattered_to_detector, Atmosphere());
							// if there's a surface between here and detector, zero out the amplitude
							for (list<segment>::iterator s(outbound_trek.Segments().begin());s!=outbound_trek.Segments().end();s++)	{	
								if (s->Zone()->issurface() && s->Start()>0.000001) {
									scattered_to_detector.amplitude(0.);
								}
							}						
							
							if (scattered_to_detector.amplitude()) {   // only keep bothering if the amplitude to be detected is non-zero
								double outbound_tau(outbound_trek(outbound_trek.Segments().back().Start()));
								if (osuppress<-2) cout << "Outbound tau calculated to be " << outbound_tau << "\n";
								
								if (osuppress<-2) cout << "Photon amplitude prior to extinction on the way out = " << scattered_to_detector.amplitude() << "\n";
								scattered_to_detector.amplitude(scattered_to_detector.amplitude()*exp(-outbound_tau));
								if (osuppress<-2) cout << "extincting by multiplying by " << exp(-outbound_tau) << " on the way out\n";
								if (osuppress<-2) cout << "Amplitude sent to detector then is " << scattered_to_detector.amplitude() << "\n";
								photon_history.push_back(scattered_to_detector);
								(*i)->detect(photon_history);
								photon_history.pop_back();
 
							}
							if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " completed step 5a\n"; cout.flush();
						}
					}				
				}
			}

			thisphoton = scatteredphoton;
		}
		
		
		if (osuppress<-1)cout << "Completed photon " << n << "/" << tophotongenerator()->size();
		if (osuppress<-1)cout << "!\n"; if (osuppress<-1)cout.flush();
		if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " completed photon, now to update total\n"; cout.flush();
		
#pragma omp atomic
		++done;
#pragma omp critical
 		if (osuppress<1) printpercent(done,tophotongenerator()->size());
		if (SRTC::ompdebug>0)  cout << "Thread " << omp_get_thread_num() << " completed totals updating.  Ending loop.\n"; cout.flush();
	
	}
	if (osuppress<1) printpercent(10,10);
	if (osuppress<1) cout << "\n";
	
// calibrate detectors
	for (int i(0);i<_detectorlist.size();i++) {
		(*_detectorlist.at(i)) /= (_tophotongenerator->photons_per_square_km() * (_detectorlist.at(i))->pixelarea_in_square_km())/(1.*Jcrap::pi); // following efy eq 11
		(*_detectorlist.at(i)).keyword("photons_per_square_km", double2str(_tophotongenerator->photons_per_square_km(),15));
	}
	
}

string SRTC::datapath()
{
	string answer("./data/");

	if (getenv("SRTCPATH")) {
		char path_cstr[1000];
		strcpy(path_cstr, getenv("SRTCPATH"));
		answer = string(path_cstr);
	}
	
	return answer;
}
