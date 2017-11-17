#include "SRTC++_testsuite.h"
#include "diffusive_reflectance.h"

cube avgcenterval(cube c)
{
	cube answer(1,1,c.N(Z));
	int searchwidth(0);
	int left((c.N(X)+1)/2 -searchwidth),right((c.N(X)+1)/2 +searchwidth),top((c.N(Y)+1)/2-searchwidth),bottom((c.N(Y)+1)/2+searchwidth);
	
	for (int z(0); z<c.N(Z); z++){	
		int count(0);
	
		for (int x(left); x<=(right); x++){
			for (int y(top); y<=(bottom);y++){
				answer(0,0,z)+=c(x,y,z);
				count++;
			
			}
		}
		
		answer(0,0,z)/=count;	
	}
	
	
	return answer;
}	
	
extern int graphno;

double test_DiffusiveReflectance(double surfalbedo) 
{
	int steps(101);
	double ssa_0(0.01),ssa_n(0.99);
	vector<double> taus;
	taus.push_back(0.1);
	taus.push_back(1.0);
	//taus.push_back(3.0);
	
	vector<cube> analyticalresults;
	vector<cube> calculatedresults;
	
	for (int t(0); t<taus.size();t++){
			// // ------------ analytical results
		diffusive_reflectance s(taus.at(t),surfalbedo);
		cube analyticalresult(s.plot(ssa_0,ssa_n,10));
		analyticalresult.graph("",outputtype);
		analyticalresults.push_back(analyticalresult);
	}
		
	for (int t(0); t<taus.size();t++){
		// // ------------ calculated results
	
		atmosphere a(atmosphere::orangerind(taus.at(t), value(100., "km"), surfalbedo));
	
		// ssalbedo spread to match the analytical
		cube ssalbedos(1,1,steps);
		for (int i(0);i<ssalbedos.N(Z);i++) {
			ssalbedos(0,0,i) = double(i)/double(ssalbedos.N(Z)-1);
			ssalbedos.Axis(Z,i) = ssalbedos(0,0,i);
		}
		a.back().at(0)->ss_albedo(ssalbedos);  
		
			// surface	
		cube albedoonevalue(1,1,1);
		albedoonevalue(0,0,0)=surfalbedo; 
		albedoonevalue.Axis(Y,0)=surfalbedo;
		atmozone *uniformsurface=new surface_atmozone(albedoonevalue);
 		a.front().addzone(uniformsurface);
		
			// sun, have to put in a wavelength, tho doesn't really matter in this case		
		vector<value> wavelengths;	
		for (int i(0);i<ssalbedos.N(Z);i++) wavelengths.push_back(value(ssalbedos(0,0,i), "um")); 	
		
		photongenerator_square hv(value(4200., "km"), 100., value(4200., "km"), wavelengths);	
		//geomvector lampposition(-2675.,0.,0.); // 100 km above the surface
		//geomvector lampdirection(1.,0.,0.); // lamp points in +x direction
		//photongenerator_diffuselamp hv(100000., wavelengths, lampposition, lampdirection);
			
		// detectors
		vector<detector*> Detectors;
		double x(-1.e9); 
		double y(0.);
		colorCCD *CCD_detector;
		CCD_detector = new colorCCD(geomvector(x,y,0.),value(3200., "km"),value(100., "km"),wavelengths);
		string name("tau_");
		name+=double2str(taus.at(t));
		CCD_detector->name(name);
		Detectors.push_back(CCD_detector);
	 
		// run SRTC++
		SRTC S(a, &hv, Detectors);
		S.run();
		
		Detectors.at(0)->CCD().write("DiffusiveReflectance.Jcube");
		cube calculatedresult(avgcenterval(Detectors.at(0)->CCD()));
		for (int z(0); z<calculatedresult.N(Z); z++){
			calculatedresult.Axis(Z,z)=Detectors.at(0)->CCD().Axis(Z,z);
		}
		calculatedresult.graph("",outputtype);
		calculatedresults.push_back(calculatedresult);
	
	}
	

	
	// // ------------  Analyze results
	
		// nice plot
	
	string command(string("graph -T "));
	command+=asstring(outputtype);
	command+=string(" -X \"Single-Scattering Albedo\" -Y \"I/F\" -L \"SRTC++ Comparison to Diffusive Reflectance (Tau=");
	for (int i(0);i<taus.size();i++) {
		command+=int2str(int(taus.at(i)));
		if (i<taus.size()-1) command+=string(",");
	}
	
	command+=string(")\" ");
	command+=string("--pen-colors \"1=red:2=purple:3=blue:4=green:5=brown\" "); 
	
	for (int i(0);i<taus.size();i++) {
		command+=string("-m 4 graph");	
		command+=int2str(graphno-2*int(taus.size())+i);
		command+=string(".out ");
	}	
	
	command+=string(" -S 3 -C -m 0 ");
	for (int i(0);i<taus.size();i++) {
		command+=string(" graph");	
		command+=int2str(graphno-int(taus.size())+i);
		command+=string(".out ");
	}
	
	command+=string(" > DiffusiveReflectance2Layer.");
	command+=asstring(outputtype);
	system(command.c_str());
	cout << command << "\n";
	
	
	
		// evaluate goodness of fit
	
	cube diff(calculatedresults.at(0));
	
	cout <<"SRTC++ \t Hapke"<<endl;
	string outname("diffusivereflectanceresults_"+double2str(surfalbedo)+".txt");
	ofstream output(outname.c_str());
	output<< "ssa\tsrtc\tdf\n";
	for (int p(0);p<diff.N(X);p++){
		if (p!=0) diff=diff.blocks(X,calculatedresults.at(p));
		for (int i(0);i<diff.N(Z);i++){
			output << diff.Axis(Z,i)<<"\t"<<diff(p,0,i) << "\t" << analyticalresults.at(p)(0,0,diff.Axis(Z,i))<<endl; 
			cout << "dividing " << diff(p,0,i) << " by " << analyticalresults.at(p)(0,0,diff.Axis(Z,i)) << "\n";
			if (analyticalresults.at(p)(0,0,diff.Axis(Z,i) != 0.))
				diff(p,0,i) /= analyticalresults.at(p)(0,0,diff.Axis(Z,i));
			else
				diff(p,0,i) = 1.; 		
		}
	}
	diff -= 1.;
	output.close(); 
	return diff.mean();
		
}
