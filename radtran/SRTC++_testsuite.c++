#include "SRTC++_testsuite.h"

const std::string textred("\033[1;31m");
const std::string textgreen("\033[1;32m");
const std::string textreset("\033[0m");

class testresult
{
	public:
		string name;
		double expected_result;
		double actual_result;
		double allowed_variance;
		string pass_fail() 
		{
			string answer;
			if (abs(actual_result-expected_result)<=allowed_variance/qualityfactor)
				answer=textgreen+string("PASS")+textreset;
			else
				answer=textred+string("FAIL")+textreset;
			return answer;
		}
	
	private:
		
};



main(int argn, char **argv)
{
	vector<testresult> allresults;

	
// Test 1:  check that A=1 Lambertian surface diskintegrates to 2/3
	testresult LambertSurface;
	LambertSurface.name="LambertianSurface";
	LambertSurface.expected_result = 2./3.;
	LambertSurface.allowed_variance = 0.01;	
	LambertSurface.actual_result = test_LambertSurface();
	allresults.push_back(LambertSurface);
	 
// Test 2:  thick atmosphere Chandrasekhar test
	testresult ChandraThickAtmo;
	ChandraThickAtmo.name="ChandraThickAtmo";
	ChandraThickAtmo.expected_result = 0.0;
	ChandraThickAtmo.allowed_variance = 0.03; 
	ChandraThickAtmo.actual_result = test_ChandraAtmo();
	allresults.push_back(ChandraThickAtmo);


// Calculate Graham Results
	testresult GrahamAnalytical;
	GrahamAnalytical.name="GrahamAnalytical";
	GrahamAnalytical.expected_result = 0.0;
	GrahamAnalytical.allowed_variance = 1.0*qualityfactor;
	GrahamAnalytical.actual_result = test_GrahamCompare();
	allresults.push_back(GrahamAnalytical);


// Test 5: DISR atmosphere compared to Sebastien's plane-parellel model
	testresult SebastienCompare;	
	SebastienCompare.name="DISR_SebastienComp";
	SebastienCompare.expected_result=0.0; 
	SebastienCompare.allowed_variance = 1.0*qualityfactor; 
	SebastienCompare.actual_result = test_SebastienCompare(); //placeholder
	allresults.push_back(SebastienCompare);
	
	
// output results report
	cout << "\n\n*******************   RESULTS  REPORT  *****************\n";
	cout << "   #    Test Name               Expected        Actual          Pass/Fail\n";
	cout << "  ---   ---------               --------        ------          ---------\n";
	for (int i(0);i<allresults.size();i++) {
		cout << "Test " << i+1 << "\t" << allresults.at(i).name << "\t";
		cout << setw(8) << allresults.at(i).expected_result << "\t" << allresults.at(i).actual_result << "\t";
		cout.flush();
		cout << "\t";
		cout << allresults.at(i).pass_fail() << "\n";
	}
	
	
}
