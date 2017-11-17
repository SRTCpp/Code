class polynomial : public Jcrap::multifunction<double>
{
	public:
		unsigned int order;
	
		polynomial(int nth) : Jcrap::multifunction<double>(nth+1){
			order = nth+1; 
			for (unsigned int i=0;i<order;i++) {
				param_names[i] = string("a_") + int2str(int(i));
			}
		}
		
		double the_function(double x, vector<double>& newa, unsigned int=0) {
			double answer(0.);
			for (unsigned int i=0;i<order;i++) {
				a[i] = newa[i];
				answer += a[i] * pow(x,int(i));
//				cout << "evaluation, answer=" << answer << ", ia[" << i << "] = "
//						<< ia[i] << "\n";
//				cout << "a[" << i << "] = " << a[i] << "\n";

			}
			return answer;
		}
		
		vector<double> the_derivitives(double x, vector<double>& aa, unsigned int=0)
		{
			vector<double> answer(aa.size());
//			vector<double> h(aa.size());

// this is numerical, and for a polynomial obviously could be modified to
// evaluate the derivitives analytically, thus saving an assload of time.
						
//			for (unsigned int i=0;i<h.size();i++) h[i] = 1.;
	
			for (unsigned int i=0;i<aa.size();i++) {
//				answer[i] = aderiv(i,x,aa,h[i]);
//				answer[i] = aa[i]*pow(x,i);
//				cout << "deriv[" << i << "] = " << answer[i] << "\n";
				answer[i] = pow(x,double(i));
			}
//			cout << "\n";
			return answer;
		}
};

class gaussian : public Jcrap::multifunction<double>
{
	public:
		gaussian() : Jcrap::multifunction<double>(4){
			param_names[0] = "Center";
			param_names[1] = "Sigma";
			param_names[2] = "Magnitude";
			param_names[3] = "Offset";
		}
		
		double the_function(double x, vector<double>& newa, unsigned int=0) {
			double answer(0.);
			for (unsigned int i=0;i<newa.size();i++) {
				a[i] = newa[i];
			}
			double &center(a[0]), &sigma(a[1]), &mag(a[2]), &offset(a[3]);
			answer = mag/(sigma*1./*sqrt(2.*3.1415926535)*/) *
					exp(-(x-center)*(x-center)/(2.*sigma*sigma)) + offset;
//				cout << "evaluation, answer=" << answer << ", ia[" << i << "] = "
//						<< ia[i] << "\n";
//				cout << "a[" << i << "] = " << a[i] << "\n";

			
			return answer;
		}
		
/*		vector<double> the_derivitives(double x, vector<double>& aa)
		{
			vector<double> answer(aa.size());
			vector<double> h(aa.size());

// this is numerical, and for a polynomial obviously could be modified to
// evaluate the derivitives analytically, thus saving an assload of time.
						
			for (int i=0;i<h.size();i++) h[i] = 1.;
	
			for (int i=0;i<aa.size();i++) {
				answer[i] = aderiv(i,x,aa,h[i]);
//				answer[i] = aa[i]*pow(x,i);
//				cout << "deriv[" << i << "] = " << answer[i] << "\n";
			}
			return answer;
		}*/
};

class sine : public Jcrap::multifunction<double>
{
	public:
		sine() : Jcrap::multifunction<double>(4){
			param_names[0] = "Amplitude";
			param_names[1] = "Wavelength";
			param_names[2] = "Phase";
			param_names[3] = "Offset";
		}
		
		double the_function(double x, vector<double>& newa, unsigned int=0) {
			double answer(0.);
			for (unsigned int i=0;i<newa.size();i++) {
				a[i] = newa[i];
			}
			double &amplitude(a[0]), &wavelength(a[1]), &phase(a[2]), offset(a[3]);
			answer = amplitude*sin(2.*Jcrap::pi*x/wavelength+phase)+offset;
			return answer;
		}
};

class powerlaw : public Jcrap::multifunction<double>
{
	public:
		powerlaw() : Jcrap::multifunction<double>(3){
			param_names[0] = "A (coefficient)";
			param_names[1] = "B (exponent)";
			param_names[2] = "C (constant offset)";
		}
		
				
		double the_function(double x, vector<double>& newa, unsigned int=0) {
			double answer(0.);
			for (unsigned int i=0;i<newa.size();i++) {
				a[i] = newa[i];
			}
			double &A(a[0]), &B(a[1]), &C(a[2]);
			answer = A*pow(x,B)+C;
			return answer;
		}
};

class errorellipseintercept : public Jcrap::multifunction<double>
{
public:
	errorellipseintercept(int indim, cube inCinv, double inDelta, 
			vector<double> inaa) :
		dim(indim), 
		Cinv(inCinv), 
		Delta(inDelta)  
		{a = inaa;  thisdim=dim;}

	double the_function(double val, vector<double>& aa, unsigned int=0) {
		a = aa;
		a.at(thisdim)=val;
/*		cout << "In errorellipseintercept.the_function with ";
		for (int j(0);j<aa.size();j++) cout << "a["<<j<<"]="<<aa[j]<<" ";
		cout << "and thisdim=" << dim << "\n";*/
		double answer(-Delta);
		for (int x(0);x<Cinv.N(X);x++)
			for (int y(0);y<Cinv.N(Y);y++)
				answer += a.at(x)*Cinv(x,y,0)*a.at(y);
		return answer;  
	}

	
	int dim;
	cube Cinv;
	double Delta;
};

class ellipsemax : public Jcrap::multifunction<double>
{
public:
	ellipsemax(errorellipseintercept E) : e(E) {
		a=e.a;
		h=e.h;
		ia=e.ia;
		param_names=e.param_names;
		thisdim=e.thisdim;
	}

	double the_function(double val, vector<double>& aa, unsigned int=0) {
		double answer(0.);
		a = aa;
//		cout << "in ellipemax.the_function with ";
//		for (unsigned int i(0);i<a.size();i++) cout << param_names.at(i) << "=" << a.at(i) << " ";
//		cout << "\n";
		
		a.at(thisdim)=val;
		e.a=a;
		pair<double, double> oneDmin(e.findminimum(e.dim,-100.*fabs(h.at(e.dim)),0.,+100.*fabs(h.at(e.dim)), 1.e-8));
//		cube testcase(e.plot(oneDmin.first-100.*fabs(h.at(e.dim)), +100.*fabs(h.at(e.dim)),1000));
//		testcase.keyword("title", "minimumfinder");
//		testcase.graph();
//		cout << "oneDmin of " << oneDmin.second << " found at " << oneDmin.first;
//		cout << " h=" << h.at(e.dim) << "\n"; 
		if (oneDmin.second<0.) {
			answer = e.findzero(oneDmin.first-100.*fabs(h.at(e.dim)), oneDmin.first, 1.e-8);
//			testcase=e.plot(oneDmin.first-100.*fabs(h.at(e.dim)), oneDmin.first, 1000);
//			testcase.keyword("title", "zerofinder");		
//			testcase.graph();
//			char c; cin>>c;
			if (answer==1.e6) answer = e.findzero(oneDmin.first, oneDmin.first-100.*e.Delta, 1.e-8);
		} else answer = oneDmin.second;
		return answer;
	}

private:
	errorellipseintercept e;
};
