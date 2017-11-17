class polynomial2d : public Jcrap::nDmultifunction<double>
{
	public:
		unsigned int order;
		vector<double> xylookup;
	
		polynomial2d(int nth) : 
				nDmultifunction<double>(2,addorial(nth+1))
		{
			cout << "Instantiating polynomial2d\n";
			order = nth; 
			
			for (int i=0;i<a.size();i++)
				a[i] = 1.;
			
			param_names[0] = string("a_0");
			for (int i=1;i<order;i++) {
				if (order >=1) {
					param_names[1] = string("a_x");
					param_names[2] = string("a_y");
				}
				if (order >=2) {
					param_names[3] = string("a_xx");
					param_names[4] = string("a_xy");
					param_names[5] = string("a_yy");
				}
				if (order >=3) {
					param_names[6] = string("a_xxx");
					param_names[7] = string("a_xxy");
					param_names[8] = string("a_xyy");
					param_names[9] = string("a_yyy");
				}
				if (order >=4) {
					param_names[10] = string("a_xxxx");
					param_names[11] = string("a_xxxy");
					param_names[12] = string("a_xxyy");
					param_names[13] = string("a_xyyy");
					param_names[14] = string("a_yyyy");
				}
			}
		}
		
		double the_function(const vector<double>& x, vector<double>& newa) {
			double answer(0.);
			for (int i=0;i<newa.size();i++) {
				a[i] = newa[i];
			}
			answer += a[0];
			if (order >= 1) {
				answer += a[1] * x[0];
				answer += a[2] * x[1];
			}
			if (order >= 2) {
				answer += a[3] * x[0] * x[0];
				answer += a[4] * x[0] * x[1];
				answer += a[5] * x[1] * x[1];
			}
			if (order >= 3) {
				answer += a[6] * x[0] * x[0] * x[0];
				answer += a[7] * x[0] * x[0] * x[1];
				answer += a[8] * x[0] * x[1] * x[1];
				answer += a[9] * x[1] * x[1] * x[1];
			}
			return answer;
		}
};
		
		
class powerlaw : public Jcrap::nDmultifunction<double>
{
	public:
		unsigned int order;
		vector<double> xylookup;
	
		powerlaw(int n) : 
				nDmultifunction<double>(2,1+6*n)
		{
			cout << "Instantiating powerlaw\n";
			order = n; 
			
			for (int i=0;i<a.size();i++)
				a[i] = 1.;
			
			param_names[0] = string("a_0");
			for (int i=0;i<order;i++) {
				param_names[1+6*i] = string("k_x")+int2str(i);
				param_names[2+6*i] = string("b_x")+int2str(i);
				param_names[3+6*i] = string("k_y")+int2str(i);
				param_names[4+6*i] = string("b_y")+int2str(i);
				param_names[5+6*i] = string("k_z")+int2str(i);
				param_names[6+6*i] = string("b_z")+int2str(i);
			}
		}
		
		double the_function(const vector<double>& x, vector<double>& newa) {
			double answer(0.);
			for (int i=0;i<newa.size();i++) {
				a[i] = newa[i];
			}
			answer += a[0];
			for (int i=0;i<(newa.size()-1)/6;i++) {
				answer += a[1+6*i]*pow(x[0],a[2+6*i]);
				answer += a[3+6*i]*pow(x[1],a[4+6*i]);
				answer += a[5+6*i]*pow(x[2],a[6+6*i]);
			}
			return answer;
		}
};
		
class exponential : public Jcrap::nDmultifunction<double>
{
	public:
		unsigned int order;
		vector<double> xylookup;
	
		exponential(int n) : 
				nDmultifunction<double>(2,1+6*n)
		{
			cout << "Instantiating exponential\n";
			order = n; 
			
			for (int i=0;i<a.size();i++)
				a[i] = 1.;
			
			param_names[0] = string("a_0");
			for (int i=0;i<order;i++) {
				param_names[1+6*i] = string("k_x")+int2str(i);
				param_names[2+6*i] = string("tau_x")+int2str(i);
				param_names[3+6*i] = string("k_y")+int2str(i);
				param_names[4+6*i] = string("tau_y")+int2str(i);
				param_names[5+6*i] = string("k_z")+int2str(i);
				param_names[6+6*i] = string("tau_z")+int2str(i);
			}
		}
		
		double the_function(const vector<double>& x, vector<double>& newa) {
			double answer(0.);
			for (int i=0;i<newa.size();i++) {
				a[i] = newa[i];
			}
			answer += a[0];
			for (int i=0;i<(newa.size()-1)/6;i++) {
				answer += a[1+6*i]*exp(x[0]*a[2+6*i]);
				answer += a[3+6*i]*exp(x[1]*a[4+6*i]);
				answer += a[5+6*i]*exp(x[2]*a[6+6*i]);
			}
			return answer;
		}
};
