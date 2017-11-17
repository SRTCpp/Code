
template <class Datatype>
Jcrap::multifunction<Datatype>::multifunction(int na)
{
//	cout << "instantiating multifunction("<<na<<")\n";
	a.resize(na);
	param_names.resize(na);
	ia.resize(na);
	quantized.resize(na);
	for (int i(0);i<na;i++) quantized[i] = 0;
	dyda.resize(na);
	h.resize(na);
	
	for (unsigned int n(0);n<h.size();n++)
		h[n] = 0.1;
	
	covar.resize(na+1);
	for (unsigned int n=0;n<covar.size();n++) covar[n].resize(na+1);
	alpha.resize(na+1);
	for (unsigned int n=0;n<alpha.size();n++) alpha[n].resize(na+1);
	
	alamda = -1.;		// initialization condition for mrqmin
	lamda_up   = 10.;
	lamda_down = 10.;
	upperstop  = 20.;
	initialalamda = 0.001;
	chidiffstop = 0.1;
	
	for (int n=0;n<na;n++) {
		char outstr[200];
		sprintf(outstr, "a%d", n);
		param_names[n] = outstr;
	}
	thisdim=-1;
	
	skiperroranalysis = false;
//	cout << "done\n";
}


template <class Datatype>
Jcrap::multifunction<Datatype>::multifunction(vector<Datatype> aa)
{
	*this=multifunction(aa.size());
	a=aa;
}


template <class Datatype>
vector<Datatype> Jcrap::multifunction<Datatype>::the_derivitives(Datatype x,
		vector<Datatype>& aa, unsigned int category)
{
	vector<Datatype> answer(aa.size());
	
	for (unsigned int i=0;i<aa.size();i++) {
		if (osuppress<0) if (ia[i]) cout << "in the_deriv[i]tives (sic), h[i]="<<h[i]<<"\n";
		if (osuppress<0) if (ia[i]) cout << "aderiv["<<i<<"]=";
		if (ia[i])
			if (quantized[i])
				answer[i] = dumbderiv(i,x,aa,h[i],category);
			else
				answer[i] = aderiv(i,x,aa,h[i],category);
		else
			answer[i] = 0.;
		if (osuppress<0) if (ia[i]) cout << answer[i]<<" ";
		if (osuppress<0) if (ia[i]) cout << " a["<<i<<"]=" << a[i] << " ";
	}
	if (osuppress<0) cout << "\n";
	return answer;
}


#define NRANSI
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::aderiv(int ii, Datatype x, vector<Datatype> av, Datatype h, unsigned int category)
{
	Datatype err;
	Datatype ans(0.);
	
	if (h == -1.) h = x/10.;
	if (h ==  0.) h = 1.;
	int i,j;
	Datatype errt,fac,hh,plus(0.),minus;
	
	vector<vector<Datatype> > aa(NTAB+1);
	for (i=0;i<=NTAB;i++) aa[i].resize(NTAB+1);
	
	hh=h;
	av[ii] += hh;
	if (osuppress<-1) cout << "Doing numerical derivative in aderiv... \n";
	if (osuppress<-1) cout << "x = " << x << "\n";
	for (unsigned int v=0;v<av.size();v++) // cout << "av[" << v << "] = " << av[v] << "\n";
	plus = the_function(x,av,category);
	av[ii] -= 2.*hh;
	minus = the_function(x,av,category);
	if (osuppress<-1)  cout << "plus= " << plus << " minus = " << minus << ", hh = " << hh << "\n";
	aa[1][1]=(plus-minus)/(2.0*hh);
	av[ii] += hh;
	err=BIG;
	for (i=2;i<=NTAB;i++) {
		if (osuppress<-1) cout << "aa[1]["<<i-1<<"]="<<aa[1][i-1]<<" ";
		if (osuppress<-1) cout << "for ii=" << ii << " over distance " << 2.*hh;
		if (osuppress<-1) cout << " slope guess = " << aa[1][i-1] <<"\n";
		hh /= CON;
		av[ii] += hh;
		plus = the_function(x,av,category);
		av[ii] -= 2.*hh;
		minus = the_function(x,av,category);
		if (osuppress<-1) cout << "plus= " << plus << " minus = " << minus << ", hh = " << hh << "\n";
		aa[1][i]=(plus-minus)/(2.0*hh);
		av[ii] += hh;
		fac=CON2;
		for (j=2;j<=i;j++) {
			aa[j][i]=(aa[j-1][i]*fac-aa[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=JMAX(jabs(aa[j][i]-aa[j-1][i]),jabs(aa[j][i]-aa[j-1][i-1]));
			
			if (errt <= err) {
				err=errt;
				ans=aa[j][i];
			}
		}
		if (jabs(aa[i][i]-aa[i-1][i-1]) >= SAFE*(err)) {
			if (i==2 && ans!=0. && osuppress < 0) {
				cout << "POTENTIAL NUMERICAL DERIVATIVE FAILURE IN aderiv!!!!!\n";
				cout << "Entire 'aa' matrix, coming up:\n";
				for (int i2(0);i2<=2;i2++) {
					for (int j2(0);j2<=2;j2++) {
						cout << "aa["<<j2<<"]["<<i2<<"]="<<aa[j2][i2]<<" ";
					}
					cout << "\n";
				}
			}
			break;
		}
	}
	return ans;
}

#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE
#undef NRANSI
		
template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::dumbderiv(int ii, Datatype x, vector<Datatype> av, Datatype h, unsigned int category)
{	
	av.size(); // just to avoid compiler warning.
	a[ii]+=h;
	Datatype plusval(the_function(x,a,category));
	a[ii]-=h*2.;
	Datatype minusval(the_function(x,a,category));
	a[ii]+=h;
	

	Datatype ans(0.);
	ans = (plusval-minusval)/(2.*h);
	
	return ans;
};

template <class Datatype>
pair<vector<Datatype>, Datatype> Jcrap::multifunction<Datatype>::fit(cube data, vector<Datatype> aa, 
		vector<int> iaa)
/* Given a cube of size (2,1,n) where x=0 contains data and x=1 contains sigmas,
	or just (1,1,n) in which case the sigmas are assumed to be uniformly 1., fit
	to the_function based on the vector aa as initial guesses for the parameters */
{	
	osuppress++;
	ia=iaa;
	if (data.N(X) < 2)
		data = data.blocksx(data*0.+1.);
	
	vector<Datatype> vector_x(0);
	vector<Datatype> vector_y(0);
	vector<Datatype> vector_sigma(0);
	vector<unsigned int> vector_category(0);
	
	if (!(osuppress-1)) cout << "Inside fit, moving from cubes to vectors\n";
	for (int i=0;i<data.N(Z);i++) {
		for (unsigned int y(0);y<data.N(Y);y++) {
			if (data(1,int(y),i)<=0.) continue;  // if the error is zero, stop considering this point
			vector_x.push_back(data.Axis(Z,i));
			vector_y.push_back(data(0,int(y),i));
			vector_sigma.push_back(data(1,int(y),i));
			vector_category.push_back(y);
//			cout << "vector_x(" << goodpoints << ")=" << vector_x.at(goodpoints) << "\n";
//			cout << "Category = " << y << " for point at " << data.Axis(Z,i) << ", " << vector_y.at(goodpoints) << "\n";
		}
	}
//	cout << "after vectors:  xsize=" << vector_x.size() << "\n";

	int degreesoffreedom=0;
	for (unsigned int i=0;i<iaa.size();i++) {
		if (iaa[i]) degreesoffreedom++;
		if (osuppress<1) cout << "iaa[" << i << "] = " << iaa[i] << "\n";
	}
		
//	cout << "First run, alamd=-1\n";
	pair<vector<Datatype>, Datatype> answer;	
	Datatype chisq=1.e98, oldchisq=1.e99, oldalamda(1.e97);
	alamda = -1.;
//	chisq = Jmrqmin(vector_x, vector_y, vector_sigma, aa, iaa);
//	cout << "new reduced chisq = " << chisq/vector_x.size() << "\n";
//	for (int i=0;i<aa.size();i++)
//		cout << " aa[" << i << "] = " << aa[i];// << "\n";
//	cout << "\t\talamda = " << alamda << "\n";

	while ((jabs(chisq-oldchisq)/vector_x.size() > chidiffstop  ||  oldalamda<alamda)  &&  alamda < pow(lamda_up,upperstop)) {
		ifstream haltcheck("halt");
		cout << "Checking for manual halt..."; cout.flush();
		if (haltcheck) {
			cout << "  halting!"; cout.flush();
			system("rm halt");
			break;
		}
		cout << "\n";
		if (osuppress<2) { cout << "Step:"; cout.flush(); }
		oldchisq = chisq;
		oldalamda = alamda;
		chisq = Jmrqmin(vector_x, vector_y, vector_sigma, aa, iaa, vector_category);
		if (osuppress<2) {
			cout << "new reduced chisq = " << chisq/(vector_x.size()-degreesoffreedom);
			cout << "; unreduced chisq = " << chisq;
			for (unsigned int i=0;i<aa.size();i++) {
				cout << " ";
				if (iaa[i]) cout << "[";
				cout << param_names[i] << "=" << aa[i];
				if (iaa[i]) cout << "]";
			}
			cout << "\n";
			cout << "alamda = " << alamda << "\n";
		}
	}
	
	alamda = 0.;
	chisq = Jmrqmin(vector_x, vector_y, vector_sigma, aa, iaa, vector_category);
	
	vector<Datatype> errors(a.size());
	if (!skiperroranalysis) errors = finderrors();
	
	cout << "\n\n*****************************  MEASURED VALUES  ********************************\n";
	for (unsigned int i(0);i<ia.size();i++) {
		if (ia.at(i)) {
			cout << "    " << setprecision(12) << param_names.at(i) << " = " << aa.at(i) << " +/- " << errors.at(i) << "\n";
		}
	}

	cout << "********************************************************************************\n\n\n";
	osuppress--;	
	return pair<vector<Datatype>, Datatype> (aa, chisq);
}



template <class Datatype>
vector<Datatype> Jcrap::multifunction<Datatype>::finderrors() const
// following formula on page 815 of 2010 NR
{
	const bool simplex(false), dumb(false), zeroversion(true), dumb2(false);  // 'zeroversion' and 'dumb2' work -- others don't
	static const double X2_limits[16]={0.,1.00,2.30,3.53,4.72,5.89,7.04,8.17636,9.30404,
			10.4235,11.5361,12.643,13.7448,14.8423,15.936,17.0262};  // faked from Chisqdist.infcdf()
	vector<Datatype> answer(ia.size());
	
	int poi(count(ia.begin(), ia.end(), true)); // parameters of interest
	if (poi > 15) cout << "FAIL -- can't handle poi=" << poi << "!!!\n";
//	double p(0.6827);  // 1-sigma errors
	Datatype Delta(X2_limits[poi]);
	
	cube C(poi, poi, 1);
	for (int ox(0),Cx(0);ox<int(ia.size());ox++) {
		if (ia.at(ox)) {
			for (int oy(0),Cy(0);oy<int(ia.size());oy++) {
				if (ia.at(oy)) {
					C(Cx,Cy,0) = covar.at(ox+1).at(oy+1);
					Cy++;
				}
			}
			Cx++;
		}
	}
	
	cube Cinv(C.matrixinvert());
	
	vector<Datatype> startinga(poi), betterstart(poi);
	vector<Datatype> relevanth;
	for (int i(0);i<poi;i++) startinga.at(i) = 0.;
	for (unsigned int i(0);i<ia.size();i++)
		if (ia.at(i)) relevanth.push_back(-h.at(i)); 

	
	for (unsigned int i(0),Ci(0);i<ia.size();i++) {
		if (ia.at(i)) {
			for (unsigned int j(0);j<startinga.size();j++) startinga.at(j)=0.;
			errorellipseintercept eei(Ci,Cinv,Delta,startinga);
			eei.param_names.resize(relevanth.size());
			for (unsigned int j(0),Ci2(0);j<ia.size();j++)
				if (ia.at(j)) {eei.param_names.at(Ci2) = param_names.at(j); Ci2++;}
				
			eei.h=relevanth;
			if (simplex || dumb) {
				for (unsigned int j(0),Ci2(0);j<ia.size();j++) {
					if (ia.at(j) && j!=Ci) {
						cout << "looking for appropriate scale for axis " << Ci2 << "\n";
						eei.dim=Ci;
						eei.a=startinga;
						double thismin(eei.findminimum(Ci2, -200.*h.at(j), 0., 200.*h.at(j), 1.e-6).first);
						eei.thisdim=Ci2;
						cube eeiplot(eei.plot(-2.*thismin, 2.*thismin, 1000));
						eeiplot.keyword("title", param_names.at(j));
						eeiplot.graph();
						relevanth.at(Ci2)=fabs(thismin);
						betterstart.at(Ci2)=thismin;
						Ci2++;					
					} else if (ia.at(j) && j==Ci && startinga.at(Ci)==0.) {
						double dimmin(eei.findminimum(Ci, -200.*h.at(i), 0., 200.*h.at(j), 1.e-6).first);
						startinga.at(Ci)=dimmin;
						betterstart.at(Ci)=dimmin;
						j=0; Ci2=0;
					}
				}
			}
			if (simplex || dumb) {
				char c; cin>>c;
				eei.h=relevanth;
				eei.thisdim=eei.dim=Ci;
				for (unsigned int j(0);j<startinga.size();j++) startinga.at(j)=0.;
				eei.a=startinga;
				startinga.at(Ci) = eei.findminimum(Ci, -200.*h.at(i), -h.at(i), 0., 1.e-6).first;
				cout << "Starting point for " << param_names.at(i) << " = " << startinga.at(Ci) << "\n";
//				startinga[2] = 0.02;
			}
			pair<vector<Datatype>, Datatype> minlocation;
			if (simplex) minlocation = eei.findminimum(betterstart, relevanth);
			if (dumb) eei.a=startinga;
			if (dumb) minlocation = eei.findminimum_dumb();
			if (dumb2) eei.a=startinga;
			if (dumb2) {
				for (unsigned int j(0);j<startinga.size();j++) startinga.at(j)=0.;
				eei.a=startinga;
				ellipsemax em(eei);
				minlocation = em.findminimum_dumb(startinga);
			}
			if (zeroversion) {
				for (unsigned int j(0);j<startinga.size();j++) startinga.at(j)=0.;
				eei.a=startinga;
				ellipsemax em(eei);
				cout << "about to findminimum, thisdim for transitfunc=" << thisdim << "\n";
				minlocation = em.findminimum(betterstart,relevanth);
				cout << "now it's " << thisdim << "\n";
			}

			answer.at(i) = -minlocation.second;
			cout << "Minimum is:  " << minlocation.second << " at " << answer.at(i) << "\n";
/*			cube eplot(1,1,1000);
			cout << "h[" << i << "]=" << h.at(i) << "\n";
			for (int j(0);j<eplot.N(Z);j++) {
				eplot.Axis(Z,j) = double(j-eplot.N(Z)/2)*h.at(i)/100.;
				eplot(0,0,j) = eei.the_function(eplot.Axis(Z,j), startinga);
				thisdim=i;
				pair<double, double> oneDmin(eei.findminimum(Ci,-h.at(i)*100.,0.,+h.at(i)*100., 1.e-8));
				cout << "minimum at " << oneDmin.first << " is " << oneDmin.second << " ? \n"; 
			}
			eplot.keyword("title", param_names.at(i));
			eplot.graph();*/
//			answer.at(i) = eei.findzero(0., h.at(i)*10., 1.e-8);  // old, wrong, non-covariant way
			startinga.at(Ci) = 0.;
			Ci++;
		} else
			answer.at(i) = 0.;
	}
	
	return answer;
}



template <class Datatype>
pair<vector<Datatype>, Datatype> Jcrap::multifunction<Datatype>::fit(cube data, vector<Datatype> aa)
// If user doesn't specify what to fit for, fit for everything.
{
	vector<int> iaa(aa.size());
	for (int i=0;i<aa.size();i++)
		iaa[i] = 1;
	ia=iaa;
	if (!osuppress) cout << "Instantiating fit(,,)\n";
	return fit(data, aa, iaa);		
}

template <class Datatype>
pair<vector<Datatype>, Datatype> Jcrap::multifunction<Datatype>::fit(cube data)
{
	return fit(data,a);
}







#define NRANSI

template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::Jmrqmin(vector<Datatype> &x, vector<Datatype> &y, 
	vector<Datatype> &sig, vector<Datatype>& aa, vector<int> &ia, vector<unsigned int> &vector_category)
{
//	int ndata(x.size());					//	(blank comment = added by JB)
	int ma(aa.size());						//
	static Datatype chisq;						//
	int j,k,l;
	static int mfit;
	static float ochisq;
	static vector<vector<Datatype> > oneda;		//
	static vector<Datatype> atry, beta, da;		//
	int ianum(0);                                   //
	for (unsigned int i(0);i<ia.size();i++) ianum+=ia[i];
			

	if (alamda < 0.0) {
		atry.resize(ma);
		beta.resize(ma+1);
		da.resize(ma+1);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j-1]) mfit++;
		oneda.resize(mfit+1);
		for (int i=0;i<=mfit;i++) oneda[i].resize(mfit+1);
		alamda= initialalamda;
		chisq = Jmrqcof(x,y,sig,aa,ia,ma,alpha,beta,vector_category);
/*		for (int newx(0);newx<ianum;newx++)
			for (int newy(0);newy<ianum;newy++)
				cout << "alpha("<<newx<<","<<newy<<")=" << alpha[newx][newy] << "\n";*/
		ochisq=chisq;
		if (osuppress < 1) cout << "First run:  chisq=";
		if (osuppress < 1) cout << chisq << " ";
		for (j=0;j<ma;j++) {
			atry[j]=aa[j];
			if (osuppress < 1) cout << " ";
			if (ia[j]) if (osuppress < 1) cout << "[";
			if (osuppress < 1) cout << param_names[j] << "=" << aa[j];
			if (ia[j]) if (osuppress < 1) cout << "]";
		}
		if (osuppress < 1) cout << "\n";
	}
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(alamda));
		oneda[j][1]=beta[j];
	}
//	if (!osuppress) cout << "before JgaussJ\n";
	Jgaussj(covar,mfit,oneda,1);
//	if (!osuppress) cout << "after JgaussJ\n";
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (alamda == 0.0) {
		Jcovsrt(covar,ma,ia,mfit);
		Jcovsrt(alpha,ma,ia,mfit);
		return chisq;
	}
//	if (!osuppress) cout << "after covsrts\n";
	for (j=0,l=1;l<=ma;l++)
		if (ia[l-1]) atry[l-1]=aa[l-1]+da[++j];
// added July 30, 2003 JB to account for other aa[] variables changing
	for (l=0;l<ma;l++) if (!ia[l]) atry[l]=aconvert(aa, ia, l);
	

//	cout << "\n\nDebug:  trying values -- ";
	
//	for (j=0;j<ma;j++) {
//		cout << " atry[" << j << "] = " << atry[j];
//	}
//	cout << "\n";
	
	chisq = Jmrqcof(x,y,sig,atry,ia,ma,covar,da,vector_category);
	if (chisq < ochisq) {
		cout << "Good step:  "; cout.flush();
		alamda /= lamda_down;
		ochisq=chisq;
		for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
		}
		for (l=0;l<ma;l++) {
			if (ia[l]) aa[l]=atry[l];
			else aa[l]=aconvert(atry, ia, l);
		}
		cout << " out of Jmrqcof\n";
	} else {
		if (osuppress < 1) {	
			int degreesoffreedom=0;
			for (unsigned int i=0;i<ia.size();i++)
				if (ia[i]) degreesoffreedom++;
			
			cout << "Bad step values: badchisq=" << 
					chisq/(x.size()-degreesoffreedom) << "  ---  ";
			for (unsigned int i=0;i<aa.size();i++) {
				cout << " ";
				if (ia[i]) cout << "[";
				cout << param_names[i] << "=" << atry[i];
				if (ia[i]) cout << "]";
			}
		}
		if (!osuppress) cout << "\n";
		alamda *= lamda_up;
		chisq=ochisq;
	}
	
	return chisq;
}
#undef NRANSI


#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

template <class Datatype>
void Jcrap::multifunction<Datatype>::Jcovsrt(vector<vector<Datatype> >& covar, int ma, vector<int> &ia, int mfit)
{
	int i,j,k;
	float swap;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j-1]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}
#undef SWAP
#include <math.h>


#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

template <class Datatype>
void Jcrap::multifunction<Datatype>::Jgaussj(vector<vector<Datatype> > &a, int n, vector<vector<Datatype> > &b, int m)
{
	vector<int> indxc(n+1), indxr(n+1), ipiv(n+1);
	int i,icol=1,irow=1,j,k,l,ll;
	float big,dum,pivinv,temp;

	for (j=1;j<=n;j++) ipiv[j]=0;
//	cout << "in Jgaussj\n";
	for (i=1;i<=n;i++) {
//		cout << "in big loop i=" << i << "\n";
		big=0.0;
		for (j=1;j<=n;j++) {
//			cout << "j = " << j << "\n";
			if (ipiv[j] != 1) {
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
//							cout << "ipiv[k]=0, irow=" << irow << ", icol=" << icol
//									<< ", j=" << j << ", k=" << k << " a[j][k]=" <<
//									a[j][k] << ", big=" << big << "\n";
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) {
						cout << "Jgaussj(Jmultifunction): Singular Matrix-1\n";
						for (unsigned int y(1);y<a[0].size();y++) {
							for (unsigned int x(1);x<a.size();x++)
								if (!(a[x][y]<=0 || a[x][y]>=0))
									cout << "("<<x<<","<<y<<")=" << a[x][y] << " ";
						}
					}
				}
			}
		}
		++(ipiv[icol]);
//		cout << "m=" << m << "; n=" << n << "\n";
//		cout << "a.size=" << a.size() << "; b.size=" << b.size() << "\n";
//		cout << "a[0].size=" << a[0].size() << "; b[0].size=" << b[0].size() << "\n";
//		cout << "irow=" << irow << "; icol=" << icol << "\n";
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
//		cout << "mid Jgaussj\n";
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) {
			cout << "Jgaussj(Jmultifunction): Singular Matrix-2\n";
			cout << "a["<<icol<<"]["<<icol<<"] = " << a[icol][icol] << "\n";
			a[icol][icol]=1.;  // fake it to keep the fit going.  might as well
		}		
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
//	cout << "nearing end of Jgaussj\n";
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
}
#undef SWAP
#undef NRANSI


#define NRANSI

template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::Jmrqcof(vector<Datatype> &x, vector<Datatype> &y,
		vector<Datatype> &sig, vector<Datatype> &aa, vector<int> &ia, int ma, 
		vector<vector<Datatype> > &alpha, vector<Datatype> &beta, vector<unsigned int> &vector_category)
{
	Datatype chisq;
	int ndata = x.size()-1;
	int mfit=0;

	for (int j(1);j<=ma;j++)
		if (ia[j-1]) mfit++;
	for (int j(1);j<=mfit;j++) {
		for (int k(1);k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	chisq=0.0;
	if ((osuppress-1)<=0) cout << "Jmrqcof progress:  00%";
//#pragma omp parallel for ordered schedule(dynamic)
	for (int i=1;i<=ndata;i++) {
		vector<Datatype> dyda(ma);
		vector<Datatype> aaa(aa);
//		cout << "i=" << i << " running in parallel\n";
		if ((osuppress-1)<=0)  {
			#pragma omp critical
			printpercent(i-1,ndata-1);
		}
		float ymod(funcs(x[i-1],aaa,dyda,vector_category.at(i-1)));
//		if (i<int(ia.size()) && ia[i-1]) cout << "after funcs, aa["<<i-1<<"] = " << aa[i-1] << " dyda=" << dyda[i-1] << "\n";
		if (!(dyda[1]<1. && dyda[1]>-1.)) {
			cout << "Potentially screwed up tau derivative of " << dyda[1];
			cout << " for point number " << i-1 << "\n";
			dyda[1] = 0.;
		}/* else {
			cout << "Tau derivative of " << dyda[1];
			cout << " for point number " << i-1 << "\n";
		}*/
		float sig2i(1.0/(sig[i-1]*sig[i-1]));
		float dy(y[i-1]-ymod);
		for (int j(0),l(1);l<=ma;l++) {
			if (ia[l-1]) {
				float wt(dyda[l-1]*sig2i);
				j++;
				for (int k(0),m(1);m<=l;m++) {
					if (ia[m-1]) {
						#pragma omp atomic
						alpha[j][++k] += wt*dyda[m-1];
					}
				}
				#pragma omp atomic
				beta[j] += dy*wt;
			}
		}
//		cout << "chisq before #=" << i << " " << chisq << "\n";
//		cout << "dy = " << dy << ", sig2i= " << sig2i;
//		cout << ", sig["<<i-1<<"]=" << sig[i-1] << "\n";
		#pragma omp atomic
		chisq += dy*dy*sig2i;
//		cout << "chisq after = " << chisq << "\n";
//		cout << "i=" << i << " finished\n";
	}
	if ((osuppress-1)<=0) printpercent(ndata, ndata);
	if ((osuppress-1)<=0) cout << "\n";
	for (int j(2);j<=mfit;j++)
		for (int k(1);k<j;k++) alpha[k][j]=alpha[j][k];
	
	cout << "Jmrqcof returning chisq=" << chisq << "\n";
	return chisq;
}
#undef NRANSI


template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::funcs(Datatype x, vector<Datatype> &a, vector<Datatype> &dyda, unsigned int category)
{
//	cout << "funcs\n";
	dyda = the_derivitives(x, a, category);
//	cout << "dyda:  ";
//	for (unsigned int i(0);i<dyda.size();i++) cout << dyda[i] << " ";
//	cout << "\n";
	return the_function(x, a, category);
}

template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::aconvert(vector<Datatype>& a, vector<int>& ia, int i)
// Override this if your a[] are functions of one another.
{
	if (0)	cout << "in multifunction::aconvert ia[0] = " << ia[0] << "\n";
	return a[i];
}


template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::chisquared(cube& data)
{
	Datatype answer(0);
	for (int i(0);i<data.N(Z);i++) {
		double thischisq((data(0,0,i)-the_function(data.Axis(Z,i)))/data(1,0,i));
		answer += thischisq*thischisq;
	}
	return answer;
} 

template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::reducedchisquared(cube& data)
{
	Datatype factor(data.N(Z));
	for (unsigned int i(0);i<ia.size();i++)
		factor += Datatype(ia.at(i));
	return chisquared(data)/factor;
} 


// some methods for MCMC
template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::Pstate(cube& data) 
{ 
	double totalprob(0.);
	for (int i(0);i<data.N(Z);i++) {
		double thisprob;
		thisprob = (data(0,0,i)-the_function(data.Axis(Z,i)))/data(1,0,i);
		thisprob *= thisprob;
		thisprob /= -2.;
//		cout << "N=" << i << "; data=" << data(0,0,i) << " This probability:  " << thisprob << "\n";
		totalprob += thisprob;
	}
//	cout << "Total probability:  " << totalprob << "\n";
	return totalprob;
	
}



template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::Pstate(cube& data, vector<Datatype> aa) 
{ 
	a = aa;
	return Pstate(data); 
}



#include "../../nr/nr3/ran.h"
#include "../../nr/nr3/deviates.h"

template <class Datatype>
vector<Datatype> Jcrap::multifunction<Datatype>::Proposal() 
{
	static Normaldev gau(0.,1.,0);
	static const Datatype stepsize(0.01);
	
	vector<Datatype> proposed_a(a);
	
	for (unsigned int i(0);i<proposed_a.size();i++) {
//		cout << param_names.at(i) << ": " << ia[i] << "\n";
		if (ia[i]) {
//			cout << "In proposal, old " << param_names.at(i) << "=";
//			cout << a.at(i) << "; proposed=";
			if (gau.doub() >= 0.5) 
				proposed_a.at(i) += gau.dev()*h[i]*stepsize;
//			cout << proposed_a.at(i) << "\n";
		}
	}
	
	return proposed_a;
}

template <class Datatype>
Datatype Jcrap::multifunction<Datatype>::MCMCstep(cube &data) 
{
	static Normaldev gau(0.,1.,0);
	
	static vector<Datatype> original_a(a);
	vector<Datatype> proposed_a(Proposal());
	
	static double currentprobability(-1.);
	if (currentprobability == -1.)
		currentprobability = Pstate(data, original_a);
	double proposedprobability(Pstate(data, proposed_a));
//	cout << "currentprob="<<currentprobability<<", prpoposed="<<proposedprobability<<"\n";
	double alpha(min(1., exp(-currentprobability+proposedprobability)));
//	cout << "alpha=" << alpha << " -- " << proposedprobability-currentprobability <<
//			"\n";
	if (gau.doub() < alpha) {
/*		cout << "Proposed prob=" << proposedprobability << ", current=";
		cout << currentprobability;
		for (unsigned int i(0);i<a.size();i++)
			if (original_a.at(i) != proposed_a.at(i))
				cout << param_names.at(i) << ":" << proposed_a.at(i)-original_a.at(i) << " ";
		cout << "\n";*/
		a = proposed_a;
		original_a = proposed_a;
		currentprobability = proposedprobability;
		return 1.;
	} else {
		a = original_a;
		return 0.;	
	}
}

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findminimum()
{
	return findminimum(a, h, 1.e-8);
} 

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findmaximum()
{
	return findmaximum(a, h, 1.e-8);
} 

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findminimum(vector<Datatype> guesspoint, Datatype offset,
		Datatype in_tolerance)
{
	 vector<Datatype> deltas(guesspoint.size(), offset);
	 return findminimum(guesspoint, deltas, in_tolerance);
}

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findminimum(vector<Datatype>
		guesspoint, vector<Datatype> dels, Datatype in_tolerance)
// provide this with an initial guess vector along with typical length scales (dels)
// to have it find an initial simplex to apply to the amoeba minimizer.
{
	int simplexstyle(0); // 0 for offsets in each dimension, 1 for attempted diagonals to follow contours
	Int ndim=guesspoint.size();
	MatDoub pp(ndim+1,ndim);
	
	if (simplexstyle==0) {
		for (Int i=0;i<ndim+1;i++) {
			for (Int j=0;j<ndim;j++)
				pp[i][j]=guesspoint[j];
			if (i != 0) {
				if (dels[i-1]!=0.)
					pp[i][i-1] += dels[i-1];
				else 
					pp[i][i-1] += 0.1;
			}
		}
	}
	if (simplexstyle==1) {
		for (Int i=0;i<ndim+1;i++) {
			for (Int j=0;j<ndim;j++){
				pp[i][j]=guesspoint[j];
				if (i != 0) pp[i][j]*=2.;
			}
			if (i >= 2) {
				pp[i][i-1] /= 2.;
			}
		}
	}
	
//	for (Int i=0;i<ndim+1;i++) {
//		for (Int j=0;j<ndim;j++)
//			cout << "pp["<<i<<"]["<<j<<"]="<<pp[i][j]<< " ";
//		cout << "\n";
//	}
	
//	cout << "automatic simplex created; now executing findminimum starting with it.\n";
	return findminimum(pp, in_tolerance);
}

template <class Datatype>
struct Jcrap::multifunction<Datatype>::Jamoeba 
{
	const Doub ftol;
	long long nfunc;
	Int mpts;
	Int ndim;
	Doub fmin;
	VecDoub y;
	MatDoub p;
	Jamoeba(const Doub ftoll) : ftol(ftoll) {}
	
	template <class T>
	VecDoub minimize(MatDoub_I &pp, Jcrap::multifunction<T> &func)
	{
		const long long NMAX=5000000;
		const Doub TINY=1.0e-5;
		Int ihi,ilo,inhi;
		mpts=pp.nrows();
		ndim=pp.ncols();
		VecDoub psum(ndim),pmin(ndim);
		vector<T> x(ndim);
		p=pp;
		y.resize(mpts);
		
		if (osuppress<2) cout << "Downhill simplex minimization.  Current minimum:  000000"; cout.flush();
		
		for (Int i=0;i<mpts;i++) {
			for (Int j=0;j<ndim;j++)
				x[j]=p[i][j];
			y[i]=func(x);
		}
		nfunc=0;
		get_psum(p,psum);
		for (;;) {
			ilo=0;
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (Int i=0;i<mpts;i++) {
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi]) {
					inhi=ihi;
					ihi=i;
				} else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			Doub rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
			if (rtol < ftol) {
				SWAP(y[0],y[ilo]);
				for (Int i=0;i<ndim;i++) {
					SWAP(p[0][i],p[ilo][i]);
					pmin[i]=p[0][i];
				}
				fmin=y[0];
				if (osuppress<2) cout << "\n";
				return pmin;
			}
			if (nfunc >= NMAX) throw("NMAX exceeded");
			if (nfunc%100==0 && osuppress<2) {
				for (int i(0);i<17;i++)
					cout << char(8); // backspace
				cout.width(8);
				cout << nfunc;
				cout << " ";
				cout.width(8);
				cout << y[ilo];
				cout.flush();
			}
			nfunc += 2;
			Doub ytry=amotry(p,y,psum,ihi,-1.0,func);
			if (ytry <= y[ilo])
				ytry=amotry(p,y,psum,ihi,2.0,func);
			else if (ytry >= y[inhi]) {
				Doub ysave=y[ihi];
				ytry=amotry(p,y,psum,ihi,0.5,func);
				if (ytry >= ysave) {
					for (Int i=0;i<mpts;i++) {
						if (i != ilo) {
							for (Int j=0;j<ndim;j++)
								p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
							y[i]=func(psum);
						}
					}
					nfunc += ndim;
					get_psum(p,psum);
				}
			} else --nfunc;
		}
	}
	inline void get_psum(MatDoub_I &p, VecDoub_O &psum)
	{
		for (Int j=0;j<ndim;j++) {
			Doub sum=0.0;
			for (Int i=0;i<mpts;i++)
				sum += p[i][j];
			psum[j]=sum;
		}
	}
	template <class T>
	Doub amotry(MatDoub_IO &p, VecDoub_O &y, VecDoub_IO &psum,
		const Int ihi, const Doub fac, Jcrap::multifunction<T> &func)
	{
		VecDoub ptry(ndim);
		Doub fac1=(1.0-fac)/ndim;
		Doub fac2=fac1-fac;
		for (Int j=0;j<ndim;j++)
			ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
		Doub ytry=func(ptry);
		if (ytry < y[ihi]) {
			y[ihi]=ytry;
			for (Int j=0;j<ndim;j++) {
				psum[j] += ptry[j]-p[ihi][j];
				p[ihi][j]=ptry[j];
			}
		}
		return ytry;
	}
};

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findminimum(MatDoub_I &pp,
		Datatype in_tolerance)
{
//	cout << "in findminimum (MatDoub, tolerance)\n";
	MatDoub new_pp(pp);
//	cout << "Creating amoeba\n";
	Jamoeba amoeba(in_tolerance);
	vector<Datatype> bestfit;
	for (int i(0);i<1;i++) {  // can dial this up to avoid algorithm stopping when it shouldn't
		new_pp=pp;
//		cout << "amoeba run #" << i << " nrows=" << new_pp.nrows() << " ncols=" << new_pp.ncols() << "\n";
		bestfit = amoeba.minimize(new_pp, *this);
//		cout << "generating new seed simplex ... nrows=" << new_pp.nrows() << " ncols=" << new_pp.ncols() << "\n";
//		cout << "bestfit.size()=" << bestfit.size() << "\n";
		for (int x(0);x<new_pp.nrows();x++)
			for (int y(0);y<new_pp.ncols();y++)
				new_pp.assign(x,y,pp[x][y]-pp[0][y]+bestfit.at(y));
	}
	pair<vector<Datatype>, Datatype> answer(bestfit,the_function(0.,bestfit));
//	cout << "Minimum of " << answer.second << " found at ";
//	for (int i(0);i<answer.first.size();i++)
//		cout << param_names.at(i) << "=" << answer.first.at(i) << " ";
//	cout << "\n";
	return answer;
}

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findmaximum(MatDoub_I &pp,
		Datatype in_tolerance)
{	
	cout << "FINDMAXIMUM NOT YET IMPLIMENTED!\n";
	return findminimum(pp, in_tolerance);
}

template <class Datatype>
pair<Datatype,Datatype> Jcrap::multifunction<Datatype>::findminimum(int d,Datatype ax,Datatype bx,Datatype cx,Datatype tol)
{
	thisdim=d;
	return Jcrap::function<Datatype>::findminimum(ax, bx, cx, tol);
}

template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findminimum_dumb()
{
	return findminimum_dumb(a);
}


template <class Datatype>
pair<vector<Datatype>,Datatype> Jcrap::multifunction<Datatype>::findminimum_dumb(vector<Datatype> initial)
{
	static const Datatype range(2.);
	a=initial;
	pair<vector<double>, double> answer;
	answer.first.resize(a.size());
	
//	cout << "Welcome to findminimum_dumb h.size()=" << h.size() << " a.size=" << a.size();
//	cout << " param_names.size()=" << param_names.size() << " thisdim=" << thisdim << "\n";
//	for (int i(0);i<h.size();i++) cout << "h["<<i<<"]="<<h.at(i)<<"\n";
//	cout << "\n";
	
	cube checkcube;
	for (int i(0);i<50;i++) {
		for (int d(0);d<int(a.size());d++) {
			pair<double, double> thisdanswer;
/*			thisdanswer = findminimum(thisdim,-range*h.at(thisdim), -h.at(thisdim), 0., 1.e-6);
			checkcube = (plot(-range*h.at(thisdim), 0., 1000));
			checkcube.keyword("title", param_names.at(thisdim));
			checkcube.graph();			
			cout << "iteration " << i << " for thisdim(" << param_names.at(thisdim) << ", min at " << thisdanswer.first << ", " << thisdanswer.second << "\n";
*/
			
			if (d!=thisdim) {
				thisdanswer = findminimum(d,-range*fabs(h.at(d))+a.at(d), a.at(d), +range*fabs(h.at(d))+a.at(d), 1.e-6);
//				checkcube = (plot(-range*fabs(h.at(d))+a.at(d), +range*fabs(h.at(d))+a.at(d), 1000));
//				checkcube.keyword("title", param_names.at(d));
//				checkcube.graph();
//				cout << "iteration " << i << ", " << param_names.at(d) << ", min at " << thisdanswer.first << ", " << thisdanswer.second << "\n";
		
			
				answer.first.at(d) = thisdanswer.first;
				answer.second = thisdanswer.second;
				a = answer.first;
//				for (int j(0);j<a.size();j++) {
//					cout << param_names.at(j) << "=" << a.at(j) << ", ";
//				}
//				cout << "\n";				
//				char c;
//				cin >> c;	
			}
		}
//		char c; cin>>c;
	}
	
	return answer;
}
