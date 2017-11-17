
int addorial(int d)
{
	int answer(0);
	for (int i=0;i<=d;i++)
		answer += i;
	return answer;
}

template <class Datatype>
Jcrap::nDmultifunction<Datatype>::nDmultifunction(int d, int na) :
		nDfunction<Datatype>(d)
{
//	cout << "instantiating multifunction("<<na<<")\n";
	a.resize(na);
	param_names.resize(na);
	ia.resize(na);
	dyda.resize(na);
	
	covar.resize(na+1);
	for (int n=0;n<covar.size();n++) covar[n].resize(na+1);
	alpha.resize(na+1);
	for (int n=0;n<alpha.size();n++) alpha[n].resize(na+1);
	
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
//	cout << "done\n";
}


template <class Datatype>
Jcrap::nDmultifunction<Datatype>::nDmultifunction(int d, vector<Datatype> aa) :
	nDfunction<Datatype>(d)
{
	*this=multifunction(aa.size());
	a=aa;
}


template <class Datatype>
vector<Datatype> Jcrap::nDmultifunction<Datatype>::the_derivitives(
		vector<Datatype> x,
		vector<Datatype>& aa)
{
	vector<Datatype> answer(aa.size());
	vector<Datatype> h(aa.size());
	
	
	for (int i=0;i<aa.size();i++) {
		h[i] = 0.1;		// this is handled in aderiv
//		cout << "aderiv["<<i<<"]=";
		answer[i] = aderiv(i,x,aa,h[i]);
//		cout << answer[i]<<" ";
	}
	return answer;
}



#define NRANSI
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

template <class Datatype>
Datatype Jcrap::nDmultifunction<Datatype>::aderiv(int ii, vector<Datatype> x, 
	vector<Datatype> av, Datatype h)
{
	Datatype err;
	Datatype ans(0.);
	
	if (h == -1.) {
		h = a[ii]/10.;
	}
	if (h ==  0.) h = 1.;
	int i,j;
	Datatype errt,fac,hh,plus,minus;
	
	vector<vector<Datatype> > aa(NTAB+1);
	for (i=0;i<=NTAB;i++) aa[i].resize(NTAB+1);
	
	hh=h;
	av[ii] += hh;
//	cout << "before plus\n";
//	cout << "x = " << x << "\n";
//	for (int v=0;v<av.size();v++) cout << "av[" << v << "] = " << av[v] << "\n";
	plus = the_function(x,av);
	av[ii] -= 2.*hh;
//	cout << "after minus\n";
	minus = the_function(x,av);
//	cout << "plus= " << plus << " minus = " << minus << ", hh = " << hh << "\n";
	aa[1][1]=(plus-minus)/(2.0*hh);
	av[ii] += hh;
	err=BIG;
	for (i=2;i<=NTAB;i++) {
//		cout << "aa[1]["<<i-1<<"]="<<aa[1][i-1]<<" ";
//		cout << "for ii=" << ii << "slope guess = " << aa[1][i-1] <<"\n";
		hh /= CON;
		av[ii] += hh;
		plus = the_function(x,av);
		av[ii] -= 2.*hh;
		minus = the_function(x,av);
//		cout << "plus= " << plus << " minus = " << minus << ", hh = " << hh << "\n";
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
		if (jabs(aa[i][i]-aa[i-1][i-1]) >= SAFE*(err)) break;
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
pair<vector<Datatype>, Datatype> Jcrap::nDmultifunction<Datatype>::fit(
		const cube& data, const cube& sigma, vector<Datatype> aa, 
		vector<int> iaa)
/* Given a cube of size (2,1,n) where x=0 contains data and x=1 contains sigmas,
	or just (1,1,n) in which case the sigmas are assumed to be uniformly 1., fit
	to the_function based on the vector aa as initial guesses for the parameters */
{
	ia=iaa;
/*	if (sigma.N() <= 1)
		sigma = data.pow(0.5);*/

	int degreesoffreedom=0;
	for (int i=0;i<iaa.size();i++) {
		if (iaa[i]) degreesoffreedom++;
		if (!osuppress) cout << "iaa[" << i << "] = " << iaa[i] << "\n";
	}
		
//	cout << "First run, alamd=-1\n";
	pair<vector<Datatype>, Datatype> answer;	
	Datatype chisq=1.e98, oldchisq=1.e99, oldalamda;
	alamda = -1.;
//	chisq = Jmrqmin(vector_x, vector_y, vector_sigma, aa, iaa);
//	cout << "new reduced chisq = " << chisq/vector_x.size() << "\n";
//	for (int i=0;i<aa.size();i++)
//		cout << " aa[" << i << "] = " << aa[i];// << "\n";
//	cout << "\t\talamda = " << alamda << "\n";

	while ((jabs(chisq-oldchisq)/data.N() > chidiffstop  ||  oldalamda<alamda)  &&  alamda < pow(lamda_up,upperstop)) {
		if (!osuppress) { cout << "Step:"; cout.flush(); }
		oldchisq = chisq;
		oldalamda = alamda;
		chisq = Jmrqmin(data, sigma, aa, iaa);
		if (!osuppress) {
			cout << "new reduced chisq = " << chisq/(data.N()-degreesoffreedom);
			for (int i=0;i<aa.size();i++) {
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
	chisq = Jmrqmin(data, sigma, aa, iaa);
	
	return pair<vector<Datatype>, Datatype> (aa, chisq);
}

template <class Datatype>
pair<vector<Datatype>, Datatype> Jcrap::nDmultifunction<Datatype>::fit(
		const cube& data, const cube& sigma, vector<Datatype> aa)
{
	vector<int> iaa(aa.size());
	for (int i=0;i<aa.size();i++)
		iaa[i] = 1;
	ia=iaa;
	if (!osuppress) cout << "Instantiating fit(,,)\n";
	return fit(data, sigma, aa, iaa);		
}

template <class Datatype>
pair<vector<Datatype>, Datatype> Jcrap::nDmultifunction<Datatype>::fit(
		const cube& data, const cube& sigma)
{
	return fit(data,sigma,a);
}




template <class Datatype>
pair<vector<Datatype>, Datatype> Jcrap::nDmultifunction<Datatype>::fit(
		const cube& data)
{
	return fit(data,cube(1,1,1),a);
}







#define NRANSI

template <class Datatype>
Datatype Jcrap::nDmultifunction<Datatype>::Jmrqmin( 
	const cube& data, const cube& sigma,
	vector<Datatype>& aa, vector<int> &ia)
{
	int ndata(data.N());					//	(blank comment = added by JB)
	int ma(aa.size());						//
	static Datatype chisq;						//
	int j,k,l;
	static int mfit;
	static float ochisq;
	static vector<vector<Datatype> > oneda;		//
	static vector<Datatype> atry, beta, da;		//

	if (alamda < 0.0) {
		atry.resize(ma);
		beta.resize(ma+1);
		da.resize(ma+1);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j-1]) mfit++;
		oneda.resize(mfit+1);
		for (int i=0;i<=mfit;i++) oneda[i].resize(mfit+1);
		alamda= initialalamda;
		chisq = Jmrqcof(data,sigma,aa,ia,ma,alpha,beta);
		ochisq=chisq;
		if (!osuppress) cout << "First run:  ";
		for (j=0;j<ma;j++) {
			atry[j]=aa[j];
			if (!osuppress) cout << " ";
			if (ia[j]) if (!osuppress) cout << "[";
			if (!osuppress) cout << param_names[j] << "=" << aa[j];
			if (ia[j]) if (!osuppress) cout << "]";
		}
		if (!osuppress) cout << "\n";
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
	
	chisq = Jmrqcof(data,sigma,atry,ia,ma,covar,da);
	if (chisq < ochisq) {
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
	} else {
		alamda *= lamda_up;
		chisq=ochisq;
	}
	
	return chisq;
}
#undef NRANSI

		

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

template <class Datatype>
void Jcrap::nDmultifunction<Datatype>::Jcovsrt(
		vector<vector<Datatype> >& covar, int ma, vector<int> &ia, int mfit)
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
void Jcrap::nDmultifunction<Datatype>::Jgaussj(vector<vector<Datatype> > &a, int n, vector<vector<Datatype> > &b, int m)
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
					} else if (ipiv[k] > 1) cout << "Jgaussj: Singular Matrix-1\n";
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
		if (a[icol][icol] == 0.0) cout << "Jgaussj: Singular Matrix-2\n";
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
Datatype Jcrap::nDmultifunction<Datatype>::Jmrqcof(const cube& data,
		const cube &sigma, vector<Datatype> &aa, vector<int> &ia, int ma, 
		vector<vector<Datatype> > &alpha, vector<Datatype> &beta)
{
	Datatype chisq;
	int ndata = data.N();
	int i,j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy;

	vector<Datatype> dyda(ma);
	for (j=1;j<=ma;j++)
		if (ia[j-1]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	chisq=0.0;
	if (!osuppress) cout << "Jmrqcof progress:  00%";
	for (int x=0;x<data.N(X);x++) {
		for (int y=0;y<data.N(Y);y++) {
			for (int z=0;z<data.N(Z);z++) {
				int i(z*data.N(Y)*data.N(X) + y*data.N(X) + x);
				if (!osuppress) printpercent(i-1,ndata-1);
				vector<double> xvec;
				if (data.N(X)>1  &&  data.N(Y)>1  &&  data.N(Z)>1) {
					xvec.resize(3);
					xvec[0] = data.Axis(X,x);
					xvec[1] = data.Axis(Y,y);
					xvec[3] = data.Axis(Z,z);
				} else if (data.N(X)>1  &&  data.N(Y)>1) {
					xvec.resize(2);
					xvec[0] = data.Axis(X,x);
					xvec[1] = data.Axis(Y,y);
				} else {
					xvec.resize(1);
					xvec[0] = data.Axis(X,x);
				}
				ymod = funcs(xvec,aa,dyda);
//				cout << "after funcs, aa[8] = " << aa[8] << "\n";
				sig2i=1.0/(sigma(x,y,z)*sigma(x,y,z));
				dy=data(x,y,z)-ymod;
				for (j=0,l=1;l<=ma;l++) {
					if (ia[l-1]) {
						wt=dyda[l-1]*sig2i;
						for (j++,k=0,m=1;m<=l;m++)
							if (ia[m-1]) alpha[j][++k] += wt*dyda[m-1];
						beta[j] += dy*wt;
					}
				}
//				cout << "chisq before = " << chisq << "\n";
//				cout << "dy = " << dy << ", sig2i= " << sig2i << "\n";
				chisq += dy*dy*sig2i;
//				cout << "chisq after = " << chisq << "\n";
			}
		}
	}
	if (!osuppress) printpercent(ndata, ndata);
	if (!osuppress) cout << "\n";
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	
	return chisq;
}
#undef NRANSI


		
template <class Datatype>
Datatype Jcrap::nDmultifunction<Datatype>::funcs(vector<Datatype> x, vector<Datatype> &a, vector<Datatype> &dyda)
{
//	cout << "funcs\n";
	dyda = the_derivitives(x, a);
	return the_function(x, a);
}

template <class Datatype>
Datatype Jcrap::nDmultifunction<Datatype>::aconvert(vector<Datatype>& a, vector<int>& ia, 
		int i)
// Override this if your a[] are functions of one another.
{
//	cout << "in multifunction::aconvert\n";
	return a[i];
}
