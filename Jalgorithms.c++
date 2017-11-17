#include "Jcrap.h"


// as found at:  http://groups.yahoo.com/group/cpptips/message/259
template <typename T>
 struct LessByPointedToValue
 : std::binary_function<T *, T *, bool>
 {
 bool operator()(T * x, T * y)
 { return (*x) < (*y); }
 };
 
pair<cube, vector<double> > cube::findprincipalcomponents() const
{	
	cube nonzero(N(X)*N(Y),1,N(Z));  // pull out just the pixels that are nonzero
	cube nzsum(N(X)*N(Y),1,1);       // for faster dupe ID
	
	cout << "Findprincipalcomponents, weeding out zeroes & dupes:    00%";
	int i=0;
	for (int x=0;x<N(X);x++) {
		printpercent(x,N(X)-1);
		for (int y=0;y<N(Y);y++) {
			// is this spectrum a repeat?
			bool dupe(0);
			float thissum((*this).skewer(Z,x,y).sum());
			for (int nzx=0;nzx<i;nzx++)
				if (thissum == nzsum(nzx,0,0)) {
					dupe = 1;
				//	cout << "Dupe"; cout.flush();
					break;
				}
			
			if (!dupe)
				for (int z=0;z<N(Z);z++) {
					if ((*this)(x,y,z) != 0) {  // copy it in only if its not all zero
						for (int j=0;j<N(Z);j++)
							nonzero(i,0,j) = (*this)(x,y,j);
						nzsum(i,0,0) = nonzero.skewer(Z,i,0).sum();
						i++;
						z=N(Z);
					}
				}
		}
	}
	printpercent(1,1);
	cout << "\n";
	nonzero = nonzero.chunkx(0,i-1);
	
// following steps from:
// http://en.wikipedia.org/wiki/Principal_components_analysis#Algorithm_details
 	 
// Step 1:  Organize your data into column vectors, so you end up with a d x n
//          matrix, D.
// Okay, well, instead, I'm going to leave things inside the cube.  This is
// effectively, then, a 1 x n x d matrix, where n = N(X)*N(Y).
	
// Step 2:  Find the empirical mean along each dimension, so you end up with a 
//          d x 1 empirical mean vector, M.
// I'll make it a 1 x 1 x d cube instead

	cube average(1,1,N(Z));
	osuppress++;
	for (int i=0;i<N(Z);i++)
		average(0,0,i) = nonzero.plane(Z,i).mean();
	osuppress--;
	
// Step 3:  Subtract the empirical mean vector M from each column of the data
//				matrix D. Store mean-subtracted data d x n matrix in S.
	
	cube subbed(nonzero-average);
	
// Step 4:  Find the empirical covariance d x d matrix C of S. C = S . S^T
// Maybe this should be a separate function, but I'll leave it embedded for now.
	
	cube covariance(N(Z), N(Z), 1, 0.);
	cout << " Calculating covariance:   00%";
	for (int cx=0;cx<covariance.N(X);cx++) {
		printpercent(cx,covariance.N(X)-1);
		for (int cy=cx;cy<covariance.N(Y);cy++) {
			double thiscov=0.;
			for (int ix=0;ix<subbed.N(X);ix++)
				for (int iy=0;iy<subbed.N(Y);iy++)
					thiscov += subbed(ix,iy,cx)*subbed(ix,iy,cy);
			thiscov /= N(X)*N(Y)-1.;
			covariance(cy,cx,0) = covariance(cx,cy,0) = thiscov;
		}
	}
	printpercent(1,1);
	cout << "\n";
		
//	cout << "\nCovariance Matrix:  " << covariance;
	
// Step 5:  Compute and sort by decreasing eigenvalue, the eigenvectors V of C.
// okay, now it gets tricky.
	
	pair<cube, vector<double> > eigen;
	eigen = covariance.Jjacobi();
	for (int z=0;z<eigen.first.N(X);z++)
		eigen.first.Axis(X,z) = Axis(Z,z);
//	cout << "Eigenvectors:  " << eigen.first;
//	cout << "Eigenvalues:  ";
//	for (int i=0;i<eigen.second.size();i++) cout << "  " << eigen.second[i];
//	cout << "\n";
	
	
	return eigen;
}


cube cube::principalcomponentprojection(cube eigenvectors) const
{
	cube zeromask(N(X),N(Y),1,0.);
	for (int x=0;x<N(X);x++)
		for (int y=0;y<N(Y);y++)
			for (int z=0;z<N(Z) && zeromask(x,y,0)==0.;z++)
				if ((*this)(x,y,z) != 0.) zeromask(x,y,0) = 1.;
	
//	cout << "Eigenvectors cube in pcprojection:  " << eigenvectors << "\n";
	cube answer(N(X),N(Y),N(Z),0.);	
	
	cube average(1,1,N(Z));
	for (int i=0;i<N(Z);i++)
		average(0,0,i) = plane(Z,i).sum()/zeromask.sum();

	cube subbed=(*this)-average;
	
	cout << "Principal Components Projecting:   00%";
	for (int x=0;x<N(X);x++) {
		printpercent(x,N(X)-1);
		for (int y=0;y<N(Y);y++)
			for (int z=0;z<N(Z);z++)  // subscript for data vector
				for (int n=0;n<N(Z);n++)  // second subscript for eigenmatrix
					answer(x,y,z) += eigenvectors(n,z,0)*subbed(x,y,n);
	}
	printpercent(1, 1); cout << "\n";
	
	return answer*zeromask;
}

cube cube::principalcomponentdeprojection(cube eigenvectors) const
{
	
//	cout << "Eigenvectors cube in pcprojection:  " << eigenvectors << "\n";
	cube answer((*this)*0.);
	
	cube eigeninverse(eigenvectors.matrixinvert());
	
	cout << "Principal Components Deprojecting:   00%";
	for (int x=0;x<N(X);x++) {
		printpercent(x,N(X)-1);
		for (int y=0;y<N(Y);y++)
			for (int z=0;z<N(Z);z++)  // subscript for data vector
				for (int n=0;n<N(Z);n++)  // second subscript for eigenmatrix
					answer(x,y,z) += eigeninverse(n,z,0)*(*this)(x,y,n);
	}
	printpercent(1, 1); cout << "\n";
	
	return answer;
}


// heirarchical cluster analysis algorithm
Jcrap::pixelnode cube::classify(methodtype mT)
// created 2005 Sept 29
// based on http://en.wikipedia.org/wiki/Cluster_analysis
{
	list<pixelnode*> ends;
	cout << "Classifying a cube with info:  " << info();
	cout << "Creating pixelnodes --   00%";
	for (int x=0;x<N(X);x++) {
		printpercent(x,N(X)-1);
		for (int y=0;y<N(Y);y++) 
			for (int z=0;z<N(Z);z++) 
				if ((*this)(x,y,z) != 0.) {
					pixelnode *a;
					cube spec(skewer(Z,x,y));
					if (mT == mFAST)
						a = new pixelnode(spec,x,y);
					if (mT == mSLOW)
						a = new pixelnode(*this,x,y,mT);
					ends.push_back(a);
					z=N(Z);
				}
	} 
	printpercent(1,1);
	cout << "\nsorting . . .    00%"; cout.flush();
	ends.sort(LessByPointedToValue<pixelnode>());
	long progress=0, hiprogress=0;
 	cout << " done\n";
		
// consolidate all indentical nodes
	progress=1;
	if (!osuppress) cout << "Consolidating identical nodes:   00%";
	cout.flush();
	long endssize(ends.size()-2);
	for (list<pixelnode*>::iterator i=++(ends.begin());i!=ends.end();i++,progress++) {
		if (!osuppress) printpercent(progress/100,endssize/100);
		list<pixelnode*>::iterator j=i;
		j++;
		if (j == ends.end()) continue;
		while ((**j) == (**i)) {
			list<pixelnode*>::iterator oldj=j;
			j++;
			if (j == ends.end()) break;
			list<pixelnode*>::iterator oldi=i;
			pixelnode* fusion=new pixelnode(*oldi,*oldj);
			i = ends.insert(i, fusion);
			ends.erase(oldi);
			ends.erase(oldj);
			
			progress++;
		}
		if (j == ends.end()) break;
		if (j == ends.end()) break;
	}
	if (!osuppress) printpercent(1,1);
	if (!osuppress) cout << "\n";
	if (!osuppress) cout << "Left with " << ends.size() << " unique nodes.\n";
		
// now that the number of nodes is sufficiently decreased, go back to mFAST
	if (mT == mSLOW) {
		cout << "Converting . . .    00%";
		int n=0,tot=ends.size()-1;
		for (list<pixelnode*>::iterator i=ends.begin();i!=ends.end();i++,n++) {
			printpercent(n,tot);
			(**i).convert();
		}
		cout << "  done.\n";
	}
	
// calculate all initial distances
	progress=0;
	cout << "assigning initial distances:   00%";
	cout.flush();
	list<pixelnode*>::iterator penultimate = ends.begin();
	for (int i=0;i!=ends.size()-1;i++) penultimate++;
	for (list<pixelnode*>::iterator i=ends.begin();i!=penultimate;i++,progress++) {
		if (!osuppress) printpercent(progress,ends.size()-2);
		
		list<pixelnode*>::iterator j;
		int maxn(int(sqrt(ends.size())*20));
		for (int n=0;n<maxn&&j!=ends.end();n++,j++) {
			if (n==0)
				j=i;
			else {
				double thisdist((*i)->distance(**j));
				if (thisdist < (*i)->closest.first) {
					(*i)->closest.first = thisdist;
					(*i)->closest.second= *j;
					if (thisdist < (*j)->closest.first) {
						(*j)->closest.first = thisdist;
						(*j)->closest.second= *i;
					}
				}
			}					
		}
	}	
	if (!osuppress) printpercent(1,1);
	if (!osuppress) cout << "\n";
	
// for debugging only:  output distances and nearest neighbors
//	int count=0;
//	for (list<pixelnode*>::iterator i=ends.begin();i!=ends.end();i++,count++) {
//		cout << count << "  : d=" << (*i)->closest.first << "\n";
//	}
	
// cluster until just one node remains:  this must be the root.
	cout << "Now clustering.\n";
	osuppress++;
	while (ends.size() != 1) {
		double closestdist;
		list<pixelnode*>::iterator closest1;
		list<pixelnode*>::iterator closest2;
		closest1=list<pixelnode*>::iterator(0);
		closest2=list<pixelnode*>::iterator(0);
		closestdist=1.e99;
		int n=0,m=0;
		for (list<pixelnode*>::iterator i=ends.begin();i!=ends.end();i++,n++) {
//			cout << "Checking for lowest -- n=" << n << "\n";
			if ((*i)->closest.first <= closestdist) {
//				cout << "  replacing " << m << ", " << closestdist << "  with ";
//				cout << n << ", " << (*i)->closest.first << "\n";
				m=n;
				if ((*i)->closest.first != closestdist) {
//					cout << "assigning closest1\n";
					closest1 = i;
				}
				else {
//					cout << "assigning closest2\n";
					closest2 = i;
				}
				
				closestdist = (*i)->closest.first;
			}
		}
		if (closest1!=list<pixelnode*>::iterator(0)  &&  
		    closest2==list<pixelnode*>::iterator(0)) {
			for (list<pixelnode*>::iterator i=ends.begin();i!=ends.end();i++) {
				if ((*i) == (*closest1)->closest.second) {
//					cout << "assigning closest2 -- backup\n";
					closest2 = i;
				}
			}
		}
//		cout << "Lowest must be #" << m << "\n";
//		if (closest1 == 0) cout << "closest1 == 0\n";
//		cout << "closest1 points to " << *closest1 << "\n";
//		cout << "closest2 points to " << *closest2 << "\n";
		
		cout << "Merging a " <<(**closest1).children();
		cout << " (#" << m << ") with a ";
		cout << (**closest2).children() << " -- "<< ends.size() << " left.\n";
		pixelnode *newparent = new pixelnode(*closest1, *closest2);
		
		int b=0;
		for (list<pixelnode*>::iterator j(ends.begin());j!=ends.end();j++,b++)
			if (closest1 == j) break;
//		cout << "Erasing #" << i << "\n";
		ends.erase(closest1);
		
		b=0;
		for (list<pixelnode*>::iterator j(ends.begin());j!=ends.end();j++,b++)
			if (closest2 == j) break;
//		cout << "Erasing #" << i << "\n";		
		ends.erase(closest2);
		
// uhoh -- now redo the closest search for those whose closest was just deleted
		int o=0;
		for (list<pixelnode*>::iterator i=ends.begin();i!=ends.end();i++,o++) {
			if ((*i)->closest.second == *closest1   ||
				 (*i)->closest.second == *closest2) {
//				cout << "Will have to recalculate closest for #"<<o<<"\n";
				(*i)->closest.first = 1.e99;
				for (list<pixelnode*>::iterator j=ends.begin();j!=ends.end();j++) {
					if (j == i) continue;
					double thisdist((*i)->distance(**j));
					if (thisdist < (*i)->closest.first) {
						(*i)->closest.first = thisdist;
						(*i)->closest.second= *j;
					}
				}
//				cout << "new closest for #" << o << " low = " << (*i)->closest.first
//						<< "\n";
			}
		}

//		place new parent in order
//		ends.push_back(newparent);
		if (ends.size()==0)
			ends.push_front(newparent);
		else if ((*newparent).sum > (*(ends.rbegin()))->sum)
			ends.push_back(newparent);
		else if ((*newparent).sum < (*(ends.begin()))->sum)
			ends.push_back(newparent);
		else {
			list<pixelnode*>::iterator k=ends.begin();
			while ((**k).sum < (*newparent).sum) k++;
			ends.insert(k, newparent);
		}
		list<pixelnode*>::iterator newparentposition;
		for (list<pixelnode*>::iterator i(ends.begin());i!=ends.end();i++) {
			if (*i == newparent) {
				newparentposition = i;
				break;
			}
		}
	
// calculate the new closest distance for the new fusion pixelnode
		for (list<pixelnode*>::iterator i(ends.begin());i!=ends.end();i++) {
			if (i == newparentposition) continue;
			double thisdist((**newparentposition).distance(**i));
			if (thisdist < (**newparentposition).closest.first) {
				(**newparentposition).closest.first = thisdist;
				(**newparentposition).closest.second= *i;
			}
			if (thisdist < (*i)->closest.first) {
				(*i)->closest.first = thisdist;
				(*i)->closest.second= *newparentposition;
			}
		}
	}
	osuppress--;
	
	return *(ends.front());
}

pair<cube,cube>& cube::pixelclassify(Jcrap::pixelnode p, cube palette, 
		unsigned int level)
{
	static pair<cube, cube> answers;
	static cube &answer(answers.first);
	static cube &spectra(answers.second);
	static cube colormap;
	if (answer.N(X)!=N(X) || answer.N(Y)!=N(Y) || answer.N(Z)!=3) {
		answer = cube(N(X),N(Y),3,0.).totif();
		cout << "Blanking pixelclassify answer cube\n";
		spectra = cube(0,1,p.centroid().N(Z));
		for (int z=0;z<spectra.N(Z);z++)
			spectra.Axis(Z,z) = p.centroid().Axis(Z,z);
		colormap = cube(0,1,3).totif();
	}
	
//  if this is a leaf node, color it in
	if (p.children() == 1) {
		int x(p.x()), y(p.y());
		answer(x,y,0) = palette(0,0,0);
		answer(x,y,1) = palette(0,0,1);
		answer(x,y,2) = palette(0,0,2);
	} else if (!(p.unispectral())) {
			pixelclassify(p.left(), palette, level);
			pixelclassify(p.right(), palette, level);
	} else { 
// otherwise, recurse, keeping track of the palette.
		if (level /*&& !(p.left().unispectral()) && !(p.right().unispectral())*/) {
//			cout << "Going left to a size " << p.left().children() << "\n";
			if (level == 0)  {
				spectra=spectra.blocksx(p.left().centroid());
				colormap = colormap.blocksx(palette.skewer(Z,0,0));
				colormap.totif().write("colormap.tif", 0., 1.);
			}
			pixelclassify(p.left(), palette(0,palette.N(X)/2-1,-1,-1,-1,-1), level);
			palette = palette(palette.N(X)/2,-1,-1,-1,-1,-1);
//			cout << "Going right to a size " << p.right().children() << "\n";
			if (level == 0) {
				spectra=spectra.blocksx(p.right().centroid());
				colormap = colormap.blocksx(palette.skewer(Z,0,0));
				colormap.totif().write("colormap.tif", 0., 1.);
			}
			pixelclassify(p.right(),palette, level);
		} else {
			pixelclassify(p.left(), palette, level);
			pixelclassify(p.right(), palette, level);
		}
	}
	return answers;
}

pair<cube,cube>& cube::pixelclassify(Jcrap::pixelnode root, cube& colorcube)
{
	static pair<cube, cube> answers;
	cube &answer(answers.first);
	if (answer.N(X) == 1)
		answer = cube(N(X), N(Y), 3, 0.).totif();
	cube &spectra(answers.second);
	if (spectra.N(Z) == 1)
		spectra = cube(colorcube.N(X),colorcube.N(Y),root.centroid().N(Z), 0.);
	
	
	for (int d=0;d<colorcube.N(Y);d++) {
		for (int e=0;e<powf(2,d);e++) {
			cout << "Pixelclassifying " << d << ", " << e << ": ";
			if (colorcube(e,d,0)==1.0 &&
				 colorcube(e,d,1)==1.0 &&
				 colorcube(e,d,2)==1.0) {
				cout << " white\n";
				continue;
			}
			if (colorcube(e,d,0)==0.25 &&
				 colorcube(e,d,1)==0.25 &&
				 colorcube(e,d,2)==0.25) {
				cout << " gray\n";
				continue;
			}
		
		// looks like we have to color it in, then.  Figure out the uppermost
		// branching point to color in this color.
			Jcrap::pixelnode *p=&root;
			for (int j=int(powf(2,d-1));j!=0;j/=2){
				cout << "Checking " << e << "&" << j;
				if (p->children()==1) {cout << " 0 "; p=p;}
				else if (e&j) { cout << " r "; p=&(p->right()); }
				else {cout << " l "; p=&(p->left()); }
			}
			cout << "\n";
			
		// have the original pixelclassify color it in for us.
			static pair<cube, cube> *partial;
			partial = &(pixelclassify(*p, colorcube.skewer(Z,e,d), 1));
			for (int x=0;x<answer.N(X);x++)
				for (int y=0;y<answer.N(Y);y++) {
					if (partial->first(x,y,0) != 0.  ||
						 partial->first(x,y,1) != 0.  ||
						 partial->first(x,y,2) != 0.)
						for (int z=0;z<answer.N(Z);z++)
							answer(x,y,z) = partial->first(x,y,z);
				}
			for (int z=0;z<spectra.N(Z);z++)
				spectra(e,d,z) = p->centroid()(0,0,z);
		}
	}
	
	return answers;
}
	
bool pixelnode::unispectral()
{
	return (left().centroid() == right().centroid());
}
