#include "Jcrap.h"

 cube pixelnode::return1;
 cube pixelnode::return2;
 int pixelnode::returncount;


Jcrap::pixelnode::pixelnode(cube& data, int x, int y, methodtype mT)
{
	Left = 0;
	Right = 0;
	Parent = 0;
	Children = 1;
	xcoord = x;
	ycoord = y;
	method = mT;
	if (mT == mFAST) {
		Centroid = new cube;
		*Centroid = data;
		sum = data.sum();
	}
	else {
		Centroid = &data;
		sum = data.skewer(Z,x,y).sum();
	}
	
	closest = pair<double, pixelnode*>(1.e99, 0);
}

Jcrap::pixelnode::pixelnode(pixelnode *l, pixelnode *r)
{
	Left = l;
	Right = r;
	l->assignparent(this);
	r->assignparent(this);
	Parent = 0;
	Children = l->children() + r->children();
	xcoord = -1;
	ycoord = -1;
	method = l->method;
	
	if (*l == *r) {
		if (method == mFAST)
			Centroid = &(l->centroid());
		if (method == mSLOW) {
			xcoord=l->xcoord;
			ycoord=l->ycoord;
			Centroid = l->Centroid;
		}
	} else {
		Centroid = new cube();
		*Centroid = l->centroid()*float(l->children())
				+    r->centroid()*float(r->children());
		*Centroid = *Centroid / float(Children);	
		Centroid->headerwipe();
	}
	sum = centroid().sum();
	
	closest = pair<double, pixelnode*>(1.e99, 0);
}

Jcrap::pixelnode::~pixelnode()
{
//	if (method == mFAST  ||  Left != 0)
//		delete Centroid;
}

bool pixelnode::operator==(pixelnode& r)
{
	if (sum==r.sum) {
		return centroid()==r.centroid();
	}
	else return 0;
}

bool pixelnode::operator<(pixelnode& r)
{
	return sum<r.sum;
}

void Jcrap::pixelnode::assignparent(pixelnode* n)
{
	Parent = n;
}

int Jcrap::pixelnode::children() { return Children; }

cube& Jcrap::pixelnode::centroid() 
{ 
	if (method==mFAST)
		return *Centroid; 
	if (returncount == 0) {
		return1 = Centroid->skewer(Z,xcoord,ycoord);
		returncount = 1;
		return1.headerwipe();
		return return1;
	}
	if (returncount == 1) {
		return2 = Centroid->skewer(Z,xcoord,ycoord);
		returncount = 0;
		return2.headerwipe();
		return return2;
	}
}

double Jcrap::pixelnode::distance(pixelnode& other)
{
	double answer(0.);
	for (int z=0;z<Centroid->N(Z);z++) {
		double diff(centroid()(0,0,z)-other.centroid()(0,0,z));
		answer += diff*diff;
	}
	return sqrt(answer);
}

void Jcrap::pixelnode::convert() 
// converts pixelnodes from mSLOW to mFAST
{
	if (Centroid->N(X) > 1) {
		method = mFAST;
		cube* transfer = new cube(Centroid->skewer(Z,xcoord,ycoord));
		Centroid = transfer;
		Centroid->headerwipe();
	}
}

int Jcrap::pixelnode::x() { return xcoord; }
int Jcrap::pixelnode::y() { return ycoord; }

pixelnode& Jcrap::pixelnode::left() { return *Left; }
pixelnode& Jcrap::pixelnode::right(){ return *Right;}
