#include "Jcoord.h"

coord::coord()
{
	x=0; y=0; z=0;
}

coord::coord(int a, int b, int c)
{
	x=a; y=b; z=c;
}

int& coord::A()
{ 
	if (a==X) return x; 
	if (a==Y) return y;
	if (a==Z) return z;
}
