#include <stdio.h>
#include "../Jcube.h"

#define pi 3.14159265
#define AMU 1
#define Z 1

main()
{
	cube P(1,1,100,0.0);
	float Re, rho;
	
	Ma = 1.67*pow(10,-24) * AMU;
	
	
	for (int z=0;z<100;z++){
		rho = .01 * expf(z/10.);
		Zaxis(z) = rho;
		Re = pow(Ma*3/(4*pi*rho, 1./3.);
		P(1,1,z) = 51.5/pow(Re,5)*(1 - Re*(.208 + .409*pow(Z, 2./3.)));
	}
	P.graph();
}
