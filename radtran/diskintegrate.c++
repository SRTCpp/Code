#include "../Jcrap.h"

main(int argn, char **argv)
{
	cube incube(argv[1]);
	
	cube output(incube.skewer(Z,0,0));
	output *= 0.;
	
	for (int z(0);z<output.N(Z);z++) {
		int n(0);
		for (int x(0);x<incube.N(X);x++) {
			for (int y(0);y<incube.N(Y);y++) {
				double d(incube.Axis(X,x)*incube.Axis(X,x) + 
					incube.Axis(Y,y)*incube.Axis(Y,y));
				d = sqrt(d);
				if (d < 2575.) {
					n++;
					output(0,0,z) += incube(x,y,z);
				}
			}
		}
		output(0,0,z) /= double(n);
	}
	
	output.write("diskaverages.Jcube");
	cout << output;
}
