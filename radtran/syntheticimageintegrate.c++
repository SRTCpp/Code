#include "../Jcrap.h"

main (int argn, char **argv)
{
	
	cube image(argv[1]);
	
	int x0(0),xf(image.N(X));
	int y0(0),yf(image.N(Y));
	

	if (argv[2]) {
		x0=1;
		xf=image.N(X)-1;
		y0=1;
		yf=image.N(Y)-1;	
	}
	
	for (int z(0); z<image.N(Z);z++){
		double sum(0.);
		
		for (int x(x0); x<xf; x++){
			for (int y(x0); y<yf; y++){
				
				sum+= image(x,y,z);
							
			}	
		}
		cout<<image.Axis(Z,z)<<"\t"<<sum<<endl;
		
	}
	
}
