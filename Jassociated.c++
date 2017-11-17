#include "Jassociated.h"

namespace Jcrap {

cube addlatlonlines(const cube& incube, float spacing, float width, cube color)
// returns a cube with latitude/longitude lines on it.  Assumes input is a
// cylindrical map.  spacing is the distance between lines, width is the
// width of the lines, and color is what to replace those areas with.
{
	cube answer(incube);
	for (int x(0);x<answer.N(X);x++)
		for (double v(-360.);v<=360.;v+=spacing)
			if (fabs(answer.Axis(X,x)-v) < width/2.)
				for (int y(0);y<answer.N(Y);y++)
					for (int z(0);z<color.N(Z);z++)
						answer(x,y,z) = color(0,0,z);
		
	for (int y(0);y<answer.N(Y);y++)
		for (double v(-360.);v<=360.;v+=spacing)
			if (fabs(answer.Axis(Y,y)-v) < width/2.)
				for (int x(0);x<answer.N(X);x++)
					for (int z(0);z<color.N(Z);z++)
						answer(x,y,z) = color(0,0,z);
	
	return answer;
		
}

}
