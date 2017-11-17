#include "SRTC++.h"

#ifndef SRTCppTESTSUITE_h
#define SRTCppTESTSUITE_h 1

const double qualityfactor(0.3);
const otype outputtype(oFIG);

cube diskintegrate(cube, double=2575.); 

double test_LambertSurface();
double test_ChandraAtmo();
double test_SebastienCompare();
double test_GrahamCompare();
double test_DiffusiveReflectance(double); 
cube avgcenterval(cube);

#endif
 
