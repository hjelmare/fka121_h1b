#ifndef _MD_functions_h
#define _MD_functions_h

extern double getDistanceSquared(double *, double *);
extern double ScalarProduct(double *, double *);
extern double GetPressure(double, double, double, int);
extern double GetInstantTemperature(double [][3], int, double, int); 
extern double GetKineticEnergy(double [][3], double, int);
extern double GetAlphaT(double, double, double, double);
extern double GetAlphaP(double, double, double, double);


#endif 
