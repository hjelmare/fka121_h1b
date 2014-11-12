#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double GetInstantTemperature(double vel [][3], int nParticles, double mass)
{
	double k = 1;		// boltzmann's constant (in metal units!!)
	double temperature = 0;
	int i,j;

	double factor = 2/(3*nParticles*k);

	for (i=0;i<nParticles;i++) {
		for (j=0;j<nParticles;j++) {
			temperature += factor * pow(vel[i][j],2)/mass;
		}
	}	

	return temperature;
}

// Calculates the kinetic energy
double GetKineticEnergy(double vel[][3], double mass, int nParticles)
{
	int i,j;
	double sum;
	double energy = 0;

	for(i = 0; i<nParticles; i++) {
		sum = 0;
		for(j = 0; j<3; j++) {
			sum +=vel[i][j]*vel[i][j];
		}
		energy += sum*mass/2;
	}
	return energy;
}

// Calculate alpha (the modifier to get the right temperature)
double GetAlpha(double wantedTemp, double currentTemp, double timestep, double timeConstant)
{
	double alpha;
	alpha = 1 + timestep/timeConstant * (wantedTemp-currentTemp)/currentTemp;

	return alpha;
}

