#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define BOLTZMANN 8.6173324*0.00001	// in eV/K

double GetPressure(double temperature, double volume, double virial, int nParticles)
{
	double pressure = 0;

	pressure = (nParticles * BOLTZMANN * temperature + virial)/volume;

	return pressure;
}

double GetInstantTemperature(double vel [][3], int nParticles, double mass)
{
	double temperature = 0;
	int i,j;

	double factor = 2/(3*nParticles*BOLTZMANN);

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

