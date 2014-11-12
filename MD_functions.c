#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Calculates the kinetic energy
double GetKineticEnergy(double vel[][3], double mass, int nParticles)
{
	int i,j;
	double sum;
	double energy = 0;
	int dim = 3;

	for(i = 0; i<nParticles; i++) {
		sum = 0;
		for(j = 0; j<dim; j++) {
			sum += vel[i][j]*vel[i][j];
		}
		energy += sqrt(sum)*mass/2;
	}
	printf("e %e \n", energy);
	return energy;
}
