#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Calculates the kinetic energy
double CalculateEnergy(void vel[][3], m, int nParticles)
{
int i,j;
double sum;
double kinEnergy;
for(i = 0; i<nParticles; i++){
	sum = 0;
	for(j = 0; j<3; j++){
		sum += vel[i][j]*vel[i][j];
	}
	kinEnergy[i] = sqrt(sum)*m/2;
}


return(kinEnergy)
}
