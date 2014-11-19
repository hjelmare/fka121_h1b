/*
 MD_main.c
 
 Created by AL on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"
#include "MD_functions.h"

/* Main program */
int main()
{
	// REMEMBER: \m/ METAL UNITS \m/
	// simulation settings
	double totalTime = 2000;
	double timestep = 0.01;
		
	// physical parameters
	int dim = 3;		// there's a few functions that are hardcoded for dim=3, not sure how to fix...
	int nCells = 4;
	int nParticles = 4*pow(nCells,dim);
	double wantedTemp = 500+273;   // The temperature that we want the system to stabilize around.
	double wantedPressure = 6.32420934 * 0.0000001;	// The pressure that we want the system to stabilize around.
//	wantedPressure = 5*0.00001;
	double timeConstantT = 0.02;  // used in determining alpha. It's the constant that determines how fast our temperature will move towards the prefered temperature
	double timeConstantP = 0.05;
	double mass = 0.00279636;  // 26.9815 u
	double latticeParameter = 4.05;
	double maxDeviation = 0.05;

	// storage of physical quantities
	double pos[nParticles][dim];
	double vel[nParticles][dim];
	double force[nParticles][dim];
	double energy, potentialEnergy, kineticEnergy;
	double square_root_of_alpha;
	double alphaT = 1;
	double alphaP = 1;
	double currentTemp;
	double currentPressure;
	double sqrtAlphaT;
	double curtAlphaP;
	double virial;
	double volume;

	// derived quantities
	double nSteps = totalTime/timestep;

	// other stuff
	int i,j,k;   

	srand(time(NULL));
	double random_value;

	init_fcc(pos, nCells, latticeParameter);

	for(i = 0; i<nParticles; i++){
		for (j = 0; j<3; j++){
			random_value = (double) rand() / (double) RAND_MAX;
			pos[i][j] = pos[i][j] + (random_value-0.5) * maxDeviation * latticeParameter;
 		}
	}

	// initialize velocities and forces to zero
	for (i = 0; i<nParticles; i++){
		for (j = 0; j<3; j++){
			vel[i][j] = 0;
			force[i][j] = 0;
		}
	}
	
	potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, nParticles);
	kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
	energy = potentialEnergy + kineticEnergy;


	//Saving initial data:
	FILE *energyFile;
	energyFile = fopen("energy.data","w");
	fprintf(energyFile, "%e \t %e \t %e \t %e \n", 0.0, energy, potentialEnergy, kineticEnergy);

	FILE *ptFile;
	ptFile = fopen("pt.data","w");

printf("lattpar innan loopen: %e \n", latticeParameter);
	for (i=1;i<nSteps;i++) {
		// Update velocities and positions
		for (j=0; j<nParticles; j++) {
			for (k = 0; k<dim; k++) {
				vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
				pos[j][k] = pos[j][k] + timestep*vel[j][k];
			}
		}
		
		// Calculate forces so that we can take the next half step
		get_forces_AL(force, pos, nCells*latticeParameter, nParticles);

		// Update velocities again
		for (j=0; j<nParticles; j++) {
			for (k=0; k<dim; k++) {
				vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
			}
		}

		potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, nParticles);
		kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
		energy = potentialEnergy + kineticEnergy;

		currentTemp = GetInstantTemperature(vel, nParticles, mass, dim);
		volume = pow(nCells*latticeParameter, 3);
		virial = get_virial_AL(pos, nCells*latticeParameter, nParticles);
		currentPressure = GetPressure(currentTemp, volume, virial, nParticles);

		//Calculate alpha (the velocity and position scaling parameters)
		alphaT = GetAlphaT(wantedTemp, currentTemp, timestep, timeConstantT);
		alphaP = GetAlphaP(wantedPressure, currentPressure, timestep, timeConstantP);
		sqrtAlphaT = sqrt(alphaT);
		curtAlphaP = pow(alphaP, 1.0 / 3);

		// Rescale velocity and position (and the size of the bounding volume, right?)
		latticeParameter = latticeParameter * curtAlphaP;
		for (j=0; j<nParticles; j++) {
			for (k=0; k<dim; k++) {
				vel[j][k] = vel[j][k] * sqrtAlphaT;
				pos[j][k] = pos[j][k] * curtAlphaP;
			}
		}
		
if(i % 5 == 0){
	printf("it %d T %e \t alphaT %e \t P %e \t alphaP %e \n", i, currentTemp, alphaT, currentPressure, alphaP);
}

		fprintf(energyFile, "%e \t %e \t %e \t %e\n", i*timestep, energy, potentialEnergy, kineticEnergy );
		fprintf(ptFile, "%e \t %e \t %e \n", i*timestep, currentPressure, currentTemp);
	}

printf("lattpar efter loopen: %e \n", latticeParameter);

FILE *donefile;
donefile = fopen("done.data", "w");
fprintf(donefile, "done");
close(donefile);


}
