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
	double totalTime = 10;
	double timestep = 0.01;
		
	// physical parameters
	int dim = 3;		// there's a few functions that are hardcoded for dim=3, not sure how to fix...
	int nCells = 4;
	int nParticles = 4*pow(nCells,dim);
	int wantedTemp = 500+273;   // The temperature that we want the system to stabilize around.
	double wantedPressure = 6.32420934 * 0.0000001;	// The pressure that we want the system to stabilize around.
	double timeConstantT = 0.01;  // used in determining alpha. It's the constant that determines how fast our temperature will move towards the prefered temperature
	double timeConstantP = 0.002;
	double mass = 0.00279636;  // 26.9815 u
	double latticeParameter = 4.05;
	double maxDeviation = 0.05;
	double curtAlphaP;
	double virial;

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
	double sqrtAlphaP;
	double volume;

	// derived quantities
	int nSteps = (int) totalTime/timestep;
	double supercellLength = nCells * latticeParameter;

	// other stuff
	int i,j,k;   

	srand(time(NULL));
	double random_value;

	init_fcc(pos, 4, latticeParameter);

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
	
	potentialEnergy = get_energy_AL(pos, supercellLength, nParticles);
	kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
	energy = potentialEnergy + kineticEnergy;


	//Saving initial data:
	FILE *kineticEnergyFile;
	kineticEnergyFile = fopen("kineticEnergy.data","w");
	fprintf(kineticEnergyFile, "%e \t %e \n", 0.0, kineticEnergy );

	FILE *potentialEnergyFile;
	potentialEnergyFile = fopen("potentialEnergy.data","w");
	fprintf(potentialEnergyFile, "%e \t %e \n", 0.0, potentialEnergy);
	
	FILE *totEnergyFile;
	totEnergyFile = fopen("totEnergy.data","w");
	fprintf(totEnergyFile, "%e \t %e \n", 0.0, energy);
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
		get_forces_AL(force, pos, supercellLength, nParticles);

		// Update velocities again
		for (j=0; j<nParticles; j++) {
			for (k=0; k<dim; k++) {
				vel[j][k] = vel[j][k] + 0.5*timestep*force[j][k] / mass;
			}
		}

		potentialEnergy = get_energy_AL(pos, supercellLength, nParticles);
		kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
		energy = potentialEnergy + kineticEnergy;

		currentTemp = GetInstantTemperature(vel, nParticles, mass, dim);
		volume = pow(nCells*latticeParameter, 3);
		virial = get_virial_AL(pos, nCells*latticeParameter, nParticles);
		currentPressure = GetPressure(currentTemp, volume, virial, nParticles);

		//Calculate alpha (the velocity and position scaling parameters)
		alphaT = GetAlphaT(wantedTemp, currentTemp, timestep, timeConstantT); // This function calculates alpha that rescales our velocity at each timestep.
		alphaP = GetAlphaP(wantedPressure, currentPressure, timestep, timeConstantP);
		sqrtAlphaT = sqrt(alphaT);
		curtAlphaP = pow(alphaP, 1.0 / 3);

		// Rescale velocity and position (and the size of the bounding volume, right?)
		latticeParameter = latticeParameter * curtAlphaP;
//printf("lattpar: %e \n", latticeParameter);
		for (j=0; j<nParticles; j++) {
			for (k=0; k<dim; k++) {
				vel[j][k] = vel[j][k] * sqrtAlphaT;
				pos[j][k] = pos[j][k] * curtAlphaP;
			}
		}
if(i%3==0) {
  printf("alphaP = %e, lattpar = %e, P = %e \n", alphaP, latticeParameter, currentPressure);
}

if(i % 10 == 0){
//		printf("alphaT = %e \t alphaP= %e \t temp= %e \t pressure= %e \t vel[0][0]=  %e \n", alphaT, alphaP, currentTemp, currentPressure, vel[0][0]);
//		printf("vol= %e \t lattpar= %e \n", volume, latticeParameter);	// latticeParameter får helt fel värde på nåt jäkla vänster....
}

		fprintf(kineticEnergyFile, "%e \t %e \n", i*timestep, kineticEnergy );
		fprintf(potentialEnergyFile, "%e \t %e \n", i*timestep, potentialEnergy);
		fprintf(totEnergyFile, "%e \t %e \n", i*timestep, energy);
	}
printf("lattpar efter loopen: %e \n", latticeParameter);

    /*
     Descriptions of the different functions in the files initfcc.c and alpotential.c are listed below.
     */
    
    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of primitive cells in each direction and a0 is the lattice parameter. The positions of all the atoms are stored in pos which should be a matrix of the size N x 3, where N is the number of atoms. The first, second and third column correspond to the x,y and z coordinate respectively.
     */
    /*
     init_fcc(pos, Nc, a0);
    */
    
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be a matrix containing the positions of all the atoms, L is the length of the supercell and N is the number of atoms.
     */
    /*
     double energy;
     energy = get_energy_AL(pos, L, N);
     */
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix containing the positions of all the atoms, L is the length of the supercell and N is the number of atoms.
     */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
     */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the forces are stored in f which should be a matrix of size N x 3, where N is the number of atoms and column 1,2 and 3 correspond to the x,y and z component of the force resepctively . pos should be a matrix containing the positions of all the atoms, L is the length of the supercell and N is the number of atoms.
     */
    /*
     get_forces_AL(f,pos, L, N);
     */
}
