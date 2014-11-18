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
	double equilibrationTime = 5;
		
	// physical parameters
	int dim = 3;		// there's a few functions that are hardcoded for dim=3, not sure how to fix...
	int nCells = 4;
	int correlationDistance = 50;	// this constant determines how sepperated the points should be in the corr.function (for temperature)
	int nParticles = 4*pow(nCells,dim);
	int equilibrationSteps = equilibrationTime/timestep;
	double wantedTemp = 500+273;   // The temperature that we want the system to stabilize around.
	double wantedPressure = 6.32420934 * 0.0000001;	// The pressure that we want the system to stabilize around.
	wantedPressure = 5*0.00001;
	double timeConstantT = 0.02;  // used in determining alpha. It's the constant that determines how fast our temperature will move towards the prefered temperature
	double timeConstantP = 0.05;
	double mass = 0.00279636;  // 26.9815 u
	double latticeParameter = 4.05;
	double maxDeviation = 0.05;

	// storage of physical quantities
	double pos[nParticles][dim];
	double vel[nParticles][dim];
	double force[nParticles][dim];
	double savedValuesT[correlationDistance];
	double meanValuesT[correlationDistance];
	double savedValuesP[correlationDistance];
	double meanValuesP[correlationDistance];
	double energy, potentialEnergy, kineticEnergy;
	double siquare_root_of_alpha;
	double alphaT = 1;
	double alphaP = 1;
	double currentTemp;
	double currentPressure;
	double sqrtAlphaT;
	double curtAlphaP;
	double virial;
	double volume;
	double temp;
	double test=0;
	double meanTemp, meanPressure, meanSquareTemp, meanSquarePressure, meanTemp_ik[correlationDistance-1], meanPressure_ik[correlationDistance-1];
	double phiTemp[correlationDistance], phiPressure[correlationDistance], meanTempSquare, meanPressureSquare;
	

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



	for (i=1;i<equilibrationSteps;i++) {
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
	}

	for (i=equilibrationSteps;i<nSteps;i++) {
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


		meanTemp += currentTemp;
		meanPressure += currentPressure;
		meanSquareTemp += currentTemp * currentTemp;
		meanSquarePressure += currentPressure * currentPressure;
		//Saves temp and preassure values in order to calculate s.
		if (i < equilibrationSteps + correlationDistance){
		savedValuesT[i % correlationDistance] = currentTemp;
		savedValuesP[i % correlationDistance] = currentPressure;
		printf("\n i=%d \t sparat T = %e \n", i, savedValuesT[0]);		
		}else{  
			test += 1;
			// updates the saved values.
			for (j = 0; j<correlationDistance - 1; j++){	
			savedValuesT[j] = savedValuesT[j+1];
			savedValuesP[j] = savedValuesP[j+1];	
			}
			savedValuesT[correlationDistance - 1] = currentTemp;
			savedValuesP[correlationDistance - 1] = currentPressure;
			
			// calculates and saves f_i*f_k
			for (j = 1; j<correlationDistance; j++){
			meanTemp_ik[j-1] += savedValuesT[0] * savedValuesT[j];
			meanPressure_ik[j-1] += savedValuesP[0]* savedValuesP[j];
			}
		}
	}

meanTemp = meanTemp/(nSteps-equilibrationSteps);
meanSquareTemp = meanSquareTemp/(nSteps-equilibrationSteps);
meanPressure = meanPressure/(nSteps-equilibrationSteps);
meanSquarePressure = meanSquarePressure/(nSteps-equilibrationSteps);

temp = nSteps-correlationDistance-equilibrationSteps;
//The last step to calculate the mean values of meanTemp_ik & meanPressure_ik.
for(i = 0; i<correlationDistance; i++){
	meanTemp_ik[i] = meanTemp_ik[i]/temp;
	meanPressure_ik[i] = meanPressure_ik[i]/temp;
}

//Calculates phi
meanTempSquare = meanTemp*meanTemp;
meanPressureSquare = meanPressure*meanPressure;
for(i = 0; i<correlationDistance; i++){
	phiTemp[i] = (meanTemp_ik[i] - meanTempSquare)/(meanSquareTemp - meanTempSquare);
	phiPressure[i] = (meanPressure_ik[i] - meanPressureSquare)/(meanSquarePressure - meanPressureSquare);

	printf("k=%d   -->   phi T = %e   phi P = %e \n", i, phiTemp[i], phiPressure[i]); 
}


		// Savesthevalues needed to calculate the statistical inefficiency (s), for T and P.
		
/*		
if(i % 5 == 0){
	printf("it %d T %e \t alphaT %e \t P %e \t alphaP %e \n", i, currentTemp, alphaT, currentPressure, alphaP);
}

		fprintf(energyFile, "%e \t %e \t %e \t %e\n", i*timestep, energy, potentialEnergy, kineticEnergy );
		fprintf(ptFile, "%e \t %e \t %e \n", i*timestep, currentPressure, currentTemp);
	}

printf("lattpar efter loopen: %e \n", latticeParameter);
*/

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
// >>>>>>> 5f87f93470cfda3cd81ecf11fe38acdcd47701cd
}
