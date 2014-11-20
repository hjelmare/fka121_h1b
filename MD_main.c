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

#define PI 3.14159265368979

/* Main program */
int main()
{
	// REMEMBER: \m/ METAL UNITS \m/
	// simulation settings
	double productionTime = 5;
	double equilibrationTime = 5;
	double timestep = 0.01;
//	timestep = 0.1;
//	timestep = 0.001;
		
	
	int nSpectrumPoints = 1000;
	double spectrumInterval = PI;
	
	double maxCorrelationTime = 1;
	
	// physical parameters
	int dim = 3;		// do not change. some functions are hardcoded for dim = 3.
	int nCells = 4;
	int nParticles = 4*pow(nCells,dim);
	int MTemp, MPressure;    // The values at which the points of Temp and Pressure no longer are correlated.
	double wantedTemp = 500+273;   // The temperature that we want the system to stabilize around.
	wantedTemp = 700+273;
//	wantedTemp = 900+273;
	double wantedPressure =  6.32420934* 0.0000001;	// The pressure that we want the system to stabilize around.
	double timeConstantT = 0.02;
	double timeConstantP = 0.05;
	//double timeConstantT = 1;
	//double timeConstantP = 1;
	double mass = 0.00279636;  // 26.9815 u
	double latticeParameter = 4.05;
	double maxDeviation = 0.05;

	// derived quantities
	double nEquilibrationSteps = equilibrationTime/timestep;
	double nProductionSteps = productionTime/timestep;
	int equilibrationSteps = equilibrationTime/timestep;
	int maxCorrelationSteps = maxCorrelationTime/timestep;


	// storage of physical quantities
	double pos[nParticles][dim];
	double vel[nParticles][dim];
	double force[nParticles][dim];
	double savedValuesT[maxCorrelationSteps];
	double meanValuesT[maxCorrelationSteps];
	double savedValuesP[maxCorrelationSteps];
	double meanValuesP[maxCorrelationSteps];
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
	double nAverageSteps;
	double meanTemp, meanSquareTemp, meanTemp_ik[maxCorrelationSteps-1];
	double meanPressure, meanSquarePressure, meanPressure_ik[maxCorrelationSteps-1];
	double phiTemp[maxCorrelationSteps], phiPressure[maxCorrelationSteps], meanTempSquare, meanPressureSquare;
	double varTemp, varPressure;
	double sTemp, sPressure;
	double savedPos[maxCorrelationSteps][nParticles][dim];
	double msd;
	double savedVelocities[maxCorrelationSteps][nParticles][dim];
	double meanVelocityScalar[maxCorrelationSteps][nParticles];
	double meanVelocityAverage[maxCorrelationSteps];
	double spectrum[nSpectrumPoints];
	double diffusionCoefficient;
	double meanDistanceScalar[maxCorrelationSteps][nParticles];
	double meanDistanceAverage[maxCorrelationSteps];

	// other stuff
	int i,j,k,m;   

	srand(time(NULL));
	double random_value;

	// End of variable declarations, start of actual code

	printf("Initializing");
	init_fcc(pos, nCells, latticeParameter);

	// applying random perturbations from strict fcc positions
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
	
	// calculating initial energies
	potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, nParticles);
	kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
	energy = potentialEnergy + kineticEnergy;


	// saving initial energies on the first line in the file
	FILE *energyFile;
	energyFile = fopen("energy.data","w");
	//energyFile = fopen("energyT1.data","w");
	//energyFile = fopen("energyT2.data","w");
	//energyFile = fopen("energyT3.data","w");
	fprintf(energyFile, "%e \t %e \t %e \t %e \n", 0.0, energy, potentialEnergy, kineticEnergy);

	FILE *ptFile;
	ptFile = fopen("pt.data","w");

	FILE *positionFile;
	positionFile = fopen("position.data","w");
//	positionFile = fopen("position1.data","w");
//	positionFile = fopen("position2.data","w");
//	positionFile = fopen("position3.data","w");


// Initialization is done, moving on to equilibration

	printf("\t");	// Progress indicator stuff
	for (i = 0; i < equilibrationSteps/100; i++) {
		printf(".");
	}
	printf("\tDone!\nEquilibration\t");

	for (i=0;i<nEquilibrationSteps;i++) {		// Start of equilibration loop
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

		// Calculating energies
		potentialEnergy = get_energy_AL(pos, nCells*latticeParameter, nParticles);
		kineticEnergy = GetKineticEnergy(vel, mass, nParticles);
		energy = potentialEnergy + kineticEnergy;
		fprintf(energyFile,"%e \t %e \t %e \t %e \n", (i+1)*timestep, energy, potentialEnergy, kineticEnergy);		

		// Get temperature and pressure
		currentTemp = GetInstantTemperature(vel, nParticles, mass, dim);
		volume = pow(nCells*latticeParameter, 3);
		virial = get_virial_AL(pos, nCells*latticeParameter, nParticles);
		currentPressure = GetPressure(currentTemp, volume, virial, nParticles);
		fprintf(ptFile, "%e \t %e \t %e \n", (i+1)*timestep, currentTemp, currentPressure);

		// Calculate scaling parameters for equilibration
		alphaT = GetAlphaT(wantedTemp, currentTemp, timestep, timeConstantT);
		alphaP = GetAlphaP(wantedPressure, currentPressure, timestep, timeConstantP);
		sqrtAlphaT = sqrt(alphaT);
		curtAlphaP = pow(alphaP, 1.0 / 3);

		// Rescale velocity and position
		latticeParameter = latticeParameter * curtAlphaP;
		for (j=0; j<nParticles; j++) {
			for (k=0; k<dim; k++) {
				vel[j][k] = vel[j][k] * sqrtAlphaT;
				pos[j][k] = pos[j][k] * curtAlphaP;
			}
		}
		if( i % 100 == 0) { 	// Progress indicator
			printf("I");
			fflush(stdout);
		}
	}	// End of equilibration loop

	printf("\tDone!\nProduction\t");
// DEBUGGING STUFF
		FILE *pt2File;
		pt2File = fopen("pt2.data","w");
// ---------
	for (i=0;i<nProductionSteps;i++) {
		//Save positions of one particle to check whether the aluminium is solid or liquid
		fprintf(positionFile, "%d \t %e \t %e \t %e \n", i, pos[0][0], pos[0][1], pos[0][2]);

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

		fprintf(pt2File, "%d\t%e\t%e\n", i, currentPressure, currentTemp);

		//Saves temp and pressure values in order to calculate s.
		if (i < maxCorrelationSteps) {
			savedValuesT[i] = currentTemp;
			savedValuesP[i] = currentPressure;		
		} else {  
			// updates the saved values.
			for (j = 0; j < maxCorrelationSteps - 1; j++){	
				savedValuesT[j] = savedValuesT[j+1];
				savedValuesP[j] = savedValuesP[j+1];	
			}
			savedValuesT[maxCorrelationSteps - 1] = currentTemp;
			savedValuesP[maxCorrelationSteps - 1] = currentPressure;
			
			// calculates and saves f_i*f_k
			for (j = 0; j<maxCorrelationSteps; j++){
				meanTemp_ik[j] += savedValuesT[maxCorrelationSteps - 1] * savedValuesT[maxCorrelationSteps - 1 - j];
				meanPressure_ik[j] += savedValuesP[maxCorrelationSteps - 1]* savedValuesP[maxCorrelationSteps - 1 - j];
			//	printf("%e \t",meanTemp_ik[j]);
			}
			//printf("\n");
			meanTemp += currentTemp;
			meanPressure += currentPressure;
			meanSquareTemp += currentTemp * currentTemp;
			meanSquarePressure += currentPressure * currentPressure;
		}

		// Saving data for the velocity correlation function
		if (i < maxCorrelationSteps) {
			for (m = 0; m < nParticles ; m++) {
				savedVelocities[i][m][0] = vel[m][0];
				savedVelocities[i][m][1] = vel[m][1];
				savedVelocities[i][m][2] = vel[m][2];
			}
		} else {
			for (j = 0; j< maxCorrelationSteps- 1; j++) {
				for (m = 0; m<nParticles;  m++) {
					savedVelocities[j][m][0] = savedVelocities[j+1][m][0];
					savedVelocities[j][m][1] = savedVelocities[j+1][m][1];
					savedVelocities[j][m][2] = savedVelocities[j+1][m][2];
				}
			}
			for (m = 0; m<nParticles; m++) {
				savedVelocities[maxCorrelationSteps - 1][m][0] = vel[m][0];
				savedVelocities[maxCorrelationSteps - 1][m][1] = vel[m][1];
				savedVelocities[maxCorrelationSteps - 1][m][2] = vel[m][2];
			}

			for (j = 0; j<maxCorrelationSteps; j++) {
				for (m = 0; m < nParticles; m++) {
					meanVelocityScalar[j][m] += ScalarProduct(savedVelocities[0][m],savedVelocities[j][m]);
				}
			}
		}

		if(i <  maxCorrelationSteps) {
			for(j = 0; j<nParticles; j++) {
				for (k = 0; k<dim; k++) {
					savedPos[i][j][k] = pos[j][k];
				} 
			}
		}else {
			for (j = 0; j<maxCorrelationSteps-1; j++) {
				for (k = 0; k<nParticles; k++) {
					for (m = 0; m<dim; m++) {
						savedPos[j][k][m] = savedPos[j+1][k][m];
					}
				}
			}
			for (j = 0; j<nParticles; j++) {
				for (k=0; k<dim; k++) {
					savedPos[maxCorrelationSteps-1][j][k] = pos[j][k];
				}
			}
			for(j = 0; j<maxCorrelationSteps; j++){
				for(k = 0; k<nParticles; k++){
					meanDistanceScalar[j][k] += getDistanceSquared(savedPos[0][k], savedPos[j][k]); 
				}
			}
		} 


		if(i%100 == 0) {	// Progress indicator
			printf("I");
			fflush(stdout);
		}

	}	// End of production loop

	printf("\tDone!\nCleaning up\n");

	// The production part of the simulation is now done.
	// The code below forms some averages and other stuff that
	// is best done with all data on hand

	nAverageSteps = nProductionSteps-maxCorrelationSteps;	// number of steps over which to take averages
	
	meanTemp = meanTemp/nAverageSteps;
	meanSquareTemp = meanSquareTemp/nAverageSteps;
	meanPressure = meanPressure/nAverageSteps;
	meanSquarePressure = meanSquarePressure/nAverageSteps;
printf("%e \t %e \t %e \t %e \n", meanTemp, meanSquareTemp, meanPressure, meanSquarePressure);
//printf("test = %e \t calculated= %e \n", test, temp);

	//The last step to calculate the mean values of meanTemp_ik & meanPressure_ik.
	for(i = 0; i < maxCorrelationSteps+1; i++){
		meanTemp_ik[i] = meanTemp_ik[i]/nAverageSteps;
		meanPressure_ik[i] = meanPressure_ik[i]/nAverageSteps;
	}
	
	msd = msd/nAverageSteps;
	printf("msd= %e DETTA ÄR INTE KLART\n", msd);

	FILE *msdFile;
	msdFile = fopen("msd3.data","w");

	for( i = 0; i < maxCorrelationSteps; i++) {
		meanDistanceAverage[i] = 0;
		for (j = 0; j<nParticles; j++) {
			meanDistanceScalar[i][j] = meanDistanceScalar[i][j]/nAverageSteps;
			meanDistanceAverage[i] += meanDistanceScalar[i][j];
		}
		meanDistanceAverage[i] = meanDistanceAverage[i]/nParticles;
		fprintf(msdFile, "%e\t%e\n", i*timestep, meanDistanceAverage[i]);
	}

	// Final processing and saving of velocity correlation function
	FILE *velcorFile;
	velcorFile = fopen("velcor.data","w");

	diffusionCoefficient = 0;

	for( i = 0; i < maxCorrelationSteps; i++) {
		meanVelocityAverage[i] = 0;
		for (j = 0; j<nParticles; j++) {
			meanVelocityScalar[i][j] = meanVelocityScalar[i][j]/nAverageSteps;
			meanVelocityAverage[i] += meanVelocityScalar[i][j];
		}
		meanVelocityAverage[i] = meanVelocityAverage[i]/nParticles;
		fprintf(velcorFile, "%e\t%e\n", i*timestep, meanVelocityAverage[i]);
		
		diffusionCoefficient += meanVelocityAverage[i];
	}

	diffusionCoefficient = (1.0/3.0) * diffusionCoefficient / maxCorrelationSteps;
//	printf("diffCoeff = %e \n", diffusionCoefficient);

	// Calculating and saving spectrum integral thingy, finding diffusion coeff from time integral
	FILE *spectrumFile;
	spectrumFile = fopen("spectrum.data","w");

	for( i = 0; i < nSpectrumPoints; i++) {
		spectrum[i] = 0;
		for( j = 0; j < maxCorrelationSteps; j++) {
			spectrum[i] += meanVelocityAverage[j] * cos(spectrumInterval * i * j / nSpectrumPoints);
		}
		spectrum[i] = 2 * spectrum[i]/maxCorrelationSteps;
		fprintf(spectrumFile,"%e \t %e \n",i*spectrumInterval/nSpectrumPoints,spectrum[i]);
	}

	//Calculates phi
	meanTempSquare = meanTemp*meanTemp;
	meanPressureSquare = meanPressure*meanPressure;

	// sparar värdena så att jag kan plotta dem och se om de blir vettiga.
	FILE *phiTempFile;
	phiTempFile = fopen("phiTemp.data","w");

	FILE *phiPFile;
	phiPFile = fopen("phiPressure.data","w");
	for(i = 0; i<maxCorrelationSteps; i++){
		phiTemp[i] = (meanTemp_ik[i] - meanTempSquare)/(meanSquareTemp - meanTempSquare);
		phiPressure[i] = (meanPressure_ik[i] - meanPressureSquare)/(meanSquarePressure - meanPressureSquare);
		fprintf(phiTempFile, "%e \t %e \n", (i+1)*timestep, phiTemp[i]); //Dessa behövs inte. De används bara för att kolla i matlab att vi får något vettigt!!!
		fprintf(phiPFile, "%e \t %e \n", (i+1)*timestep, phiPressure[i]);
	}

	// Finding M (the longest time for which correlation is above some limit)
	double limit = exp(-2.0);
	for(i = 0; i < maxCorrelationSteps; i++){
		MTemp  = i;
		if( phiTemp[i] < limit){
			i = maxCorrelationSteps;
		}
	}


	for(i = 0; i<maxCorrelationSteps; i++){
		MPressure  = i+1;
		if( phiPressure[i] < limit){
			i = maxCorrelationSteps;
		}
	}
	

	//Sum over phi to get s - Temperature
	for(i = 0; i<MTemp; i++){
		sTemp += phiTemp[i];
	}
	sTemp = sTemp*2; // phi is symmetric and we sum from 0 to M not -M to M as it should.


	//Sum over phi to get s . Pressure
	for(i = 0; i<MPressure; i++){
		sPressure += phiPressure[i];
	}
	sPressure = sPressure*2; // phi is symmetric and we sum from 0 to M not -M to M as it should.

	varTemp = meanSquareTemp - meanTempSquare;
	varPressure = meanSquarePressure - meanPressureSquare;

	varTemp = varTemp*sTemp/nProductionSteps;
	varPressure = varPressure*sPressure/nProductionSteps;
printf("sTemp= %e \t sPressure= %e \nVAR[temp]= %e \t VAR[P]= %e \n",sTemp, sPressure, varTemp, varPressure);
	
	FILE *valuesFile;
	valuesFile = fopen("values.data","w");
	fprintf(valuesFile,"Ds-int\n");
	fprintf(valuesFile,"%e \t",diffusionCoefficient);
	fprintf(valuesFile,"\n");
	close(valuesFile);

	FILE *donefile;
	donefile = fopen("done.data", "w");
	fprintf(donefile, "done");
	close(donefile);

	printf("Done!\n");
	return 0;
}
